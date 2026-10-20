from .MatchedSpectra import MatchedSpectra
from .MgfInstance import MgfInstance
from .Util import path_gnps, Parametres, get_correct_inchi, _mass_selection_by_tolerance, _get_match
from matchms.importing import load_from_mgf
from matchms import Spectrum
import json
import requests
import chardet
import sys
import pandas as pd
from pathlib import Path
import ftplib
import time
import ssl
import socket
import os


def _gnps_user_path(username, remote_path):
    """
    Build the GNPS path of a file uploaded by FTP in the user space.

    GNPS displays such files as 'f.<user>/<folder>/<file>' and its file selectors
    store them with a trailing semicolon. The prefix can be overridden with the
    environment variable MS2DECIDE_GNPS_PREFIX (default: 'f.').

    Args:
        username (str): GNPS user name.
        remote_path (str): Path returned by the FTP upload (e.g. '/folder/file.mgf').

    Returns:
        str: Path usable in 'spec_on_server' / 'quantification_table'.
    """
    prefix = os.getenv("MS2DECIDE_GNPS_PREFIX", "f.")
    rel = str(remote_path).replace("\\", "/").lstrip("/")
    # Avoid a duplicated user name if the FTP root already contains it.
    if rel.startswith(username + "/"):
        rel = rel[len(username) + 1:]
    return f"{prefix}{username}/{rel};"


def _invoke_workflow(auth, base_url, parameters):
    """
    Log in to a specific URL using username and password to launch a workflow with given parameters.

    Args:
        auth (Auth): Authentication object with username and password.
        base_url (str): Base URL to log in (e.g., 'gnps.ucsd.edu').
        parameters (dict): Parameters of the workflow.

    Returns:
        str: Task ID of the launched workflow.

    Raises:
        RuntimeError: if the GNPS task cannot be created.
    """

    username = auth.username
    password = auth.password

    verify_https = str(os.getenv("MS2DECIDE_HTTPS_VERIFY_TLS", "0")).lower() in (
        "1", "true", "yes"
    )

    if not verify_https:
        try:
            from urllib3.exceptions import InsecureRequestWarning
            requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)
        except Exception:
            pass

    s = requests.Session()

    payload = {
        'user': username,
        'password': password,
        'login': 'Sign in'
    }

    r = s.post(
        'https://' + base_url + '/ProteoSAFe/user/login.jsp',
        data=payload,
        verify=verify_https,
    )
    if r.status_code != 200:
        raise PermissionError(
            f"Connexion GNPS refusee (HTTP {r.status_code}).")

    # No automatic retry here: a retry after a late failure could create duplicate jobs.
    r = s.post('https://' + base_url + '/ProteoSAFe/InvokeTools',
               data=parameters, verify=verify_https)
    task_id = r.text.strip()

    if 4 < len(task_id) < 60 and "error" not in task_id.lower():
        print("Launched Task: : " + task_id)
        return task_id

    print(task_id, file=sys.stderr, flush=True)
    raise RuntimeError(
        "Echec de la creation de la tache GNPS (HTTP "
        + str(r.status_code) + "). Reponse du serveur :\n" + task_id)


def _upload_to_gnps(auth, input_file_mgf, input_file_quan, folder):
    username = auth.username
    password = auth.password
    url = "massive-ftp.ucsd.edu"

    # GNPS/MassIVE FTPS can expose certificate-chain issues on some Conda/Windows
    # stacks. Keep encryption but disable certificate verification by default to
    # match historical ftplib FTP_TLS behavior and maintain compatibility.
    verify_tls = str(os.getenv("MS2DECIDE_FTPS_VERIFY_TLS", "0")).lower() in (
        "1", "true", "yes"
    )
    ssl_context = ssl.create_default_context() if verify_tls else ssl._create_unverified_context()
    if hasattr(ssl, "TLSVersion"):
        ssl_context.minimum_version = ssl.TLSVersion.TLSv1_2
        ssl_context.maximum_version = ssl.TLSVersion.TLSv1_2

    connection_errors = (
        ConnectionResetError,
        ConnectionAbortedError,
        BrokenPipeError,
        TimeoutError,
        EOFError,
        ssl.SSLError,
        OSError,
        ftplib.error_temp,
        ftplib.error_proto,
    )

    folder = str(folder).strip().strip("/\\")
    if not folder:
        raise ValueError("Le nom du dossier distant ne peut pas être vide.")

    debug_ftp = str(os.getenv("MS2DECIDE_FTPS_DEBUG", "0")).lower() in (
        "1", "true", "yes"
    )

    class _FTP_TLS_SessionReuse(ftplib.FTP_TLS):
        """FTPS client that reuses control-channel TLS session for data channels.

        Some FTPS servers close data connections unless TLS session reuse is enabled.
        """

        def ntransfercmd(self, cmd, rest=None):
            conn, size = ftplib.FTP.ntransfercmd(self, cmd, rest)
            if self._prot_p:
                wrap_kwargs = {"server_hostname": self.host}
                session = getattr(self.sock, "session", None)
                if session is not None:
                    wrap_kwargs["session"] = session
                conn = self.context.wrap_socket(conn, **wrap_kwargs)
            return conn, size

    def _connect(strategy):
        ftp_cls = _FTP_TLS_SessionReuse if strategy["session_reuse"] else ftplib.FTP_TLS
        ftp_conn = ftp_cls(context=ssl_context)
        if debug_ftp:
            ftp_conn.set_debuglevel(2)
        ftp_conn.connect(url, 21, timeout=120)
        ftp_conn.af = socket.AF_INET
        ftp_conn.login(username, password)
        if strategy["private_data_channel"]:
            ftp_conn.prot_p()
        else:
            ftp_conn.prot_c()
        ftp_conn.set_pasv(strategy["passive"])
        ftp_conn.voidcmd("TYPE I")
        try:
            ftp_conn.sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
        except OSError:
            pass
        return ftp_conn

    def _safe_close(ftp_conn):
        if ftp_conn is None:
            return
        try:
            ftp_conn.quit()
        except Exception:
            try:
                ftp_conn.close()
            except Exception:
                pass

    def _prepare_remote_folder(ftp_conn):
        print("Current remote directory:", ftp_conn.pwd())
        # Try to enter the user directory when the FTP server exposes it.
        for candidate in (username, f"/{username}"):
            try:
                ftp_conn.cwd(candidate)
                print("Entered user directory:", ftp_conn.pwd())
                break
            except ftplib.error_perm:
                continue
        else:
            print("User directory not accessible directly, staying in:", ftp_conn.pwd())

        try:
            ftp_conn.mkd(folder)
            print("Folder created")
        except ftplib.error_perm as e:
            print("MKD skipped:", e)

        try:
            ftp_conn.cwd(folder)
        except ftplib.error_perm as e:
            raise Exception(f"Impossible d'accéder au dossier '{folder}': {e}")

        print("Working remote directory:", ftp_conn.pwd())

    def _upload_one_file(ftp_conn, local_path, remote_name, strategy, retries=4):
        for attempt in range(1, retries + 1):
            try:
                with open(str(local_path), "rb") as file_in:
                    ftp_conn.storbinary(
                        f"STOR {remote_name}",
                        file_in,
                        blocksize=256 * 1024,
                    )
                return ftp_conn
            except connection_errors as exc:
                if attempt == retries:
                    raise
                print(
                    f"Upload failed for {remote_name} (attempt {attempt}/{retries}): {exc}. Reconnecting..."
                )
                _safe_close(ftp_conn)
                time.sleep(min(attempt * 2, 8))
                ftp_conn = _connect(strategy)
                _prepare_remote_folder(ftp_conn)

        return ftp_conn

    def _check_remote_files(ftp_conn, names):
        """Verify that the uploaded files are listed in the remote folder."""
        try:
            listed = [n.replace("\\", "/").split("/")[-1] for n in ftp_conn.nlst()]
        except ftplib.error_perm:
            # An empty listing is reported as an error by some servers.
            listed = []
        missing = [n for n in names if n not in listed]
        if missing:
            raise ConnectionError(
                f"Fichiers absents du dossier distant apres upload : {missing}")

    upload_strategies = [
        {
            "name": "FTPS passive (PROT P + session reuse)",
            "passive": True,
            "private_data_channel": True,
            "session_reuse": True,
        },
        {
            "name": "FTPS passive (PROT P)",
            "passive": True,
            "private_data_channel": True,
            "session_reuse": False,
        },
        {
            "name": "FTPS active (PROT P)",
            "passive": False,
            "private_data_channel": True,
            "session_reuse": False,
        },
        {
            "name": "FTPS passive (PROT C)",
            "passive": True,
            "private_data_channel": False,
            "session_reuse": False,
        },
    ]

    last_error = None

    def _safe_remote_name(name):
        # Avoid FTP path parsing issues with whitespace and control chars.
        cleaned = "".join(ch if ch.isprintable() else "_" for ch in name)
        cleaned = "_".join(cleaned.split())
        return cleaned or "upload_file"

    remote_mgf = _safe_remote_name(Path(input_file_mgf).name)
    remote_csv = _safe_remote_name(Path(input_file_quan).name)

    for strategy in upload_strategies:
        ftp = None
        try:
            print(f"Trying upload strategy: {strategy['name']}")
            ftp = _connect(strategy)
            print("Login successful")
            _prepare_remote_folder(ftp)

            ftp = _upload_one_file(ftp, input_file_mgf, remote_mgf, strategy)
            path_mgf = ftp.pwd() + "/" + remote_mgf
            print("MGF uploaded:", path_mgf)

            ftp = _upload_one_file(ftp, input_file_quan, remote_csv, strategy)
            path_csv = ftp.pwd() + "/" + remote_csv
            print("CSV uploaded:", path_csv)

            _check_remote_files(ftp, [remote_mgf, remote_csv])

            print("Files uploaded")
            return (path_mgf, path_csv)
        except connection_errors as exc:
            last_error = exc
            print(f"Strategy failed: {strategy['name']} -> {exc}")
        finally:
            _safe_close(ftp)

    raise ConnectionError(
        f"Impossible de televerser les fichiers sur GNPS apres plusieurs strategies FTPS. Derniere erreur: {last_error}"
    )


def _get_networking_parameters():
    """
    Return the FBMN parameters for GNPS.

    Returns:
        dict: Molecular network parameters.
    """
    invokeParameters = {}
    invokeParameters["workflow"] = "FEATURE-BASED-MOLECULAR-NETWORKING"
    invokeParameters["protocol"] = "None"
    invokeParameters["library_on_server"] = "d.speclibs;"

    # Networking
    invokeParameters["tolerance.PM_tolerance"] = "0.02"
    invokeParameters["tolerance.Ion_tolerance"] = "0.02"
    invokeParameters["PAIRS_MIN_COSINE"] = "0.70"
    invokeParameters["MIN_MATCHED_PEAKS"] = "6"
    invokeParameters["TOPK"] = "10"
    invokeParameters["MAX_SHIFT"] = "500"

    # Network Pruning
    invokeParameters["MAXIMUM_COMPONENT_SIZE"] = "100"

    # Library Search
    invokeParameters["MIN_MATCHED_PEAKS_SEARCH"] = "6"
    invokeParameters["SCORE_THRESHOLD"] = "0.001"
    invokeParameters["TOP_K_RESULTS"] = "1"
    invokeParameters["ANALOG_SEARCH"] = "1"
    invokeParameters["MAX_SHIFT_MASS"] = "100.0"
    invokeParameters["FILTER_STDDEV_PEAK_datasetsINT"] = "0.0"
    invokeParameters["MIN_PEAK_INT"] = "0.0"
    invokeParameters["FILTER_PRECURSOR_WINDOW"] = "0"
    invokeParameters["FILTER_LIBRARY"] = "0"
    invokeParameters["WINDOW_FILTER"] = "0"

    # Quant
    invokeParameters["GROUP_COUNT_AGGREGATE_METHOD"] = "Mean"
    invokeParameters["QUANT_FILE_NORM"] = "RowSum"
    invokeParameters["QIIME2_PCOA_DISTANCE"] = 'cosine'
    invokeParameters["QUANT_TABLE_SOURCE"] = "MZMINE2"
    invokeParameters["RUN_STATS"] = "No"
    invokeParameters["RUN_DEREPLICATOR"] = "0"
    return invokeParameters


def _get_iterative_parameters():
    """
    Return iterative parameters for peak matching in GNPS.

    Returns:
        dict: Iterative parameters for peak matching.
    """
    N1 = 1
    N2 = round(N1-0.025*12, 3)
    N3 = round(N2-0.025*10, 3)

    alpha = {6: {0.02: round(N1, 3),
                 0.1: round(N1-0.025*1, 3),
                 10: round(N1-0.025*2, 3),
                 25: round(N1-0.025*4, 3),
                 50: round(N1-0.025*6, 3),
                 100: round(N1-0.025*8, 3),
                 250: round(N1-0.025*10, 3),
                 500: round(N1-0.025*12, 3),
             0: round(N1-0.025*14, 3)},
             5: {0.02: round(N2, 3),
                 0.1: round(N2-0.025*1, 3),
                 10: round(N2-0.025*2, 3),
                 25: round(N2-0.025*4, 3),
                 50: round(N2-0.025*6, 3),
                 100: round(N2-0.025*8, 3),
                 250: round(N2-0.025*10, 3),
                 500: round(N2-0.025*12, 3),
             0: round(N2-0.025*14, 3)},
             4: {0.02: N3,
                 0.1: round(N3-0.025*1, 3),
                 10: round(N3-0.025*2, 3),
                 25: round(N3-0.025*4, 3),
                 50: round(N3-0.025*6, 3),
                 100: round(N3-0.025*8, 3),
                 250: round(N3-0.025*10, 3),
                 500: round(N3-0.025*12, 3),
             0: round(N3-0.025*14, 3)}}
    return alpha


def _wait_after_upload():
    """Short pause so that GNPS can index the files uploaded by FTP."""
    delay = int(os.getenv("MS2DECIDE_GNPS_DELAY", "5"))
    if delay > 0:
        print(f"Waiting {delay}s for GNPS to index the uploaded files...")
        time.sleep(delay)


def _launch_GNPS_workflow(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description, score_threshold):
    """
    Log in to GNPS using username and password, upload file to MassIVE repository, and launch feature-based molecular networking workflow.

    Args:
        auth (Auth): Authentication object with username, password, and email.
        path_file_mgf_in_gnps (str): Path of the .mgf file uploaded by FTP.
        path_file_quan_in_gnps (str): Path of the quantification table uploaded by FTP.
        job_description (str): Title of the job.
        score_threshold: Minimum cosine score for library search.

    Returns:
        str: Task ID of the launched workflow.
    """
    invokeParameters = _get_networking_parameters()
    invokeParameters["SCORE_THRESHOLD"] = score_threshold
    invokeParameters["desc"] = job_description
    invokeParameters["MAX_SHIFT_MASS"] = str(0.02)
    invokeParameters["quantification_table"] = _gnps_user_path(
        auth.username, path_file_quan_in_gnps)
    invokeParameters["spec_on_server"] = _gnps_user_path(
        auth.username, path_file_mgf_in_gnps)
    invokeParameters["email"] = auth.mail

    print("spec_on_server:", invokeParameters["spec_on_server"])
    print("quantification_table:", invokeParameters["quantification_table"])

    task_id = _invoke_workflow(auth, "gnps.ucsd.edu", invokeParameters)

    return task_id


def _launch_GNPS_workflow_iterative(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description):
    """
    Log in to GNPS, upload files to the MassIVE repository, and launch iterative weighted analog search jobs.

    Args:
        auth (Auth): Authentication object with username, password, and email.
        path_file_mgf_in_gnps (str): The path of the .mfg file uploaded to GNPS.
        path_file_quan_in_gnps (str): The path of the quantification table uploaded to GNPS.
        job_description (str): Title of the job.

    Returns:
        dict: A dictionary containing task IDs for different peak and mass combinations.
    """
    SCORE_THRESHOLD = input(
        "INPUT A SCORE_THRESHOLD FOR ITERATIVE GNPS. DEFAULT 0.001").strip()
    if not SCORE_THRESHOLD:
        SCORE_THRESHOLD = "0.001"
    path_quan = _gnps_user_path(auth.username, path_file_quan_in_gnps)
    path_mgf = _gnps_user_path(auth.username, path_file_mgf_in_gnps)
    D_job = {}
    for peak in [6, 5, 4]:
        d = {i: 0 for i in [0.02, 0.1, 10, 25, 50, 100, 250, 500, 0]}
        for mass_diff in [0.02, 0.1, 10, 25, 50, 100, 250, 500, 0]:
            invokeParameters = _get_networking_parameters()
            invokeParameters["SCORE_THRESHOLD"] = SCORE_THRESHOLD
            invokeParameters["MIN_MATCHED_PEAKS_SEARCH"] = str(peak)
            invokeParameters["MAX_SHIFT_MASS"] = str(mass_diff)
            invokeParameters["desc"] = job_description + \
                '_'+str(peak)+'_'+str(mass_diff)
            invokeParameters["quantification_table"] = path_quan
            invokeParameters["spec_on_server"] = path_mgf
            invokeParameters["email"] = auth.mail
            task_id = _invoke_workflow(auth, "gnps.ucsd.edu", invokeParameters)
            d[mass_diff] = task_id
        D_job[peak] = d
    return D_job


def _gnps_annotations_download_results(task_id, retries=20, wait=60):
    """
    Download GNPS molecular networking job results locally.

    GNPS sometimes does not answer (or is down for a long time), so the download
    is retried (default: 20 attempts, 60 s apart).

    Args:
        task_id (str): GNPS task ID.
        retries (int): Number of attempts.
        wait (int): Seconds to wait between two attempts.

    Returns:
        pd.DataFrame: Dataframe containing job annotations.
    """

    # Base link to download GNPS job annotations
    if (task_id[:3].lower() == 'id='):
        task_id = task_id[3:]
    link = ('https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp'
            '?task={}&view=view_all_annotations_DB').format(task_id)
    verify_https = str(os.getenv("MS2DECIDE_HTTPS_VERIFY_TLS", "0")).lower() in (
        "1", "true", "yes"
    )
    if not verify_https:
        try:
            from urllib3.exceptions import InsecureRequestWarning
            requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)
        except Exception:
            pass
    headers = {"User-Agent": "Mozilla/5.0 (MS2DECIDE)"}

    data = None
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(link, headers=headers,
                             verify=verify_https, timeout=120)
            r.raise_for_status()
            raw_data = r.content
            detected_encoding = chardet.detect(raw_data)
            decoded_data = raw_data.decode(
                detected_encoding['encoding'] or 'utf-8', errors='replace')
            data = json.loads(decoded_data)
            print("JSON Loaded Successfully")
            break
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            print(
                f"GNPS results unavailable (attempt {attempt}/{retries}): {exc}")
            if attempt < retries:
                time.sleep(wait)

    if data is None:
        raise RuntimeError(
            f"CAN NOT GET THE DATA FROM THE JOB: {task_id}. Last error: {last_error}\n"
            f"Check in a browser: {link}")

    df_annotations = pd.DataFrame.from_dict(data.get('blockData', []))
    if df_annotations.empty:
        print('No spectral library annotation found in the job:' + task_id)
        return pd.DataFrame({'#Scan#': pd.Series(dtype='int64'),
                             'MQScore': pd.Series(dtype='float64')})
    print('==================')
    print(' FBMN job detected with ' +
          str(df_annotations.shape[0]-1) + ' spectral library annotations in the job:' + task_id)
    print('==================')
    df_annotations[['#Scan#', 'MQScore']] = df_annotations[[
        '#Scan#', 'MQScore']].apply(pd.to_numeric)
    return (df_annotations)


def _workflows_statues(dict_task_id, retries=5, wait=30):
    """
    Return the set of GNPS statuses of all given tasks.

    GNPS sometimes answers with an empty or non-JSON page, so each request is retried.
    Raises RuntimeError if a status cannot be read after all attempts.
    """
    state = set()
    url = 'https://gnps.ucsd.edu/ProteoSAFe/status_json.jsp?task={}'
    verify_https = str(os.getenv("MS2DECIDE_HTTPS_VERIFY_TLS", "0")).lower() in (
        "1", "true", "yes"
    )
    if not verify_https:
        try:
            from urllib3.exceptions import InsecureRequestWarning
            requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)
        except Exception:
            pass

    for peak in dict_task_id:
        for mass in dict_task_id[peak]:
            task = dict_task_id[peak][mass]
            last_error = None
            for attempt in range(1, retries + 1):
                try:
                    r = requests.get(
                        url.format(task), verify=verify_https, timeout=60)
                    r.raise_for_status()
                    state.add(json.loads(r.text)["status"])
                    break
                except (requests.RequestException, ValueError, KeyError) as exc:
                    last_error = exc
                    print(
                        f"GNPS status unavailable (attempt {attempt}/{retries}): {exc}")
                    if attempt < retries:
                        time.sleep(wait)
            else:
                raise RuntimeError(
                    f"Impossible de lire le statut du job {task}: {last_error}")
    return (state)


def _wait_for_workflow_finish(dict_task_id, poll_interval=60, max_wait_seconds=None,
                              max_outage_seconds=None):
    """
    Wait until all the given GNPS jobs are DONE.

    A GNPS outage does not stop the wait: the status is polled again every
    poll_interval seconds, and the function only gives up after max_outage_seconds
    of continuous unavailability (default 24 h, or env MS2DECIDE_GNPS_MAX_OUTAGE).
    """
    if poll_interval <= 0:
        raise ValueError("poll_interval doit etre > 0.")
    if max_outage_seconds is None:
        max_outage_seconds = int(os.getenv("MS2DECIDE_GNPS_MAX_OUTAGE", str(24 * 3600)))

    start_time = time.time()
    outage_start = None
    last_printed = None
    while True:
        try:
            statues = _workflows_statues(dict_task_id, retries=1, wait=0)
            outage_start = None
        except RuntimeError:
            if outage_start is None:
                outage_start = time.time()
            down = time.time() - outage_start
            if down >= max_outage_seconds:
                raise
            print(f"GNPS unreachable for {int(down // 60)} min, still trying...")
            time.sleep(poll_interval)
            continue

        if statues != last_printed:
            print(f"Current GNPS status: {statues}")
            last_printed = statues
        if statues == {'DONE'}:
            print('all Jobs done')
            return statues
        if statues & {'FAILED', 'SUSPENDED'}:
            raise RuntimeError(
                f"Un job GNPS a echoue ou a ete suspendu. Statut: {statues}")
        if max_wait_seconds is not None and (time.time() - start_time) >= max_wait_seconds:
            raise TimeoutError(
                f"Timeout atteint apres {max_wait_seconds} secondes. Dernier statut: {statues}"
            )
        time.sleep(poll_interval)


def closest_gnps(auth, input_file_mgf, input_file_quan, folder, score_threshold):
    """
    Returns a dictionary of MatchedSpectra objects for the closest matches in GNPS.

    Args:
        auth (Auth): Authentication object with username, password, and email.
        mgf_path (str): The local path of the .mfg file.

    Returns:
        dict: Dictionary {id: MatchedSpectra}.
    """
    mgf_instance = MgfInstance(Path(input_file_mgf))
    dict_data = mgf_instance.data
    path_file_mgf_in_gnps, path_file_quan_in_gnps = _upload_to_gnps(
        auth, input_file_mgf, input_file_quan, folder)
    _wait_after_upload()
    job_description = Path(input_file_mgf).stem
    task_id = _launch_GNPS_workflow(
        auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description, score_threshold)
    _wait_for_workflow_finish({0: {0: task_id}})
    df_anno = _gnps_annotations_download_results(task_id)
    df_anno = df_anno.rename(columns=str.lower)
    dic_anno = df_anno.set_index("#scan#").to_dict('index')
    dic_anno = {int(i): dic_anno[i] for i in dic_anno}
    matched_spectra_dict = {}
    for ID in dict_data:
        if (ID in dic_anno):
            matched_spectra_dict[ID] = MatchedSpectra(
                dict_data[ID], get_correct_inchi(dic_anno[ID]), dic_anno[ID]['mqscore'])
        else:
            matched_spectra_dict[ID] = MatchedSpectra(dict_data[ID], "#", 0)
    return (matched_spectra_dict)


def closest_gnps_local(mgf):
    """
    Returns a dictionary of MatchedSpectra objects for the closest matches in GNPS using local processing.

    Args:
        mgf : mgfInstance

    Returns:
        dict: Dictionary {id: MatchedSpectra}.
    """
    gnps_mass, gnps = _load_gnps()
    tolerance, mz_power, intensity_power, shift = Parametres()
    print('==================')
    print('==================')
    RES = {}
    for i in mgf.data:
        sp = mgf.data[i]
        selected = _mass_selection_by_tolerance(sp, gnps_mass, gnps, tolerance)
        rsp, res = _get_match(sp, selected, tolerance,
                              mz_power, intensity_power, shift)
        if (type(rsp) == Spectrum):
            RES[i] = MatchedSpectra(i, get_correct_inchi(rsp), res[0])
        else:
            RES[i] = MatchedSpectra(i, rsp, res)
    return RES


def _load_gnps():
    """
    Load GNPS data from an MGF file.

    This function retrieves the path to the GNPS data file, loads the data, and extracts
    precursor mass information along with the GNPS data objects. It returns a list of
    precursor masses and a list of GNPS data objects.

    Returns:
        tuple: A tuple containing two lists:
            - mass (list): A list of precursor m/z values.
            - gnps (list): A list of GNPS data objects.
    """
    path_data_gnps = path_gnps()
    gnps_data = load_from_mgf(path_data_gnps)
    mass = []
    gnps = []
    for i in gnps_data:
        mass.append(i.metadata['precursor_mz'])
        gnps.append(i)
    return mass, gnps


def closest_gnps_by_id(job_id, mgf_path):
    """
    Returns a dictionary of MatchedSpectra objects for the closest matches in GNPS.

    Args:
        job_id (str): GNPS task ID.
        mgf_path (str): The local path of the .mfg file.

    Returns:
        dict: Dictionary {id: MatchedSpectra}.
    """
    df_anno = _gnps_annotations_download_results(job_id)
    df_anno = df_anno.rename(columns=str.lower)
    dic_anno = df_anno.set_index("#scan#").to_dict('index')
    dic_anno = {int(i): dic_anno[i] for i in dic_anno}
    mgf_instance = MgfInstance(Path(mgf_path))
    dict_data = mgf_instance.data
    matched_spectra_dict = {}
    for ID in dict_data:
        if (ID in dic_anno):
            matched_spectra_dict[ID] = MatchedSpectra(
                dict_data[ID], get_correct_inchi(dic_anno[ID]), dic_anno[ID]['mqscore'])
        else:
            matched_spectra_dict[ID] = MatchedSpectra(dict_data[ID], "#", 0)
    return (matched_spectra_dict)


def closest_gnps_iterative(auth, input_file_mgf, input_file_quan):
    """
    Returns a dictionary of MatchedSpectra objects for the closest matches in GNPS.

    Args:
        auth (Auth): Authentication object with username, password, and email.
        input_file_mgf (str): The local path of the .mfg file.
        input_file_quan (str): The local path of the quantification table.

    Returns:
        dict: Dictionary {id: MatchedSpectra}.
    """
    job_description = Path(input_file_mgf).stem
    path_file_mgf_in_gnps, path_file_quan_in_gnps = _upload_to_gnps(
        auth, input_file_mgf, input_file_quan, job_description)
    _wait_after_upload()
    dict_task_id = _launch_GNPS_workflow_iterative(
        auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description)
    _wait_for_workflow_finish(dict_task_id)

    return closest_gnps_ietrative_by_id(dict_task_id, input_file_mgf)


def closest_gnps_ietrative_by_id(dict_task_id, input_file_mgf):
    """
    Returns a dictionary of MatchedSpectra objects for the closest matches in GNPS using job IDs.

    Args:
        dict_task_id (dict): Dictionary containing task IDs for different peak matches.
        input_file_mgf (str): The local path of the .mfg file.

    Returns:
        dict: Dictionary {id: MatchedSpectra} for the closest matches based on the provided task IDs.
    """
    res = {}
    for peak in dict_task_id:
        s = {i: 0 for i in dict_task_id[peak]}
        for mass in dict_task_id[peak]:
            job = dict_task_id[peak][mass]
            s[mass] = closest_gnps_by_id(job, input_file_mgf)
        res[peak] = s

    mgf_instance = MgfInstance(Path(input_file_mgf))
    dict_data = mgf_instance.data

    alpha = _get_iterative_parameters()
    ordre = [0.02, 0.1, 10, 25, 50, 100, 250, 500, 0]
    matched_spectra_dict = {}
    for ID in dict_data:
        for peak in [6, 5, 4]:
            for mass in ordre:
                if (ID not in matched_spectra_dict):
                    # ['#','*','?']):
                    if (res[peak][mass][ID].inchi not in ['#', '?']):
                        sp = res[peak][mass][ID].spectrum
                        inc = res[peak][mass][ID].inchi
                        score = res[peak][mass][ID].score
                        matched_spectra_dict[ID] = MatchedSpectra(
                            sp, inc, score*alpha[peak][mass])
        if (ID not in matched_spectra_dict):
            sp = res[4][0][ID].spectrum
            inc = res[4][0][ID].inchi
            score = res[4][0][ID].score*alpha[4][0]
            matched_spectra_dict[ID] = MatchedSpectra(sp, inc, score)
    return (matched_spectra_dict)

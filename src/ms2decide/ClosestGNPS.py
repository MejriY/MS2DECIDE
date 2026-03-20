from .MatchedSpectra import MatchedSpectra
from .MgfInstance import MgfInstance
from .Util import path_gnps, Parametres, get_correct_inchi, _mass_selection_by_tolerance, _get_match
from matchms.importing import load_from_mgf
from matchms import Spectrum
import urllib.request
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


def _invoke_workflow(auth, base_url, parameters):
    """
    Log in to a specific URL using username and password to launch a workflow with given parameters.

    Args:
        auth (Auth): Authentication object with username and password.
        base_url (str): Base URL to log in (e.g., 'gnps.ucsd.edu').
        parameters (dict): Parameters of the workflow.

    Returns:
        str: Task ID of the launched workflow.
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
    r = s.post('https://' + base_url + '/ProteoSAFe/InvokeTools',
               data=parameters, verify=verify_https)
    task_id = r.text

    print(r.text, file=sys.stderr, flush=True)

    if len(task_id) > 4 and len(task_id) < 60:
        print("Launched Task: : " + r.text)
        return task_id
    else:
        print(task_id)
        return None


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


def _launch_GNPS_workflow(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description, score_threshold):
    """
    Log in to GNPS using username and password, upload file to MassIVE repository, and launch feature-based molecular networking workflow.

    Args:
        auth (Auth): Authentication object with username, password, and email.
        file_path (str): The local path of the .mfg file.
        job_description (str): Title of the job.

    Returns:
        str: Task ID of the launched workflow.
    """
    invokeParameters = {}
    invokeParameters = _get_networking_parameters()
    invokeParameters["SCORE_THRESHOLD"] = score_threshold
    invokeParameters["desc"] = job_description
    invokeParameters["MAX_SHIFT_MASS"] = str(0.02)
    invokeParameters["quantification_table"] = "d./" + \
        auth.username + path_file_quan_in_gnps
    invokeParameters["spec_on_server"] = "d./" + \
        auth.username + path_file_mgf_in_gnps
    invokeParameters["email"] = auth.mail

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
        "INPUT A SCORE_THRESHOLD FOR ITERATIVE GNPS. DEFAULT 0.001")
    D_job = {}
    for peak in [6, 5, 4]:
        d = {i: 0 for i in [0.02, 0.1, 10, 25, 50, 100, 250, 500, 0]}
        for mass_diff in [0.02, 0.1, 10, 25, 50, 100, 250, 500, 0]:
            invokeParameters = {}
            invokeParameters = _get_networking_parameters()
            invokeParameters["SCORE_THRESHOLD"] = SCORE_THRESHOLD
            invokeParameters["MIN_MATCHED_PEAKS_SEARCH"] = str(peak)
            invokeParameters["MAX_SHIFT_MASS"] = str(mass_diff)
            invokeParameters["desc"] = job_description + \
                '_'+str(peak)+'_'+str(mass_diff)
            invokeParameters["quantification_table"] = "d./" + \
                auth.username + path_file_quan_in_gnps
            invokeParameters["spec_on_server"] = "d./" + \
                auth.username + path_file_mgf_in_gnps
            invokeParameters["email"] = auth.mail
            task_id = _invoke_workflow(auth, "gnps.ucsd.edu", invokeParameters)
            d[mass_diff] = task_id
        D_job[peak] = d
    return D_job


def _gnps_annotations_download_results(task_id):
    """
    Download GNPS molecular networking job results locally.

    Args:
        task_id (str): GNPS task ID.

    Returns:
        pd.DataFrame: Dataframe containing job annotations.
    """

    # Base link to download GNPS job annotations
    if (task_id[:3].lower() == 'id='):
        task_id = task_id[3:]
    gnps_download_link = 'https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp?task={}&view=view_all_annotations_DB'
    # Launch the download
    with urllib.request.urlopen(gnps_download_link.format(task_id)) as url:
        raw_data = url.read()  # Read the binary data
        # Detect encoding
        detected_encoding = chardet.detect(raw_data)
        decoded_data = raw_data.decode(
            detected_encoding['encoding'], errors='replace')
    # Parse the string as JSON
    try:
        # Convert the decoded string to a Python dictionary
        data = json.loads(decoded_data)
        print("JSON Loaded Successfully")
        df_annotations = pd.DataFrame.from_dict(data['blockData'])
        print('==================')
        print(' FBMN job detected with ' +
              str(df_annotations.shape[0]-1) + ' spectral library annotations in the job:' + task_id)
        print('==================')
        df_annotations[['#Scan#', 'MQScore']] = df_annotations[[
            '#Scan#', 'MQScore']].apply(pd.to_numeric)
        return (df_annotations)
    except json.JSONDecodeError as e:
        print(f"JSON Decode Error: {e}")
        raise Exception('CAN NOT GET THE DATA FROM THE JOB: ' + task_id)


def _workflows_statues(dict_task_id):
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
            json_obj = json.loads(
                requests.get(
                    url.format(dict_task_id[peak][mass]),
                    verify=verify_https,
                ).text
            )
            state.add(json_obj["status"])
    return (state)


def _wait_for_workflow_finish(dict_task_id, poll_interval=60, max_wait_seconds=None):
    if poll_interval <= 0:
        raise ValueError("poll_interval doit etre > 0.")

    start_time = time.time()
    statues = _workflows_statues(dict_task_id)
    print(f"Current GNPS status: {statues}")
    while statues != {'DONE'} :
        if max_wait_seconds is not None and (time.time() - start_time) >= max_wait_seconds:
            raise TimeoutError(
                f"Timeout atteint apres {max_wait_seconds} secondes. Dernier statut: {statues}"
            )
        time.sleep(poll_interval)
        statues = _workflows_statues(dict_task_id)
        print(f"Current GNPS status: {statues}")
    print('all Jobs done')
    return statues


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
    job_description = input_file_mgf.stem
    task_id = _launch_GNPS_workflow(
        auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description, score_threshold)
    _wait_for_workflow_finish({0:{0:task_id}})
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
        auth (Auth): Authentication object with username, password, and email.
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
        mgf_path (str): The local path of the .mfg file.

    Returns:
        dict: Dictionary {id: MatchedSpectra}.
    """
    mgf_instance = MgfInstance(Path(input_file_mgf))
    dict_data = mgf_instance.data
    job_description = Path(input_file_mgf).stem
    path_file_mgf_in_gnps, path_file_quan_in_gnps = _upload_to_gnps(
        auth, input_file_mgf, input_file_quan, job_description)
    dict_task_id = _launch_GNPS_workflow_iterative(
        auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description)
    _wait_for_workflow_finish(dict_task_id)

    res = {}
    for peak in dict_task_id:
        s = {i: 0 for i in dict_task_id[peak]}
        for mass in dict_task_id[peak]:
            job = dict_task_id[peak][mass]
            s[mass] = closest_gnps_by_id(job, input_file_mgf)
        res[peak] = s

    alpha = _get_iterative_parameters
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

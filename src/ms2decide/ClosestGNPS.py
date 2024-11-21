from .MatchedSpectra import MatchedSpectra
from .MgfInstance import MgfInstance
from .Util import path_gnps, Parametres, get_correct_inchi, _mass_selection_by_tolerance, _get_match
from matchms.importing import load_from_mgf
from matchms import Spectrum
import urllib.request
import json
import requests
import sys
import pandas as pd
from pathlib import Path
import ftplib
# import datetime


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

    s = requests.Session()

    payload = {
        'user': username,
        'password': password,
        'login': 'Sign in'
    }

    r = s.post('https://' + base_url +
               '/ProteoSAFe/user/login.jsp', data=payload, verify=False)
    r = s.post('https://' + base_url + '/ProteoSAFe/InvokeTools',
               data=parameters, verify=False)
    task_id = r.text

    print(r.text, file=sys.stderr, flush=True)

    if len(task_id) > 4 and len(task_id) < 60:
        print("Launched Task: : " + r.text)
        return task_id
    else:
        print(task_id)
        return None


def _upload_to_gnps(auth, input_file_mgf, input_file_quan, folder):
    """
    Log in to the MassIVE repository (GNPS library) to upload a given file.

    Args:
        auth (Auth): Authentication object with username and password.
        input_file (str): The local path of the .mfg file.

    Returns:
        str: The name of the uploaded file.
    """
    username = auth.username
    password = auth.password
    url = "ccms-ftp01.ucsd.edu"
    ftp = ftplib.FTP(url, username, password)
    print("Login successful \n")
    try:
        ftp.mkd(folder)
        ftp.cwd('/'+folder+'/')
    except:
        raise Exception('Folder already in server')
    print('Folder created')
    pp = ftp.pwd()
    file_mgf = open(str(input_file_mgf), "rb")
    input_file_mgf = Path(input_file_mgf).stem+'.mgf'
    ftpCommand = "STOR " + input_file_mgf
    ftpResponseMessage = ftp.storbinary(ftpCommand, fp=file_mgf)
    file_mgf.close()
    path_mgf = pp+'/'+input_file_mgf

    file_quan = open(str(input_file_quan), "rb")
    input_file_csv = Path(input_file_quan).stem+'.csv'
    ftpCommand = "STOR " + input_file_csv
    ftpResponseMessage = ftp.storbinary(ftpCommand, fp=file_quan)
    file_quan.close()
    ftp.close()
    path_csv = pp+'/'+input_file_csv
    print('Files uploaded')
    return (path_mgf, path_csv)


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


def _launch_GNPS_workflow(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description):
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
    invokeParameters["desc"] = job_description
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
    D_job = {}
    for peak in [6, 5, 4]:
        d = {i: 0 for i in [0.02, 0.1, 10, 25, 50, 100, 250, 500, 0]}
        for mass_diff in [0.02, 0.1, 10, 25, 50, 100, 250, 500, 0]:
            invokeParameters = {}
            invokeParameters = _get_networking_parameters()
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
        data = json.load(url)
    df_annotations = pd.DataFrame.from_dict(data['blockData'])
    print('==================')
    print(' FBMN job detected with ' +
          str(df_annotations.shape[0]-1) + ' spectral library annotations in the job:' + task_id)
    print('==================')
    return (df_annotations)


def _wait_for_workflow_finish(dict_task_id):
    state = set()
    url = 'https://gnps.ucsd.edu/ProteoSAFe/status_json.jsp?task={}'
    for peak in dict_task_id:
        for mass in dict_task_id[peak]:
            json_obj = json.loads(requests.get(url.format(
                dict_task_id[peak][mass]), verify=False).text)
            state.add(json_obj["status"])
    return (state)


def closest_gnps(auth, input_file_mgf, input_file_quan):
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
        auth, input_file_mgf, input_file_quan)
    job_description = input_file_mgf.stem
    task_id = _launch_GNPS_workflow(
        auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description)
    _wait_for_workflow_finish("gnps.ucsd.edu", task_id)
    df_anno = _gnps_annotations_download_results(task_id)
    df_anno = df_anno.rename(columns=str.lower)
    dic_anno = df_anno.set_index("#scan#").to_dict('index')
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
    _wait_for_workflow_finish("gnps.ucsd.edu", dict_task_id[4][0])

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

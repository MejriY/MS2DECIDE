from .MatchedSpectra import MatchedSpectra
from .MgfInstance import MgfInstance
from .Util import path_gnps, Parametres, get_correct_inchi, _mass_selection_by_tolerance, _get_match
from matchms.importing import load_from_mgf
from matchms import Spectrum
import os
import subprocess
import shlex
import zipfile
import json
import time
import requests
import sys
import pandas as pd
from pathlib import Path
import ftplib
import datetime
import tarfile


# def _invoke_workflow(auth, base_url, parameters):
# return None


# def _upload_to_gnps(auth, input_file_mgf, input_file_quan, folder):
# return None


# def _get_networking_parameters():
# return invokeParameters


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
                 1000: round(N1-0.025*14, 3)},
             5: {0.02: round(N2, 3),
                 0.1: round(N2-0.025*1, 3),
                 10: round(N2-0.025*2, 3),
                 25: round(N2-0.025*4, 3),
                 50: round(N2-0.025*6, 3),
                 100: round(N2-0.025*8, 3),
                 250: round(N2-0.025*10, 3),
                 500: round(N2-0.025*12, 3),
                 1000: round(N2-0.025*14, 3)},
             4: {0.02: N3,
                 0.1: round(N3-0.025*1, 3),
                 10: round(N3-0.025*2, 3),
                 25: round(N3-0.025*4, 3),
                 50: round(N3-0.025*6, 3),
                 100: round(N3-0.025*8, 3),
                 250: round(N3-0.025*10, 3),
                 500: round(N3-0.025*12, 3),
                 1000: round(N3-0.025*14, 3)}}
    return alpha


# def _launch_GNPS_workflow(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description):
    # return None


# def _launch_GNPS_workflow_iterative(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description):
    # return None


def _gnps2_annotations_download_results(task_id):
    # Base link to download GNPS job annotations
    if (task_id[:3].lower() == 'id='):
        task_id = task_id[3:]
    gnps_download_link = "https://gnps2.org/taskzip?task="+task_id
    # Launch the download
    print('==================')
    print('This is the GNPS job link:' + gnps_download_link)
    # output_folder = folder to storge gnps download but will be deleted in the end
    output_folder = os.getcwd()  # LAZEM TBADLOU
    zip_gnps_annotation = os.path.join(output_folder, task_id+'.tar')
    # Clearing local files
    if (task_id+'.tar' in os.listdir(output_folder)):
        os.remove(zip_gnps_annotation)
    r = requests.get(gnps_download_link)
    # Call the process and catch errors
    try:
        with open(zip_gnps_annotation, 'wb') as f:
            f.write(r.content)
    except:
        raise Exception('problem in file from GNPS2')

    try:
        with tarfile.open(zip_gnps_annotation) as tar:
            l = tar.getmembers()
            path = [i for i in l if i.name.startswith(
                "nf_output/library/merged_results_with_gnps.tsv")]
            path_data = tar.extractfile(path[0])
            df_annotations = pd.read_csv(path_data, sep='\t')
            print('==================')
            print('   FBMN job detected')
            print('==================')
            print('      '+str(df_annotations.shape[0]-1) +
                  ' spectral library annotations in the job.')
    except:
        raise Exception('check the job ID and/or GNPS version')
    df_annotations['InChI'] = df_annotations['INCHI']
    os.remove(zip_gnps_annotation)
    return (df_annotations)


# def _wait_for_workflow_finish(base_url, task_id):
    # return None


# def closest_gnps2(task_id,input_file_mgf):
    # return (matched_spectra_dict)


def closest_gnps2_by_id(job_id, mgf_path):
    """
    Returns a dictionary of MatchedSpectra objects for the closest matches in GNPS.

    Args:
        auth (Auth): Authentication object with username, password, and email.
        mgf_path (str): The local path of the .mfg file.

    Returns:
        dict: Dictionary {id: MatchedSpectra}.
    """
    df_anno = _gnps2_annotations_download_results(job_id)
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


# def closest_gnps2_iterative(auth, input_file_mgf, input_file_quan):
    # return (matched_spectra_dict)


def closest_gnps2_ietrative_by_id(dict_task_id, input_file_mgf):
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
            s[mass] = closest_gnps2_by_id(job, input_file_mgf)
        res[peak] = s

    mgf_instance = MgfInstance(Path(input_file_mgf))
    dict_data = mgf_instance.data

    alpha = _get_iterative_parameters()
    ordre = [0.02, 0.1, 10, 25, 50, 100, 250, 500, 1000]
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
            sp = res[4][1000][ID].spectrum
            inc = res[4][1000][ID].inchi
            score = res[4][1000][ID].score*alpha[4][1000]
            matched_spectra_dict[ID] = MatchedSpectra(sp, inc, score)
    return (matched_spectra_dict)

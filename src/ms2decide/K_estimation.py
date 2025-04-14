from .MgfInstance import MgfInstance
from .AuthMail import AuthMail
from .ClosestGNPS import _launch_GNPS_workflow_iterative, _upload_to_gnps, closest_gnps_ietrative_by_id, _wait_for_workflow_finish, _launch_GNPS_workflow, closest_gnps_by_id
from .IsdbAnnotation import get_cfm_annotation
from .SiriusAnnotation import SiriusAnnotation
from .gnps_empty_repport import gnps_empty_repport
from .MultipleSourceAnnotation import MultipleSourceAnnotation, MultipleSourceAnnotation_to_dataframe
import time
import json
import getpass
from pathlib import Path


def K_estimation():
    """
    Estimate K values for features prioritization based on MS/MS-based GNPS, ISDB-LOTUS, and Sirius annotations.

    This function prompts the user for various input file paths, including an authentication information,
    quantitative data, MGF data, and the results of a Sirius job. It then retrieves and processes
    annotation data from multiple sources and compiles the results into a dataframe filtered by specific
    precursor mass criteria.

    Returns:
        DataFrame: A dataframe containing the filtered results based on K estimation.
    """
    username = input('GNPS username: ')
    password = getpass.getpass("GNPS password: ")
    mail = input('GNPS mail: ')
    auth = AuthMail(username, password, mail)

    quan_path = input(
        'SELECT THE PATH FOR YOUR QUANTITATIVE FILE. This path needs to terminate with a *.csv at the end \n :')
    mgf_path = input(
        'SELECT THE PATH FOR YOU MGF FILE. This path needs to terminate with a *.mgf at the end \n :')
    mgf = MgfInstance(Path(mgf_path))
    job_description = input('INPUT A TITLE FOR YOUR FBMN GNPS JOB \n :')
    path_file_mgf_in_gnps, path_file_quan_in_gnps = _upload_to_gnps(
        auth, Path(mgf_path), Path(quan_path), job_description)
    reaserach_type = input(
        "PLEASE SPECIFY THE TYPE OF GNPS LIBRARY SEARCH THAT YOU WOULD LIKE TO APPLY. \n strict: for mass difference of 0.02 Da and six matched peaks with very low cosine threshold. \n iterative: for iterative weighted analog search (can take up to three hours). \n")
    if (reaserach_type.lower() == 'iterative'):
        dict_task_id = _launch_GNPS_workflow_iterative(
            auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description)
        fald_job = input(
            "DO YOU WANT TO SAVE ITERATIVE GNPS JOBS? (yes or no)")
        if (fald_job.lower() == 'yes'):
            save_path_job = input(
                'SELECT THE SAVE PATH FOR THE .txt FILE OF ITERATIVE GNPS JOBS. \n This path needs to terminate with a *.txt at the end \n :')
            f = open(save_path_job, "w")
            f.write(str(dict_task_id))
            f.close()
    else:
        score_threshold = input(
            'Please can you provide a minimum cosine score that MS/MS spectra should get in spectral matching ')
        task_id = _launch_GNPS_workflow(auth, path_file_mgf_in_gnps,
                                        path_file_quan_in_gnps, job_description, score_threshold)

    ISDBtol = float(
        input('SELECT MASS TOLERANCE FOR ISDB-LOTUS ANNOTATION (LESS THAN 0.5, DEFAULT 0.02) '))
    while (ISDBtol > 0.5):
        ISDBtol = input(
            'SELECT MASS TOLERANCE FOR ISDB-LOTUS ANNOTATION (LESS THAN 0.5, DEFAULT 0.02) ')
    isdb_res = get_cfm_annotation(mgf, ISDBtol)

    sirius_path = input(
        'SELECT THE PATH FOR YOUR SIRIUS6 ANNOTATION FILE. This path needs to terminate with structure_identifications.tsv at the end \n :')
    indexScore = input('Select index score = exact, approximate: ')
    sirius_res = SiriusAnnotation(sirius_path, mgf, indexScore)

    if (reaserach_type.lower() == 'iterative'):
        state = _wait_for_workflow_finish(dict_task_id)
        while (state != {'DONE'}):
            if ('FAILED' in state):
                raise Exception('FAILED in GNPS Job')
            time.sleep(100)
            state = _wait_for_workflow_finish(dict_task_id)
        gnps_res = closest_gnps_ietrative_by_id(dict_task_id, mgf_path)
    else:
        state = _wait_for_workflow_finish({6: {0.02: task_id}})
        while (state != {'DONE'}):
            if ('FAILED' in state):
                raise Exception('FAILED in GNPS Job')
            time.sleep(100)
            state = _wait_for_workflow_finish({6: {0.02: task_id}})
        gnps_res = closest_gnps_by_id(task_id, mgf_path)

    get_gnps, get_sirius, get_isdb = True, True, True
    results = {}
    results['K'] = MultipleSourceAnnotation(
        mgf, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res, isdb_res, 'K')
    dfw = MultipleSourceAnnotation_to_dataframe(
        mgf, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res, isdb_res, results)
    # taken = [i for i in mgf.data if mgf.data[i].metadata['precursor_mz']
    # > 200 and mgf.data[i].metadata['precursor_mz'] < 850]
    # taken = [i for i in taken if len(mgf.data[i].peaks.mz) > 3]
    # dfw_taken = dfw[dfw.ID.isin(taken)]
    dfw = dfw.sort_values(by=['K'])
    dfw['ranking by k'] = [i+1 for i in range(len(dfw))]
    save_path = input(
        'SELECT THE SAVE PATH FOR THE .TSV FILE OF MS2DECIDE OUPUT. This path needs to terminate with a *.tsv at the end \n :')
    dfw.to_csv(save_path, sep='\t', index=False)
    print('YOU CAN ANALYZE THE RANKING BY K BY MAPPING IT USING YOUR GRAPH SOFTWARE.')
    if (reaserach_type.lower() == 'iterative'):
        get_empty = input('GET EMPTY ANNOTATION MSLIB ID (yes or no): ')
        if (get_empty == 'yes'):
            save_path_empty = Path(save_path).parent.joinpath('empty.tsv')
            gnps_empty_repport('GNPS', dict_task_id, dfw, save_path_empty)
            print('EMPTY ANNOTATION FILE IS SAVED IN THE SAME DIRECTORY AS K RANKING')

    return (dfw)

from .MgfInstance import MgfInstance
from .AuthMail import AuthMail
from .ClosestGNPS import _launch_GNPS_workflow_iterative, _upload_to_gnps, closest_gnps_ietrative_by_id, _wait_for_workflow_finish
from .IsdbAnnotation import get_cfm_annotation
from .SiriusAnnotation import SiriusAnnotation
from .MultipleSourceAnnotation import MultipleSourceAnnotation, MultipleSourceAnnotation_to_dataframe
import time
import getpass
from pathlib import Path


def K_estimation():
    """
    Estimate K values for compounds based on GNPS, ISDB, and Sirius annotations.

    This function prompts the user for various input file paths, including an authentication file,
    quantitative data, MGF data, and the results of a Sirius job. It then retrieves and processes
    annotation data from multiple sources and compiles the results into a dataframe filtered by specific
    precursor mass criteria.

    Returns:
        DataFrame: A dataframe containing the filtered results based on K estimation.
    """
    username = input('GNPS username : ')
    password = getpass.getpass("GNPS passeword : ")
    mail = input('GNPS mail : ')
    auth = AuthMail(username, password, mail)

    quan_path = input('SELECT THE PATH FOR YOUR QUANTITAIVE DATA \n :')
    mgf_path = input('SELECT THE PATH OF YOU MGF DATA \n :')
    mgf = MgfInstance(Path(mgf_path))
    job_description = input('Name for the job description in GNPS \n :')
    path_file_mgf_in_gnps, path_file_quan_in_gnps = _upload_to_gnps(
        auth, Path(mgf_path), Path(quan_path), job_description)
    dict_task_id = _launch_GNPS_workflow_iterative(
        auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description)

    ISDBtol = float(
        input('SELECT MASS TOLENRENCE FOR ISDB-LOTUS LESS THEN 0.5, DEFAULT 0.02 '))
    while (ISDBtol > 0.5):
        ISDBtol = input(
            'SELECT MASS TOLENRENCE FOR ISDB-LOTUS LESS THEN 0.5, DEFAULT 0.02 ')
    isdb_res = get_cfm_annotation(mgf, ISDBtol)

    sirius_path = input(
        'SELECT THE PATH FOR THE RESULT OF SIRIUS JOB (THE structure_identifications.tsv) \n :')
    indexScore = input('Select index score =  exact, approximate')
    sirius_res = SiriusAnnotation(sirius_path, mgf, indexScore)

    flag = _wait_for_workflow_finish("gnps.ucsd.edu", dict_task_id[4][0])
    while (flag != 'DONE'):
        time.sleep(100)
        flag = _wait_for_workflow_finish("gnps.ucsd.edu", dict_task_id[4][0])

    gnps_res = closest_gnps_ietrative_by_id(dict_task_id, mgf_path)

    get_gnps, get_sirius, get_isdb = True, True, True
    results = {}
    results['K'] = MultipleSourceAnnotation(
        mgf, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res.data, isdb_res, 'K')
    dfw = MultipleSourceAnnotation_to_dataframe(
        mgf, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res.data, isdb_res, results)
    taken = [i for i in mgf.data if mgf.data[i].metadata['precursor_mz']
             > 200 and mgf.data[i].metadata['precursor_mz'] < 850]
    taken = [i for i in taken if len(mgf.data[i].peaks.mz) > 3]
    dfw_taken = dfw[dfw.ID.isin(taken)]
    dfw = dfw.sort_values(by=['K'])
    return (dfw_taken)

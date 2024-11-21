from .MgfInstance import MgfInstance
from .AuthMail import AuthMail
from .ClosestGNPS import _launch_GNPS_workflow_iterative, _upload_to_gnps, closest_gnps_ietrative_by_id, _wait_for_workflow_finish
from .IsdbAnnotation import get_cfm_annotation, get_cfm_annotation_GUI
from .SiriusAnnotation import SiriusAnnotation
from .MultipleSourceAnnotation import MultipleSourceAnnotation, MultipleSourceAnnotation_to_dataframe
import time
from pathlib import Path


def K_estimation_GUI(username, password, mail, mgf_path, quan_path, job_description, path_isdb, ion_mode, ISDBtol, sirius_path, indexScore, save_path):
    """
    Estimate K values for compounds based on GNPS, ISDB, and Sirius annotations.

    This function prompts the user for various input file paths, including an authentication file,
    quantitative data, MGF data, and the results of a Sirius job. It then retrieves and processes
    annotation data from multiple sources and compiles the results into a dataframe filtered by specific
    precursor mass criteria.

    Returns:
        DataFrame: A dataframe containing the filtered results based on K estimation.
    """

    auth = AuthMail(username, password, mail)

    mgf = MgfInstance(Path(mgf_path))
    path_file_mgf_in_gnps, path_file_quan_in_gnps = _upload_to_gnps(
        auth, Path(mgf_path), Path(quan_path), job_description)
    dict_task_id = _launch_GNPS_workflow_iterative(
        auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description)

    isdb_res = get_cfm_annotation_GUI(path_isdb, mgf, ion_mode, ISDBtol)

    sirius_res = SiriusAnnotation(sirius_path, mgf, indexScore)

    state = _wait_for_workflow_finish(dict_task_id)
    while (state != {'DONE'}):
        if ('FAILED' in state):
            raise Exception('FAILED in GNPS Job')
        time.sleep(100)
        state = _wait_for_workflow_finish(dict_task_id)

    gnps_res = closest_gnps_ietrative_by_id(dict_task_id, mgf_path)

    get_gnps, get_sirius, get_isdb = True, True, True
    results = {}
    results['K'] = MultipleSourceAnnotation(
        mgf, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res.data, isdb_res, 'K')
    dfw = MultipleSourceAnnotation_to_dataframe(
        mgf, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res.data, isdb_res, results)
    #taken = [i for i in mgf.data if mgf.data[i].metadata['precursor_mz']
             #> 200 and mgf.data[i].metadata['precursor_mz'] < 850]
    #taken = [i for i in taken if len(mgf.data[i].peaks.mz) > 3]
    #dfw_taken = dfw[dfw.ID.isin(taken)]
    dfw = dfw.sort_values(by=['K'])
    dfw.to_csv(save_path,sep='\t',inde=False)
    return (dfw)

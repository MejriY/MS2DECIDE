from .MgfInstance import MgfInstance
from .AuthMail import AuthMail
from .ClosestGNPS import _launch_GNPS_workflow_iterative, _upload_to_gnps, closest_gnps_ietrative_by_id, _wait_for_workflow_finish, _launch_GNPS_workflow, closest_gnps_by_id
from .IsdbAnnotation import get_cfm_annotation
from .SiriusAnnotation import SiriusAnnotation
from .gnps_empty_repport import gnps_empty_repport
from .MultipleSourceAnnotation import MultipleSourceAnnotation, MultipleSourceAnnotation_to_dataframe
import os
import re
import pickle
import getpass
from pathlib import Path


def _clean_path(raw):
    """Remove spaces and quotes added by copy/paste or drag-and-drop."""
    return str(raw).strip().strip('"\'').strip()


def _ask_float(prompt, default, low, high):
    """Ask for a float in ]low, high]; empty answer returns default."""
    while True:
        raw = input(prompt).strip()
        try:
            value = float(raw) if raw else default
        except ValueError:
            print("Valeur invalide.")
            continue
        if low < value <= high:
            return value
        print(f"La valeur doit etre comprise entre {low} (exclu) et {high}.")


def _ask_choice(prompt, choices):
    """Ask until the answer is one of the given choices (case-insensitive)."""
    while True:
        raw = input(prompt).strip().lower()
        if raw in choices:
            return raw
        print(f"Choix possibles : {', '.join(choices)}")


def _isdb_cache_path(mgf_path):
    """Default location of the ISDB results cache, next to the MGF file."""
    mgf_path = Path(mgf_path)
    return mgf_path.with_name(mgf_path.stem + '_isdb_results.pkl')


def _isdb_signature(mgf_path, isdb_tol):
    """Identify the inputs of an ISDB run (MGF file size/date and tolerance)."""
    stat = os.stat(str(mgf_path))
    return {
        'mgf_name': Path(mgf_path).name,
        'mgf_size': stat.st_size,
        'mgf_mtime': int(stat.st_mtime),
        'isdb_tol': float(isdb_tol),
    }


def _save_isdb_results(cache_path, signature, isdb_res):
    """Save ISDB results to disk. A failure never stops the workflow."""
    tmp_path = Path(str(cache_path) + '.tmp')
    try:
        with open(tmp_path, 'wb') as f:
            pickle.dump({'signature': signature, 'isdb_res': isdb_res},
                        f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp_path, cache_path)
        print(f"ISDB results saved in: {cache_path}")
    except Exception as exc:
        print(f"WARNING: ISDB results could not be saved ({exc}). The run continues.")
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass


def _load_isdb_results(cache_path, signature):
    """Return the saved ISDB results if they match the current inputs, else None."""
    cache_path = Path(cache_path)
    if not cache_path.exists():
        return None
    try:
        with open(cache_path, 'rb') as f:
            content = pickle.load(f)
    except Exception as exc:
        print(f"WARNING: ISDB cache unreadable ({exc}); it will be recomputed.")
        return None
    if content.get('signature') != signature:
        print("ISDB cache found but it was made with a different MGF file or tolerance; it will be recomputed.")
        return None
    return content.get('isdb_res')


def K_estimation():
    """
    Estimate K values for features prioritization based on MS/MS-based GNPS, ISDB-LOTUS, and Sirius annotations.

    This function prompts the user for various input file paths, including an authentication information,
    quantitative data, MGF data, and the results of a Sirius job. It then retrieves and processes
    annotation data from multiple sources and compiles the results into a dataframe filtered by specific
    precursor mass criteria.

    An existing GNPS FBMN job can be reused (strict mode) by providing its task ID,
    which skips the upload and the launch of a new GNPS job.

    The ISDB results are saved next to the MGF file (<mgf name>_isdb_results.pkl) and can be
    reused on the next run if the MGF file and the tolerance are unchanged.

    Returns:
        DataFrame: A dataframe containing the filtered results based on K estimation.
    """
    username = input('GNPS username: ').strip()
    password = getpass.getpass("GNPS password: ")
    mail = input('GNPS mail: ').strip()
    auth = AuthMail(username, password, mail)

    quan_path = _clean_path(input(
        'SELECT THE PATH FOR YOUR QUANTITATIVE FILE. This path needs to terminate with a *.csv at the end \n :'))
    mgf_path = _clean_path(input(
        'SELECT THE PATH FOR YOU MGF FILE. This path needs to terminate with a *.mgf at the end \n :'))
    mgf = MgfInstance(Path(mgf_path))

    existing_task = input(
        'EXISTING GNPS FBMN TASK ID TO REUSE (leave empty to launch a new GNPS job) \n :').strip()
    if existing_task.lower().startswith('task='):
        existing_task = existing_task[5:]

    dict_task_id = None
    task_id = None
    if existing_task:
        reaserach_type = 'strict'
        task_id = existing_task
        print("Reusing GNPS task:", task_id)
    else:
        job_description = input('INPUT A TITLE FOR YOUR FBMN GNPS JOB \n :').strip()
        # Safe remote folder name (no spaces / special characters).
        remote_folder = re.sub(r'[^A-Za-z0-9_.-]', '_', job_description) or 'ms2decide_job'
        path_file_mgf_in_gnps, path_file_quan_in_gnps = _upload_to_gnps(
            auth, Path(mgf_path), Path(quan_path), remote_folder)
        reaserach_type = _ask_choice(
            "PLEASE SPECIFY THE TYPE OF GNPS LIBRARY SEARCH THAT YOU WOULD LIKE TO APPLY. \n"
            " strict: for mass difference of 0.02 Da and six matched peaks with very low cosine threshold. \n"
            " iterative: for iterative weighted analog search (can take up to three hours). \n",
            ('strict', 'iterative'))
        if reaserach_type == 'iterative':
            dict_task_id = _launch_GNPS_workflow_iterative(
                auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description)
            print("GNPS iterative task IDs:", dict_task_id)
            save_jobs = _ask_choice(
                "DO YOU WANT TO SAVE ITERATIVE GNPS JOBS? (yes or no) ", ('yes', 'no'))
            if save_jobs == 'yes':
                save_path_job = _clean_path(input(
                    'SELECT THE SAVE PATH FOR THE .txt FILE OF ITERATIVE GNPS JOBS. \n This path needs to terminate with a *.txt at the end \n :'))
                with open(save_path_job, "w") as f:
                    f.write(str(dict_task_id))
        else:
            score_threshold = _ask_float(
                'Please can you provide a minimum cosine score that MS/MS spectra should get in spectral matching (0-1, default 0.1) ',
                0.1, 0.0, 1.0)
            task_id = _launch_GNPS_workflow(
                auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description, str(score_threshold))
            # Keep the task ID visible so the run can be resumed if a later step fails.
            print("GNPS task_id (keep it to resume):", task_id)

    ISDBtol = _ask_float(
        'SELECT MASS TOLERANCE FOR ISDB-LOTUS ANNOTATION (LESS THAN 0.5, DEFAULT 0.02) ',
        0.02, 0.0, 0.5)

    # ISDB annotation: reuse the saved results when the inputs are unchanged.
    cache_path = _isdb_cache_path(mgf_path)
    signature = _isdb_signature(mgf_path, ISDBtol)
    isdb_res = None
    if cache_path.exists():
        cached = _load_isdb_results(cache_path, signature)
        if cached is not None:
            reuse = _ask_choice(
                f"SAVED ISDB RESULTS FOUND ({cache_path}). REUSE THEM? (yes or no) ",
                ('yes', 'no'))
            if reuse == 'yes':
                isdb_res = cached
                print("ISDB results loaded from file.")
    if isdb_res is None:
        isdb_res = get_cfm_annotation(mgf, ISDBtol)
        _save_isdb_results(cache_path, signature, isdb_res)

    sirius_path = _clean_path(input(
        'SELECT THE PATH FOR YOUR SIRIUS6 ANNOTATION FILE. This path needs to terminate with structure_identifications.tsv at the end \n :'))
    indexScore = _ask_choice('Select index score = exact, approximate: ',
                             ('exact', 'approximate'))
    sirius_res = SiriusAnnotation(sirius_path, mgf, indexScore)

    # _wait_for_workflow_finish blocks until all jobs are DONE and raises on FAILED/SUSPENDED.
    if reaserach_type == 'iterative':
        _wait_for_workflow_finish(dict_task_id)
        gnps_res = closest_gnps_ietrative_by_id(dict_task_id, mgf_path)
    else:
        _wait_for_workflow_finish({6: {0.02: task_id}})
        gnps_res = closest_gnps_by_id(task_id, mgf_path)

    get_gnps, get_sirius, get_isdb = True, True, True
    results = {}
    results['K'] = MultipleSourceAnnotation(
        mgf, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res, isdb_res, 'K')
    dfw = MultipleSourceAnnotation_to_dataframe(
        mgf, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res, isdb_res, results)
    dfw = dfw.sort_values(by=['K'])
    dfw['ranking by k'] = range(1, len(dfw) + 1)
    save_path = _clean_path(input(
        'SELECT THE SAVE PATH FOR THE .TSV FILE OF MS2DECIDE OUPUT. This path needs to terminate with a *.tsv at the end \n :'))
    dfw.to_csv(save_path, sep='\t', index=False)
    print('YOU CAN ANALYZE THE RANKING BY K BY MAPPING IT USING YOUR GRAPH SOFTWARE.')
    if reaserach_type == 'iterative':
        get_empty = _ask_choice('GET EMPTY ANNOTATION MSLIB ID (yes or no): ',
                                ('yes', 'no'))
        if get_empty == 'yes':
            save_path_empty = Path(save_path).parent.joinpath('empty.tsv')
            gnps_empty_repport('GNPS', dict_task_id, dfw, save_path_empty)
            print('EMPTY ANNOTATION FILE IS SAVED IN THE SAME DIRECTORY AS K RANKING')

    return (dfw)

from .MgfInstance import MgfInstance
from .ClosestGNPS import closest_gnps_ietrative_by_id, closest_gnps_by_id
from .IsdbAnnotation import get_cfm_annotation
from .SiriusAnnotation import SiriusAnnotation
from .gnps_empty_repport import gnps_empty_repport
from .MultipleSourceAnnotation import MultipleSourceAnnotation, MultipleSourceAnnotation_to_dataframe
from pathlib import Path


def K_estimation_by_id():
    """
    Estimate K values for compounds based on GNPS, ISDB, and Sirius annotations.

    This function prompts the user for various input file paths, including an authentication file,
    quantitative data, MGF data, and the results of a Sirius job. It then retrieves and processes
    annotation data from multiple sources and compiles the results into a dataframe filtered by specific
    precursor mass criteria.

    Returns:
        DataFrame: A dataframe containing the filtered results based on K estimation.
    """
    mgf_path = input('SELECT THE PATH FOR YOU MGF FILE \n :')
    mgf = MgfInstance(Path(mgf_path))

    reaserach_type = input(
        "PLEASE SPECIFY THE TYPE OF GNPS RESEARCH YOU WOULD LIKE TO APPLY. \n strict: for mass difference of 0.02 Da. \n iterative: for Weighted iterative GNPS analog search. \n")
    if (reaserach_type.lower() == 'iterative'):
        dict_task_id_path = input(
            'SELECT THE PATH FOR YOUR *.txt FILE THAT ITERATIVE GNPS JOBS')
        dict_task_id = eval(open(dict_task_id_path, 'r').read())
        gnps_res = closest_gnps_ietrative_by_id(dict_task_id, mgf_path)
    else:
        job_id = input(
            'PLEASE WRITE THE TASK JOB IN GNPS')
        gnps_res = closest_gnps_by_id(job_id, mgf_path)

    ISDBtol = float(
        input('SELECT MASS TOLERANCE FOR ISDB-LOTUS ANNOTATION (LESS THAN 0.5, DEFAULT 0.02 '))
    while (ISDBtol > 0.5):
        ISDBtol = input(
            'SELECT MASS TOLERANCE FOR ISDB-LOTUS ANNOTATION LESS THAN 0.5, DEFAULT 0.02 ')
    isdb_res = get_cfm_annotation(mgf, ISDBtol)

    sirius_path = input(
        'SELECT THE PATH FOR YOUR SIRIUS6 ANNOTATION FILE (structure_identifications.tsv) \n :')
    indexScore = input('Select index score =  exact, approximate: ')
    sirius_res = SiriusAnnotation(sirius_path, mgf, indexScore)

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
        'SELECT THE SAVE PATH FOR THE .TSV FILE OF MS2DECIDE. \n This path need to terminate with a *.tsv at the end')
    dfw.to_csv(save_path, sep='\t', index=False)
    print('YOU CAN ANALYZE THE RANKING BY K BY MAPPING IT ONTO YOUR VISUALIZATION SOFTWARE.')
    if (reaserach_type.lower() == 'iterative'):
        get_empty = input('GET EMPTY ANNOTATION MSLIB ID (yes or no): ')
        if (get_empty == 'yes'):
            save_path_empty = Path(save_path).parent.joinpath('empty.tsv')
            gnps_empty_repport('GNPS', dict_task_id, dfw, save_path_empty)
            print('EMPTY ANNOTATION ARE SAVED IN THE SAME DIRECTERY AS K RANKING')
    return (dfw)

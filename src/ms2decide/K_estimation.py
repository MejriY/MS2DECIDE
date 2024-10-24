from .MgfInstance import MgfInstance
from .AuthMail import AuthMail
from .ClosestGNPS import closest_gnps_iterative
from .IsdbAnnotation import get_cfm_annotation
from .SiriusAnnotation import SiriusAnnotation
from .MultipleSourceAnnotation import MultipleSourceAnnotation, MultipleSourceAnnotation_to_dataframe


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
    auth_path = input('SELECT THE PATH FOR YOUR AUTH.TXT FILE \n :')
    auth = AuthMail.from_txt(auth_path)

    quan_path = input('SELECT THE PATH FOR YOUR QUANTITAIVE DATA \n :')
    mgf_path = input('SELECT THE PATH OF YOU MGF DATA \n :')
    mgf = MgfInstance(mgf_path)
    sirius_path = input(
        'SELECT THE PATH FOR THE RESULT OF SIRIUS JOB (THE structure_identifications.tsv) \n :')
    indexScore = input('Select index score =  exact, approximate')
    gnps_res = closest_gnps_iterative(auth, mgf_path, quan_path)        
    ISDBtol = input('SELECT MASS TOLENRENCE FOR ISDB-LOTUS, DEFAULT 0.02')
    if(type(ISDBtol)==int)or(type(ISDBtol)==float):
        isdb_res = get_cfm_annotation(mgf,ISDBtol)
    else:
        isdb_res = get_cfm_annotation(mgf,ISDBtol)
    sirius_res = SiriusAnnotation(sirius_path,mgf,indexScore)

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
    return (dfw_taken)

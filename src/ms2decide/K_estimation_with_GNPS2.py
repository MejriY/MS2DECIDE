from .MgfInstance import MgfInstance
from .ClosestGNPS2 import closest_gnps2_ietrative_by_id
from .IsdbAnnotation import get_cfm_annotation
from .SiriusAnnotation import SiriusAnnotation
from .MultipleSourceAnnotation import MultipleSourceAnnotation, MultipleSourceAnnotation_to_dataframe


def K_estimation_with_GNPS2(dict_task_id, mgf_path, sirius_path, siriusIndex):

    mgf = MgfInstance(mgf_path)
    gnps2_res = closest_gnps2_ietrative_by_id(dict_task_id, mgf_path)
    isdb_res = get_cfm_annotation(mgf)
    sirius_res = SiriusAnnotation(sirius_path, mgf, siriusIndex)
    get_gnps, get_sirius, get_isdb = True, True, True
    results = {}
    results['K'] = MultipleSourceAnnotation(
        mgf, get_gnps, get_sirius, get_isdb, gnps2_res, sirius_res.data, isdb_res, 'K')
    dfw = MultipleSourceAnnotation_to_dataframe(
        mgf, get_gnps, get_sirius, get_isdb, gnps2_res, sirius_res.data, isdb_res, results)
    # taken = [i for i in mgf.data if mgf.data[i].metadata['precursor_mz']> 200 and mgf.data[i].metadata['precursor_mz'] < 850]
    taken = [i for i in mgf.data if len(mgf.data[i].peaks.mz) > 3]
    dfw_taken = dfw[dfw.ID.isin(taken)]
    return (dfw_taken)

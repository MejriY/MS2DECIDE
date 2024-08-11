import pandas as pd
from .Tool import Tool
from .Tanimotos import Tanimotos


def MultipleSourceAnnotation_to_dataframeForSiriusFormula(mgf, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res_anno, sirius_res_calcule, isdb_res, results):
    if (type(results) != dict):
        Exception(
            'results type is no dict \n please insert a resylt : dict {estimator : results_estimator} \n estimator : Matching, K, K_old, ConvEval')
    if (len(set(results.keys()) - set(['Matching', 'K', 'K_old', 'ConvEval'])) != 0):
        Exception(
            'please insert a resylt : dict {estimator : results_estimator} \n estimator : Matching, K, K_old, ConvEval')
    tool_dict = {}

    # GNPS
    if (get_gnps == True) and (gnps_res != None):
        tool_dict['GNPS'] = gnps_res

    # sirius
    if (get_sirius == True) and (sirius_res_calcule != None):
        tool_dict['SIRIUS'] = sirius_res_calcule

    # isdb
    if (get_isdb == True) and (isdb_res != None):
        tool_dict['ISDB'] = isdb_res

    data = []
    for ID in mgf.data:
        inchi_gnps = tool_dict['GNPS'][ID].inchi
        score_gnps = tool_dict['GNPS'][ID].score
        inchi_sirius_ann = sirius_res_anno[ID].inchi
        score_sirius = tool_dict[Tool(2).name][ID].score
        inchi_isdb = tool_dict['ISDB'][ID].inchi
        score_isdb = tool_dict['ISDB'][ID].score

        dict_inchi = {i: tool_dict[i][ID].inchi for i in tool_dict}
        t = Tanimotos(dict_inchi)
        tgs, tgi, tsi = t.compute_tanimoto()
        l = [ID, inchi_gnps, score_gnps, inchi_sirius_ann,
             score_sirius, inchi_isdb, score_isdb, tgs, tgi, tsi]
        for estimator in results:
            l.append(results[estimator][ID])
        data.append(l)
    df = pd.DataFrame(data, columns=['ID', 'inchi_gnps', 'score_gnps', 'inchi_sirius',
                      'score_sirius', 'inchi_isdb', 'score_isdb', 'tgs', 'tgi', 'tsi']+list(results.keys()))
    return (df)

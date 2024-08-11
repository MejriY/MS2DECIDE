from .Tool import Tool
from .Matching import Matching
from .K import K
from .K_old import K_old
from .ConfEval import ConfEval
from .Tanimotos import Tanimotos
import pandas as pd


def MultipleSourceAnnotation(mgf_instance, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res, isdb_res, estimator):
    """
    Get combined annotations from GNPS, Sirius, and ISDB sources.

    Args:
        mgf_instance (MgfInstance): Instance of MgfInstance containing spectra data.
        get_gnps (bool): Flag to determine whether to get GNPS annotation.
        get_sirius (bool): Flag to determine whether to get Sirius annotation.
        get_isdb (bool): Flag to determine whether to get ISDB annotation.
        gnps_res (dict): Result of GNPS annotation (default=None).
        sirius_res (dict): Result of Sirius annotation (default=None).
        isdb_res (dict): Result of ISDB annotation (default=None).
        estimator (str): ['ConvEval', 'Matching', 'K', 'K_old']
    Returns:
        dict: Dictionary containing combined annotations for each spectrum ID.
    """
    if(estimator not in ['ConvEval', 'Matching', 'K', 'K_old']):
        raise TypeError("error in estimator")
    tool_dict = {}

    # GNPS
    if(get_gnps == True) and (gnps_res != None):
        tool_dict[Tool(1).name] = gnps_res

    # sirius
    if(get_sirius == True) and (sirius_res != None):
        tool_dict[Tool(2).name] = sirius_res

    # isdb
    if(get_isdb == True) and (isdb_res != None):
        tool_dict[Tool(3).name] = isdb_res

    if(len(tool_dict) == 0):
        print("ERROR")
    else:
        data = {}
        for ID in mgf_instance.data:
            dict_index = {i: 0 for i in tool_dict}
            dict_inchi = {i: '' for i in tool_dict}
            for tool in tool_dict:
                try:
                    dict_index[tool] = tool_dict[tool][ID].score
                    dict_inchi[tool] = tool_dict[tool][ID].inchi
                except:
                    raise TypeError(
                        "can't get annotaion for {ID="+str(ID)+", tool="+str(tool)+"}")
            if(estimator == 'ConvEval'):
                data[ID] = ConfEval(dict_index, Tanimotos(
                    dict_inchi)).recommendation()
            elif(estimator == 'Matching'):
                data[ID] = Matching(dict_index, Tanimotos(dict_inchi)).k()
            elif(estimator == 'K'):
                data[ID] = K(dict_index, Tanimotos(dict_inchi)).k()
            elif(estimator == 'K_old'):
                data[ID] = K_old(dict_index, Tanimotos(dict_inchi)).k()
        return(data)


def MultipleSourceAnnotation_to_dataframe(mgf_instance, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res, isdb_res, results):
    if(type(results) != dict):
        Exception(
            'results type is no dict \n please insert a resylt : dict {estimator : results_estimator} \n estimator : Matching, K, K_old, ConvEval')
    if(len(set(results.keys()) - set(['Matching', 'K', 'K_old', 'ConvEval'])) != 0):
        Exception(
            'please insert a resylt : dict {estimator : results_estimator} \n estimator : Matching, K, K_old, ConvEval')
    tool_dict = {}

    # GNPS
    if(get_gnps == True) and (gnps_res != None):
        tool_dict[Tool(1).name] = gnps_res

    # sirius
    if(get_sirius == True) and (sirius_res != None):
        tool_dict[Tool(2).name] = sirius_res

    # isdb
    if(get_isdb == True) and (isdb_res != None):
        tool_dict[Tool(3).name] = isdb_res

    data = []
    for ID in mgf_instance.data:
        inchi_gnps = tool_dict[Tool(1).name][ID].inchi
        score_gnps = tool_dict[Tool(1).name][ID].score
        inchi_sirius = tool_dict[Tool(2).name][ID].inchi
        score_sirius = tool_dict[Tool(2).name][ID].score
        inchi_isdb = tool_dict[Tool(3).name][ID].inchi
        score_isdb = tool_dict[Tool(3).name][ID].score

        dict_inchi = {i: tool_dict[i][ID].inchi for i in tool_dict}
        t = Tanimotos(dict_inchi)
        tgs, tgi, tsi = t.compute_tanimoto()
        l = [ID, inchi_gnps, score_gnps, inchi_sirius,
             score_sirius, inchi_isdb, score_isdb, tgs, tgi, tsi]
        for estimator in results:
            l.append(results[estimator][ID])
        data.append(l)
    df = pd.DataFrame(data, columns=['ID', 'inchi_gnps', 'score_gnps', 'inchi_sirius',
                      'score_sirius', 'inchi_isdb', 'score_isdb', 'tgs', 'tgi', 'tsi']+list(results.keys()))
    return(df)

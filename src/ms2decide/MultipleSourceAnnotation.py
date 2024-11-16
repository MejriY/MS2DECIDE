"""
This module provides functionality to annotate spectra data using multiple sources: GNPS, Sirius, and ISDB-LOTUS. 
It offers methods for combining annotations from different tools and organizing the results in a structured format.

Functions:
    1. MultipleSourceAnnotation: 
       Combines annotation results from multiple sources (GNPS, Sirius, ISDB-LOTUS) and provides a recommendation based on an estimator.
    
    2. MultipleSourceAnnotation_to_dataframe:
       Converts the combined annotation results into a Pandas DataFrame, making it easier to view and analyze the data.

Classes and Libraries Used:
    - Tool: Represents a specific annotation tool.
    - Matching, K, K_old: Estimators used for the recommendation process.
    - ConfEval: Evaluates the consensus of multiple sources.
    - Tanimotos: Computes similarity scores between annotations based on InChI.
    - pandas (pd): Used for data manipulation and organization.

Typical usage example:
    results = MultipleSourceAnnotation(mgf_instance, get_gnps=True, get_sirius=True, get_isdb=False, gnps_res, sirius_res, isdb_res, estimator='Matching')
    df = MultipleSourceAnnotation_to_dataframe(mgf_instance, get_gnps=True, get_sirius=True, get_isdb=False, gnps_res, sirius_res, isdb_res, results)
"""

from .Tool import Tool
from .Matching import Matching
from .K import K
from .K_old import K_old
from .ConfEval import ConfEval
from .Tanimotos import Tanimotos
import pandas as pd


def MultipleSourceAnnotation(mgf_instance, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res, isdb_res, estimator):
    """
    Get aggregated annotations from GNPS, Sirius, and ISDB-LOTUS sources.

    Args:
        mgf_instance (MgfInstance): Instance of MgfInstance containing spectra data.
        get_gnps (bool): Flag to determine whether to get GNPS annotation.
        get_sirius (bool): Flag to determine whether to get Sirius annotation.
        get_isdb (bool): Flag to determine whether to get ISDB-Lotus annotation.
        gnps_res (dict): Result of GNPS annotation (default=None).
        sirius_res (SiriusAnnotation): Result of Sirius annotation (default=None). (sirius_res.data is the dict)
        isdb_res (dict): Result of ISDB-Lotus annotation (default=None).
        estimator (str): ['ConvEval', 'Matching', 'K', 'K_old']
    Returns:
        dict: Dictionary containing combined annotations for each spectrum ID.
    """
    if (estimator not in ['ConvEval', 'Matching', 'K', 'K_old']):
        raise TypeError("error in estimator")
    tool_dict = {}

    # GNPS
    if (get_gnps == True) and (gnps_res != None):
        tool_dict[Tool(1).name] = gnps_res

    # sirius
    if (get_sirius == True) and (sirius_res != None):
        tool_dict[Tool(2).name] = sirius_res.data

    # isdb
    if (get_isdb == True) and (isdb_res != None):
        tool_dict[Tool(3).name] = isdb_res

    if (len(tool_dict) == 0):
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
            if (estimator == 'ConvEval'):
                data[ID] = ConfEval(dict_index, Tanimotos(
                    dict_inchi)).recommendation()
            elif (estimator == 'Matching'):
                data[ID] = Matching(dict_index, Tanimotos(dict_inchi)).k()
            elif (estimator == 'K'):
                data[ID] = K(dict_index, Tanimotos(dict_inchi)).k()
            elif (estimator == 'K_old'):
                data[ID] = K_old(dict_index, Tanimotos(dict_inchi)).k()
        return (data)


def MultipleSourceAnnotation_to_dataframe(mgf_instance, get_gnps, get_sirius, get_isdb, gnps_res, sirius_res, isdb_res, results):
    """
    Converts the combined annotation results into a Pandas DataFrame for easy viewing and analysis.

    Args:
        mgf_instance (MgfInstance): Instance of MgfInstance containing spectra data.
        get_gnps (bool): Flag to determine whether to get GNPS annotation.
        get_sirius (bool): Flag to determine whether to get Sirius annotation.
        get_isdb (bool): Flag to determine whether to get ISDB-LOTUS annotation.
        gnps_res (dict): Result of GNPS annotation.
        sirius_res (SiriusAnnotation): Result of Sirius annotation.
        isdb_res (dict): Result of ISDB-LOTUS annotation.
        results (dict): Combined results of different estimators in the form {estimator: results_estimator}.

    Returns:
        pd.DataFrame: DataFrame containing combined annotations and scores for each spectrum ID.
    """
    if (type(results) != dict):
        Exception(
            'results type is no dict \n please insert a resylt : dict {estimator : results_estimator} \n estimator : Matching, K, K_old, ConvEval')
    if (len(set(results.keys()) - set(['Matching', 'K', 'K_old', 'ConvEval'])) != 0):
        Exception(
            'please insert a resylt : dict {estimator : results_estimator} \n estimator : Matching, K, K_old, ConvEval')
    tool_dict = {}

    # GNPS
    if (get_gnps == True) and (gnps_res != None):
        tool_dict[Tool(1).name] = gnps_res

    # sirius
    if (get_sirius == True) and (sirius_res != None):
        tool_dict[Tool(2).name] = sirius_res.data

    # isdb
    if (get_isdb == True) and (isdb_res != None):
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
        if(inchi_sirius=='#'):
            tgs+=0.25
            tsi+=0.25
            score_sirius=0.5
        if(inchi_gnps=='*'):
            tgs+=0.25
            tgi+=0.25
        l = [ID, inchi_gnps, score_gnps, inchi_sirius,
             score_sirius, inchi_isdb, score_isdb, tgs, tgi, tsi]
        for estimator in results:
            l.append(results[estimator][ID])
        data.append(l)
    df = pd.DataFrame(data, columns=['ID', 'inchi_gnps', 'score_gnps', 'inchi_sirius',
                      'score_sirius', 'inchi_isdb', 'score_isdb', 'tgs', 'tgi', 'tsi']+list(results.keys()))
    return (df)

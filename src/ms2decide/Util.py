import os
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from pathlib import Path
from urllib import request
from matchms.Spectrum import Spectrum
from decimal import Decimal
from matchms.similarity import ModifiedCosine
from matchms import calculate_scores
import warnings


# tool need to be gnps=> for gnps / isdb_pos => for isdb pos / isdb_neg => for isdb neg


def _get_download_link(tool):
    """
    Get the download link for the specified tool, GNPS or ISDB.

    Args:
        tool (str): The tool name.

    Returns:
        str: The download link for the specified tool.
    """
    if (tool == 'gnps'):
        return ('https://gnps-external.ucsd.edu/gnpslibrary/ALL_GNPS.mgf')
    elif (tool == 'isdb_pos'):
        return ('https://zenodo.org/record/7534250/files/isdb_pos.mgf?download=1')
    elif (tool == 'isdb_neg'):
        return ('https://zenodo.org/record/7534250/files/isdb_neg.mgf?download=1')
    else:
        raise Exception('check download link')
        return (None)


def _downolad_file(path, tool):
    """
    Download the file for the specified tool.

    Args:
        tool (str): The tool name.
    """
    remote_url = _get_download_link(tool)
    try:
        local_file = path+tool+'.mgf'
        request.urlretrieve(remote_url, local_file)
    except:
        raise Exception('check the tool should be isdb_pos/ isdb_neg/ gnps ')


def _check_file_and_download(tool):
    """
    Check if the file for the specified tool exists, and download if not.

    Args:
        tool (str): The tool name.

    Returns:
        Path: The path to the downloaded file.
    """
    isdb_path = input(
        'SELECT DIRECTORY OF ISDB DATA \n IF MGF FILE NOT THERE, FILE WILL BE DONWLOADED : \n')
    # if (isdb_path[-2:] != "\\"):
    #     isdb_path += '\\'

    if (tool+'.mgf' in os.listdir(isdb_path)):
        print('==========> FILE ALREADY IN DIRECTORY')
        s = Path(isdb_path+tool+'.mgf')
    else:
        print('==========> DIRECTORY EXIST WITHOUT FILES')
        print('==================')
        print('BEGIN DOWNLOAD')
        _downolad_file(isdb_path, tool)
        print('==========> FILE DOWNLOAD IN DIRECTORY')
        s = Path(isdb_path+tool+'.mgf')
    return (s)


def path_isdb(ion_mod):
    """
    Get the path to the ISDB file based on the ion mode.

    Args:
        ion_mod (str): The ion mode ('pos' or 'neg').

    Returns:
        Path: The path to the ISDB file.
    """
    if (ion_mod == 'pos'):
        return (_check_file_and_download('isdb_'+ion_mod))
    elif (ion_mod == 'neg'):
        return (_check_file_and_download('isdb_'+ion_mod))
    else:
        raise Exception('check ion_mod should be pos or neg')


def path_gnps():
    """
    Get the path to the GNPS file.

    Returns:
        Path: The path to the GNPS file.
    """
    return (_check_file_and_download('gnps'))


def tanimoto(inc1, inc2):
    """
    Compute Tanimoto similarity between two InChI representations.

    Args:
        inc1 (str): The first InChI representation.
        inc2 (str): The second InChI representation.

    Returns:
        float: Tanimoto similarity score.
    """
    if (inc1 == '?') or (inc2 == '?'):
        return (0)
    else:
        mol1 = Chem.MolFromInchi(inc1)
        mol2 = Chem.MolFromInchi(inc2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2,  nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
        s = DataStructs.TanimotoSimilarity(fp1, fp2)
        return s


def _correct_inchi_string(inchi):
    """
    Correct InChI string by adding 'InChI=' if missing.

    Args:
        inchi (str): The InChI string.

    Returns:
        str: The corrected InChI string.
    """
    inchi = str(inchi).strip()
    if (inchi == 'nan') or (inchi == ''):
        return ('*')
    else:
        if (inchi[0] == '"'):
            inchi = inchi[1:]
        if (inchi[-1] == '"'):
            inchi = inchi[:-1]
        if ('InChI=' not in inchi):
            inchi = 'InChI='+inchi
        return (inchi)


def _is_mol(inc):
    """
    Check if the input is a valid molecular structure and return its InChI representation.

    Args:
        inc (str): The InChI or SMILES representation.

    Returns:
        str: The InChI representation of the molecule, or '?' if invalid.
    """
    try:
        mol = Chem.MolFromInchi(inc)
    except:
        try:
            mol = Chem.MolFromSmiles(inc)
        except:
            mol = ''
    if (type(mol) == Chem.rdchem.Mol):
        return (Chem.MolToInchi(mol))
    else:
        return ('?')


def get_correct_inchi(data):
    """
    Get the correct InChI representation from Spectrum data.

    Args:
        data (Spectrum or dict): Spectrum data.

    Returns:
        str: The correct InChI representation.
    """
    if (type(data) == str):
        return ('#')
    elif (type(data) == Spectrum):
        inchi = _correct_inchi_string(data.metadata['inchi'])
        if (str(inchi) == 'nan') or (len(inchi) == 0) or (str(inchi) == '*'):
            try:
                inchi = Chem.MolToInchi(
                    Chem.MolFromSmiles(data.metadata['smiles']))
            except:
                warnings.warn('no valide inchi found.')
                return ('*')
    elif (type(data) == dict):
        inchi = _correct_inchi_string(data['inchi'])
        if (str(inchi) == 'nan') or (len(inchi) == 0) or (str(inchi) == '*'):
            try:
                inchi = Chem.MolToInchi(
                    Chem.MolFromSmiles(data['smiles']))
            except:
                warnings.warn('no valide inchi found.')
                return ('*')
    else:
        raise Exception('no valide valide inpute data')
    return _is_mol(inchi)


def _sirius_score_calcule(score, d):
    """
    Calculate the score based on the specified criterion.

    Args:
        score (str): The scoring criterion ('confidencescore', 'zodiacscore', 'confzodiac', etc.).
        d (dict): A dictionary containing relevant scores.

    Returns:
        float: The calculated score based on the input parameters.
    """
    if (score == 'confidencescore'):
        if (str(d['confidencescore']) == 'nan'):
            if (d['zodiacscore'] == 1):
                return (1)
            else:
                return (0)
        else:
            return (d['confidencescore'])
    elif (score == 'zodiacscore'):
        return (d['zodiacscore'])
    elif (score == 'confzodiac'):
        if (str(d['confidencescore']) == 'nan'):
            if (d['zodiacscore'] == 1):
                return (1)
            else:
                return (0)
        else:
            return (d['confidencescore']*d['zodiacscore'])
    else:
        x = d['csi:fingeridscore']/(d['#predictedfps']*d['formularank'])
        x = 1/(1+np.exp(x))
        con = d['confidencescore']
        if (str(d['confidencescore']) == 'nan'):
            if (d['zodiacscore'] == 1):
                con = 1
        return x*con


def Weights():
    """
    Get default weights and parameters for scoring.

    Returns:
        Tuple: A tuple containing default weights and parameters.
    """
    a_priori_confidences = {'GNPS': 0.6, 'SIRIUS': 0.3, 'ISDB': 0.1}
    weight_tanimoto = 0.5
    weights = [0.7, 0.2, 0.1]
    default_cos = 0.8
    return (a_priori_confidences, weight_tanimoto, weights, default_cos)


def Parametres():
    tolerance, mz_power, intensity_power, shift = 0.02, 0.0, 0.5, 0.0
    return tolerance, mz_power, intensity_power, shift


def owa(weights, values):
    """
    Compute the Ordered Weighted Averaging (OWA) score.

    Args:
        weights (list): List of weights.
        values (list): List of values.

    Returns:
        float: OWA score.
    """
    sorted_values = np.nan_to_num(np.array(sorted(values, reverse=True)))
    # return(sum([weights[i]*sorted_values[i] for i in range(len(weights))]))
    return (float(Decimal(np.dot(np.array(weights), sorted_values))))


def f_T(t, c1, c2, delta):
    """
    Compute a score based on Tanimoto similarity and a confidence function.

    Args:
        t (float): Tanimoto similarity score.
        c1 (float): Confidence for the first tool.
        c2 (float): Confidence for the second tool.
        delta (float): Delta parameter.

    Returns:
        float: Computed score.
    """
    s = (delta/0.2)*((max(c1, 0.8)-0.8) + (max(c2, 0.8)-0.8))*(t >= 0.7)
    return ((delta/0.2) * (min(max(t, 0.5), 0.7)-0.5) + (0.5/0.3)*(max(t, 0.7)-0.7) + s)

########


def Cutoff():
    """
    Returns:
        float: Cutoff of tools.
    """
    plg = 0.4
    pmg = 0.6
    phg = 0.8

    pls = 0.5
    pms = 0.7
    phs = 0.8

    pli = 0.5
    pmi = 0.8
    phi = 0.9
    return plg, pmg, phg, pls, pms, phs, pli, pmi, phi


def Answers():
    """
    Returns:
        float: aswer of our DM
    """
    x1 = 0.5
    x2 = 0.6
    x3 = 0.66
    x4 = 0.6
    x5 = 0.65
    x6 = 0.45
    x7 = 0.61
    x8 = 0.81

    tb_gs = 0.8
    tb_gi = 1
    tb_si = 0.8

    x_gs = 1
    x_gi = 0.9
    x_si = 0.6
    return x1, x2, x3, x4, x5, x6, x7, x8, tb_gs, tb_gi, tb_si, x_gs, x_gi, x_si


def _ratios():
    """
    Compute the ratios of slopes of each function and the ratios of slope between functions

    Returns:
        float: Computed ratios.
    """
    x1, x2, x3, x4, x5, x6, x7, x8, tb_gs, tb_gi, tb_si, x_gs, x_gi, x_si = Answers()
    plg, pmg, phg, pls, pms, phs, pli, pmi, phi = Cutoff()
    dg_ml = (x1+x2-pmg-plg)/(phg-pmg)
    dg_hm = ((x3-pmg)/(1-phg))*(1+((pms-x4)/(pms-pls))*(((phg-pmg) *
                                                         (x1-plg) - (x3-pmg)*(x1+x2-pmg-plg))/((x3-pmg)*(x1+x2-pmg-plg))))

    ds_ml = ((pms-pls)/(phs-pms))*(((x3-pmg)*(x1+x2-pmg-plg)) /
                                   ((phg-pmg)*(x1-plg) - (x3-pmg)*(x1+x2-pmg-plg)))
    ds_hm = ((x5-pms)/(1-phs))+((pms-x4)/(1-phs))*(1/ds_ml)

    di_ml = ((pmi-pli)*(x7-pmg)*(x1+x2-pmg-plg))/((phi-pmi)*(x6-plg)*(phg-pmg))
    di_hm = ((phi-pmi)*(x8-phg)*(x3-pmg)*(phs-x4)) / \
        ((1-phi)*(x7 - pmg)*(1-phg)*(phs-pms))

    alpha_s = (x1+x2-2*plg)/(x1-plg)
    alpha_i = ((phg-pmg)*(x1+x2-2*plg)) / \
        ((x7-pmg)*(x1+x2-pmg-plg)+(phg-pmg)*(x6-plg))

    # this expression without gamma
    deta_g_l = (alpha_s * alpha_i) / \
        ((alpha_s+alpha_i+alpha_s*alpha_i)*(pmg-plg+(phg-pmg)*dg_ml))
    deta_s_l = alpha_i/((alpha_s+alpha_i+alpha_s*alpha_i)
                        * (pms-pls+(phs-pms)*ds_ml))
    deta_i_l = alpha_s/((alpha_s+alpha_i+alpha_s*alpha_i)
                        * (pmi-pli+(phi-pmi)*di_ml))
    ratios_g = [deta_g_l, dg_ml, dg_hm]
    ratios_s = [deta_s_l, ds_ml, ds_hm]
    ratios_i = [deta_i_l, di_ml, di_hm]
    return ratios_g, ratios_s, ratios_i


def _calcule(values):
    c = [0, 0.4, 0.6, 0.8, 1]
    a = [round((values[i]-values[i-1])/((c[i]-c[i-1])), 5)
         for i in range(1, len(values))]
    b = [values[i+1]-a[i]*c[i+1] for i in range(0, len(c)-1)]
    return (a, b)


def _eta(a, b, x):
    if (x <= 0.4):
        return (a[0]*x+b[0])
    if (x <= 0.6):
        return (a[1]*x+b[1])
    if (x <= 0.8):
        return (a[2]*x+b[2])
    else:
        return (a[3]*x+b[3])


def _calcule_matching(gamma, cutoff, ratios):
    deta_l = ratios[0]
    ratios = ratios[1:]
    pentes = [deta_l/8, deta_l, deta_l*ratios[0], deta_l*ratios[0]*ratios[1]]
    a = [i*(1-gamma) for i in pentes]
    b = [-1*a[0]*cutoff[0], -1*a[1]*cutoff[0]]
    b.append(a[1]*cutoff[1]+b[1]-a[2]*cutoff[1])
    b.append(a[2]*cutoff[2]+b[2]-a[3]*cutoff[2])
    return (a, b)


def _eta_matching(a, b, cutoff, x):
    if (x <= cutoff[0]):
        return (a[0]*x+b[0])
    if (x <= cutoff[1]):
        return (a[1]*x+b[1])
    if (x <= cutoff[2]):
        return (a[2]*x+b[2])
    else:
        return (a[3]*x+b[3])


######## local librery search ##

def _find_matches(spec1_mz, spec2_mz, tolerance, shift):
    """
    Find matching peaks between two mass spectra within a specified tolerance.

    Args:
        spec1_mz (numpy.ndarray): Array of m/z values from the first spectrum.
        spec2_mz (numpy.ndarray): Array of m/z values from the second spectrum.
        tolerance (float): The allowed deviation in m/z values for peaks to be considered a match.
        shift (float): An additional value to shift the m/z values of the second spectrum.

    Returns:
        list: A list of tuples, where each tuple contains the indices of matching peaks 
              from the first and second spectra.
    """
    lowest_idx = 0
    matches = []
    for peak1_idx in range(spec1_mz.shape[0]):
        mz = spec1_mz[peak1_idx]
        low_bound = mz - tolerance
        high_bound = mz + tolerance
        for peak2_idx in range(lowest_idx, spec2_mz.shape[0]):
            mz2 = spec2_mz[peak2_idx] + shift
            if mz2 > high_bound:
                break
            if mz2 < low_bound:
                lowest_idx = peak2_idx
            else:
                matches.append((peak1_idx, peak2_idx))
    return matches


def _mass_selection_by_tolerance(sp, query_mass, query, tolerance):
    """
    Select query elements based on their proximity to a specified mass within a given tolerance.

    Args:
        sp (object): An object containing metadata, particularly 'precursor_mz'.
        query_mass (list or numpy.ndarray): Array of mass values to be matched against the precursor mass.
        query (list): List of query objects corresponding to the masses.
        tolerance (float): The maximum allowable difference between the query mass and the precursor mass for a match.

    Returns:
        list: A list of selected query objects that match the specified mass within the tolerance.
    """
    m = sp.metadata['precursor_mz']
    idx = np.argwhere(np.abs(np.array(query_mass)-m) <= tolerance).flatten()
    selected = [query[i] for i in idx]
    return selected


def _get_match(sp, query, tolerance, mz_power, intensity_power, shift):
    """
    Retrieve the best match from a spectrum based on a given query and scoring method.

    Args:
        sp (object): The spectrum object containing mass peaks and metadata.
        query (list): A list of query objects to match against the spectrum.
        tolerance (float): The allowable mass difference for matching peaks.
        mz_power (float): Power applied to mass values for scoring (not directly used here).
        intensity_power (float): Power applied to intensity values for scoring.
        shift (float): Value to shift the mass values during matching.

    Returns:
        list: A list containing the best matching query object and its score.
    """
    score = ModifiedCosine(tolerance=tolerance,
                           intensity_power=intensity_power)
    if (len(query) == 0):
        rsp, res = '#', 0
    else:
        try:
            SCORES = calculate_scores(query, [sp], score)
            rsp, res = SCORES.scores_by_query(
                sp, 'ModifiedCosine_score', sort=True)[0]
        except:
            rsp, res = '#', 0
            for i in query:
                match_mz = _find_matches(
                    sp.peaks.mz, i.peaks.mz, tolerance, shift)
                if (len(match_mz) != 0):
                    res_score = score.pair(i, sp)['score']
                    print(i)
                    if (res_score > res):
                        res = res_score
                        rsp = i
    return [rsp, res]

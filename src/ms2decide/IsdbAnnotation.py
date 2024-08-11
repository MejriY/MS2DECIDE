from .MatchedSpectra import MatchedSpectra
from .Util import path_isdb, Parametres, get_correct_inchi, _mass_selection_by_tolerance, _get_match
from matchms.importing import load_from_mgf
from matchms import Spectrum
import numpy as np
import os
import sys
current = os.path.dirname(os.path.realpath('__file__'))
parent = os.path.dirname(current)
sys.path.append(parent)


def get_cfm_annotation(mgf_instance):
    """
    Returns ISDB annotations based on ISDB subdivision.

    Args:
        mgf_instance (MgfInstance): An instance of MgfInstance containing Mass Spectrometry data.

    Returns:
        dict: A dictionary where the key is the spectrum ID and the value is a MatchedSpectra object with ISDB annotation.
    """
    return(_get_cfm_annotation(mgf_instance))


def _get_cfm_annotation(mgf, ion_mode='pos'):
    isdb_mass, isdb = _load_isdb(ion_mode)
    tolerance, mz_power, intensity_power, shift = Parametres()
    print('==================')
    print('cfm file is loaded')
    print('==================')
    RES = {}
    for i in mgf.data:
        sp = mgf.data[i]
        selected = _mass_selection_by_tolerance(sp, isdb_mass, isdb, tolerance)
        rsp, res = _get_match(sp, selected, tolerance,
                              mz_power, intensity_power, shift)
        if(type(rsp) == Spectrum):
            if(type(res)==np.ndarray):  
                res=float(res['score'])
            else:
                res=res[0]
            RES[i] = MatchedSpectra(i, get_correct_inchi(rsp), res)
        else:
            RES[i] = MatchedSpectra(i, rsp, res)
    return RES


def _load_isdb(ion_mode):
    path_isdb_data = path_isdb(ion_mode)
    isdb_data = load_from_mgf(str(path_isdb_data))
    mass = []
    isdb = []
    for i in isdb_data:
        mass.append(i.metadata['precursor_mz'])
        isdb.append(i)
    return mass, isdb

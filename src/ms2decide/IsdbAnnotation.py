"""
This module is designed to provide ISDB-Lotus (In Silico DataBase) annotations for spectra data using the CFM (Competitive Fragmentation Modeling) method.
It handles mass spectrometry data, identifies matches based on mass tolerance, and returns structured annotations.

Functions:
    1. get_cfm_annotation:
       Public function that returns ISDB-Lotus annotations for a given mgf_instance.

    2. _get_cfm_annotation:
       Internal function that processes the mass spectrometry data from an MGF instance and matches it against the ISDB database. 
       It returns the best matching spectra based on the parameters and tolerance.

    3. _load_isdb:
       Loads the ISDB-Lotus data from the specified file based on the ion mode (positive or negative) and organizes it by precursor mass for easy lookup.

Classes and Libraries Used:
    - MatchedSpectra: Represents the matched spectra for a given spectrum ID.
    - Parametres: Provides parameters such as mass tolerance, mz_power, intensity_power, and shift used in the matching process.
    - matchms.Spectrum: Represents a mass spectrum from the `matchms` library.
    - numpy (np): Used for handling numerical data like arrays.

Typical usage example:
    annotations = get_cfm_annotation(mgf_instance)
"""

from .MatchedSpectra import MatchedSpectra
from .Util import path_isdb, Parametres, get_correct_inchi, _mass_selection_by_tolerance, _get_match, _downolad_file
from matchms.importing import load_from_mgf
from matchms import Spectrum
import numpy as np
import os
import sys
current = os.path.dirname(os.path.realpath('__file__'))
parent = os.path.dirname(current)
sys.path.append(parent)


def get_cfm_annotation(mgf_instance, tol):
    """
    Returns ISDB-Lotus annotations based on ISDB-Lotus subdivision.

    Args:
        mgf_instance (MgfInstance): An instance of MgfInstance containing Mass Spectrometry data.

    Returns:
        dict: A dictionary where the key is the spectrum ID and the value is a MatchedSpectra object with ISDB-Lotus annotation.
    """
    ion_mode = input(
        'SELECT ION MODE FOR ISDB-LOTUS ANNOTATION (POS for positive and NEG for negative)')
    ion_mode = ion_mode.lower()
    return (_get_cfm_annotation(mgf_instance, ion_mode, tol))


def _get_cfm_annotation(mgf, ion_mode, tol):
    """
    Internal function to get ISDB-Lotus annotations for the given MGF instance.

    Args:
        mgf (MgfInstance): An instance of MgfInstance containing spectra data.
        ion_mode (str, optional): Specifies the ionization mode ('pos' for positive, 'neg' for negative). Defaults to 'pos'.

    Returns:
        dict: A dictionary where the key is the spectrum ID and the value is a MatchedSpectra object with ISDB-Lotus annotation.
    """
    isdb_mass, isdb = _load_isdb(ion_mode)
    tolerance, mz_power, intensity_power, shift = Parametres()
    tolerance = tol

    print('==================')
    print('cfm file is loaded')
    print('==================')
    RES = {}
    for i in mgf.data:
        sp = mgf.data[i]
        selected = _mass_selection_by_tolerance(sp, isdb_mass, isdb, tolerance)
        rsp, res = _get_match(sp, selected, tolerance,
                              mz_power, intensity_power, shift)
        if (type(rsp) == Spectrum):
            if (type(res) == np.ndarray):
                res = float(res['score'])
            else:
                res = res[0]
            RES[i] = MatchedSpectra(i, get_correct_inchi(rsp), res)
        else:
            RES[i] = MatchedSpectra(i, rsp, res)
    return RES


def _load_isdb(ion_mode):
    """
    Loads ISDB-Lotus data for the specified ionization mode ('pos' or 'neg').

    Args:
        ion_mode (str): The ionization mode ('pos' for positive, 'neg' for negative).

    Returns:
        tuple: A tuple containing two lists:
            - mass (list): Precursor masses from the ISDB-Lotus data.
            - isdb (list): List of ISDB-Lotus spectra.
    """
    path_isdb_data = path_isdb(ion_mode)
    isdb_data = load_from_mgf(str(path_isdb_data))
    mass = []
    isdb = []
    for i in isdb_data:
        mass.append(i.metadata['precursor_mz'])
        isdb.append(i)
    return mass, isdb


def get_cfm_annotation_GUI(path_isdb, mgf, ion_mode='pos', tol=0.02):
    """
    Internal function to get ISDB-Lotus annotations for the given MGF instance.

    Args:
        mgf (MgfInstance): An instance of MgfInstance containing spectra data.
        ion_mode (str, optional): Specifies the ionization mode ('pos' for positive, 'neg' for negative). Defaults to 'pos'.

    Returns:
        dict: A dictionary where the key is the spectrum ID and the value is a MatchedSpectra object with ISDB-Lotus annotation.
    """
    if (str(path_isdb)[-4:] != '.mgf'):
        tool = 'isdb_'+ion_mode
        _downolad_file(path_isdb, tool)
        path_isdb = str(path_isdb)+'\\'+tool+'.mgf'
    isdb_data = load_from_mgf(str(path_isdb))
    isdb_mass = []
    isdb = []
    for i in isdb_data:
        isdb_mass.append(i.metadata['precursor_mz'])
        isdb.append(i)

    tolerance, mz_power, intensity_power, shift = Parametres()
    tolerance = tol

    print('==================')
    print('cfm file is loaded')
    print('==================')
    RES = {}
    for i in mgf.data:
        sp = mgf.data[i]
        selected = _mass_selection_by_tolerance(sp, isdb_mass, isdb, tolerance)
        rsp, res = _get_match(sp, selected, tolerance,
                              mz_power, intensity_power, shift)
        if (type(rsp) == Spectrum):
            if (type(res) == np.ndarray):
                res = float(res['score'])
            else:
                res = res[0]
            RES[i] = MatchedSpectra(i, get_correct_inchi(rsp), res)
        else:
            RES[i] = MatchedSpectra(i, rsp, res)
    return RES

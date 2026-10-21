from .MatchedSpectra import MatchedSpectra
from .MgfInstance import MgfInstance
from .Util import get_correct_inchi
from pathlib import Path
import pandas as pd


_EMPTY_VALUES = {'', 'nan', 'none', 'n/a', 'na', 'null', '*', '#', '?'}


def _clean_text(value):
    """Return '' for empty-like values (NaN, 'N/A', ...), the stripped text otherwise."""
    if value is None:
        return ''
    text = str(value).strip()
    return '' if text.lower() in _EMPTY_VALUES else text


def read_gnps2_annotations(annotation_file):
    """
    Read the library-search table of a GNPS2 FBMN job.

    The expected file is the one downloaded from the finished GNPS2 task
    (<task id>-merged_results_with_gnps.tsv). It needs at least the columns
    '#Scan#' (feature ID = SCANS of the MGF file) and 'MQScore' (cosine score).
    The columns 'INCHI' and 'Smiles' are used when present.

    Args:
        annotation_file (str or Path): Path of the .tsv file.

    Returns:
        dict: {scan (int): {'inchi': str, 'smiles': str, 'mqscore': float}}
    """
    annotation_file = Path(str(annotation_file).strip().strip('"\''))
    if not annotation_file.exists():
        raise FileNotFoundError(f"GNPS2 results file not found: {annotation_file}")

    df = pd.read_csv(annotation_file, sep='\t', low_memory=False)
    by_lower = {c.strip().lower(): c for c in df.columns}
    missing = [c for c in ('#scan#', 'mqscore') if c not in by_lower]
    if missing:
        raise ValueError(
            f"Column(s) {missing} not found in {annotation_file.name}. "
            f"Columns found: {list(df.columns)}")

    scan_col = by_lower['#scan#']
    score_col = by_lower['mqscore']
    inchi_col = by_lower.get('inchi')
    smiles_col = by_lower.get('smiles')

    df = df.copy()
    df[scan_col] = pd.to_numeric(df[scan_col], errors='coerce')
    df[score_col] = pd.to_numeric(df[score_col], errors='coerce')
    df = df.dropna(subset=[scan_col, score_col])
    # One annotation per spectrum: keep the best cosine score.
    df = df.sort_values(score_col, ascending=False).drop_duplicates(
        subset=[scan_col], keep='first')

    annotations = {}
    for _, row in df.iterrows():
        annotations[int(row[scan_col])] = {
            'inchi': _clean_text(row[inchi_col]) if inchi_col else '',
            'smiles': _clean_text(row[smiles_col]) if smiles_col else '',
            'mqscore': float(row[score_col]),
        }
    return annotations


def closest_gnps2_from_file(annotation_file, mgf_path):
    """
    Returns a dictionary of MatchedSpectra objects built from the results of a
    Feature-Based Molecular Networking job run on GNPS2.

    The GNPS2 job is run manually on https://gnps2.org with the same MGF file
    as the one given to MS2DECIDE; its library-search table is then downloaded
    and used here instead of the results of a classic GNPS job.

    Args:
        annotation_file (str or Path): GNPS2 '...merged_results_with_gnps.tsv' file.
        mgf_path (str or Path): The local path of the .mgf file.

    Returns:
        dict: Dictionary {id: MatchedSpectra}.
    """
    annotations = read_gnps2_annotations(annotation_file)
    dict_data = MgfInstance(Path(mgf_path)).data

    matched_spectra_dict = {}
    n_matched = 0
    for ID in dict_data:
        if ID in annotations:
            matched_spectra_dict[ID] = MatchedSpectra(
                dict_data[ID], get_correct_inchi(annotations[ID]), annotations[ID]['mqscore'])
            n_matched += 1
        else:
            matched_spectra_dict[ID] = MatchedSpectra(dict_data[ID], "#", 0)

    not_in_mgf = [i for i in annotations if i not in dict_data]
    print('==================')
    print(f" GNPS2 results: {len(annotations)} library annotations read, "
          f"{n_matched} matched to the {len(dict_data)} spectra of the MGF file")
    if not_in_mgf:
        print(f" WARNING: {len(not_in_mgf)} annotated scans are not in the MGF file "
              f"(first ones: {not_in_mgf[:5]}). Check that the GNPS2 job used this MGF file.")
    if annotations and n_matched == 0:
        raise ValueError(
            "None of the GNPS2 annotations matches the scans of the MGF file. "
            "The GNPS2 job was probably run with another MGF file.")
    print('==================')
    return matched_spectra_dict

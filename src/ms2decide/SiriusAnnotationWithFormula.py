import pandas as pd
import os
from .MatchedSpectra import MatchedSpectra
from .Util import get_correct_inchi


def get_sirius_by_id(ID, racine_path, Infolder):
    racine_path = racine_path+'\\'+Infolder[ID]
    if ('structure_candidates.tsv' in os.listdir(racine_path)):
        df = pd.read_table(racine_path+'\\structure_candidates.tsv')
        if (len(df) != 0):
            df2 = df[df.adduct == '[M + H]+']
            if (len(df2) != 0):
                try:
                    data_ann = df2[df2.ConfidenceScore == df2.ConfidenceScore.max()].rename(
                        columns=str.lower)
                    sc = data_ann.confidencescore.values.tolist()[0]
                    data_ann = data_ann.to_dict('index')
                    data_ann = data_ann[list(data_ann.keys())[0]]
                except:
                    sc = 0
                    data_ann = df2.iloc[0].rename(str.lower).to_dict()
            else:
                try:
                    data_ann = df[df.ConfidenceScore == df.ConfidenceScore.max()].rename(
                        columns=str.lower)
                    sc = data_ann.confidencescore.values.tolist()[0]*0.75
                    data_ann = data_ann.to_dict('index')
                    data_ann = data_ann[list(data_ann.keys())[0]]
                except:
                    data_ann = df.iloc[0].rename(str.lower).to_dict()
                    sc = 0
        else:
            df = pd.read_table(racine_path+'\\formula_candidates.tsv')
            try:
                df = df[df.explainedIntensity == df.explainedIntensity.max()]
                sc, data_ann = df[['explainedIntensity', 'molecularFormula']].values.tolist()[
                    0]
                if (df.adduct.values.tolist()[0] != '[M + H]+'):
                    sc = sc*0.75
            except:
                sc = 0
                data_ann = '#'
    else:
        df = pd.read_table(racine_path+'\\formula_candidates.tsv')
        try:
            df = df[df.explainedIntensity == df.explainedIntensity.max()]
            sc, data_ann = df[['explainedIntensity', 'molecularFormula']].values.tolist()[
                0]
            if (df.adduct.values.tolist()[0] != '[M + H]+'):
                sc = sc*0.75
        except:
            sc = 0
            data_ann = '#'
    return (sc, data_ann)


def SiriusAnnotationWithFormula(racine_path, mgf_instance):
    data_calcule = {}
    data_anno = {}
    l = [i for i in os.listdir(racine_path) if (
        ('.tsv' not in i) and ('.mztab' not in i))]
    d = {int(i.split('_')[3]): i for i in l}

    for ID in mgf_instance.data:
        if (ID in d):
            sc, data_ann = get_sirius_by_id(ID, racine_path, d)
            if (type(data_ann) != str):
                inchi_calcule = get_correct_inchi(data_ann)
                inchi_anno = inchi_calcule
            else:
                inchi_calcule = get_correct_inchi(data_ann)
                inchi_anno = data_ann
            if (inchi_calcule != '#'):
                data_calcule[ID] = MatchedSpectra(
                    mgf_instance.data[ID], inchi_calcule, sc)
                data_anno[ID] = MatchedSpectra(
                    mgf_instance.data[ID], inchi_anno, sc)
            else:
                data_calcule[ID] = MatchedSpectra(
                    mgf_instance.data[ID], '#', sc)
                data_anno[ID] = MatchedSpectra(
                    mgf_instance.data[ID], inchi_anno, sc)
        else:
            data_calcule[ID] = MatchedSpectra(mgf_instance.data[ID], '#', sc)
            data_anno[ID] = MatchedSpectra(
                mgf_instance.data[ID], inchi_anno, sc)
    return (data_anno, data_calcule)

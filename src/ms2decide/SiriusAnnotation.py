import pandas as pd
from .MatchedSpectra import MatchedSpectra
from .Util import get_correct_inchi, _sirius_score_calcule


class SiriusAnnotation():
    """
    Represents Sirius annotations for Mass Spectrometry data.

    Attributes:
        path (str): The path to the Sirius annotation file.
        data (dict): A dictionary containing MatchedSpectra objects with Sirius annotations.

    Methods:
        __init__(path, mgf_instance): Initializes SiriusAnnotation object.
    """

    def __init__(self, path, mgf_instance, indexScore):
        """
        Initializes a SiriusAnnotation object.

        Args:
            path (str): The path to the Sirius annotation project.
            mgf_instance (MgfInstance): An instance of MgfInstance containing Mass Spectrometry data.
        """
        self.path = path
        indexScore = str(indexScore).lower()
        while (indexScore not in ['exact', 'approximate']):
            indexScore = input('Select index score =  exact, approximate')
            indexScore = str(indexScore).lower()
        if (indexScore == 'exact'):
            indexScore = 'confidencescoreexact'
        else:
            indexScore = 'confidencescoreapproximate'

        data_sirius = pd.read_table(path)
        data_sirius = data_sirius.rename(columns=str.lower)
        # data_sirius['#scan#'] = data_sirius['id'].apply(lambda x: int(x.split('_')[0]))
        data_sirius = data_sirius.set_index(
            "mappingfeatureid").to_dict('index')

        data = {}
        for ID in mgf_instance.data:
            if (ID in data_sirius):
                try:
                    sc = data_sirius[ID][indexScore]
                except:
                    sc = 0
                data[ID] = MatchedSpectra(
                    mgf_instance.data[ID], get_correct_inchi(data_sirius[ID]), sc)
            else:
                data[ID] = MatchedSpectra(mgf_instance.data[ID], "#", 0)
        self.data = data

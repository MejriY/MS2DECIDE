from matchms.importing import load_from_mgf
from matchms.filtering import normalize_intensities
from pathlib import Path


class MgfInstance():
    """
    Represents an instance for handling Mass Spectrometry data in MGF format.

    This class provides methods to initialize the instance and subdivide the loaded spectra.

    Attributes:
        data (dict): A dictionary containing spectra data with spectrum IDs as keys.

    Methods:
        __init__(data_input):
            Initializes an MgfInstance with data from either a file or a list of spectra.

        subdiv_file(step: int) -> List[List[Spectrum]]:
            Subdivides the loaded spectra into smaller chunks.

    Example:
        Initializing an MgfInstance from a file:

        >>> mgf_instance = MgfInstance(Path('path/to/your/file.mgf'))

        Initializing an MgfInstance from a list of spectra:

        >>> spectra_list = [...]  # Replace with actual list of spectra
        >>> mgf_instance = MgfInstance(spectra_list)

        Subdividing the loaded spectra into chunks of 100:

        >>> subdivided_data = mgf_instance.subdiv_file(100)
    """

    def __init__(self, data_input):
        """
        Initializes an MgfInstance with data from either a file or a list of spectra.

        Args:
            data_input: A Path object pointing to an MGF file or a list of Spectrum objects.
        """
        self.data = {}
        if (isinstance(data_input, Path) == True):
            data = load_from_mgf(str(data_input))
            # set_spectra = [add_precursor_mz(i) for i in data]
            # set_spectra = [normalize_intensities(i) for i in set_spectra]
            # set_spectra = [set_spectra[i].set('id', i) for i in range(len(set_spectra))]
            # dict_spectra = {i.metadata['id']: i for i in set_spectra}
            try:
                # dict_spectra = {int(i.metadata['scans']): normalize_intensities(i) for i in data}
                dict_spectra = {int(i.metadata['scans']): i for i in data}
                self.data = dict_spectra
            except:
                raise Exception('mgf file not form cytoscape')
        elif (type(data_input) == list):
            # set_spectra = [add_precursor_mz(i) for i in data_input]
            # set_spectra = [set_spectra[i].set('id', i) for i in range(len(set_spectra))]
            try:
                dict_spectra = {
                    int(i.metadata['scans']): normalize_intensities(i) for i in data}
                self.data = dict_spectra
            except:
                raise Exception('mgf file not form cytoscape')

    def subdiv_file(self, step):
        """
        Subdivides the loaded spectra into smaller chunks.

        Args:
            step (int): The size of each chunk.

        Returns:
            List[List[Spectrum]]: A list of lists containing subdivided spectra data.
        """
        db_div = []
        k = 0
        data = list(self.data.values())
        while ((k+1)*step < len(data)):
            db_div.append(data[k*step:(k+1)*step])
            k += 1
        db_div.append(data[k*step:])
        return (db_div)

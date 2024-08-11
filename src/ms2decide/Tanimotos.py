from rdkit import Chem
from itertools import combinations
from .Tool import Tool
from .Util import tanimoto


class Tanimotos():
    """
    Computes Tanimoto similarity between different tools based on their InChI annotations.

    This class takes a dictionary of tool annotations (InChI representations) and computes Tanimoto
    similarity scores between pairs of tools.

    Attributes:
        data (dict): A dictionary containing tool annotations with tool names as keys.

    Methods:
        __init__(dict_annotations: dict):
            Initializes a Tanimotos instance with a dictionary of tool annotations.

        compute_tanimoto() -> Tuple[float, float, float]:
            Computes Tanimoto similarity scores between pairs of tools and returns the results.

    Example:
        Initializing a Tanimotos instance with tool annotations:

        >>> tool_annotations = {'GNPS': 'InChI=1S/C17H14N2/c1-10-12-7-8-18-11(2)14(12)9-15-13-5-3-4-6-16(13)19-17(10)15/h3-9,19H,1-2H3',
        >>>                     'SIRIUS': 'InChI=1S/C12H20N2O/c1-5-14(6-2)10-7-8-11(13(3)4)12(15)9-10/h7-9,15H,5-6H2,1-4H3',
        >>>                     'ISDB': 'InChI=1S/C10H11N/c1-7-4-3-5-10-9(7)6-8(2)11-10/h3-6,11H,1-2H3'}
        >>> tanimoto_instance = Tanimotos(tool_annotations)

        Computing Tanimoto similarity scores:

        >>> tgs, tgi, tsi = tanimoto_instance.compute_tanimoto()
    """

    def __init__(self,  dict_annotations):
        """
        Initializes a Tanimotos instance with a dictionary of tool annotations.

        Args:
            dict_annotations (dict): A dictionary containing tool annotations with tool names as keys.
        """
        d = {i: 0 for i in dict_annotations}
        for i in d:
            # TODO Je ne comprends pas ce code, mais il faut éviter de coder des absences de valeur par le string 'nan' (bien qu'au fromage, ce soit bon).
            # TODO the values '#' represent miss match in isdb
            if(dict_annotations[i] == '#') or (dict_annotations[i] == '?') or (dict_annotations[i] == '*'):
                d[i] = '?'
            elif(str(dict_annotations[i]) != 'nan'):
                try:
                    mol = Chem.MolFromInchi(dict_annotations[i])
                    inchi = Chem.MolToInchi(mol)
                    d[i] = inchi
                except:
                    raise Exception('check inchi for tool '+i)
        self.data = d
        self.averages = None
        self.tgs = None
        self.tgi = None
        self.tsi = None

    def _compute_tanimoto(self):
        """
        Computes Tanimoto similarity scores between pairs of tools.
        """
        tools = list(self.data.keys())
        couples = list(combinations(tools, 2))
        data = self.data
        for i in couples:
            if(i == ('GNPS', 'SIRIUS')) or (i == ('SIRIUS', 'GNPS')):
                self.tgs = tanimoto(data[i[0]], data[i[1]])
            if(i == ('GNPS', 'ISDB')) or (i == ('ISDB', 'GNPS')):
                self.tgi = tanimoto(data[i[0]], data[i[1]])
            if(i == ('ISDB', 'SIRIUS')) or (i == ('SIRIUS', 'ISDB')):
                self.tsi = tanimoto(data[i[0]], data[i[1]])

    def compute_tanimoto(self):
        """
        Computes Tanimoto similarity scores between pairs of tools and returns the results.

        Returns:
            Tuple[float, float, float]: Tanimoto similarity scores for ('GNPS', 'SIRIUS'), ('GNPS', 'ISDB'),
            and ('ISDB', 'SIRIUS') pairs.
        """
        if(self.tgs == None) or (self.tgi == None) or (self.tsi == None):
            self._compute_tanimoto()
        return(self.tgs, self.tgi, self.tsi)

    def _averages(self):
        d = {}
        self._compute_tanimoto()
        for i in self.data:
            if i == Tool(1).name:
                d[i] = (self.tgs+self.tgi)/2
            elif i == Tool(2).name:
                d[i] = (self.tgs+self.tsi)/2
            else:
                d[i] = (self.tgi+self.tsi)/2
        self.averages = d

    def average_tool(self, tool):
        if(self.averages == None):
            self._averages()
        return(self.averages[tool])

from .Tool import Tool
from .Util import _calcule, _eta, f_T


class K_old():
    """
    Estimator_k class calculates an aggregated score (k score) based on
    similarities scores and Tanimoto values for different tools.

    Args:
        similarities (dict): Dictionary of similarity scores values for each tool.
        tanimotos (Tanimotos): Instance of Tanimotos class.
    """

    def __init__(self, similarities, tanimotos):
        """
        Initialize Estimator_k with similarity scores and Tanimoto values.

        Args:
            similarities (dict): Dictionary of similarity scores values for each tool.
            tanimotos (Tanimotos): Instance of Tanimotos class.
        """
        if (set(similarities.keys()) == set(tanimotos.data.keys())):
            self.similarities = similarities
            self.tanimotos = tanimotos
        else:
            raise Exception(
                'similarity dict and tanimoto dict should have the same keys')

        Ag, Bg = _calcule([0, 0.5, 2.5, 7.5, 8.9])
        As, Bs = _calcule([0, 0.5, 1.5, 3.5, 4.1])
        Ai, Bi = _calcule([0, 0.5, 0.5, 3, 3.5])
        self.Ag = Ag
        self.Bg = Bg
        self.As = As
        self.Bs = Bs
        self.Ai = Ai
        self.Bi = Bi

    def eta_g(self, x):
        return _eta(self.Ag, self.Bg, x)

    def eta_s(self, x):
        return _eta(self.As, self.Bs, x)

    def eta_i(self, x):
        return _eta(self.Ai, self.Bi, x)

    def k(self):
        """
        Calculate the aggregated k score based on similarities scores and Tanimotos.

        Returns:
            float: Aggregated k score.
        """
        tgs, tgi, tsi = self.tanimotos.compute_tanimoto()
        cg, cs, ci = self.similarities[Tool(1).name], self.similarities[Tool(
            2).name], self.similarities[Tool(3).name]

        g = self.eta_g(cg)
        s = self.eta_s(cs)
        i = self.eta_i(ci)

        gs = f_T(tgs, cg, cs, 3)
        gi = f_T(tgi, cg, ci, 2)
        si = f_T(tsi, cs, ci, 1.5)

        return (g+s+i+gs+gi+si)

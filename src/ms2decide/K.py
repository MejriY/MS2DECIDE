import os
import sys
from .Tool import Tool
from .Util import _calcule, _eta


class K():
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

        Ag,Bg = _calcule([-4,0,8,24,50])
        As,Bs = _calcule([-1,0,2,6,12.5])
        Ai,Bi = _calcule([-0.5,0,1,3,6.25])
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

    def _psi_g(self, c):
        return ( self.eta_g(c)-self.eta_g(0) ) / ( self.eta_g(1)-self.eta_g(0) )

    def _psi_s(self, c):
        return ( self.eta_s(c)-self.eta_s(0) ) / ( self.eta_s(1)-self.eta_s(0) )

    def _psi_i(self, c):
        return ( self.eta_i(c)-self.eta_i(0) ) / ( self.eta_i(1)-self.eta_i(0) )

    def phi_gs(self, t, cg, cs):
        if (t <= 0.5):
            return (0)
        else:
            return ( (t-0.5) / 0.2 ) * ( (self._psi_g(cg)+self._psi_s(cs))/(self._psi_g(0.8)+self._psi_s(0.8)) ) * 10

    def phi_gi(self, t, cg, ci):
        if (t <= 0.5):
            return (0)
        else:
            return ( (t-0.5) / 0.2 ) * ( (self._psi_g(cg)+self._psi_i(ci))/(self._psi_g(0.8)+self._psi_i(0.8)) ) * 6

    def phi_si(self, t, cs, ci):
        if (t <= 0.5) :
            return (0)
        else:
            return ( (t-0.5) / 0.2 ) * ( (self._psi_s(cs)+self._psi_i(ci))/(self._psi_s(0.8)+self._psi_i(0.8)) ) * 15

    def k(self):
        """
        Calculate the aggregated k score based on similarities scores and Tanimotos.

        Returns:
            float: Aggregated k score.
        """
        tgs, tgi, tsi = self.tanimotos.compute_tanimoto()
        cg, cs, ci = self.similarities[Tool(1).name], self.similarities[Tool(2).name], self.similarities[Tool(3).name]
        if(self.tanimotos.data['SIRIUS']=='#'):
            tgs+=0.25
            tsi+=0.25
            cs=0.5
        if(self.tanimotos.data['GNPS']=='*'):
            tgs+=0.25
            tgi+=0.25
        g = self.eta_g(cg)
        s = self.eta_s(cs)
        i = self.eta_i(ci)

        gs = self.phi_gs(tgs, cg, cs)
        gi = self.phi_gi(tgi, cg, ci)
        si = self.phi_si(tsi, cs, ci)

        return (g+s+i+gs+gi+si)

from .Tool import Tool
from .Util import _calcule_matching, _ratios, Cutoff, _eta_matching, Answers


class Matching():
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
        ratios_g, ratios_s, ratios_i = _ratios()
        plg, pmg, phg, pls, pms, phs, pli, pmi, phi = Cutoff()
        x1, x2, x3, x4, x5, x6, x7, x8, tb_gs, tb_gi, tb_si, x_gs, x_gi, x_si = Answers()
        Ag, Bg = _calcule_matching(0.33, [plg, pmg, phg], ratios_g)
        As, Bs = _calcule_matching(0.33, [pls, pms, phs], ratios_s)
        Ai, Bi = _calcule_matching(0.33, [pli, pmi, phi], ratios_i)
        self.Ag = Ag
        self.Bg = Bg
        self.As = As
        self.Bs = Bs
        self.Ai = Ai
        self.Bi = Bi

        self.tb_gs = tb_gs
        self.tb_gi = tb_gi
        self.tb_si = tb_si
        self.x_gs = x_gs
        self.x_gi = x_gi
        self.x_si = x_si

        self.plg = plg
        self.pmg = pmg
        self.phg = phg
        self.pls = pls
        self.pms = pms
        self.phs = phs
        self.pli = pli
        self.pmi = pmi
        self.phi = phi

    def eta_g(self, x):
        return _eta_matching(self.Ag, self.Bg, [self.plg, self.pmg, self.phg], x)

    def eta_s(self, x):
        return _eta_matching(self.As, self.Bs, [self.pls, self.pms, self.phs], x)

    def eta_i(self, x):
        return _eta_matching(self.Ai, self.Bi, [self.pli, self.pmi, self.phi], x)

    def a_gs(self, x, y):
        s = (self.eta_g(x)-self.eta_g(self.plg))/(
            self.eta_g(self.phg)-self.eta_g(self.plg))
        s += (self.eta_s(y)-self.eta_s(self.pls))/(
            self.eta_s(self.phs)-self.eta_s(self.pls))
        return s/2

    def a_gi(self, x, y):
        s = (self.eta_g(x)-self.eta_g(self.plg))/(
            self.eta_g(self.phg)-self.eta_g(self.plg))
        s += (self.eta_i(y)-self.eta_i(self.pli))/(
            self.eta_i(self.phi)-self.eta_i(self.pli))
        return s/2

    def a_si(self, x, y):
        s = (self.eta_s(x)-self.eta_s(self.pls))/(
            self.eta_s(self.phs)-self.eta_s(self.pls))
        s += (self.eta_i(y)-self.eta_i(self.pli))/(
            self.eta_i(self.phi)-self.eta_i(self.pli))
        return s/2

    def _psi_gs(self, t, x):
        return (((t-0.5)*x)/((self.tb_gs-0.5)*self.a_gs(self.phg, self.phs)))*(self.eta_g(self.x_gs)-self.eta_g(self.pmg)-self.eta_s(self.pms))

    def _psi_gi(self, t, x):
        return (((t-0.5)*x)/((self.tb_gi-0.5)*self.a_gi(self.phg, self.phi)))*(self.eta_g(self.x_gi)-self.eta_g(self.pmg)-self.eta_i(self.pmi))

    def _psi_si(self, t, x):
        return (((t-0.5)*x)/((self.tb_si-0.5)*self.a_si(self.phs, self.phi)))*(self.eta_g(self.x_si)-self.eta_s(self.phs)-self.eta_i(self.phi))

    def phi_gs(self, t, c_g, c_s):
        if (t <= 0.5) or (c_g < self.plg) or (c_s < self.pls):
            return (0)
        else:
            c = self.a_gs(c_g, c_s)
            return self._psi_gs(t, c)

    def phi_gi(self, t, c_g, c_i):
        if (t <= 0.5) or (c_g < self.plg) or (c_i < self.pli):
            return (0)
        else:
            c = self.a_gi(c_g, c_i)
            return self._psi_gi(t, c)

    def phi_si(self, t, c_s, c_i):
        if (t <= 0.5) or (c_s < self.pls) or (c_i < self.pli):
            return (0)
        else:
            c = self.a_si(c_s, c_i)
        return self._psi_si(t, c)

    def k(self):
        """
        Calculate the aggregated k score based on similarities scores and Tanimotos.

        Returns:
            float: Aggregated k score.
        """
        g = self.eta_g(self.similarities[Tool(1).name])
        s = self.eta_s(self.similarities[Tool(2).name])
        i = self.eta_i(self.similarities[Tool(3).name])

        tgs, tgi, tsi = self.tanimotos.compute_tanimoto()

        gs = self.phi_gs(tgs, self.similarities[Tool(
            1).name], self.similarities[Tool(2).name])
        gi = self.phi_gi(tgi, self.similarities[Tool(
            1).name], self.similarities[Tool(3).name])
        si = self.phi_si(tsi, self.similarities[Tool(
            2).name], self.similarities[Tool(3).name])

        return (g+s+i+gs+gi+si)

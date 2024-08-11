from unittest.mock import MagicMock
import unittest
import os
import sys
from pathlib import Path
from K import K
from Tanimotos import Tanimotos


class TestMultipleSourceAnnotation(unittest.TestCase):

    def test_k(self):
        dict_index = {'GNPS': 0.200693,
                      'ISDB': 0.6242559705119061, 'SIRIUS': 0.3395112597307259}
        dict_inchi = {'GNPS': 'InChI=1S/C43H54N4O6/c1-7-23-17-25-20-43(22-53-52-6)39-28(15-16-47(40(23)43)41(25)48)27-13-14-34(50-4)36(38(27)45-39)31-18-29-24(8-2)21-46(3)33(35(29)42(49)51-5)19-30-26-11-9-10-12-32(26)44-37(30)31/h9-14,23-25,29,31,33,35,40-41,44-45,48H,7-8,15-21H2,1-6H3/t23-,24+,25-,29-,31-,33-,35?,40-,41?,43?/m0/s1',
                      'ISDB': 'InChI=1S/C23H27N2O3/c1-25-8-7-23-16-4-2-3-5-17(16)24-20(26)11-18-21(22(23)24)15(10-19(23)28-13-25)14(12-25)6-9-27-18/h2-6,15,18-19,21-22H,7-13H2,1H3/q+1',
                      'SIRIUS': 'InChI=1S/C43H52N4O5/c1-7-24-15-23-20-43(42(49)52-6)39-27(13-14-47(21-23)40(24)43)29-17-30(36(50-4)19-34(29)45-39)31-16-28-25(8-2)22-46(3)35(37(28)41(48)51-5)18-32-26-11-9-10-12-33(26)44-38(31)32/h8-12,17,19,23-24,28,31,35,37,40,44-45H,7,13-16,18,20-22H2,1-6H3'}

        tanimotos = Tanimotos(dict_inchi)
        data = K(dict_index, tanimotos)
        g = data.eta_g(data.similarities['GNPS'])
        s = data.eta_s(data.similarities['SIRIUS'])
        i = data.eta_i(data.similarities['ISDB'])
        k = data.k()

        self.assertEqual(g, -1.99307)
        self.assertEqual(s, -0.15122185067318528)
        self.assertEqual(i, 1.2425597051190609)
        self.assertEqual(k, -0.9017321455541243)

    def test_k_max(self):
        dict_index = {'GNPS': 1,
                      'ISDB': 1, 'SIRIUS': 1}
        dict_inchi = {'GNPS':   'InChI=1S/C43H54N4O6/c1-7-23-17-25-20-43(22-53-52-6)39-28(15-16-47(40(23)43)41(25)48)27-13-14-34(50-4)36(38(27)45-39)31-18-29-24(8-2)21-46(3)33(35(29)42(49)51-5)19-30-26-11-9-10-12-32(26)44-37(30)31/h9-14,23-25,29,31,33,35,40-41,44-45,48H,7-8,15-21H2,1-6H3/t23-,24+,25-,29-,31-,33-,35?,40-,41?,43?/m0/s1',
                      'ISDB':   'InChI=1S/C43H54N4O6/c1-7-23-17-25-20-43(22-53-52-6)39-28(15-16-47(40(23)43)41(25)48)27-13-14-34(50-4)36(38(27)45-39)31-18-29-24(8-2)21-46(3)33(35(29)42(49)51-5)19-30-26-11-9-10-12-32(26)44-37(30)31/h9-14,23-25,29,31,33,35,40-41,44-45,48H,7-8,15-21H2,1-6H3/t23-,24+,25-,29-,31-,33-,35?,40-,41?,43?/m0/s1',
                      'SIRIUS': 'InChI=1S/C43H54N4O6/c1-7-23-17-25-20-43(22-53-52-6)39-28(15-16-47(40(23)43)41(25)48)27-13-14-34(50-4)36(38(27)45-39)31-18-29-24(8-2)21-46(3)33(35(29)42(49)51-5)19-30-26-11-9-10-12-32(26)44-37(30)31/h9-14,23-25,29,31,33,35,40-41,44-45,48H,7-8,15-21H2,1-6H3/t23-,24+,25-,29-,31-,33-,35?,40-,41?,43?/m0/s1'}

        tanimotos = Tanimotos(dict_inchi)

        data = K(dict_index, tanimotos)

        g = data.eta_g(data.similarities['GNPS'])
        s = data.eta_s(data.similarities['SIRIUS'])
        i = data.eta_i(data.similarities['ISDB'])
        k = data.k()

        self.assertEqual(g, 50.0)
        self.assertEqual(s, 12.5)
        self.assertEqual(i, 6.25)
        self.assertEqual(k, 218.21428571428572)

    def test_k_with_two_eqaul_tools(self):
        dict_index = {'GNPS': 1.0,
                      'ISDB': 0.1, 'SIRIUS': 0.1}
        dict_inchi = {'GNPS':   'InChI=1S/C43H54N4O6/c1-7-23-17-25-20-43(22-53-52-6)39-28(15-16-47(40(23)43)41(25)48)27-13-14-34(50-4)36(38(27)45-39)31-18-29-24(8-2)21-46(3)33(35(29)42(49)51-5)19-30-26-11-9-10-12-32(26)44-37(30)31/h9-14,23-25,29,31,33,35,40-41,44-45,48H,7-8,15-21H2,1-6H3/t23-,24+,25-,29-,31-,33-,35?,40-,41?,43?/m0/s1',
                      'ISDB':   'InChI=1S/C23H27N2O3/c1-25-8-7-23-16-4-2-3-5-17(16)24-20(26)11-18-21(22(23)24)15(10-19(23)28-13-25)14(12-25)6-9-27-18/h2-6,15,18-19,21-22H,7-13H2,1H3/q+1',
                      'SIRIUS': 'InChI=1S/C23H27N2O3/c1-25-8-7-23-16-4-2-3-5-17(16)24-20(26)11-18-21(22(23)24)15(10-19(23)28-13-25)14(12-25)6-9-27-18/h2-6,15,18-19,21-22H,7-13H2,1H3/q+1'}

        tanimotos = Tanimotos(dict_inchi)
        data = K(dict_index, tanimotos)
        data.tanimotos.compute_tanimoto()

        t = data.tanimotos.tsi
        g = data.eta_g(data.similarities['GNPS'])
        s = data.eta_s(data.similarities['SIRIUS'])
        i = data.eta_i(data.similarities['ISDB'])
        k = data.k()

        self.assertEqual(g, 50.0)
        self.assertEqual(s, -0.75)
        self.assertEqual(i, -0.375)
        self.assertEqual(k, 50.214285714285715)
        self.assertEqual(t, 1.0)

    def test_k_empty_case(self):
        dict_index = {'GNPS': 0, 'ISDB': 0, 'SIRIUS': 0}
        dict_inchi = {'GNPS': '*', 'ISDB': '*', 'SIRIUS': '*'}

        tanimotos = Tanimotos(dict_inchi)

        data = K(dict_index, tanimotos)

        g = data.eta_g(data.similarities['GNPS'])
        s = data.eta_s(data.similarities['SIRIUS'])
        i = data.eta_i(data.similarities['ISDB'])
        k = data.k()

        self.assertEqual(g, -4.0)
        self.assertEqual(s, -1.0)
        self.assertEqual(i, -0.5)
        self.assertEqual(k, -5.5)

    def test_k_explicites(self):
        dict_index = {'GNPS': 0.4, 'ISDB': 0.4, 'SIRIUS': 0.4}
        dict_inchi = {'GNPS': '*', 'ISDB': '*', 'SIRIUS': '*'}

        tanimotos = Tanimotos(dict_inchi)
        data = K(dict_index, tanimotos)
        data.tanimotos.tgs=0.7
        data.tanimotos.tgi=0.7
        data.tanimotos.tsi=0.7
        
        g = data.eta_g(data.similarities['GNPS'])
        s = data.eta_s(data.similarities['SIRIUS'])
        i = data.eta_i(data.similarities['ISDB'])
        gs=data.phi_gs(data.tanimotos.tgs,data.similarities['GNPS'],data.similarities['SIRIUS'])
        gi=data.phi_gi(data.tanimotos.tgi,data.similarities['GNPS'],data.similarities['ISDB'])
        si=data.phi_si(data.tanimotos.tsi,data.similarities['SIRIUS'],data.similarities['ISDB'])
        
        k = data.k()

        self.assertEqual(g, 0)
        self.assertEqual(s, 0)
        self.assertEqual(i, 0)
        self.assertEqual(gs, 1.4285714285714282)
        self.assertEqual(gi, 0.857142857142857)
        self.assertEqual(si, 2.1428571428571423)
        self.assertEqual(k, 4.428571428571427)
        
        self.assertEqual(k, g+s+i+gs+gi+si)
if __name__ == '__main__':
    unittest.main()

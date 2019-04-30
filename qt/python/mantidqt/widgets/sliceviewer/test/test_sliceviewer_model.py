# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantid workbench.
#
#
from __future__ import (absolute_import, division, print_function)

from mantid.simpleapi import CreateMDHistoWorkspace
from mantidqt.widgets.sliceviewer.model import SliceViewerModel
from numpy.testing import assert_equal
import numpy as np
import unittest


class SliceViewerModelTest(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.ws_MD_3D = CreateMDHistoWorkspace(Dimensionality=3,
                                               Extents='-3,3,-10,10,-1,1',
                                               SignalInput=range(100),
                                               ErrorInput=range(100),
                                               NumberOfBins='5,5,4',
                                               Names='Dim1,Dim2,Dim3',
                                               Units='MomentumTransfer,EnergyTransfer,Angstrom',
                                               OutputWorkspace='ws_MD_2d')

    def test_model_MDH(self):

        model = SliceViewerModel(self.ws_MD_3D)

        self.assertEqual(model.get_ws(), self.ws_MD_3D)

        assert_equal(model.get_data((None, 2, 2)), range(90,95))
        assert_equal(model.get_data((1, 2, None)), range(18,118,25))
        assert_equal(model.get_data((None, None, 0)), np.reshape(range(50,75), (5,5)).T)

        dim_info = model.get_dim_info(0)
        self.assertEqual(dim_info['minimum'], -3)
        self.assertEqual(dim_info['maximum'], 3)
        self.assertEqual(dim_info['number_of_bins'], 5)
        self.assertAlmostEqual(dim_info['width'], 1.2)
        self.assertEqual(dim_info['name'], 'Dim1')
        self.assertEqual(dim_info['units'], 'MomentumTransfer')

        dim_infos = model.get_dimensions_info()
        self.assertEqual(len(dim_infos), 3)

        dim_info = dim_infos[2]
        self.assertEqual(dim_info['minimum'], -1)
        self.assertEqual(dim_info['maximum'], 1)
        self.assertEqual(dim_info['number_of_bins'], 4)
        self.assertAlmostEqual(dim_info['width'], 0.5)
        self.assertEqual(dim_info['name'], 'Dim3')
        self.assertEqual(dim_info['units'], 'Angstrom')


if __name__ == '__main__':
    unittest.main()

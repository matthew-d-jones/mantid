from __future__ import (absolute_import, division, print_function)
import unittest
import sys
from sans.test_helper.mock_objects import create_mock_diagnostics_tab
from sans.gui_logic.models.beam_centre_model import BeamCentreModel
from sans.common.enums import FindDirectionEnum,SANSInstrument
from sans.gui_logic.presenter.beam_centre_presenter import BeamCentrePresenter, find_beam_centre
from sans.test_helper.mock_objects import (create_run_tab_presenter_mock)
from sans.common.enums import IntegralEnum
if sys.version_info.major == 3:
    from unittest import mock
else:
    import mock

class DiagnosticsPagePresenterTest(unittest.TestCase):
    def setUp(self):
        self.parent_presenter = create_run_tab_presenter_mock(use_fake_state = False)
        self.view = create_mock_diagnostics_tab()
        self.WorkHandler = mock.MagicMock()
        self.run_integral = mock.MagicMock()
        self.presenter = DiagnosticsPagePresenter(self.parent_presenter, self.WorkHandler, self.run_integral)

    def test_that_calls_run_integral_with_correct_parameters(self):
        self.on_det1_horizontal_clicked()

        self.assertEqual(self.presenter._work_handler.process, 1)
        self.assertEqual(self.presenter._work_handler.process.call_args[0][1], self.run_integral)
        self.assertEqual(self.presenter._work_handler.process.call_args[0][2], '')
        self.assertEqual(self.presenter._work_handler.process.call_args[0][3], '')
        self.assertEqual(self.presenter._work_handler.process.call_args[0][4], '')
        self.assertEqual(self.presenter._work_handler.process.call_args[0][5], '')
        self.assertEqual(self.presenter._work_handler.process.call_args[0][6], IntegralEnum.Det1Horizontal)



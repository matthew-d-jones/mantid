// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2012 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

// Mantid Coding standars <http://www.mantidproject.org/Coding_Standards>
// Mantid Headers from the same project
#include "ui_DPDFBackgroundRemover.h"
// Mantid headers from other projects
#include "DllConfig.h"
#include "MantidQtWidgets/Common/UserSubWindow.h"
// 3rd party library headers
#include <QDialog>
#include <boost/shared_ptr.hpp>

// Class forward declarations
namespace MantidQt {
namespace CustomInterfaces {
namespace DynamicPDF {
class InputDataControl;
class SliceSelector;
class DisplayControl;
class FitControl;
class FourierTrasform;
} // namespace DynamicPDF
} // namespace CustomInterfaces
} // namespace MantidQt

namespace MantidQt {
namespace CustomInterfaces {
namespace DynamicPDF {
/** An interface to remove the multiphonon background from a
   S(Q,E) structure factor.

  @date 2016-03-11
*/
class MANTIDQT_DYNAMICPDF_DLL BackgroundRemover
    : public MantidQt::API::UserSubWindow {
  Q_OBJECT

public:
  /// The name of the interface as registered into the factory
  static std::string name() { return "Dynamic PDF Background Remover"; }
  // This interface's categories.
  static QString categoryInfo() { return "DynamicPDF"; }
  BackgroundRemover(QWidget *parent = nullptr);
  ~BackgroundRemover() override;

private slots:
  void showHelp();
  void summonSliceSelector();

private:
  void initLayout() override;
  /// The form generated by Qt Designer
  Ui::BackgroundRemover m_uiForm;
  /// GUI to load slices and select slice for fitting
  std::unique_ptr<SliceSelector> m_sliceSelector;
  /// Manage data of the loaded slice
  std::unique_ptr<InputDataControl> m_inputDataControl;
  /// Format the curves to be displayed
  std::unique_ptr<DisplayControl> m_displayControl;
  /// Managing the model and browser
  FitControl *m_fitControl;
  FourierTransform *m_fourierTransform;
};
} // namespace DynamicPDF
} // namespace CustomInterfaces
} // namespace MantidQt

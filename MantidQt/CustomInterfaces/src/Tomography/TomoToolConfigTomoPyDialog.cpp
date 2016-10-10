#include "MantidQtCustomInterfaces/Tomography/TomoToolConfigTomoPyDialog.h"
#include "MantidQtCustomInterfaces/Tomography/ToolConfigTomoPy.h"
#include "MantidQtCustomInterfaces/Tomography/TomoReconToolsUserSettings.h"

namespace MantidQt {
namespace CustomInterfaces {

void TomoToolConfigTomoPyDialog::setupDialogUi() {
  m_tomoPyUi.setupUi(this);
  m_tomoPyUi.comboBox_method->clear();

  const auto methods = ToolConfigTomoPy::methods();
  for (size_t i = 0; i < methods.size(); i++) {
    m_tomoPyUi.comboBox_method->addItem(
        QString::fromStdString(methods[i].second));
  }
}

/** Calls the execute of the QDialog
*/
int TomoToolConfigTomoPyDialog::executeQt() { return this->exec(); }

void TomoToolConfigTomoPyDialog::setupToolConfig() {
  // move to member/global variable and use more space OR keep here
  const auto methods = ToolConfigTomoPy::methods();

  int mi = m_tomoPyUi.comboBox_method->currentIndex();
  // TODO: for the output path, probably better to take the sample path,
  // then up one level
  m_toolSettings.tomoPy = ToolConfigTomoPy(
      m_runPath, m_pathOut + m_localOutNameAppendix, m_paths.pathDarks(),
      m_paths.pathOpenBeam(), m_paths.pathSamples());
  // maybe comboBox_method->currentText
  m_toolMethod = methods[mi].first;
}

} // Custominterfaces
} // MantidQt

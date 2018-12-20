// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef MANTIDQTCUSTOMDIALOGS_STARTLIVEDATADIALOG_H_
#define MANTIDQTCUSTOMDIALOGS_STARTLIVEDATADIALOG_H_

//----------------------
// Includes
//----------------------
#include "MantidAPI/Algorithm.h"
#include "MantidAPI/IAlgorithm.h"
#include "MantidQtWidgets/Common/AlgorithmDialog.h"
#include "MantidQtWidgets/Common/AlgorithmInputHistory.h"
#include "MantidQtWidgets/Common/WidgetScrollbarDecorator.h"
#include "ui_StartLiveDataDialog.h"

namespace MantidQt {
namespace CustomDialogs {
class StartLiveDataDialog : public MantidQt::API::AlgorithmDialog {
  Q_OBJECT

public:
  /// Default Constructor
  StartLiveDataDialog(QWidget *parent = nullptr);
  ~StartLiveDataDialog() override;

public slots:
  void radioProcessClicked();
  void radioPostProcessClicked();
  void changeProcessingAlgorithm();
  void changePostProcessingAlgorithm();
  void radioTimeClicked();
  void chkPreserveEventsToggled();

private slots:
  void setDefaultAccumulationMethod(const QString &);
  void updateUiElements(const QString &);
  void accept() override;
  void initListenerPropLayout(const QString &);
  void updateConnectionChoices(const QString &inst_name);
  void updateConnectionDetails(const QString &connection);

private:
  /// Initialize the layout
  void initLayout() override;

  /// Parse the input from the dialog when it has been accepted
  void parseInput() override;

  Mantid::API::Algorithm_sptr
  changeAlgorithm(MantidQt::MantidWidgets::AlgorithmSelectorWidget *selector,
                  MantidQt::API::AlgorithmPropertiesWidget *propWidget);

private:
  /// The form generated by Qt Designer
  Ui::StartLiveDataDialog ui;

  /// This dialog needs scrollbars on small screens
  API::WidgetScrollbarDecorator m_scrollbars;

  bool m_useProcessAlgo;
  bool m_useProcessScript;
  bool m_usePostProcessAlgo;
  bool m_usePostProcessScript;

  /// The algorithm for processing chunks
  Mantid::API::Algorithm_sptr m_processingAlg;

  /// The algorithm for processing the accumulated workspace
  Mantid::API::Algorithm_sptr m_postProcessingAlg;

  /// Constant used for custom listener connection setups
  static const QString CUSTOM_CONNECTION;
};
} // namespace CustomDialogs
} // namespace MantidQt

#endif // MANTIDQTCUSTOMDIALOGS_STARTLIVEDATADIALOG_H_

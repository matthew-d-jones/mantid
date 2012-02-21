#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/IAlgorithm.h"
#include "MantidAPI/IMDHistoWorkspace.h"
#include "MantidGeometry/MDGeometry/IMDDimension.h"
#include "MantidGeometry/MDGeometry/MDTypes.h"
#include "MantidKernel/VMD.h"
#include "MantidQtSliceViewer/LineViewer.h"
#include <QIntValidator>
#include <qwt_plot_curve.h>
#include <qwt_plot.h>
#include "MantidQtAPI/MantidQwtIMDWorkspaceData.h"
#include "MantidAPI/NullCoordTransform.h"
#include "MantidQtSliceViewer/LinePlotOptions.h"

using namespace Mantid;
using namespace Mantid::API;
using namespace Mantid::Kernel;
using Mantid::Geometry::IMDDimension_const_sptr;

namespace MantidQt
{
namespace SliceViewer
{


LineViewer::LineViewer(QWidget *parent)
 : QWidget(parent),
   m_planeWidth(0),
   m_numBins(100),
   m_allDimsFree(false), m_freeDimX(0), m_freeDimY(1),
   m_fixedBinWidthMode(false), m_fixedBinWidth(0.1), m_binWidth(0.1)
{
	ui.setupUi(this);

	// Other setup
  ui.textBinWidth->setValidator(new QDoubleValidator(ui.textBinWidth));

	// --------- Create the plot -----------------
  m_plotLayout = new QVBoxLayout(ui.frmPlot);
  m_plot = new QwtPlot();
  m_plot->autoRefresh();
  m_plot->setBackgroundColor(QColor(255,255,255)); // White background
  m_plotLayout->addWidget(m_plot, 1);

  // Make the 2 curves
  m_previewCurve = new QwtPlotCurve("Preview");
  m_fullCurve = new QwtPlotCurve("Integrated");
  m_previewCurve->attach(m_plot);
  m_fullCurve->attach(m_plot);
  m_previewCurve->setVisible(false);
  m_fullCurve->setVisible(false);

  // The plotOptions
  m_lineOptions = new LinePlotOptions(this);
  m_plotLayout->addWidget(m_lineOptions, 0);


  // Make the splitter use the minimum size for the controls and not stretch out
  ui.splitter->setStretchFactor(0, 0);
  ui.splitter->setStretchFactor(1, 1);

  //----------- Connect signals -------------
  QObject::connect(ui.btnApply, SIGNAL(clicked()), this, SLOT(apply()));
  QObject::connect(ui.chkAdaptiveBins, SIGNAL(  stateChanged(int)), this, SLOT(adaptiveBinsChanged()));
  QObject::connect(ui.spinNumBins, SIGNAL(valueChanged(int)), this, SLOT(numBinsChanged()));
  QObject::connect(ui.textPlaneWidth, SIGNAL(textEdited(QString)), this, SLOT(thicknessTextEdited()));
  QObject::connect(ui.radNumBins, SIGNAL(toggled(bool)), this, SLOT(on_radNumBins_toggled()));
  QObject::connect(ui.textBinWidth, SIGNAL(editingFinished()), this, SLOT(textBinWidth_changed()));

  QObject::connect(m_lineOptions, SIGNAL(changedPlotAxis()), this, SLOT(refreshPlot()));
  QObject::connect(m_lineOptions, SIGNAL(changedNormalization()), this, SLOT(refreshPlot()));
}

LineViewer::~LineViewer()
{

}

//-----------------------------------------------------------------------------------------------
/** With the workspace set, create the dimension text boxes */
void LineViewer::createDimensionWidgets()
{
  // Create all necessary widgets
  if (m_startText.size() < int(m_ws->getNumDims()))
  {
    for (size_t d=m_startText.size(); d<m_ws->getNumDims(); d++)
    {
      QLabel * dimLabel = new QLabel(this);
      dimLabel->setAlignment(Qt::AlignHCenter);
      ui.gridLayout->addWidget(dimLabel, 0, int(d)+1);
      m_dimensionLabel.push_back(dimLabel);

      QLineEdit * startText = new QLineEdit(this);
      QLineEdit * endText = new QLineEdit(this);
      QLineEdit * thicknessText = new QLineEdit(this);
      startText->setMaximumWidth(100);
      endText->setMaximumWidth(100);
      thicknessText->setMaximumWidth(100);
      startText->setToolTip("Start point of the line in this dimension");
      endText->setToolTip("End point of the line in this dimension");
      thicknessText->setToolTip("Integration thickness (above and below plane) in this dimension");
      startText->setValidator(new QDoubleValidator(startText));
      endText->setValidator(new QDoubleValidator(endText));
      thicknessText->setValidator(new QDoubleValidator(thicknessText));
      ui.gridLayout->addWidget(startText, 1, int(d)+1);
      ui.gridLayout->addWidget(endText, 2, int(d)+1);
      ui.gridLayout->addWidget(thicknessText, 3, int(d)+1);
      m_startText.push_back(startText);
      m_endText.push_back(endText);
      m_thicknessText.push_back(thicknessText);
      // Signals that don't change
      QObject::connect(thicknessText, SIGNAL(textEdited(QString)), this, SLOT(thicknessTextEdited()));
    }
  }

  // ------ Update the widgets -------------------------
  for (int d=0; d<int(m_ws->getNumDims()); d++)
  {
    m_dimensionLabel[d]->setText( QString::fromStdString(m_ws->getDimension( size_t(d))->getName() ) );
  }
}


//-----------------------------------------------------------------------------------------------
/** Disable any controls relating to dimensions that are not "free"
 * e.g. if you are in the X-Y plane, the Z position cannot be changed.
 * Also updates the radio buttons for the choice of X axis.
 */
void LineViewer::updateFreeDimensions()
{
  for (int d=0; d<int(m_ws->getNumDims()); d++)
  {
    // Can always change the start value
    m_startText[d]->setEnabled(true);

    // This dimension is free to move if b == true
    bool b = (m_allDimsFree || d == m_freeDimX || d == m_freeDimY);
    m_endText[d]->setEnabled(b);
    // If all dims are free, width makes little sense. Only allow one (circular) width
    if (m_allDimsFree)
      m_thicknessText[d]->setVisible(d != 0);
    else
      m_thicknessText[d]->setVisible(!b);

    // --- Adjust the signals ---
    m_startText[d]->disconnect();
    m_endText[d]->disconnect();

    if (d == m_freeDimX || d == m_freeDimY)
    {
      // Free dimension - update the preview
      QObject::connect(m_startText[d], SIGNAL(textEdited(QString)), this, SLOT(startEndTextEdited()));
      QObject::connect(m_endText[d], SIGNAL(textEdited(QString)), this, SLOT(startEndTextEdited()));
    }
    else
    {
      // Non-Free dimension - link start to end
      QObject::connect(m_startText[d], SIGNAL(textEdited(QString)), this, SLOT(startLinkedToEndText()));
    }
  }
  if (!m_allDimsFree)
  {
    std::string s = "(in " + m_ws->getDimension(m_freeDimX)->getName() + "-" +  m_ws->getDimension(m_freeDimY)->getName()
        + " plane)";
    ui.lblPlaneWidth->setText(QString::fromStdString(s));
  }

}

//-----------------------------------------------------------------------------------------------
/** Show the start/end/width points in the GUI */
void LineViewer::updateStartEnd()
{
  for (int d=0; d<int(m_ws->getNumDims()); d++)
  {
    m_startText[d]->setText(QString::number(m_start[d]));
    m_endText[d]->setText(QString::number(m_end[d]));
    m_thicknessText[d]->setText(QString::number(m_thickness[d]));
  }
  ui.textPlaneWidth->setText(QString::number(m_planeWidth));

  // Now show the width
  this->updateBinWidth();
}

//-----------------------------------------------------------------------------------------------
/** Calculate the number of bins (fixed-bin-width mode)
 * or show the bin width (fixed-#-bins mode) */
void LineViewer::updateBinWidth()
{
  // If partially initialized, vectors might be wrong
  if (m_start.getNumDims() != m_end.getNumDims())
    return;
  double length = (m_start - m_end).norm();
  if (m_fixedBinWidthMode)
  {
    // Fixed bin width. Find the number of bins.
    m_numBins = size_t(length / m_fixedBinWidth + 0.5);
    if (m_numBins < 1) m_numBins = 1;
    // Show the # of bins
    ui.spinNumBins->blockSignals(true);
    ui.spinNumBins->setValue(int(m_numBins));
    ui.spinNumBins->blockSignals(false);
    // Show the fixed bin width
    m_binWidth = length / double(m_numBins);
    ui.textBinWidth->setText(QString::number(m_fixedBinWidth));
  }
  else
  {
    // Fixed number of bins mode
    m_binWidth = length / double(m_numBins);
    ui.textBinWidth->setText(QString::number(m_binWidth));
  }
}

//-----------------------------------------------------------------------------------------------
/** Read all the text boxes and interpret their values.
 * Does not refresh.
 */
void LineViewer::readTextboxes()
{
  VMD start = m_start;
  VMD end = m_start;
  VMD width = m_thickness;
  bool allOk = true;
  bool ok;
  for (int d=0; d<int(m_ws->getNumDims()); d++)
  {
    start[d] = VMD_t(m_startText[d]->text().toDouble(&ok));
    allOk = allOk && ok;

    end[d] = VMD_t(m_endText[d]->text().toDouble(&ok));
    allOk = allOk && ok;

    width[d] = VMD_t(m_thicknessText[d]->text().toDouble(&ok));
    allOk = allOk && ok;
  }
  // Now the planar width
  double tempPlaneWidth = ui.textPlaneWidth->text().toDouble(&ok);
  allOk = allOk && ok;

  // Only continue if all values typed were valid numbers.
  if (!allOk) return;
  m_start = start;
  m_end = end;
  m_thickness = width;
  m_planeWidth = tempPlaneWidth;
}

//-----------------------------------------------------------------------------------------------
/** Perform the 1D integration using the current parameters.
 * @throw std::runtime_error if an error occurs.
 * */
void LineViewer::apply()
{
  if (m_allDimsFree)
    throw std::runtime_error("Not currently supported with all dimensions free!");

  // BinMD fails on MDHisto.
  IMDHistoWorkspace_sptr mdhws = boost::dynamic_pointer_cast<IMDHistoWorkspace>(m_ws);

  std::string outWsName = m_ws->getName() + "_line" ;
  bool adaptive = ui.chkAdaptiveBins->isChecked();

  // (half-width in the plane)
  double planeWidth = this->getPlanarWidth();
  // Length of the line
  double length = (m_end - m_start).norm();
  double dx = m_end[m_freeDimX] - m_start[m_freeDimX];
  double dy = m_end[m_freeDimY] - m_start[m_freeDimY];
  // Angle of the line
  double angle = atan2(dy, dx);
  double perpAngle = angle + M_PI / 2.0;

  // Build the basis vectors using the angles
  VMD basisX = m_start * 0;
  basisX[m_freeDimX] = VMD_t(cos(angle));
  basisX[m_freeDimY] = VMD_t(sin(angle));
  VMD basisY = m_start * 0;
  basisY[m_freeDimX] = VMD_t(cos(perpAngle));
  basisY[m_freeDimY] = VMD_t(sin(perpAngle));

  // This is the origin = Translation parameter
  VMD origin = m_start;

  IAlgorithm * alg = NULL;
  size_t numBins = m_numBins;
  if (adaptive)
  {
    alg = FrameworkManager::Instance().createAlgorithm("SliceMD");
    // "SplitInto" parameter
    numBins = 2;
  }
  else
    alg = FrameworkManager::Instance().createAlgorithm("BinMD");

  alg->setProperty("InputWorkspace", m_ws);
  alg->setPropertyValue("OutputWorkspace", outWsName);
  alg->setProperty("AxisAligned", false);

  std::vector<int> OutputBins;
  std::vector<double> OutputExtents;

  // The X basis vector
  alg->setPropertyValue("BasisVector0", "X,units," + basisX.toString(",") );
  OutputExtents.push_back(0);
  OutputExtents.push_back(length);
  OutputBins.push_back(int(numBins));

  // The Y basis vector, with one bin
  alg->setPropertyValue("BasisVector1", "Y,units," + basisY.toString(","));
  OutputExtents.push_back(-planeWidth);
  OutputExtents.push_back(+planeWidth);
  OutputBins.push_back(1);

  // Now each remaining dimension
  std::string dimChars = "012345"; // SlicingAlgorithm::getDimensionChars();
  size_t propNum = 2;
  for (int d=0; d<int(m_ws->getNumDims()); d++)
  {
    if ((d != m_freeDimX) && (d != m_freeDimY))
    {
      // Letter of the dimension
      std::string dim(" "); dim[0] = dimChars[propNum];
      // Simple basis vector going only in this direction
      VMD basis = m_start * 0;
      basis[d] = 1.0;
      // Set the basis vector with the width *2 and 1 bin
      alg->setPropertyValue("BasisVector" + dim, dim +",units," + basis.toString(",") );
      OutputExtents.push_back(-m_thickness[d]);
      OutputExtents.push_back(+m_thickness[d]);
      OutputBins.push_back(1);

      propNum++;
      if (propNum > dimChars.size())
        throw std::runtime_error("LineViewer::apply(): too many dimensions!");

    }
  }

  alg->setPropertyValue("Translation", origin.toString(",") );
  alg->setProperty("OutputBins", OutputBins );
  alg->setProperty("OutputExtents", OutputExtents );
  if (!adaptive)
  {
    alg->setProperty("IterateEvents", true);
  }
  alg->execute();

  if (alg->isExecuted())
  {
    //m_sliceWS = alg->getProperty("OutputWorkspace");
    m_sliceWS = boost::dynamic_pointer_cast<IMDWorkspace>(AnalysisDataService::Instance().retrieve(outWsName));
    this->showFull();
  }
  else
  {
    // Unspecified error in algorithm
    this->showPreview();
    m_plot->setTitle("Error integrating workspace - see log.");
  }


}

// ==============================================================================================
// ================================== SLOTS =====================================================
// ==============================================================================================

//-------------------------------------------------------------------------------------------------
/** Slot called when the start text of a non-free dimensions is changed.
 * Changes the end text correspondingly
 */
void LineViewer::startLinkedToEndText()
{
  for (int d=0; d<int(m_ws->getNumDims()); d++)
  {
    if (d != m_freeDimX && d != m_freeDimY)
    {
      // Copy the start text to the end text
      m_endText[d]->setText( m_startText[d]->text() );
    }
  }
  // Call the slot to update the preview
  startEndTextEdited();
}


//-------------------------------------------------------------------------------------------------
/** Slot called when any of the start/end text boxes are edited
 * in GUI. Only changes the values if they are all valid.
 */
void LineViewer::startEndTextEdited()
{
  this->readTextboxes();
  this->showPreview();
  // Send the signal that the positions changed
  emit changedStartOrEnd(m_start, m_end);
}

/** Slot called when the width text box is edited */
void LineViewer::thicknessTextEdited()
{
  this->readTextboxes();
  //TODO: Don't always auto-apply
  this->apply();
  // Send the signal that the width changed
  emit changedPlanarWidth(this->getPlanarWidth());
}

/** Slot called when the number of bins changes */
void LineViewer::numBinsChanged()
{
  m_numBins = ui.spinNumBins->value();
  // Show the bin width
  this->updateBinWidth();
  //TODO: Don't always auto-apply
  this->apply();
}

/** Slot called when checking the adaptive box */
void LineViewer::adaptiveBinsChanged()
{
  //TODO: Don't always auto-apply
  this->apply();
}

/** Slot called when the num bins/bin width radio choice changes */
void LineViewer::on_radNumBins_toggled()
{
  setFixedBinWidthMode(ui.radNumBins->isChecked(), m_fixedBinWidth);
}


/** Slot called when the desired fixed bin width text box
 * is edited and the user pressed Return or lost focus.
 */
void LineViewer::textBinWidth_changed()
{
  if (m_fixedBinWidthMode)
  {
    bool ok;
    double width = ui.textBinWidth->text().toDouble(&ok);
    if (ok && width > 0)
    {
      // Change the desired bin size and update necessary GUI
      this->setFixedBinWidthMode(m_fixedBinWidthMode, width);
    }
    else
    {
      // Bad number! Reset to the old value
      this->updateBinWidth();
    }
  }
}

// ==============================================================================================
// ================================== External Getters ==========================================
// ==============================================================================================
/** @return the width in the plane, or the width in dimension 0 if not restricted to a plane */
double LineViewer::getPlanarWidth() const
{
  return m_planeWidth;
}

/// @return the full width vector in each dimensions. The values in the X-Y dimensions should be ignored
Mantid::Kernel::VMD LineViewer::getWidth() const
{
  return m_thickness;
}

/** For fixed-bin-width mode, get the desired fixed bin width.
 * @return the desired fixed bin width
 */
double LineViewer::getFixedBinWidth() const
{
  return m_fixedBinWidth;
}

/** Is the LineViewer in fixed-bin-width mode?
 * @return True if in fixed bin width mode.
 */
bool LineViewer::getFixedBinWidthMode() const
{
  return m_fixedBinWidthMode;
}

// ==============================================================================================
// ================================== External Setters ==========================================
// ==============================================================================================
//-----------------------------------------------------------------------------------------------
/** Set the workspace being sliced
 *
 * @param ws :: IMDWorkspace */
void LineViewer::setWorkspace(Mantid::API::IMDWorkspace_sptr ws)
{
  m_ws = ws;
  m_thickness = VMD(ws->getNumDims());
  createDimensionWidgets();
  // Update the dimensions shown in the original workspace
  m_lineOptions->setOriginalWorkspace(m_ws);
}


/** Set the start point of the line to integrate
 * @param start :: vector for the start point */
void LineViewer::setStart(Mantid::Kernel::VMD start)
{
  if (m_ws && start.getNumDims() != m_ws->getNumDims())
    throw std::runtime_error("LineViewer::setStart(): Invalid number of dimensions in the start vector.");
  m_start = start;
  updateStartEnd();
}

/** Set the end point of the line to integrate
 * @param end :: vector for the end point */
void LineViewer::setEnd(Mantid::Kernel::VMD end)
{
  if (m_ws && end.getNumDims() != m_ws->getNumDims())
    throw std::runtime_error("LineViewer::setEnd(): Invalid number of dimensions in the end vector.");
  m_end = end;
  updateStartEnd();
}


/** Set the width of the line in each dimensions
 * @param width :: vector for the width in each dimension. X dimension stands in for the XY plane width */
void LineViewer::setThickness(Mantid::Kernel::VMD width)
{
  if (m_ws && width.getNumDims() != m_ws->getNumDims())
    throw std::runtime_error("LineViewer::setThickness(): Invalid number of dimensions in the width vector.");
  m_thickness = width;
  updateStartEnd();
}

/** Set the width of the line in the planar dimension only.
 * Other dimensions' widths will follow unless they were manually changed
 * @param width :: width in the plane. */
void LineViewer::setPlanarWidth(double width)
{
  if (m_allDimsFree)
  {
    for (size_t d=0; d<m_thickness.getNumDims(); d++)
      m_thickness[d] = VMD_t(width);
  }
  else
  {
    double oldPlanarWidth = this->getPlanarWidth();
    for (size_t d=0; d<m_thickness.getNumDims(); d++)
    {
      // Only modify the locked onese
      if (m_thickness[d] == oldPlanarWidth)
        m_thickness[d] = VMD_t(width);
    }
    // And always set the planar one
    m_planeWidth = width;
  }
  updateStartEnd();
  // Emit the signal that the width changed
  emit changedPlanarWidth(width);
}

/** Set the number of bins in the line.
 * @param numBins :: # of bins
 * @throw std::invalid_argument if numBins < 1
 * */
void LineViewer::setNumBins(int numBins)
{
  if (numBins < 1)
    throw std::invalid_argument("LineViewer::setNumBins(): must be > 0");
  m_numBins = size_t(numBins);
  //ui.spinNumBins->blockSignals(true);
  ui.spinNumBins->setValue( int(numBins) );
  //ui.spinNumBins->blockSignals(false);
}

/** Set the free dimensions - dimensions that are allowed to change
 *
 * @param all :: Flag that is true when all dimensions are allowed to change
 * @param dimX :: Index of the X dimension in the 2D slice
 * @param dimY :: Index of the Y dimension in the 2D slice
 */
void LineViewer::setFreeDimensions(bool all, int dimX, int dimY)
{
  int nd = int(m_ws->getNumDims());
  if (dimX < 0 || dimX >= nd)
    throw std::runtime_error("LineViewer::setFreeDimensions(): Free X dimension index is out of range.");
  if (dimY < 0 || dimY >= nd)
    throw std::runtime_error("LineViewer::setFreeDimensions(): Free Y dimension index is out of range.");
  m_allDimsFree = all;
  m_freeDimX = dimX;
  m_freeDimY = dimY;
  this->updateFreeDimensions();
}

/** Slot called to set the free dimensions (called from the SliceViewer widget)
 *
 * @param dimX :: index of the X-dimension of the plane
 * @param dimY :: index of the Y-dimension of the plane
 */
void LineViewer::setFreeDimensions(size_t dimX, size_t dimY)
{
  m_allDimsFree = false;
  m_freeDimX = int(dimX);
  m_freeDimY = int(dimY);
  this->updateFreeDimensions();
}

/** Sets the fixed bin width mode on or off.
 *
 * In fixed bin width mode, the width of each bin along the line length
 * is constant, and the number of bins is adjusted to as the line
 * gets longer.
 * If off, then you use a fixed number of bins, and the bin width is
 * then simply: width = length / number_of_bins.
 *
 * @param fixedWidth :: if True, then keep the bin width fixed.
 * @param binWidth :: for fixed bin width mode, this specified the desired
 *        bin width. Must be > 0. Ignored for non-fixed-bin-width mode.
 * @throw std::invalid_argument if binWidth <= 0
 */
void LineViewer::setFixedBinWidthMode(bool fixedWidth, double binWidth)
{
  if (binWidth <= 0)
    throw std::invalid_argument("LineViewer::setFixedBinWidthMode(): binWidth must be > 0");

  m_fixedBinWidthMode = fixedWidth;
  if (m_fixedBinWidthMode)
  {
    m_fixedBinWidth = binWidth;
    ui.textBinWidth->setReadOnly(false);
    ui.textBinWidth->setToolTip("Desired bin width (will adjust the number of bins).");
    ui.spinNumBins->setReadOnly(true);
    ui.spinNumBins->setToolTip("Current number of bins (calculated from the fixed bin width)");
  }
  else
  {
    ui.textBinWidth->setReadOnly(true);
    ui.textBinWidth->setToolTip("Current bin width, given the number of bins.");
    ui.spinNumBins->setReadOnly(false);
    ui.spinNumBins->setToolTip("Desired number of bins.");
  }

  // Show in GUI
  ui.radNumBins->blockSignals(true);
  ui.radNumBins->setChecked(!m_fixedBinWidthMode);
  ui.radBinWidth->setChecked(m_fixedBinWidthMode);
  ui.radNumBins->blockSignals(false);

  // Signal the LineOverlay to used fixed bin width
  emit changedFixedBinWidth(fixedWidth, binWidth);
  // Show the start/end, and update # of bins
  this->updateStartEnd();

  //TODO: Don't always auto-apply
  this->apply();
}

// ==============================================================================================
// ================================== Methods for Python ==========================================
// ==============================================================================================
/** Set the start point of the line to integrate
 *
 * @param x :: position of the start in the "X" dimension
 *        (as shown in the SliceViewer).
 * @param y :: position of the start in the "Y" dimension
 *        (as shown in the SliceViewer).
 */
void LineViewer::setStartXY(double x, double y)
{
  if (m_allDimsFree)
    throw std::runtime_error("LineViewer::setStartXY(): cannot use with all dimensions free.");
  m_start[m_freeDimX] = VMD_t(x);
  m_start[m_freeDimY] = VMD_t(y);
  updateStartEnd();
  // Send the signal that the positions changed
  emit changedStartOrEnd(m_start, m_end);
}

/** Set the start point of the line to integrate
 *
 * @param x :: position of the start in the "X" dimension
 *        (as shown in the SliceViewer).
 * @param y :: position of the start in the "Y" dimension
 *        (as shown in the SliceViewer).
 */
void LineViewer::setEndXY(double x, double y)
{
  if (m_allDimsFree)
    throw std::runtime_error("LineViewer::setEndXY(): cannot use with all dimensions free.");
  m_end[m_freeDimX] = VMD_t(x);
  m_end[m_freeDimY] = VMD_t(y);
  updateStartEnd();
  // Send the signal that the positions changed
  emit changedStartOrEnd(m_start, m_end);
}

//------------------------------------------------------------------------------
/** Return the start of the line to integrate,
 * in the X/Y coordinates as shown in the SliceViewer
 * @return [X,Y] coordinates
 */
QPointF LineViewer::getStartXY() const
{
  if (m_allDimsFree)
    throw std::runtime_error("LineViewer::getStartXY(): cannot use with all dimensions free.");
  return QPointF(m_start[m_freeDimX], m_start[m_freeDimY]);
}

/** Return the end of the line to integrate,
 * in the X/Y coordinates as shown in the SliceViewer
 * @return [X,Y] coordinates
 */
QPointF LineViewer::getEndXY() const
{
  if (m_allDimsFree)
    throw std::runtime_error("LineViewer::getEndXY(): cannot use with all dimensions free.");
  return QPointF(m_end[m_freeDimX], m_end[m_freeDimY]);
}

//------------------------------------------------------------------------------
/** Set the thickness to integrate to be the same in all dimensions
 *
 * This sets the planar width and all the other dimensions' thicknesses
 * to the same value.
 *
 * @param width :: width of integration, in the units of all dimensions
 */
void LineViewer::setThickness(double width)
{
  if (!m_ws) return;
  for (int i=0; i<int(m_ws->getNumDims()); i++)
    m_thickness[i] = VMD_t(width);
  this->setPlanarWidth(width);
}

/** Set the thickness to integrate in a particular dimension.
 *
 * Integration is performed perpendicular to the XY plane,
 * from -thickness below to +thickness above the center.
 *
 * Use setPlanarWidth() to set the width along the XY plane.
 *
 * @param dim :: index of the dimension to change
 * @param width :: width of integration, in the units of the dimension.
 * @throw std::invalid_argument if the index is invalid
 */
void LineViewer::setThickness(int dim, double width)
{
  if (!m_ws) return;
  if (dim >= int(m_ws->getNumDims()) || dim < 0)
    throw std::invalid_argument("There is no dimension # " + Strings::toString(dim) + " in the workspace.");
  m_thickness[dim] = VMD_t(width);
  updateStartEnd();
}

/** Set the thickness to integrate in a particular dimension.
 *
 * Integration is performed perpendicular to the XY plane,
 * from -thickness below to +thickness above the center.
 *
 * Use setPlanarWidth() to set the width along the XY plane.
 *
 * @param dim :: name of the dimension to change
 * @param width :: thickness of integration, in the units of the dimension.
 * @throw std::runtime_error if the name is not found in the workspace
 */
void LineViewer::setThickness(const QString & dim, double width)
{
  if (!m_ws) return;
  int index = int(m_ws->getDimensionIndexByName(dim.toStdString()));
  return this->setThickness(index, width);
}

/** Get the number of bins
 *
 * @return the number of bins in the line to integrate (int)
 */
int LineViewer::getNumBins() const
{
  return int(m_numBins);
}

/** Get the width of each bin
 *
 * @return the width of each bin (double)
 */
double LineViewer::getBinWidth() const
{
  return m_binWidth;
}

/** Choose which coordinates to use as the X axis to plot in the line view.
 *
 * @param choice :: PlotAxisChoice, either Auto, X, Y or Distance.
 */
void LineViewer::setPlotAxis(int choice)
{
  m_lineOptions->setPlotAxis(choice);
}

/** Return which coordinates to use as the X axis to plot in the line view.
 *
 * @return PlotAxisChoice, either Auto, X, Y or Distance.
 */
int LineViewer::getPlotAxis() const
{
  return m_lineOptions->getPlotAxis();
}


// ==============================================================================================
// ================================== Rendering =================================================
// ==============================================================================================


//-----------------------------------------------------------------------------
/** Calculate and show the preview (non-integrated) line,
 * using the current parameters. */
void LineViewer::showPreview()
{
  MantidQwtIMDWorkspaceData curveData(m_ws, false,
      m_start, m_end, m_lineOptions->getNormalization());
  curveData.setPreviewMode(true);
  curveData.setPlotAxisChoice(m_lineOptions->getPlotAxis());
  m_previewCurve->setData(curveData);

  if (m_fullCurve->isVisible())
  {
    m_fullCurve->setVisible(false);
    m_fullCurve->detach();
    m_previewCurve->attach(m_plot);
  }
  m_previewCurve->setVisible(true);
  m_plot->replot();
  m_plot->setTitle("Preview Plot");

  m_plot->setAxisTitle( QwtPlot::xBottom, QString::fromStdString( curveData.getXAxisLabel() ));;
  m_plot->setAxisTitle( QwtPlot::yLeft, QString::fromStdString( curveData.getYAxisLabel() ));;
}


//-----------------------------------------------------------------------------
/** Calculate and show the full (integrated) line, using the latest
 * integrated workspace. The apply() method must have been called
 * before calling this. */
void LineViewer::showFull()
{
  if (!m_sliceWS) return;
  MantidQwtIMDWorkspaceData curveData(m_sliceWS, false,
      VMD(), VMD(), m_lineOptions->getNormalization());
  curveData.setPreviewMode(false);
  curveData.setPlotAxisChoice(m_lineOptions->getPlotAxis());
  m_fullCurve->setData(curveData);

  if (m_previewCurve->isVisible())
  {
    m_previewCurve->setVisible(false);
    m_previewCurve->detach();
    m_fullCurve->attach(m_plot);
  }
  m_fullCurve->setVisible(true);
  m_plot->replot();
  m_plot->setTitle("Integrated Line Plot");
  m_plot->setAxisTitle( QwtPlot::xBottom, QString::fromStdString( curveData.getXAxisLabel() ));;
  m_plot->setAxisTitle( QwtPlot::yLeft, QString::fromStdString( curveData.getYAxisLabel() ));;
}

//-----------------------------------------------------------------------------
/** Slot called when the options of the plot display change (normalization
 * or plot axis.
 * Refreshes the preview or full plot, whichever is visible.
 */
void LineViewer::refreshPlot()
{
  if (m_previewCurve->isVisible())
    showPreview();
  else
    showFull();
}


} // namespace
}

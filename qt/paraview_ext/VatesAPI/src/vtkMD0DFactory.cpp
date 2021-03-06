// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidVatesAPI/vtkMD0DFactory.h"
#include "MantidAPI/IMDWorkspace.h"
#include "MantidKernel/Logger.h"
#include "MantidVatesAPI/ProgressAction.h"
#include "MantidVatesAPI/vtkNullUnstructuredGrid.h"

using namespace Mantid::API;

namespace {
Mantid::Kernel::Logger g_log("vtkMD0DFactory");
}

namespace Mantid {
namespace VATES {
/**
Constructor
*/
vtkMD0DFactory::vtkMD0DFactory() {}

/// Destructor
vtkMD0DFactory::~vtkMD0DFactory() {}

/**
Create the vtkStructuredGrid from the provided workspace
@param progressUpdating: Reporting object to pass progress information up the
stack.
@return fully constructed vtkDataSet.
*/
vtkSmartPointer<vtkDataSet> vtkMD0DFactory::create(ProgressAction &) const {
  g_log.warning() << "Factory " << this->getFactoryTypeName()
                  << " is being used. You are viewing data with less than "
                     "three dimensions in the VSI. \n";
  vtkNullUnstructuredGrid nullGrid;
  auto visualDataSet =
      vtkSmartPointer<vtkDataSet>::Take(nullGrid.createNullData());
  return visualDataSet;
}

/// Initalize with a target workspace.
void vtkMD0DFactory::initialize(
    const Mantid::API::Workspace_sptr & /*workspace*/) {}

/// Validate the workspace
void vtkMD0DFactory::validate() const {}
} // namespace VATES
} // namespace Mantid

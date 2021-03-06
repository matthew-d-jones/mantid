// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2015 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include <memory>

namespace Mantid {
namespace API {

/**
  This file provides forward declarations for Mantid::API::IMDEventWorkspace
*/

/// forward declare of Mantid::API::IMDEventWorkspace
class IMDEventWorkspace;
/// Shared pointer to Mantid::API::IMDEventWorkspace
using IMDEventWorkspace_sptr = std::shared_ptr<IMDEventWorkspace>;
/// Shared pointer to Mantid::API::IMDEventWorkspace (const version)
using IMDEventWorkspace_const_sptr = std::shared_ptr<const IMDEventWorkspace>;
/// unique pointer to Mantid::API::IMDEventWorkspace
using IMDEventWorkspace_uptr = std::unique_ptr<IMDEventWorkspace>;
/// unique pointer to Mantid::API::IMDEventWorkspace (const version)
using IMDEventWorkspace_const_uptr = std::unique_ptr<const IMDEventWorkspace>;

} // namespace API
} // namespace Mantid

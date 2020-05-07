// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2016 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/ISpectrum.h"
#include "MantidAlgorithms/DllConfig.h"
#include "MantidAlgorithms/SampleCorrections/MCInteractionVolume.h"
#include "MantidHistogramData/Histogram.h"
#include <tuple>

namespace Mantid {
namespace API {
class Sample;
}
namespace Kernel {
class PseudoRandomNumberGenerator;
class V3D;
class Logger;
} // namespace Kernel
namespace Algorithms {
class IBeamProfile;

/**
  Implements the algorithm for calculating absorption weighted path lengths
  using a Monte Carlo. A single instance has a fixed nominal source position,
  nominal sample position & sample + containers shapes.
*/
class MANTID_ALGORITHMS_DLL MCAbsorptionWeightedPathStrategy {
public:
  MCAbsorptionWeightedPathStrategy(
      const IBeamProfile &beamProfile,
                       const API::Sample &sample,
                       const size_t nevents,
                       const size_t maxScatterPtAttempts,
                       Kernel::Logger &logger);
  void calculate(Kernel::PseudoRandomNumberGenerator &rng,
                 const Kernel::V3D &finalPos,
                 const Mantid::HistogramData::Points &lambdas,
                 double lambdaFixed,
                 Mantid::API::ISpectrum &attenuationFactorsSpectrum);

private:
  const IBeamProfile &m_beamProfile;
  MCInteractionVolume m_scatterVol;
  const size_t m_nevents;
  const size_t m_maxScatterAttempts;
  const double m_error;
};

} // namespace Algorithms
} // namespace Mantid

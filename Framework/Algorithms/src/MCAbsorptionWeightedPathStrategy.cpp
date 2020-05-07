// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidAlgorithms/MCAbsorptionWeightedPathStrategy.h"
#include "MantidAlgorithms/SampleCorrections/IBeamProfile.h"
#include "MantidKernel/PseudoRandomNumberGenerator.h"
#include "MantidKernel/V3D.h"

#include "MantidAlgorithms/SampleCorrections/RectangularBeamProfile.h"
#include "MantidGeometry/Objects/CSGObject.h"

namespace Mantid {
using Kernel::PseudoRandomNumberGenerator;

namespace Algorithms {

/**
 * Constructor
 * @param beamProfile A reference to the object the beam profile
 * @param sample A reference to the object defining details of the sample
 * @param EMode The energy mode of the instrument
 * @param nevents The number of Monte Carlo events used in the simulation
 * @param nlambda The number of steps to use across the wavelength range
 * @param maxScatterPtAttempts The maximum number of tries to generate a random
 * @param useSparseInstrument Whether a sparse instrument representation is
 * being used
 * @param interpolateOpt Method to be used to interpolate between wavelength
 * points
 * @param regenerateTracksForEachLambda Whether to resimulate tracks for each
 * wavelength point or not
 * @param logger Logger from parent algorithm to write logging info
 * point within the object.
 */
MCAbsorptionWeightedPathStrategy::MCAbsorptionWeightedPathStrategy(
    const IBeamProfile &beamProfile, const API::Sample &sample,
    const size_t nevents, const size_t maxScatterPtAttempts,
    Kernel::Logger &logger)
    : m_beamProfile(beamProfile),
      m_scatterVol(MCInteractionVolume(
          sample, beamProfile.defineActiveRegion(sample), logger)),
      m_nevents(nevents), m_maxScatterAttempts(maxScatterPtAttempts),
      m_error(1.0 / std::sqrt(m_nevents)) {}

/**
 * Compute the weighted path length for a final position of the neutron and
 * wavelengths before and after scattering
 * @param rng A reference to a PseudoRandomNumberGenerator
 * @param finalPos Defines the final position of the neutron, assumed to be
 * where it is detected
 */
void MCAbsorptionWeightedPathStrategy::calculate(
    Kernel::PseudoRandomNumberGenerator &rng, const Kernel::V3D &finalPos,
    const Mantid::HistogramData::Points &lambdas, double lambdaFixed,
    Mantid::API::ISpectrum &attenuationFactorsSpectrum) {
  const auto scatterBounds = m_scatterVol.getBoundingBox();
  auto &attenuationFactors = attenuationFactorsSpectrum.mutableY();

  for (size_t i = 0; i < m_nevents; ++i) {
    Geometry::Track beforeScatter;
    Geometry::Track afterScatter;
  }

  m_scatterVol.generateScatterPointStats();

  attenuationFactors = attenuationFactors / static_cast<double>(m_nevents);

  // Interpolate through points not simulated. Simulation WS only has
  // reduced X values if using sparse instrument so no interpolation required
  auto attenuationFactorsHist = attenuationFactorsSpectrum.histogram();

  std::fill(attenuationFactorsHist.mutableE().begin(),
            attenuationFactorsHist.mutableE().end(), m_error);
  // ISpectrum::histogram() returns a value rather than reference so need to
  // reapply
  attenuationFactorsSpectrum.setHistogram(attenuationFactorsHist);
}

} // namespace Algorithms
} // namespace Mantid

// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2014 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef MANTID_GEOMETRY_PRODUCTGROUP_H_
#define MANTID_GEOMETRY_PRODUCTGROUP_H_

#include "MantidGeometry/Crystal/Group.h"
#include "MantidGeometry/DllConfig.h"

namespace Mantid {
namespace Geometry {

/**
    @class ProductOfCyclicGroups

    ProductOfCyclicGroups expands a bit on the explanations given in
    CyclicGroup. As shown for example in [1], some point groups
    cannot be expressed solely as a cyclic group. Instead it's
    necessary to multiply two or three cyclic groups to obtain all
    symmetry operations of that group.

    For this purpose, ProductOfCyclicGroups was created. It takes a set of
    n symmetry operations, each of which is seen as a generator of
    a cyclic group C_i. The resulting n groups ("factor groups") are
    multiplied to form a product group G:

        G = C_1 * C_2 * ... * C_n

    Where C_i is generated by the symmetry operation S_i. The notation
    in code to generate even large groups from a few generators
    becomes very short using this class:

        Group_const_sptr pointGroup422 =
            GroupFactory::create<ProductOfCyclicGroups>("-y,x,z; x,-y,-z");

    This is for example used in SpaceGroupFactory to create space groups
    from a small set of generators supplied in the International Tables
    for Crystallography A.

    [1] Shmueli, U. Acta Crystallogr. A 40, 559–567 (1984).
        http://dx.doi.org/10.1107/S0108767384001161

      @author Michael Wedel, Paul Scherrer Institut - SINQ
      @date 08/10/2014
  */
class MANTID_GEOMETRY_DLL ProductOfCyclicGroups : public Group {
public:
  ProductOfCyclicGroups(const std::string &generators);
  ProductOfCyclicGroups(const std::vector<Group_const_sptr> &factorGroups);

protected:
  Group_const_sptr getGeneratedGroup(const std::string &generators) const;
  std::vector<Group_const_sptr> getFactorGroups(
      const std::vector<SymmetryOperation> &symmetryOperations) const;
  Group_const_sptr getProductOfCyclicGroups(
      const std::vector<Group_const_sptr> &factorGroups) const;
};

} // namespace Geometry
} // namespace Mantid

#endif /* MANTID_GEOMETRY_PRODUCTGROUP_H_ */

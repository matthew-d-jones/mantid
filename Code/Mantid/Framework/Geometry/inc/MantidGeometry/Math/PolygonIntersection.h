#ifndef MANTID_GEOMETRY_POLYGONINTERSECTION_H_
#define MANTID_GEOMETRY_POLYGONINTERSECTION_H_

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include "MantidGeometry/DllConfig.h"
#include "MantidGeometry/Math/ConvexPolygon.h"

namespace Mantid
{
  namespace Geometry
  {
    /** 
    This header defines various algorithms for defining polygon intersection

    @author Martyn Gigg, Tessella plc

    Copyright &copy; 2011 ISIS Rutherford Appleton Laboratory & NScD Oak Ridge National Laboratory

    This file is part of Mantid.

    Mantid is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Mantid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    File change history is stored at: <https://svn.mantidproject.org/mantid/trunk/Code/Mantid>
    Code Documentation is available at: <http://doxygen.mantidproject.org>
    */

    /// Define an enum for the intersection algorithm to use
    namespace PolygonIntersection
    {
      enum Method{Laszlo, ORourke};
    };

    /// Compute the intersection of two polygons using a chosen method
    MANTID_GEOMETRY_DLL 
    ConvexPolygon intersection(const ConvexPolygon &p, const ConvexPolygon &q, 
      PolygonIntersection::Method method = PolygonIntersection::ORourke);

  } //namespace Geometry
} //namespace Mantid

#endif //MANTID_GEOMETRY_POLYGONINTERSECTION_H_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "MantidKernel/Logger.h"
#include "MantidKernel/Exception.h"
#include "MantidGeometry/Matrix.h"
#include "MantidGeometry/V3D.h"
#include "MantidGeometry/Track.h"
#include "MantidGeometry/Surface.h"


//const double surfaceTolerance(1e-5);       ///< Below this two point touch.

namespace Mantid
{

  namespace Geometry
  {

    Kernel::Logger& Track::PLog(Kernel::Logger::get("Track"));

    TUnit::TUnit(const Geometry::V3D& A,const Geometry::V3D& B,
      const double D,const int ID) :
    PtA(A),PtB(B),Dist(D),Length(A.distance(B)),ObjID(ID)
      /*!
      Constuctor
      \param A :: V3D point to start
      \param B :: V3D point to end track
      \param D :: Total distance from start of track
      \param ID :: ID number
      */
    {}

    TPartial::TPartial(const int ID,const int flag,
      const Geometry::V3D& PVec,const double D) :
    ObjID(ID),Direction(flag),PtA(PVec),Dist(D)
      /*!
      Constuctor
      \param PVec :: V3D point to end track
      \param ID :: ID number
      \param flag :: Indicated the direction
      \param D :: Total distance from start of track
      */
    {}

    int
      TPartial::operator<(const TPartial& A) const
      /*!
      This calculates < based on the length AND
      the flag state. Such that is the length are
      closer than +/- Tolerance then the flag state
      take precidence. 
      \param A :: TPartial object to compare
      \return A is < this
      */
    {
      const double dL=fabs(Dist-A.Dist);
      return (dL>Surface::getSurfaceTolerance()) ? 
        Dist<A.Dist : Direction<A.Direction;
    }

    //---------------------------------------------
    //           TRACK
    // --------------------------------------------

    Track::Track(const Geometry::V3D& StartPt,
      const Geometry::V3D& UV,const int initObj) : 
    iPt(StartPt),uVec(UV),iObj(initObj)
      /*!
      Constructor
      \param StartPt :: Initial Point
      \param UV :: unit vector of direction
      \param initObj :: inital object identifier
      */ 
    {}

    Track::Track(const Track& A) :
    iPt(A.iPt),uVec(A.uVec),iObj(A.iObj),
      Link(A.Link),surfPoints(A.surfPoints)
      /*!
      Copy Constructor
      \param A :: Track to copy
      */ 
    {}

    Track&
      Track::operator=(const Track& A)
      /*!
      Assignment operator
      \param A :: Track to copy
      \return *this
      */ 
    {
      if (this != &A)
      {
        iPt=A.iPt;
        uVec=A.uVec;
        iObj=A.iObj;
        Link=A.Link;
        surfPoints=A.surfPoints;
      }
      return *this;
    }

    Track::~Track()
      /*!
      Destructor
      */
    {}

    void 
      Track::setFirst(const Geometry::V3D& StartPoint,
      const Geometry::V3D& UV)
      /*!
      Sets the first Point
      \param StartPoint :: First Point
      \param UV :: Unit vector
      */
    {
      iPt=StartPoint;
      uVec=UV;
      return;
    }

    int
      Track::nonComplete() const
      /*!
      Determines if the track is complete
      \retval 0 :: Complete from Init to end without gaps
      \retval +ve :: Index number of incomplete segment +1
      */
    {
      const double TrackTolerance(1e-6);
      if (Link.size()<2)
        return 0;

      LType::const_iterator ac=Link.begin();
      if (iPt.distance(ac->PtA)>TrackTolerance)
        return 1;
      LType::const_iterator bc=ac;
      bc++;

      while(bc!=Link.end())
      {
        if ((ac->PtB).distance(bc->PtA)>TrackTolerance)
          return distance(Link.begin(),bc)+1;
        ac++;
        bc++;
      }
      // success
      return 0;
    }

    void
      Track::removeCoJoins()
      /*!
      Remove touching TUnits that have identical
      components
      */
    {
      // No work to do:
      if (Link.empty())
        return; 
      // ac == previous : bc = next node.
      LType::iterator ac=Link.begin();
      LType::iterator bc=Link.begin();
      bc++;
      while(bc!=Link.end())
      {
        if (ac->ObjID==bc->ObjID)
        {
          ac->PtB=bc->PtB;
          ac->Dist=(ac->PtA).distance(ac->PtB);
          ac->Length=bc->Length;
          Link.erase(bc);
          bc=ac;
          bc++;
        }
        else
        {
          ac++;
          bc++;
        }
      }
      return;
    }

    void Track::addPoint(const int ID,const int Direct,
      const Geometry::V3D& Pt) 
      /*!
      Objective is to merge in partial information
      about the beginning and end of the tracks.
      We do not need to keep surfPoints in order
      because that will be done when they are converted into
      TUnits.  

      \param ID :: Id number of the object
      \param Direct :: direction of travel
      \param Pt :: Point to go
      */
    {
      surfPoints.push_back(TPartial(ID,Direct,Pt,Pt.distance(iPt)));
      return;
    }

    int Track::addTUnit(const int ID,const Geometry::V3D& Apt,
      const Geometry::V3D& Bpt,const double D)
      /*!
      This adds a whole segment to the track : This currently assumes that links are added in order
      \param ID :: Id number of the object
      \param Apt :: first Point
      \param Bpt :: second Point
      \param D :: Distance along track
      \retval Index point 
      */
    {
      // Process First Point
      TUnit newTUnit(Apt,Bpt,D,ID);
      if (Link.empty())
      {
        Link.push_back(newTUnit);
        return 0;
      }
      std::vector<TUnit>::iterator xV = 
        lower_bound(Link.begin(),Link.end(),newTUnit);

      //must extract the distance before you insert otherwise the iterators are incompatible
      int index = distance(Link.begin(),xV);

      Link.insert(xV,newTUnit);

      return index;
    }

    void Track::buildLink()
      /*!
      Builds a set of linking track components.
      This version deals with touching surfaces 
      */
    {
      if (surfPoints.empty())
        return;

      // First sort surfPoints
      sort(surfPoints.begin(),surfPoints.end());
      PType::const_iterator ac=surfPoints.begin();
      PType::const_iterator bc=ac;
      bc++;
      Geometry::V3D workPt=iPt;            // last good point
      // First point is not necessarily in an object
      // Process first point:
      while(ac!=surfPoints.end() && ac->Direction != 1)    // stepping from an object.
      {
        if (ac->Direction==-1)
        {
          addTUnit(ac->ObjID,iPt,ac->PtA,ac->Dist);  // from the void
          workPt=ac->PtA;
        }
        ac++;
        if (bc!=surfPoints.end())
          bc++;
      } 

      //have we now passed over all of the potential intersections without actually hitting the object
      if (ac==surfPoints.end())
      {
        //yes
        surfPoints.clear();
        return;
      }

      workPt=ac->PtA;      

      while(bc!=surfPoints.end())      // Since bc > ac
      {
        if (ac->Direction==1 && bc->Direction==-1)
        {
          // Touching surface / identical surface
          if (fabs(ac->Dist-bc->Dist)>Surface::getSurfaceTolerance())
          {
            // track leave ac into bc.
            addTUnit(ac->ObjID,ac->PtA,bc->PtA,bc->Dist);
          }
          // Points with intermediate void
          else
          {
            addTUnit(ac->ObjID,workPt,ac->PtA,ac->Dist);
          }
          workPt=bc->PtA;

          // ADDING to ac twice: since processing pairs
          ac++;
          ac++;
          bc++;    // can I do this past the end ? 
          if (bc!=surfPoints.end())
            bc++;
        }
        else         // Test for glacing point / or void edges
        {          // These all can be skipped
          ac++;
          bc++;
        }
      }	

      surfPoints.clear();        // While vector 
      return;
    }

  } // NAMESPACE Geometry

}  // NAMESPACE Mantid

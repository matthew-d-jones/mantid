#ifndef LOADEMPTYINSTRUMENTTEST_H_
#define LOADEMPTYINSTRUMENTTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidDataHandling/LoadEmptyInstrument.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/Instrument.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidKernel/Exception.h"
#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/Workspace.h"
#include "MantidAPI/Algorithm.h"
#include "MantidAPI/SpectraDetectorMap.h"
#include "MantidGeometry/Instrument/Component.h"
#include <vector>

using namespace Mantid::API;
using namespace Mantid::Kernel;
using namespace Mantid::Geometry;
using namespace Mantid::DataHandling;
using namespace Mantid::DataObjects;

class LoadEmptyInstrumentTest : public CxxTest::TestSuite
{
public:

  void testExecSLS()
  {
    LoadEmptyInstrument loaderSLS;

    TS_ASSERT_THROWS_NOTHING(loaderSLS.initialize());
    TS_ASSERT( loaderSLS.isInitialized() );
    loaderSLS.setPropertyValue("Filename", "../../../../Test/Instrument/SANDALS_Definition.xml");
    inputFile = loaderSLS.getPropertyValue("Filename");
    wsName = "LoadEmptyInstrumentTestSLS";
    loaderSLS.setPropertyValue("OutputWorkspace", wsName);

    std::string result;
    TS_ASSERT_THROWS_NOTHING( result = loaderSLS.getPropertyValue("Filename") )
    TS_ASSERT_EQUALS( result, inputFile);

    TS_ASSERT_THROWS_NOTHING( result = loaderSLS.getPropertyValue("OutputWorkspace") )
    TS_ASSERT( ! result.compare(wsName));

    TS_ASSERT_THROWS_NOTHING(loaderSLS.execute());

    TS_ASSERT( loaderSLS.isExecuted() );


    MatrixWorkspace_sptr output;
    output = boost::dynamic_pointer_cast<MatrixWorkspace>(AnalysisDataService::Instance().retrieve(wsName));
    
    // Check the total number of elements in the map for SLS
    TS_ASSERT_EQUALS(output->spectraMap().nElements(),683);
  }


  void testParameterTags()
  {
    LoadEmptyInstrument loader;

    TS_ASSERT_THROWS_NOTHING(loader.initialize());
    TS_ASSERT( loader.isInitialized() );
    loader.setPropertyValue("Filename", "../../../../Test/Instrument/IDF_for_unit_testing2.xml");
    inputFile = loader.getPropertyValue("Filename");
    wsName = "LoadEmptyInstrumentParamTest";
    loader.setPropertyValue("OutputWorkspace", wsName);

    TS_ASSERT_THROWS_NOTHING(loader.execute());
    TS_ASSERT( loader.isExecuted() );


    MatrixWorkspace_sptr ws;
    ws = boost::dynamic_pointer_cast<MatrixWorkspace>(AnalysisDataService::Instance().retrieve(wsName));
    ws->populateInstrumentParameters();

    // get parameter map
    ParameterMap& paramMap = ws->instrumentParameters();

    // check that parameter have been read into the instrument parameter map
    std::vector<V3D> ret1 = paramMap.getV3D("monitors", "pos");
    TS_ASSERT_DELTA( ret1[0].X(), 5.0, 0.0001);
    TS_ASSERT_DELTA( ret1[0].Y(), 0.0, 0.0001);
    TS_ASSERT_DELTA( ret1[0].Z(), 0.0, 0.0001);

    // get detector corresponding to workspace index 0
    IDetector_sptr det = ws->getDetector(0);  

    TS_ASSERT_EQUALS( det->getID(), 1001);
    TS_ASSERT_EQUALS( det->getName(), "upstream_monitor_det");

    Parameter_sptr param = paramMap.get(&(*det), "boevs2");
    TS_ASSERT_DELTA( param->value<double>(), 16.0, 0.0001);

    param = paramMap.get(&(*det), "boevs3");
    TS_ASSERT_DELTA( param->value<double>(), 32.0, 0.0001);

    param = paramMap.get(&(*det), "boevs");
    TS_ASSERT( param == NULL );

    param = paramMap.getRecursive(&(*det), "boevs", "spade");
    TS_ASSERT_DELTA( param->value<double>(), 8.0, 0.0001);


    // check reserved keywords
    std::vector<double> dummy = paramMap.getDouble("nickel-holder", "klovn");
    TS_ASSERT_DELTA( dummy[0], 2.0, 0.0001);
    dummy = paramMap.getDouble("nickel-holder", "pos");
    TS_ASSERT_EQUALS (dummy.size(), 0);
    dummy = paramMap.getDouble("nickel-holder", "rot");
    TS_ASSERT_EQUALS (dummy.size(), 0);
    dummy = paramMap.getDouble("nickel-holder", "taabe");
    TS_ASSERT_DELTA (dummy[0], 200.0, 0.0001);
    dummy = paramMap.getDouble("nickel-holder", "mistake");
    TS_ASSERT_EQUALS (dummy.size(), 0);


    AnalysisDataService::Instance().remove(wsName);
  }

  void testToscaParameterTags()
  {
    LoadEmptyInstrument loader;

    TS_ASSERT_THROWS_NOTHING(loader.initialize());
    TS_ASSERT( loader.isInitialized() );
    loader.setPropertyValue("Filename", "../../../../Test/Instrument/TOSCA_Definition.xml");
    inputFile = loader.getPropertyValue("Filename");
    wsName = "LoadEmptyInstrumentParamToscaTest";
    loader.setPropertyValue("OutputWorkspace", wsName);

    TS_ASSERT_THROWS_NOTHING(loader.execute());
    TS_ASSERT( loader.isExecuted() );


    MatrixWorkspace_sptr ws;
    ws = boost::dynamic_pointer_cast<MatrixWorkspace>(AnalysisDataService::Instance().retrieve(wsName));
    ws->populateInstrumentParameters();

    // get parameter map
    ParameterMap& paramMap = ws->instrumentParameters();

    // get detector corresponding to workspace index 0
    IDetector_sptr det = ws->getDetector(69);  

    TS_ASSERT_EQUALS( det->getID(), 78);
    TS_ASSERT_EQUALS( det->getName(), "Detector #70");

    Parameter_sptr param = paramMap.get(&(*det), "Efixed");
    TS_ASSERT_DELTA( param->value<double>(), 4.000, 0.0001);

    AnalysisDataService::Instance().remove(wsName);
  }

void testCheckIfVariousInstrumentsLoad()
  {
    LoadEmptyInstrument loader;

    TS_ASSERT_THROWS_NOTHING(loader.initialize());
    TS_ASSERT( loader.isInitialized() );
    loader.setPropertyValue("Filename", "../../../../Test/Instrument/SANS2D_Definition.xml");
    inputFile = loader.getPropertyValue("Filename");
    wsName = "LoadEmptyInstrumentParaSans2dTest";
    loader.setPropertyValue("OutputWorkspace", wsName);

    TS_ASSERT_THROWS_NOTHING(loader.execute());
    TS_ASSERT( loader.isExecuted() );

    AnalysisDataService::Instance().remove(wsName);

    LoadEmptyInstrument loaderPolRef;
    loaderPolRef.initialize();
    loaderPolRef.setPropertyValue("Filename", "../../../../Test/Instrument/POLREF_Definition.xml");
    inputFile = loaderPolRef.getPropertyValue("Filename");
    wsName = "LoadEmptyInstrumentParamPOLREFTest";
    loaderPolRef.setPropertyValue("OutputWorkspace", wsName);

    TS_ASSERT_THROWS_NOTHING(loaderPolRef.execute());
    TS_ASSERT( loaderPolRef.isExecuted() );

    AnalysisDataService::Instance().remove(wsName);
  }


private:
  std::string inputFile;
  std::string wsName;

};

#endif /*LOADEMPTYINSTRUMENTTEST_H_*/

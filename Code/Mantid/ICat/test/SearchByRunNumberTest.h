#ifndef SEARCHBYADVANCED_H_
#define SEARCHBYADVANCED_H_

#include <cxxtest/TestSuite.h>
#include "MantidICat/SearchByRunNumber.h"
#include "MantidICat/Session.h"
#include "MantidICat/Login.h"
#include "MantidDataObjects/WorkspaceSingleValue.h"

using namespace Mantid;
using namespace Mantid::ICat;
class CSearchByAdvancedTest: public CxxTest::TestSuite
{
public:
	
	void testInit()
	{
		CSearchByRunNumber searchobj;
		Login loginobj;
		TS_ASSERT_THROWS_NOTHING( searchobj.initialize());
		TS_ASSERT( searchobj.isInitialized() );
	}
	void testSearchByRunNumberandInstrument()
	{
		/*std::string s;
		std::getline(std::cin,s);*/

		CSearchByRunNumber searchobj;
		Login loginobj;
		Session::Instance();
		if ( !loginobj.isInitialized() ) loginobj.initialize();

		loginobj.setPropertyValue("Username", "mantid_test");
		loginobj.setPropertyValue("Password", "mantidtestuser");
		//loginobj.setPropertyValue("DBServer", "");
		
		TS_ASSERT_THROWS_NOTHING(loginobj.execute());
		TS_ASSERT( loginobj.isExecuted() );

		if ( !searchobj.isInitialized() ) searchobj.initialize();
		
		searchobj.setPropertyValue("StartRun", "100.0");
		searchobj.setPropertyValue("EndRun", "109.0");
		searchobj.setPropertyValue("Instrument","LOQ");
		searchobj.setPropertyValue("OutputWorkspace","SearchBy_RunNumber");
				
		TS_ASSERT_THROWS_NOTHING(searchobj.execute());
		TS_ASSERT( searchobj.isExecuted() );

	}
	void testSearchByKeywords()
	{
		
		/*std::string s;
		std::getline(std::cin,s);*/
		
		CSearchByRunNumber searchobj;
		Login loginobj;
		Session::Instance();
		if ( !loginobj.isInitialized() ) loginobj.initialize();

		loginobj.setPropertyValue("Username", "mantid_test");
		loginobj.setPropertyValue("Password", "mantidtestuser");
		//loginobj.setPropertyValue("DBServer", "");
		
		TS_ASSERT_THROWS_NOTHING(loginobj.execute());
		TS_ASSERT( loginobj.isExecuted() );

		if ( !searchobj.isInitialized() ) searchobj.initialize();
				
		searchobj.setPropertyValue("Keywords","000117");
		searchobj.setPropertyValue("Instrument","HRPD");
		searchobj.setPropertyValue("OutputWorkspace","SearchBy_RunNumber");
				
		TS_ASSERT_THROWS_NOTHING(searchobj.execute());
		TS_ASSERT( searchobj.isExecuted() );

	}
	void testSearchBybyStartDate()
	{
		/*std::string s;
		std::getline(std::cin,s);*/

		CSearchByRunNumber searchobj;
		Login loginobj;
		Session::Instance();
		if ( !loginobj.isInitialized() ) loginobj.initialize();

		loginobj.setPropertyValue("Username", "mantid_test");
		loginobj.setPropertyValue("Password", "mantidtestuser");
		//loginobj.setPropertyValue("DBServer", "");
		
		TS_ASSERT_THROWS_NOTHING(loginobj.execute());
		TS_ASSERT( loginobj.isExecuted() );

		if ( !searchobj.isInitialized() ) searchobj.initialize();
		
		//searchobj.setPropertyValue("StartRun", "100.0");
		//searchobj.setPropertyValue("EndRun", "109.0");
		searchobj.setPropertyValue("Instrument","MERLIN");
		searchobj.setPropertyValue("StartDate","10/08/2008");
		searchobj.setPropertyValue("EndDate","22/08/2008");
		searchobj.setPropertyValue("OutputWorkspace","SearchBy_RunNumber");
				
		TS_ASSERT_THROWS_NOTHING(searchobj.execute());
		TS_ASSERT( searchobj.isExecuted() );

	}
	void testSearchByRunNumberInvalidInput()
	{		
		Login loginobj;
		Session::Instance();
		if ( !loginobj.isInitialized() ) loginobj.initialize();

		// Now set it...
		loginobj.setPropertyValue("Username", "mantid_test");
		loginobj.setPropertyValue("Password", "mantidtestuser");
		//loginobj.setPropertyValue("DBServer", "");
		
		TS_ASSERT_THROWS_NOTHING(loginobj.execute());
		TS_ASSERT( loginobj.isExecuted() );
		CSearchByRunNumber searchobj;
		if ( !searchobj.isInitialized() ) searchobj.initialize();
		
		// start run number > end run number
		searchobj.setPropertyValue("StartRun", "150.0");
		searchobj.setPropertyValue("EndRun", "102.0");
		searchobj.setPropertyValue("Instrument","LOQ");
				
    	searchobj.setPropertyValue("OutputWorkspace","SearchBy_RunNumber");
		TS_ASSERT_THROWS_NOTHING(searchobj.execute());
		//should fail
		TS_ASSERT( !searchobj.isExecuted() );

	}
		
};
#endif
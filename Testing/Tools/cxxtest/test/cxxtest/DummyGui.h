// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

//
// The DummyGui is a "GUI" that prints messages to cout
// It is used for testing CxxTest
//

#include <cxxtest/Gui.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/ErrorPrinter.h>

namespace CxxTest
{
    class DummyGui : public GuiListener
    {
    public:
        void guiEnterWorld( unsigned numTotalTests )
        {
            (CXXTEST_STD(cout) << " {Start " << numTotalTests << " tests} ").flush();
        }
        
        void guiEnterTest( const char *suiteName, const char *testName )
        {
            (CXXTEST_STD(cout) << " {" << suiteName << "::" << testName << "()} ").flush();
        }
        
        void yellowBar()
        {
            (CXXTEST_STD(cout) << " {Yellow} ").flush();
        }

        void redBar()
        {
            (CXXTEST_STD(cout) << " {Red} ").flush();
        }

        void leaveWorld( const WorldDescription & )
        {
            (CXXTEST_STD(cout) << " {Stop} ").flush();
        }
    };
}

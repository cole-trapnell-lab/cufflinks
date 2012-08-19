/*
 *  test_abundances.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/23/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

//class AutoFail : public CppUnit::TestFixture
//{
//	CPPUNIT_TEST_SUITE( AutoFail );
//	CPPUNIT_TEST( testFalse );
//	CPPUNIT_TEST( testNotEqual );
//	CPPUNIT_TEST_SUITE_END();
//	
//public:	
//	void testFalse();
//	void testNotEqual();
//};
//
//void AutoFail::testFalse()
//{
//	CPPUNIT_ASSERT(false);
//}
//
//void AutoFail::testNotEqual()
//{
//	CPPUNIT_ASSERT_EQUAL(1, 2);
//}
//
//// For testing Xcode integration and CppUnit installation
//CPPUNIT_TEST_SUITE_REGISTRATION( AutoFail );
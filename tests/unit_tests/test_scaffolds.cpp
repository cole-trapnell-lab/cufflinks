/*
 *  test_scaffolds.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 1/27/11.
 *  Copyright 2011 Cole Trapnell. All rights reserved.
 *
 */

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "common.h"
#include "scaffolds.h"

class AugmentedCuffOpTests : public CppUnit::TestFixture
{
    // Nomenclature:
    // I = Intron
    // M = Match
    // U = Unknown
    // d = disjoint (no olap)
    // s = olap <= bowtie_overhang_tolerance
    // l = olap > bowtie_overhang_tolerance
	CPPUNIT_TEST_SUITE( AugmentedCuffOpTests );

    CPPUNIT_TEST( testMMs );
    CPPUNIT_TEST( testMMl );

    CPPUNIT_TEST( testMUs );
    CPPUNIT_TEST( testMUl );

    CPPUNIT_TEST( testIMs );
    CPPUNIT_TEST( testIMl );

    CPPUNIT_TEST( testIUs );
    CPPUNIT_TEST( testIUl );

    CPPUNIT_TEST( testUUs );
    CPPUNIT_TEST( testUUl );

    CPPUNIT_TEST( testIIs );
    CPPUNIT_TEST( testIIl );
    
	CPPUNIT_TEST_SUITE_END();
	
public:	
    void testMMs();
    void testMMl();
    
    void testMUs();
    void testMUl();
    
    void testIMs();
    void testIMl();
    
    void testIUs();
    void testIUl();
    
    void testUUs();
    void testUUl();
    
    void testIIs();
    void testIIl(); 
};

void AugmentedCuffOpTests::testMMs()
{
    AugmentedCuffOp lhs(CUFF_MATCH, 0, 50);
    AugmentedCuffOp rhs(CUFF_MATCH, 50 - bowtie_overhang_tolerance, 50);
	CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs));
}

void AugmentedCuffOpTests::testMMl()
{
    AugmentedCuffOp lhs(CUFF_MATCH, 0, 50);
    AugmentedCuffOp rhs(CUFF_MATCH, 50 - 2 * bowtie_overhang_tolerance, 50);
	CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs));
}

void AugmentedCuffOpTests::testMUs()
{
    AugmentedCuffOp lhs(CUFF_MATCH, 0, 50);
    AugmentedCuffOp rhs(CUFF_UNKNOWN, 50 - bowtie_overhang_tolerance, 50);
	CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs));
}

void AugmentedCuffOpTests::testMUl()
{
    AugmentedCuffOp lhs(CUFF_MATCH, 0, 50);
    AugmentedCuffOp rhs(CUFF_UNKNOWN, 50 - 2 * bowtie_overhang_tolerance, 50);
	CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs));
}

void AugmentedCuffOpTests::testIMs()
{
    AugmentedCuffOp lhs(CUFF_INTRON, 0, 50);
    AugmentedCuffOp rhs(CUFF_MATCH, 50 - bowtie_overhang_tolerance, 50);
	CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs));
}

void AugmentedCuffOpTests::testIMl()
{
    AugmentedCuffOp lhs(CUFF_INTRON, 0, 50);
    AugmentedCuffOp rhs(CUFF_MATCH, 50 - 2 * bowtie_overhang_tolerance, 50);
	CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs) == false);
}

// Note: the cases below will introduce transitivity hazards if unknowns overlapping
// non-constitutive regions are not ablated before assembly.  We need this
// compatibility rule for constitutive regions and all quantification
void AugmentedCuffOpTests::testIUs()
{
    AugmentedCuffOp lhs(CUFF_INTRON, 0, 50);
    AugmentedCuffOp rhs(CUFF_UNKNOWN, 50 - bowtie_overhang_tolerance, 50);
	CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs));
}

void AugmentedCuffOpTests::testIUl()
{
    AugmentedCuffOp lhs(CUFF_INTRON, 0, 50);
    AugmentedCuffOp rhs(CUFF_UNKNOWN, 50 - 2 * bowtie_overhang_tolerance, 50);
    CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs));
}


void AugmentedCuffOpTests::testUUs()
{
    AugmentedCuffOp lhs(CUFF_UNKNOWN, 0, 50);
    AugmentedCuffOp rhs(CUFF_UNKNOWN, 50 - bowtie_overhang_tolerance, 50);
	CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs));
}

void AugmentedCuffOpTests::testUUl()
{
    AugmentedCuffOp lhs(CUFF_UNKNOWN, 0, 50);
    AugmentedCuffOp rhs(CUFF_UNKNOWN, 50 - 2 * bowtie_overhang_tolerance, 50);
    CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs));
}

void AugmentedCuffOpTests::testIIs()
{
    AugmentedCuffOp lhs(CUFF_INTRON, 0, 50);
    AugmentedCuffOp rhs(CUFF_INTRON, 50 - 2 * bowtie_overhang_tolerance, 50);
    CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs) == false);
}

void AugmentedCuffOpTests::testIIl()
{
    AugmentedCuffOp lhs(CUFF_INTRON, 0, 50);
    AugmentedCuffOp rhs(CUFF_INTRON, 50 - 2 * bowtie_overhang_tolerance, 50);
    CPPUNIT_ASSERT(AugmentedCuffOp::compatible(lhs, rhs) == false);
}


// For testing Xcode integration and CppUnit installation
CPPUNIT_TEST_SUITE_REGISTRATION( AugmentedCuffOpTests );

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
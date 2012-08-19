/*
 *  test_main.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/23/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <cppunit/Portability.h>
#include <cppunit/Outputter.h>
#include <cppunit/portability/Stream.h>
#include <cppunit/config/SourcePrefix.h>
#include <cppunit/Exception.h>
#include <cppunit/SourceLine.h>
#include <cppunit/TestFailure.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/CompilerOutputter.h>
#include <algorithm>
#include <cppunit/tools/StringTools.h>

#include <string>

using namespace std;
using namespace CppUnit;

class CustomFormatter : public Outputter
{
public:
	/*! \brief Constructs a CompilerOutputter object.
	 * \param result Result of the test run.
	 * \param stream Stream used to output test result.
	 * \param locationFormat Error location format used by your compiler. Default
	 *                       to \c CPPUNIT_COMPILER_LOCATION_FORMAT which is defined
	 *                       in the configuration file. See setLocationFormat() for detail.
	 * \see setLocationFormat().
	 */
	CustomFormatter( TestResultCollector *result,
					  OStream &stream,
					  const std::string &locationFormat = "%p:%l: error: " );
	
	/// Destructor.
	virtual ~CustomFormatter();
	
	/*! \brief Sets the error location format.
	 * 
	 * Indicates the format used to report location of failed assertion. This format should
	 * match the one used by your compiler.
	 *
	 * The location format is a string in which the occurence of the following character
	 * sequence are replaced:
	 *
	 * - "%l" => replaced by the line number
	 * - "%p" => replaced by the full path name of the file ("G:\prg\vc\cppunit\MyTest.cpp")
	 * - "%f" => replaced by the base name of the file ("MyTest.cpp")
	 *
	 * Some examples:
	 *
	 * - VC++ error location format: "%p(%l):" => produce "G:\prg\MyTest.cpp(43):"
	 * - GCC error location format: "%f:%l:" => produce "MyTest.cpp:43:"
	 * 
	 * Thoses are the two compilers currently <em>supported</em> (gcc format is used if
	 * VC++ is not detected). If you want your compiler to be automatically supported by
	 * CppUnit, send a mail to the mailing list (preferred), or submit a feature request
	 * that indicates how to detect your compiler with the preprocessor (\#ifdef...) and
	 * your compiler location format.
	 */
	void setLocationFormat( const std::string &locationFormat );
	
	/*! \brief Creates an instance of an outputter that matches your current compiler.
	 * \deprecated This class is specialized through parameterization instead of subclassing...
	 *             Use CompilerOutputter::CompilerOutputter instead.
	 */
	static CustomFormatter *defaultOutputter(TestResultCollector *result,
											 OStream &stream );
	
	void write();
	
	void setNoWrap();
	
	void setWrapColumn( int wrapColumn );
	
	int wrapColumn() const;
	
	virtual void printSuccess();
	virtual void printFailureReport();
	virtual void printFailuresList();
	virtual void printStatistics();
	virtual void printFailureDetail( TestFailure *failure );
	virtual void printFailureLocation( SourceLine sourceLine );
	virtual void printFailureType( TestFailure *failure );
	virtual void printFailedTestName( TestFailure *failure );
	virtual void printFailureMessage( TestFailure *failure );
	
private:
	/// Prevents the use of the copy constructor.
	CustomFormatter( const CustomFormatter &copy );
	
	/// Prevents the use of the copy operator.
	void operator =( const CustomFormatter &copy );
	
	virtual bool processLocationFormatCommand( char command, 
											  const SourceLine &sourceLine );
	
	virtual std::string extractBaseName( const std::string &fileName ) const;
	
private:
	TestResultCollector *m_result;
	OStream &m_stream;
	std::string m_locationFormat;
	int m_wrapColumn;
};

CustomFormatter::CustomFormatter( TestResultCollector *result,
									 OStream &stream,
									 const std::string &locationFormat )
: m_result( result )
, m_stream( stream )
, m_locationFormat( locationFormat )
, m_wrapColumn( CPPUNIT_WRAP_COLUMN )
{
}


CustomFormatter::~CustomFormatter()
{
}


void CustomFormatter::setLocationFormat( const std::string &locationFormat )
{
	m_locationFormat = locationFormat;
}


CustomFormatter* CustomFormatter::defaultOutputter( TestResultCollector *result,
									OStream &stream )
{
	return new CustomFormatter( result, stream );
}


void CustomFormatter::write()
{
	if ( m_result->wasSuccessful() )
		printSuccess();
	else
		printFailureReport();
}


void CustomFormatter::printSuccess()
{
	m_stream  << "OK (" << m_result->runTests()  << ")\n";
}


void CustomFormatter::printFailureReport()
{
	printFailuresList();
	printStatistics();
}


void CustomFormatter::printFailuresList()
{
	for ( int index =0; index < m_result->testFailuresTotal(); ++index)
	{
		printFailureDetail( m_result->failures()[ index ] );
	}
}


void CustomFormatter::printFailureDetail( TestFailure *failure )
{
	printFailureLocation( failure->sourceLine() );
	//printFailureType( failure );
	printFailedTestName( failure ); 
	m_stream << " ";
	printFailureMessage( failure );
	m_stream << endl;
}


void CustomFormatter::printFailureLocation( SourceLine sourceLine )
{
	if ( !sourceLine.isValid() )
	{
		m_stream  <<  "##Failure Location unknown## : ";
		return;
	}
	
	std::string location;
	for ( unsigned int index = 0; index < m_locationFormat.length(); ++index )
	{
		char c = m_locationFormat[ index ];
		if ( c == '%'  &&  ( index+1 < m_locationFormat.length() ) )
		{
			char command = m_locationFormat[index+1];
			if ( processLocationFormatCommand( command, sourceLine ) )
			{
				++index;
				continue;
			}
		}
		
		m_stream  << c;
	}
}


bool CustomFormatter::processLocationFormatCommand( char command, 
												const SourceLine &sourceLine )
{
	switch ( command )
	{
		case 'p':
			m_stream  <<  sourceLine.fileName();
			return true;
		case 'l':
			m_stream  <<  sourceLine.lineNumber();
			return true;
		case 'f':
			m_stream  <<  extractBaseName( sourceLine.fileName() );
			return true;
	}
	
	return false;
}


std::string CustomFormatter::extractBaseName( const std::string &fileName ) const
{
	int indexLastDirectorySeparator = fileName.find_last_of( '/' );
	
	if ( indexLastDirectorySeparator < 0 )
		indexLastDirectorySeparator = fileName.find_last_of( '\\' );
	
	if ( indexLastDirectorySeparator < 0 )
		return fileName;
	
	return fileName.substr( indexLastDirectorySeparator +1 );
}


void CustomFormatter::printFailureType( TestFailure *failure )
{
	m_stream  <<  (failure->isError() ? "Error" : "Assertion");
}


void CustomFormatter::printFailedTestName( TestFailure *failure )
{
	m_stream << failure->failedTestName();
}


void CustomFormatter::printFailureMessage( TestFailure *failure )
{
	
	Exception *thrownException = failure->thrownException();
	m_stream  << thrownException->message().shortDescription()  <<  " ";
	
	std::string message = thrownException->message().details();
	//if ( m_wrapColumn > 0 )
	//	message = StringTools::wrap( message, m_wrapColumn );
	
	m_stream  <<  message;
}


void CustomFormatter::printStatistics()
{
	m_stream  <<  "Failures !!!\n";
	m_stream  <<  "Run: "  <<  m_result->runTests()  << "   "
	<<  "Failure total: "  <<  m_result->testFailuresTotal()  << "   "
	<<  "Failures: "  <<  m_result->testFailures()  << "   "
	<<  "Errors: "  <<  m_result->testErrors()
	<<  "\n";
}


void CustomFormatter::setWrapColumn( int wrapColumn )
{
	m_wrapColumn = wrapColumn;
}


void CustomFormatter::setNoWrap()
{
	m_wrapColumn = 0;
}


int CustomFormatter::wrapColumn() const
{
	return m_wrapColumn;
}

// Harness
int main(int argc, char* argv[])
{
	// Get the top level suite from the registry
	CppUnit::Test *suite = CppUnit::TestFactoryRegistry::getRegistry().makeTest();
	
	// Adds the test to the list of test to run
	CppUnit::TextUi::TestRunner runner;
	runner.addTest( suite );
	
	CustomFormatter* outputter = new CustomFormatter( &runner.result(), std::cerr );
	//outputter->setLocationFormat("%p:%l: error: ");
	//outputter->setNoWrap();
	// Change the default outputter to a compiler error format outputter
	runner.setOutputter( outputter );
	// Run the tests.
	bool wasSucessful = runner.run("",false);
	
	// Return error code 1 if the one of test failed.
	return wasSucessful ? 0 : 1;
}

/* -*- C++ -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library
 *
 * Copyright (C) 2003-2008
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#ifndef LEMON_ERROR_H
#define LEMON_ERROR_H

//! \ingroup exceptions
//! \file
//! \brief Basic exception classes and error handling.

#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <memory>

namespace lemon {

  /// \addtogroup exceptions
  /// @{
  
  /// \brief Exception safe wrapper class.
  ///
  /// Exception safe wrapper class to implement the members of exceptions.
  template <typename _Type>
  class ExceptionMember {
  public:
    typedef _Type Type;

    ExceptionMember() throw () {
      try {
	ptr.reset(new Type());
      } catch (...) {}
    }

    ExceptionMember(const Type& type) throw () {
      try {
	ptr.reset(new Type());
	if (ptr.get() == 0) return;
	*ptr = type;
      } catch (...) {}
    }

    ExceptionMember(const ExceptionMember& copy) throw() {
      try {
	if (!copy.valid()) return;
	ptr.reset(new Type());
	if (ptr.get() == 0) return;
	*ptr = copy.get();
      } catch (...) {}
    }

    ExceptionMember& operator=(const ExceptionMember& copy) {
      if (ptr.get() == 0) return;
      try {
	if (!copy.valid()) return;
 	*ptr = copy.get();
      } catch (...) {}
    }

    void set(const Type& type) {
      if (ptr.get() == 0) return;
      try {
	*ptr = type;
      } catch (...) {}
    }

    const Type& get() const {
      return *ptr;
    }

    bool valid() const {
      return ptr.get() != 0;
    }
    
  private:
    std::auto_ptr<_Type> ptr;
  };

  /// Exception-safe convenient "error message" class.

  /// Helper class which provides a convenient ostream-like (operator <<
  /// based) interface to create a string message. Mostly useful in
  /// exception classes (therefore the name).
  class ErrorMessage {
  protected:
    ///\e 

    ///\todo The good solution is boost::shared_ptr...
    ///
    mutable
    std::auto_ptr<std::ostringstream> buf;
    
    ///\e 
    bool init() throw() {
      try {
	buf.reset(new std::ostringstream);
      }
      catch(...) {
	buf.reset();
      }
      return buf.get();
    }

  public:

    ///\e 
    ErrorMessage() throw() { init(); }

    ErrorMessage(const ErrorMessage& em) throw() : buf(em.buf) { }

    ///\e 
    ErrorMessage(const char *msg) throw() {
      init();
      *this << msg;
    }

    ///\e 
    ErrorMessage(const std::string &msg) throw() {
      init();
      *this << msg;
    }

    ///\e 
    template <typename T>
    ErrorMessage& operator<<(const T &t) throw() {
      if( ! buf.get() ) return *this;

      try {
	*buf << t;
      }
      catch(...) {
	buf.reset();
      }
      return *this;
    }

    ///\e 
    const char* message() throw() {
      if( ! buf.get() ) return 0;

      const char* mes = 0;
      try {
	mes = buf->str().c_str();
      }
      catch(...) {}
      return mes;
    }
    
  };

  /**
   * \brief Generic exception class.
   *
   * Base class for exceptions used in LEMON.
   */
  class Exception : public std::exception {
  public:
    ///\e 
    Exception() {}
    ///\e 
    virtual ~Exception() throw() {}
    ///\e 
    virtual const char* what() const throw() {
      return "lemon::Exception";
    }
  };

  /**
   * \brief One of the two main subclasses of \ref Exception.
   *
   * Logic errors represent problems in the internal logic of a program;
   * in theory, these are preventable, and even detectable before the
   * program runs (e.g., violations of class invariants).
   *
   * A typical example for this is \ref UninitializedParameter.
   */
  class LogicError : public Exception {
  public:
    virtual const char* what() const throw() {
      return "lemon::LogicError";
    }
  };

  /**
   * \brief \ref Exception for uninitialized parameters.
   *
   * This error represents problems in the initialization
   * of the parameters of the algorithms.
   */
  class UninitializedParameter : public LogicError {
  public:
    virtual const char* what() const throw() {
      return "lemon::UninitializedParameter";
    }
  };

  
  /**
   * \brief One of the two main subclasses of \ref Exception.
   *
   * Runtime errors represent problems outside the scope of a program;
   * they cannot be easily predicted and can generally only be caught as
   * the program executes.
   */
  class RuntimeError : public Exception {
  public:
    virtual const char* what() const throw() {
      return "lemon::RuntimeError";
    }
  };

  ///\e
  class RangeError : public RuntimeError {
  public:
    virtual const char* what() const throw() {
      return "lemon::RangeError";
    }
  };

  ///\e 
  class IoError : public RuntimeError {
  public:
    virtual const char* what() const throw() {
      return "lemon::IoError";
    }
  };

  ///\e 
  class DataFormatError : public IoError {
  protected:
    ExceptionMember<std::string> _message;
    ExceptionMember<std::string> _file;
    int _line;

    mutable ExceptionMember<std::string> _message_holder;
  public:

    DataFormatError(const DataFormatError &dfe) : 
      IoError(dfe), _message(dfe._message), _file(dfe._file),
      _line(dfe._line) {}

    ///\e 
    explicit DataFormatError(const char *the_message)
      : _message(the_message), _line(0) {}

    ///\e 
    DataFormatError(const std::string &file_name, int line_num,
		    const char *the_message)
      : _message(the_message), _line(line_num) { file(file_name); }

    ///\e 
    void line(int ln) { _line = ln; }
    ///\e 
    void message(const std::string& msg) { _message.set(msg); }
    ///\e 
    void file(const std::string &fl) { _file.set(fl); }
 
    ///\e
    int line() const { return _line; }
    ///\e
    const char* message() const { 
      if (_message.valid() && !_message.get().empty()) {
	return _message.get().c_str();
      } else {
	return 0;
      }
    }

    /// \brief Returns the filename.
    ///
    /// Returns \e null if the filename was not specified.
    const char* file() const {
      if (_file.valid() && !_file.get().empty()) {
	return _file.get().c_str();
      } else {
	return 0;
      }
    }

    ///\e 
    virtual const char* what() const throw() {
      try {
	std::ostringstream ostr;
	ostr << "lemon:DataFormatError" << ": ";
	if (message()) ostr << message();
	if( file() || line() != 0 ) {
	  ostr << " (";
	  if( file() ) ostr << "in file '" << file() << "'";
	  if( file() && line() != 0 ) ostr << " ";
	  if( line() != 0 ) ostr << "at line " << line();
	  ostr << ")";
	}
	_message_holder.set(ostr.str());
      }
      catch (...) {}
      if( _message_holder.valid()) return _message_holder.get().c_str();
      return "lemon:DataFormatError";
    }

    virtual ~DataFormatError() throw() {}
  };

  ///\e 
  class FileOpenError : public IoError {
  protected:
    ExceptionMember<std::string> _file;

    mutable ExceptionMember<std::string> _message_holder;
  public:

    FileOpenError(const FileOpenError &foe) : 
      IoError(foe), _file(foe._file) {}

    ///\e 
    explicit FileOpenError(const std::string& fl)
      : _file(fl) {}


    ///\e 
    void file(const std::string &fl) { _file.set(fl); }
 
    /// \brief Returns the filename.
    ///
    /// Returns \e null if the filename was not specified.
    const char* file() const {
      if (_file.valid() && !_file.get().empty()) {
	return _file.get().c_str();
      } else {
	return 0;
      }
    }

    ///\e 
    virtual const char* what() const throw() {
      try {
	std::ostringstream ostr;
	ostr << "lemon::FileOpenError" << ": ";
	ostr << "Cannot open file - " << file();
	_message_holder.set(ostr.str());
      }
      catch (...) {}
      if( _message_holder.valid()) return _message_holder.get().c_str();
      return "lemon::FileOpenError";
    }
    virtual ~FileOpenError() throw() {}
  };

  class IoParameterError : public IoError {
  protected:
    ExceptionMember<std::string> _message;
    ExceptionMember<std::string> _file;

    mutable ExceptionMember<std::string> _message_holder;
  public:

    IoParameterError(const IoParameterError &ile) : 
      IoError(ile), _message(ile._message), _file(ile._file) {}

    ///\e 
    explicit IoParameterError(const char *the_message)
      : _message(the_message) {}

    ///\e 
    IoParameterError(const char *file_name, const char *the_message)
      : _message(the_message), _file(file_name) {}

     ///\e 
    void message(const std::string& msg) { _message.set(msg); }
    ///\e 
    void file(const std::string &fl) { _file.set(fl); }
 
     ///\e
    const char* message() const { 
      if (_message.valid()) {
	return _message.get().c_str();
      } else {
	return 0;
      }
    }

    /// \brief Returns the filename.
    ///
    /// Returns \e null if the filename was not specified.
    const char* file() const {
      if (_file.valid()) {
	return _file.get().c_str();
      } else {
	return 0;
      }
    }

    ///\e 
    virtual const char* what() const throw() {
      try {
	std::ostringstream ostr;
	if (message()) ostr << message();
	if (file()) ostr << "(when reading file '" << file() << "')";
	_message_holder.set(ostr.str());
      }
      catch (...) {}
      if( _message_holder.valid() ) return _message_holder.get().c_str();
      return "lemon:IoParameterError";
    }
    virtual ~IoParameterError() throw() {}
  };


  ///\e
  class AssertionFailedError : public LogicError {
  protected:
    const char *assertion;
    const char *file;
    int line;
    const char *function;
    const char *message;

    mutable ExceptionMember<std::string> _message_holder;
  public:
    ///\e
    AssertionFailedError(const char *_file, int _line, const char *func,
			 const char *msg, const char *_assertion = 0) :
      assertion(_assertion), file(_file), line(_line), function(func),
      message(msg) {}

    ///\e
    const char* get_assertion() const { return assertion; }
    ///\e
    const char* get_message() const { return message; }
    ///\e
    const char* get_file() const { return file; }
    ///\e
    const char* get_function() const { return function; }
    ///\e
    int get_line() const { return line; }


    virtual const char* what() const throw() {
      try {
	std::ostringstream ostr;
	ostr << file << ":" << line << ": ";
	if( function )
	  ostr << function << ": ";
	ostr << message;
	if( assertion )
	   ostr << " (assertion '" << assertion << "' failed)";
	_message_holder.set(ostr.str());
	return ostr.str().c_str();
      }
      catch(...) {}
      if( _message_holder.valid() ) return _message_holder.get().c_str();
      return "lemon::AssertionFailedError";
    }
   virtual ~AssertionFailedError() throw() {}
  };


  /****************  Macros  ****************/


  template <typename Exception>
  inline void assert_fail(const char *file, int line, 
                          const char *func,
                          Exception exception, 
                          const char *assertion = 0,
                          bool do_abort=true)
  {
    using namespace std;
    cerr << file << ":" << line << ": ";
    if( func )
      cerr << func << ": ";
    cerr << exception.what();
    if( assertion )
      cerr << " (assertion '" << assertion << "' failed)";
    cerr << endl;
    if(do_abort)
      abort();
  }

  template <>
  inline void assert_fail<const char *>(const char *file, int line, 
                                        const char *func,
                                        const char *message, 
                                        const char *assertion,
                                        bool do_abort)
  {
    using namespace std;
    cerr << file << ":" << line << ": ";
    if( func )
      cerr << func << ": ";
    cerr << message;
    if( assertion )
      cerr << " (assertion '" << assertion << "' failed)";
    cerr << endl;
    if(do_abort)
      abort();
  }

  template <>
  inline void assert_fail<std::string>(const char *file, int line, 
                                       const char *func,
                                       std::string message, 
                                       const char *assertion,
                                       bool do_abort)
  {
    assert_fail(file, line, func, message.c_str(), assertion, do_abort);
  }

  template <typename Exception>
  inline void assert_fail_failure(const char *file, int line, const char *func,
			   Exception exception, 
			   const char *assertion = 0,
			   bool = true)
  {
    throw AssertionFailedError(file, line, func, exception.what(), assertion);
  }

  template <>
  inline void assert_fail_failure<const char *>(const char *file, int line, 
                                                const char *func,
                                                const char *message, 
                                                const char *assertion,
                                                bool)
  {
    throw AssertionFailedError(file, line, func, message, assertion);
  }

  template <>
  inline void assert_fail_failure<std::string>(const char *file, int line, 
                                               const char *func,
                                               std::string message, 
                                               const char *assertion,
                                               bool)
  {
    assert_fail_failure(file, line, func, message.c_str(), assertion, true);
  }

  template <typename Exception> 
  inline void assert_fail_exception(const char *file, int line, const char *func,
			     Exception exception, 
			     const char *assertion = 0, bool = true)
  {
    throw exception;
  }

  template <> 
  inline void assert_fail_exception<const char *>(const char *file, int line, 
					   const char *func,
					   const char *message, 
					   const char *assertion,
					   bool)
  {
    throw AssertionFailedError(file, line, func, message, assertion);
  }

  template <>
  inline void assert_fail_exception<std::string>(const char *file, int line, 
                                                 const char *func,
                                                 std::string message, 
                                                 const char *assertion,
                                                 bool)
  {
    assert_fail_exception(file, line, func, message.c_str(), assertion, true);    
  }

/// @}

}
#endif // LEMON_ERROR_H

#undef LEMON_ASSERT
#undef LEMON_FIXME

#ifdef LEMON_ENABLE_ASSERTS
#  define LEMON_ASSERT_ABORT
#endif

#ifndef LEMON_ASSERT_DO_ABORT
#  define LEMON_ASSERT_DO_ABORT 1
#endif

#ifndef LEMON_ASSERT_HANDLER
#  if defined LEMON_ASSERT_EXCEPTION
#    define LEMON_ASSERT_HANDLER ::lemon::assert_fail_exception
#  elif defined LEMON_ASSERT_FAILURE
#    define LEMON_ASSERT_HANDLER ::lemon::assert_fail_failure
#  elif defined LEMON_ASSERT_ABORT
#    define LEMON_ASSERT_HANDLER ::lemon::assert_fail
#  else
#    define LEMON_DISABLE_ASSERTS
#  endif
#endif

#ifdef DOXYGEN

/// \brief Macro for assertions with customizable message
///
/// Macro for assertions with customizable message.
///
/// The assertions are disabled in the default behaviour. You can
/// enable the assertions with the
/// \code
/// #define LEMON_ENABLE_ASSERTS
/// \endcode
/// Then an assert
/// provides a log on the standard error about the assertion and aborts
/// the program if LEMON_ASSERT_DO_ABORT is also defined (otherwise the
/// program keeps on running).
/// By defining LEMON_ASSERT_FAILURE or
/// LEMON_ASSERT_EXCEPTION, you can set other behaviour to the
/// assertions. In case LEMON_ASSERT_FAILURE is given, LEMON_ASSERT
/// will always throw an \c AssertionFailedError exception with
/// the \c msg error message. By using
/// LEMON_ASSERT_EXCEPTION, one can define an arbitrary exception to be thrown.
///
/// The LEMON_ASSERT macro should be called with the \c exp parameter
/// which should be an expression convertible to bool. If the given
/// parameter is false the assertion is raised and one of the assertion
/// behaviour will be activated. The \c msg should be either a const
/// char* message or an exception. When the \c msg is an exception the
/// \ref lemon::Exception::what() "what()" function is called to retrieve and
/// display the error message.
///
/// \todo We should provide some way to reset to the default behaviour,
/// shouldn't we?
///
/// \todo This whole 'assert' business should be placed in a separate
/// include file. The boost assert is not guarded by header sentries
/// which may help to change the behaviour of the assertions in 
/// the files.
///
/// \todo __PRETTY_FUNCTION__ should be replaced by something
/// compiler-independent, like BOOST_CURRENT_FUNCTION

#  define LEMON_ASSERT(exp, msg)                 \
     (static_cast<void> (!!(exp) ? 0 : (         \
       LEMON_ASSERT_HANDLER(__FILE__, __LINE__,  \
                            __PRETTY_FUNCTION__, \
                            msg, #exp, LEMON_ASSERT_DO_ABORT), 0)))

#else 
#  if defined LEMON_DISABLE_ASSERTS

#    define LEMON_ASSERT(exp, msg)  (static_cast<void> (0))

#  else
#    define LEMON_ASSERT(exp, msg)                 \
       (static_cast<void> (!!(exp) ? 0 : (         \
         LEMON_ASSERT_HANDLER(__FILE__, __LINE__,  \
                              __PRETTY_FUNCTION__, \
                              msg, #exp, LEMON_ASSERT_DO_ABORT), 0)))
#  endif
#endif

/**
 * \brief Macro for mark not yet implemented features.
 *
 * \todo Is this the right place for this? It should be used only in
 * modules under development.
 *
 * \todo __PRETTY_FUNCTION__ should be replaced by something
 * compiler-independent, like BOOST_CURRENT_FUNCTION
 */

# define LEMON_FIXME(msg) \
    (LEMON_ASSERT_HANDLER(__FILE__, __LINE__, __PRETTY_FUNCTION__, \
			  "FIXME: " msg))

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

#ifndef LEMON_TOLERANCE_H
#define LEMON_TOLERANCE_H

///\ingroup misc
///\file
///\brief A basic tool to handle the anomalies of calculation with
///floating point numbers.
///
///\todo It should be in a module like "Basic tools"


namespace lemon {

  /// \addtogroup misc
  /// @{
  
  ///\brief A class to provide a basic way to
  ///handle the comparison of numbers that are obtained
  ///as a result of a probably inexact computation.
  ///
  ///Tolerance is a class to provide a basic way to
  ///handle the comparison of numbers that are obtained
  ///as a result of a probably inexact computation.
  ///
  ///This is an abstract class, it should be specialized for all numerical
  ///data types. These specialized classes like \ref Tolerance\<double\>
  ///may offer additional tuning parameters.
  ///
  ///\sa Tolerance<float>
  ///\sa Tolerance<double>
  ///\sa Tolerance<long double>
  ///\sa Tolerance<int>
  ///\sa Tolerance<long long int>
  ///\sa Tolerance<unsigned int>
  ///\sa Tolerance<unsigned long long int>

  template<class T>
  class Tolerance
  {
  public:
    typedef T Value;

    ///\name Comparisons
    ///The concept is that these bool functions return with \c true only if
    ///the related comparisons hold even if some numerical error appeared
    ///during the computations.

    ///@{

    ///Returns \c true if \c a is \e surely strictly less than \c b
    static bool less(Value a,Value b) {return false;}
    ///Returns \c true if \c a is \e surely different from \c b
    static bool different(Value a,Value b) {return false;}
    ///Returns \c true if \c a is \e surely positive
    static bool positive(Value a) {return false;}
    ///Returns \c true if \c a is \e surely negative
    static bool negative(Value a) {return false;}
    ///Returns \c true if \c a is \e surely non-zero
    static bool nonZero(Value a) {return false;}

    ///@}

    ///Returns the zero value.
    static Value zero() {return T();}

    //   static bool finite(Value a) {}
    //   static Value big() {}
    //   static Value negativeBig() {}
  };


  ///Float specialization of \ref Tolerance.

  ///Float specialization of \ref Tolerance.
  ///\sa Tolerance
  ///\relates Tolerance
  template<>
  class Tolerance<float>
  {
    static float def_epsilon;
    float _epsilon;
  public:
    ///\e
    typedef float Value;

    ///Constructor setting the epsilon tolerance to the default value.
    Tolerance() : _epsilon(def_epsilon) {}
    ///Constructor setting the epsilon tolerance.
    Tolerance(float e) : _epsilon(e) {}

    ///Return the epsilon value.
    Value epsilon() const {return _epsilon;}
    ///Set the epsilon value.
    void epsilon(Value e) {_epsilon=e;}

    ///Return the default epsilon value.
    static Value defaultEpsilon() {return def_epsilon;}
    ///Set the default epsilon value.
    static void defaultEpsilon(Value e) {def_epsilon=e;}

    ///\name Comparisons
    ///See class Tolerance for more details.

    ///@{

    ///Returns \c true if \c a is \e surely strictly less than \c b
    bool less(Value a,Value b) const {return a+_epsilon<b;}
    ///Returns \c true if \c a is \e surely different from \c b
    bool different(Value a,Value b) const { return less(a,b)||less(b,a); }
    ///Returns \c true if \c a is \e surely positive
    bool positive(Value a) const { return _epsilon<a; }
    ///Returns \c true if \c a is \e surely negative
    bool negative(Value a) const { return -_epsilon>a; }
    ///Returns \c true if \c a is \e surely non-zero
    bool nonZero(Value a) const { return positive(a)||negative(a); };

    ///@}

    ///Returns zero
    static Value zero() {return 0;}
  };

  ///Double specialization of \ref Tolerance.

  ///Double specialization of \ref Tolerance.
  ///\sa Tolerance
  ///\relates Tolerance
  template<>
  class Tolerance<double>
  {
    static double def_epsilon;
    double _epsilon;
  public:
    ///\e
    typedef double Value;

    ///Constructor setting the epsilon tolerance to the default value.
    Tolerance() : _epsilon(def_epsilon) {}
    ///Constructor setting the epsilon tolerance.
    Tolerance(double e) : _epsilon(e) {}

    ///Return the epsilon value.
    Value epsilon() const {return _epsilon;}
    ///Set the epsilon value.
    void epsilon(Value e) {_epsilon=e;}

    ///Return the default epsilon value.
    static Value defaultEpsilon() {return def_epsilon;}
    ///Set the default epsilon value.
    static void defaultEpsilon(Value e) {def_epsilon=e;}

    ///\name Comparisons
    ///See class Tolerance for more details.

    ///@{

    ///Returns \c true if \c a is \e surely strictly less than \c b
    bool less(Value a,Value b) const {return a+_epsilon<b;}
    ///Returns \c true if \c a is \e surely different from \c b
    bool different(Value a,Value b) const { return less(a,b)||less(b,a); }
    ///Returns \c true if \c a is \e surely positive
    bool positive(Value a) const { return _epsilon<a; }
    ///Returns \c true if \c a is \e surely negative
    bool negative(Value a) const { return -_epsilon>a; }
    ///Returns \c true if \c a is \e surely non-zero
    bool nonZero(Value a) const { return positive(a)||negative(a); };

    ///@}

    ///Returns zero
    static Value zero() {return 0;}
  };

  ///Long double specialization of \ref Tolerance.

  ///Long double specialization of \ref Tolerance.
  ///\sa Tolerance
  ///\relates Tolerance
  template<>
  class Tolerance<long double>
  {
    static long double def_epsilon;
    long double _epsilon;
  public:
    ///\e
    typedef long double Value;

    ///Constructor setting the epsilon tolerance to the default value.
    Tolerance() : _epsilon(def_epsilon) {}
    ///Constructor setting the epsilon tolerance.
    Tolerance(long double e) : _epsilon(e) {}

    ///Return the epsilon value.
    Value epsilon() const {return _epsilon;}
    ///Set the epsilon value.
    void epsilon(Value e) {_epsilon=e;}

    ///Return the default epsilon value.
    static Value defaultEpsilon() {return def_epsilon;}
    ///Set the default epsilon value.
    static void defaultEpsilon(Value e) {def_epsilon=e;}

    ///\name Comparisons
    ///See class Tolerance for more details.

    ///@{

    ///Returns \c true if \c a is \e surely strictly less than \c b
    bool less(Value a,Value b) const {return a+_epsilon<b;}
    ///Returns \c true if \c a is \e surely different from \c b
    bool different(Value a,Value b) const { return less(a,b)||less(b,a); }
    ///Returns \c true if \c a is \e surely positive
    bool positive(Value a) const { return _epsilon<a; }
    ///Returns \c true if \c a is \e surely negative
    bool negative(Value a) const { return -_epsilon>a; }
    ///Returns \c true if \c a is \e surely non-zero
    bool nonZero(Value a) const { return positive(a)||negative(a); };

    ///@}

    ///Returns zero
    static Value zero() {return 0;}
  };

  ///Integer specialization of \ref Tolerance.

  ///Integer specialization of \ref Tolerance.
  ///\sa Tolerance
  template<>
  class Tolerance<int>
  {
  public:
    ///\e
    typedef int Value;

    ///\name Comparisons
    ///See \ref Tolerance for more details.

    ///@{

    ///Returns \c true if \c a is \e surely strictly less than \c b
    static bool less(Value a,Value b) { return a<b;}
    ///Returns \c true if \c a is \e surely different from \c b
    static bool different(Value a,Value b) { return a!=b; }
    ///Returns \c true if \c a is \e surely positive
    static bool positive(Value a) { return 0<a; }
    ///Returns \c true if \c a is \e surely negative
    static bool negative(Value a) { return 0>a; }
    ///Returns \c true if \c a is \e surely non-zero
    static bool nonZero(Value a) { return a!=0; };

    ///@}

    ///Returns zero
    static Value zero() {return 0;}
  };

  ///Unsigned integer specialization of \ref Tolerance.

  ///Unsigned integer specialization of \ref Tolerance.
  ///\sa Tolerance
  template<>
  class Tolerance<unsigned int>
  {
  public:
    ///\e
    typedef unsigned int Value;

    ///\name Comparisons
    ///See \ref Tolerance for more details.

    ///@{

    ///Returns \c true if \c a is \e surely strictly less than \c b
    static bool less(Value a,Value b) { return a<b;}
    ///Returns \c true if \c a is \e surely different from \c b
    static bool different(Value a,Value b) { return a!=b; }
    ///Returns \c true if \c a is \e surely positive
    static bool positive(Value a) { return 0<a; }
    ///Returns \c true if \c a is \e surely negative
    static bool negative(Value) { return false; }
    ///Returns \c true if \c a is \e surely non-zero
    static bool nonZero(Value a) { return a!=0; };

    ///@}

    ///Returns zero
    static Value zero() {return 0;}
  };
  

  ///Long integer specialization of \ref Tolerance.

  ///Long integer specialization of \ref Tolerance.
  ///\sa Tolerance
  template<>
  class Tolerance<long int>
  {
  public:
    ///\e
    typedef long int Value;

    ///\name Comparisons
    ///See \ref Tolerance for more details.

    ///@{

    ///Returns \c true if \c a is \e surely strictly less than \c b
    static bool less(Value a,Value b) { return a<b;}
    ///Returns \c true if \c a is \e surely different from \c b
    static bool different(Value a,Value b) { return a!=b; }
    ///Returns \c true if \c a is \e surely positive
    static bool positive(Value a) { return 0<a; }
    ///Returns \c true if \c a is \e surely negative
    static bool negative(Value a) { return 0>a; }
    ///Returns \c true if \c a is \e surely non-zero
    static bool nonZero(Value a) { return a!=0;};

    ///@}

    ///Returns zero
    static Value zero() {return 0;}
  };

  ///Unsigned long integer specialization of \ref Tolerance.

  ///Unsigned long integer specialization of \ref Tolerance.
  ///\sa Tolerance
  template<>
  class Tolerance<unsigned long int>
  {
  public:
    ///\e
    typedef unsigned long int Value;

    ///\name Comparisons
    ///See \ref Tolerance for more details.

    ///@{

    ///Returns \c true if \c a is \e surely strictly less than \c b
    static bool less(Value a,Value b) { return a<b;}
    ///Returns \c true if \c a is \e surely different from \c b
    static bool different(Value a,Value b) { return a!=b; }
    ///Returns \c true if \c a is \e surely positive
    static bool positive(Value a) { return 0<a; }
    ///Returns \c true if \c a is \e surely negative
    static bool negative(Value) { return false; }
    ///Returns \c true if \c a is \e surely non-zero
    static bool nonZero(Value a) { return a!=0;};

    ///@}

    ///Returns zero
    static Value zero() {return 0;}
  };

#if defined __GNUC__ && !defined __STRICT_ANSI__

  ///Long long integer specialization of \ref Tolerance.

  ///Long long integer specialization of \ref Tolerance.
  ///\warning This class (more exactly, type <tt>long long</tt>)
  ///is not ansi compatible.
  ///\sa Tolerance
  template<>
  class Tolerance<long long int>
  {
  public:
    ///\e
    typedef long long int Value;

    ///\name Comparisons
    ///See \ref Tolerance for more details.

    ///@{

    ///Returns \c true if \c a is \e surely strictly less than \c b
    static bool less(Value a,Value b) { return a<b;}
    ///Returns \c true if \c a is \e surely different from \c b
    static bool different(Value a,Value b) { return a!=b; }
    ///Returns \c true if \c a is \e surely positive
    static bool positive(Value a) { return 0<a; }
    ///Returns \c true if \c a is \e surely negative
    static bool negative(Value a) { return 0>a; }
    ///Returns \c true if \c a is \e surely non-zero
    static bool nonZero(Value a) { return a!=0;};

    ///@}

    ///Returns zero
    static Value zero() {return 0;}
  };

  ///Unsigned long long integer specialization of \ref Tolerance.

  ///Unsigned long long integer specialization of \ref Tolerance.
  ///\warning This class (more exactly, type <tt>unsigned long long</tt>)
  ///is not ansi compatible.
  ///\sa Tolerance
  template<>
  class Tolerance<unsigned long long int>
  {
  public:
    ///\e
    typedef unsigned long long int Value;

    ///\name Comparisons
    ///See \ref Tolerance for more details.

    ///@{

    ///Returns \c true if \c a is \e surely strictly less than \c b
    static bool less(Value a,Value b) { return a<b;}
    ///Returns \c true if \c a is \e surely different from \c b
    static bool different(Value a,Value b) { return a!=b; }
    ///Returns \c true if \c a is \e surely positive
    static bool positive(Value a) { return 0<a; }
    ///Returns \c true if \c a is \e surely negative
    static bool negative(Value) { return false; }
    ///Returns \c true if \c a is \e surely non-zero
    static bool nonZero(Value a) { return a!=0;};

    ///@}

    ///Returns zero
    static Value zero() {return 0;}
  };

#endif

  /// @}

} //namespace lemon

#endif //LEMON_TOLERANCE_H

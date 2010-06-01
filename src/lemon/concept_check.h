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

// Modified for use in LEMON.
// We should really consider using Boost...

//
// (C) Copyright Jeremy Siek 2000.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// Revision History:
//   05 May   2001: Workarounds for HP aCC from Thomas Matelich. (Jeremy Siek)
//   02 April 2001: Removed limits header altogether. (Jeremy Siek)
//   01 April 2001: Modified to use new <boost/limits.hpp> header. (JMaddock)
//

// See http://www.boost.org/libs/concept_check for documentation.

#ifndef LEMON_BOOST_CONCEPT_CHECKS_HPP
#define LEMON_BOOST_CONCEPT_CHECKS_HPP

namespace lemon {

  /*
    "inline" is used for ignore_unused_variable_warning()
    and function_requires() to make sure there is no
    overtarget with g++.
  */

  template <class T> inline void ignore_unused_variable_warning(const T&) { }

  template <class Concept>
  inline void function_requires()
  {
#if !defined(NDEBUG)
    void (Concept::*x)() = & Concept::constraints;
    ignore_unused_variable_warning(x);
#endif
  }

  template <typename Concept, typename Type>
  inline void checkConcept() {
#if !defined(NDEBUG)
    typedef typename Concept::template Constraints<Type> ConceptCheck;
    void (ConceptCheck::*x)() = & ConceptCheck::constraints;
    ignore_unused_variable_warning(x);
#endif
  }
#if 0
#define BOOST_CLASS_REQUIRE(type_var, ns, concept) \
  typedef void (ns::concept <type_var>::* func##type_var##concept)(); \
  template <func##type_var##concept Tp1_> \
  struct concept_checking_##type_var##concept { }; \
  typedef concept_checking_##type_var##concept< \
    BOOST_FPTR ns::concept<type_var>::constraints> \
    concept_checking_typedef_##type_var##concept

#define BOOST_CLASS_REQUIRE2(type_var1, type_var2, ns, concept) \
  typedef void (ns::concept <type_var1,type_var2>::* \
     func##type_var1##type_var2##concept)(); \
  template <func##type_var1##type_var2##concept Tp1_> \
  struct concept_checking_##type_var1##type_var2##concept { }; \
  typedef concept_checking_##type_var1##type_var2##concept< \
    BOOST_FPTR ns::concept<type_var1,type_var2>::constraints> \
    concept_checking_typedef_##type_var1##type_var2##concept

#define BOOST_CLASS_REQUIRE3(tv1, tv2, tv3, ns, concept) \
  typedef void (ns::concept <tv1,tv2,tv3>::* \
     func##tv1##tv2##tv3##concept)(); \
  template <func##tv1##tv2##tv3##concept Tp1_> \
  struct concept_checking_##tv1##tv2##tv3##concept { }; \
  typedef concept_checking_##tv1##tv2##tv3##concept< \
    BOOST_FPTR ns::concept<tv1,tv2,tv3>::constraints> \
    concept_checking_typedef_##tv1##tv2##tv3##concept

#define BOOST_CLASS_REQUIRE4(tv1, tv2, tv3, tv4, ns, concept) \
  typedef void (ns::concept <tv1,tv2,tv3,tv4>::* \
     func##tv1##tv2##tv3##tv4##concept)(); \
  template <func##tv1##tv2##tv3##tv4##concept Tp1_> \
  struct concept_checking_##tv1##tv2##tv3##tv4##concept { }; \
  typedef concept_checking_##tv1##tv2##tv3##tv4##concept< \
    BOOST_FPTR ns::concept<tv1,tv2,tv3,tv4>::constraints> \
    concept_checking_typedef_##tv1##tv2##tv3##tv4##concept
#endif

} // namespace lemon

#endif // LEMON_BOOST_CONCEPT_CHECKS_HPP

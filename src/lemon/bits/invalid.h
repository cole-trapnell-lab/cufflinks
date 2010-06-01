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

#ifndef LEMON_BITS_INVALID_H
#define LEMON_BITS_INVALID_H

///\file
///\brief Definition of INVALID.

namespace lemon {

  /// \brief Dummy type to make it easier to make invalid iterators.
  ///
  /// See \ref INVALID for the usage.
  struct Invalid {
  public:
    bool operator==(Invalid) { return true;  }
    bool operator!=(Invalid) { return false; }
    bool operator< (Invalid) { return false; }
  };
  
  /// Invalid iterators.
  
  /// \ref Invalid is a global type that converts to each iterator
  /// in such a way that the value of the target iterator will be invalid.

  //Some people didn't like this:
  //const Invalid &INVALID = *(Invalid *)0;

#ifdef LEMON_ONLY_TEMPLATES
  const Invalid INVALID = Invalid();
#else
  extern const Invalid INVALID;
#endif

} //namespace lemon

#endif
  

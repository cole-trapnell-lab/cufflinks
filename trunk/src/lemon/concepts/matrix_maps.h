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

#ifndef LEMON_CONCEPT_MATRIX_MAPS_H
#define LEMON_CONCEPT_MATRIX_MAPS_H

#include <lemon/bits/utility.h>
#include <lemon/concept_check.h>

///\ingroup concept
///\file
///\brief MatrixMap concepts checking classes for testing and documenting.

namespace lemon {

  namespace concepts {
  
    /// \addtogroup concept
    /// @{

    /// Readable matrix map concept
    template <typename K1, typename K2, typename V>
    class ReadMatrixMap
    {
    public:
      /// Map's first key type.
      typedef K1 FirstKey;    
      /// Map's second key type.
      typedef K2 SecondKey;    
      /// \brief Map's value type. 
      /// (The type of objects associated with the pairs of keys).
      typedef V Value;

      // \bug Value don't need to be default constructible.
      /// Returns the value associated with a key.
      Value operator()(const FirstKey&, const SecondKey&) const {
	return Value();
      }

      template <typename _ReadMatrixMap>
      struct Constraints {

	void constraints() {
	  Value val = m(first_key, second_key);
	  val = m(first_key, second_key);
	  typename _ReadMatrixMap::Value own_val = 
	    m(own_first_key, own_second_key); 
	  own_val = m(own_first_key, own_second_key);
	  ignore_unused_variable_warning(val);
	  ignore_unused_variable_warning(own_val);
	}

	FirstKey& first_key;
	SecondKey& second_key;	
	typename _ReadMatrixMap::FirstKey& own_first_key;
	typename _ReadMatrixMap::SecondKey& own_second_key;
	_ReadMatrixMap& m;
      };
      
    };


    /// Writable map concept
    template <typename K1, typename K2, typename V>
    class WriteMatrixMap {
    public:
      /// Map's first key type.
      typedef K1 FirstKey;    
      /// Map's second key type.
      typedef K2 SecondKey;    
      /// \brief Map's value type. 
      /// (The type of objects associated with the pairs of keys).
      typedef V Value;

      /// Sets the value associated with the pair of keys.
      void set(const FirstKey&, const SecondKey& ,const Value&) {}

      template <typename _WriteMatrixMap>
      struct Constraints {
	void constraints() {
	  // No constraints for constructor.
	  m.set(first_key, second_key, val);
	  m.set(own_first_key, own_second_key, own_val);
	}

	Value& val;
	typename _WriteMatrixMap::Value own_val;
	FirstKey& first_key;
	SecondKey& second_key;
	typename _WriteMatrixMap::FirstKey& own_first_key;
	typename _WriteMatrixMap::SecondKey& own_second_key;
	_WriteMatrixMap& m;

      };
    };

    ///Read/Writable map concept
    template<typename K1, typename K2, typename V>
    class ReadWriteMatrixMap 
      : public ReadMatrixMap<K1, K2, V>, public WriteMatrixMap<K1, K2, V> {
    public:
      /// Map's first key type.
      typedef K1 FirstKey;    
      /// Map's second key type.
      typedef K2 SecondKey;    
      /// \brief Map's value type. 
      /// (The type of objects associated with the pairs of keys).
      typedef V Value;

      /// Returns the value associated with a pair of keys.
      Value operator()(const FirstKey&, const SecondKey&) const { 
	return Value(); 
      }
      /// Sets the value associated with the pair of keys.
      void set(const FirstKey&, const SecondKey& ,const Value&) {}

      template<typename _ReadWriteMatrixMap>
      struct Constraints {
	void constraints() {
	  checkConcept<ReadMatrixMap<K1, K2, V>, _ReadWriteMatrixMap >();
	  checkConcept<WriteMatrixMap<K1, K2, V>, _ReadWriteMatrixMap >();
	}
      };
    };
  
  
    ///Dereferable matrix map concept
    template<typename K1, typename K2, typename V, typename R, typename CR>
    class ReferenceMatrixMap : public ReadWriteMatrixMap<K1, K2, V>
    {
    public:
      /// Tag for reference maps.
      typedef True ReferenceMapTag;
      /// Map's first key type.
      typedef K1 FirstKey;    
      /// Map's second key type.
      typedef K1 SecondKey;    
      /// Map's value type. (The type of objects associated with the keys).
      typedef V Value;
      /// Map's reference type.
      typedef R Reference;
      /// Map's const reference type.
      typedef CR ConstReference;

    protected:
      Value tmp;
    public:

      ///Returns a reference to the value associated to a pair of keys.
      Reference operator()(const FirstKey&, const SecondKey&) { 
	return tmp; 
      }
      ///Returns a const reference to the value associated to a pair of keys.
      ConstReference operator()(const FirstKey&, const SecondKey&) const { 
	return tmp; 
      }
      /// Sets the value associated with the pair of keys.
      void set(const FirstKey&, const SecondKey& ,const Value&) {}

      // \todo rethink this concept
      template<typename _ReferenceMatrixMap>
      struct ReferenceMapConcept {

	void constraints() {
	  checkConcept<ReadWriteMatrixMap, _ReferenceMatrixMap >();
	  m(first_key, second_key) = val;
	  val  = m(first_key, second_key);
	  m(first_key, second_key) = ref;
	  ref = m(first_key, second_key);
	  m(own_first_key, own_second_key) = own_val;
	  own_val  = m(own_first_key, own_second_key);
	  m(own_first_key, own_second_key) = own_ref;
	  own_ref = m(own_first_key, own_second_key); 
	}

	typename _ReferenceMatrixMap::Key& own_first_key;
	typename _ReferenceMatrixMap::Key& own_second_key;
	typename _ReferenceMatrixMap::Value& own_val;
	typename _ReferenceMatrixMap::Reference& own_ref;
	FirstKey& first_key;
	SecondKey& second_key;
	Value& val;
	Reference& ref;
	_ReferenceMatrixMap& m;
      };
    };

    // @}

  } //namespace concepts
} //namespace lemon
#endif // LEMON_CONCEPT_MATRIX_MAPS_H

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

#ifndef LEMON_CONCEPT_MAPS_H
#define LEMON_CONCEPT_MAPS_H

#include <lemon/bits/utility.h>
#include <lemon/concept_check.h>

///\ingroup concept
///\file
///\brief Map concepts checking classes for testing and documenting.

namespace lemon {
	
	namespace concepts {
		
		/// \addtogroup concept
		/// @{
		
		/// Readable map concept
		
		/// Readable map concept.
		///
		template<typename K, typename T>
		class ReadMap
			{
			public:
				/// The key type of the map.
				typedef K Key;    
				/// The value type of the map. (The type of objects associated with the keys).
				typedef T Value;
				
				/// Returns the value associated with a key.
				
				/// Returns the value associated with a key.
				/// \bug Value shouldn't need to be default constructible.
				///
				Value operator[](const Key &) const {return Value();}
				
				template<typename _ReadMap>
				struct Constraints {
					
					void constraints() {
						Value val = m[key];
						val = m[key];
						typename _ReadMap::Value own_val = m[own_key]; 
						own_val = m[own_key]; 
						
						ignore_unused_variable_warning(val);
						ignore_unused_variable_warning(own_val);
						ignore_unused_variable_warning(key);
					}
					Key& key;
					typename _ReadMap::Key& own_key;
					_ReadMap& m;
				};
				
			};
		
		
		/// Writable map concept
		
		/// Writable map concept.
		///
		template<typename K, typename T>
		class WriteMap
			{
			public:
				/// The key type of the map.
				typedef K Key;    
				/// The value type of the map. (The type of objects associated with the keys).
				typedef T Value;
				
				/// Sets the value associated with a key.
				void set(const Key &,const Value &) {}
				
				///Default constructor
				WriteMap() {}
				
				template <typename _WriteMap>
				struct Constraints {
					void constraints() {
						// No constraints for constructor.
						m.set(key, val);
						m.set(own_key, own_val);
						ignore_unused_variable_warning(key);
						ignore_unused_variable_warning(val);
						ignore_unused_variable_warning(own_key);
						ignore_unused_variable_warning(own_val);
					}
					
					Value& val;
					typename _WriteMap::Value own_val;
					Key& key;
					typename _WriteMap::Key& own_key;
					_WriteMap& m;
					
				};
			};
		
		/// Read/writable map concept
		
		/// Read/writable map concept.
		///
		template<typename K, typename T>
		class ReadWriteMap : public ReadMap<K,T>,
		public WriteMap<K,T>
		{
		public:
			/// The key type of the map.
			typedef K Key;    
			/// The value type of the map. (The type of objects associated with the keys).
			typedef T Value;
			
			/// Returns the value associated with a key.
			Value operator[](const Key &) const {return Value();}
			/// Sets the value associated with a key.
			void set(const Key & ,const Value &) {}
			
			template<typename _ReadWriteMap>
			struct Constraints {
				void constraints() {
					checkConcept<ReadMap<K, T>, _ReadWriteMap >();
					checkConcept<WriteMap<K, T>, _ReadWriteMap >();
				}
			};
		};
		
		
		///Dereferable map concept
		
		/// Dereferable map concept.
		///
		/// \todo Rethink this concept.
		template<typename K, typename T, typename R, typename CR>
		class ReferenceMap : public ReadWriteMap<K,T>
		{
		public:
			/// Tag for reference maps.
			typedef True ReferenceMapTag;
			/// The key type of the map.
			typedef K Key;    
			/// The value type of the map. (The type of objects associated with the keys).
			typedef T Value;
			/// The reference type of the map.
			typedef R Reference;
			/// The const reference type of the map.
			typedef CR ConstReference;
			
		protected:
			Value tmp;
		public:
			
			///Returns a reference to the value associated with a key.
			Reference operator[](const Key &) { return tmp; }
			///Returns a const reference to the value associated with a key.
			ConstReference operator[](const Key &) const { return tmp; }
			/// Sets the value associated with a key.
			void set(const Key &k,const Value &t) { operator[](k)=t; }
			
			template<typename _ReferenceMap>
			struct Constraints {
				
				void constraints() {
					checkConcept<ReadWriteMap<K, T>, _ReferenceMap >();
					m[key] = val;
					val  = m[key];
					m[key] = ref;
					ref = m[key];
					m[own_key] = own_val;
					own_val  = m[own_key];
					m[own_key] = own_ref;
					own_ref = m[own_key];	  	  
				}
				
				typename _ReferenceMap::Key& own_key;
				typename _ReferenceMap::Value& own_val;
				typename _ReferenceMap::Reference own_ref;
				Key& key;
				Value& val;
				Reference ref;
				_ReferenceMap& m;
			};
		};
		
		// @}
		
	} //namespace concepts
	
} //namespace lemon

#endif // LEMON_CONCEPT_MAPS_H

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

#ifndef LEMON_BITS_VECTOR_MAP_H
#define LEMON_BITS_VECTOR_MAP_H

#include <vector>
#include <algorithm>

#include <lemon/bits/traits.h>
#include <lemon/bits/utility.h>

#include <lemon/bits/alteration_notifier.h>

#include <lemon/concept_check.h>
#include <lemon/concepts/maps.h>

///\ingroup graphbits
///
///\file
///\brief Vector based graph maps.
namespace lemon {

  /// \ingroup graphbits
  ///
  /// \brief Graph map based on the std::vector storage.
  ///
  /// The VectorMap template class is graph map structure what
  /// automatically updates the map when a key is added to or erased from
  /// the map. This map type uses the std::vector to store the values.
  ///
  /// \param Notifier The AlterationNotifier that will notify this map.
  /// \param Item The item type of the graph items.
  /// \param Value The value type of the map.
  /// 
  /// \author Balazs Dezso  	
  template <typename _Graph, typename _Item, typename _Value>
  class VectorMap 
    : public ItemSetTraits<_Graph, _Item>::ItemNotifier::ObserverBase {
  private:
		
    /// The container type of the map.
    typedef std::vector<_Value> Container;	

  public:

    /// The graph type of the map. 
    typedef _Graph Graph;
    /// The item type of the map.
    typedef _Item Item;
    /// The reference map tag.
    typedef True ReferenceMapTag;

    /// The key type of the map.
    typedef _Item Key;
    /// The value type of the map.
    typedef _Value Value;

    /// The notifier type.
    typedef typename ItemSetTraits<_Graph, _Item>::ItemNotifier Notifier;

    /// The map type.
    typedef VectorMap Map;
    /// The base class of the map.
    typedef typename Notifier::ObserverBase Parent;

    /// The reference type of the map;
    typedef typename Container::reference Reference;
    /// The const reference type of the map;
    typedef typename Container::const_reference ConstReference;


    /// \brief Constructor to attach the new map into the notifier.
    ///
    /// It constructs a map and attachs it into the notifier.
    /// It adds all the items of the graph to the map.
    VectorMap(const Graph& graph) {
      Parent::attach(graph.notifier(Item()));
      container.resize(Parent::notifier()->maxId() + 1);
    }

    /// \brief Constructor uses given value to initialize the map. 
    ///
    /// It constructs a map uses a given value to initialize the map. 
    /// It adds all the items of the graph to the map.
    VectorMap(const Graph& graph, const Value& value) {
      Parent::attach(graph.notifier(Item()));
      container.resize(Parent::notifier()->maxId() + 1, value);
    }

    /// \brief Copy constructor
    ///
    /// Copy constructor.
    VectorMap(const VectorMap& _copy) : Parent() {
      if (_copy.attached()) {
	Parent::attach(*_copy.notifier());
	container = _copy.container;
      }
    }

    /// \brief Assign operator.
    ///
    /// This operator assigns for each item in the map the
    /// value mapped to the same item in the copied map.  
    /// The parameter map should be indiced with the same
    /// itemset because this assign operator does not change
    /// the container of the map. 
    VectorMap& operator=(const VectorMap& cmap) {
      return operator=<VectorMap>(cmap);
    }


    /// \brief Template assign operator.
    ///
    /// The given parameter should be conform to the ReadMap
    /// concecpt and could be indiced by the current item set of
    /// the NodeMap. In this case the value for each item
    /// is assigned by the value of the given ReadMap. 
    template <typename CMap>
    VectorMap& operator=(const CMap& cmap) {
      checkConcept<concepts::ReadMap<Key, _Value>, CMap>();
      const typename Parent::Notifier* nf = Parent::notifier();
      Item it;
      for (nf->first(it); it != INVALID; nf->next(it)) {
        set(it, cmap[it]);
      }
      return *this;
    }
    
  public:

    /// \brief The subcript operator.
    ///
    /// The subscript operator. The map can be subscripted by the
    /// actual items of the graph.      
    Reference operator[](const Key& key) {
      return container[Parent::notifier()->id(key)];
    } 
		
    /// \brief The const subcript operator.
    ///
    /// The const subscript operator. The map can be subscripted by the
    /// actual items of the graph. 
    ConstReference operator[](const Key& key) const {
      return container[Parent::notifier()->id(key)];
    }


    /// \brief The setter function of the map.
    ///
    /// It the same as operator[](key) = value expression.
    void set(const Key& key, const Value& value) {
      (*this)[key] = value;
    }

  protected:

    /// \brief Adds a new key to the map.
    ///		
    /// It adds a new key to the map. It called by the observer notifier
    /// and it overrides the add() member function of the observer base.     
    virtual void add(const Key& key) {
      int id = Parent::notifier()->id(key);
      if (id >= int(container.size())) {
	container.resize(id + 1);
      }
    }

    /// \brief Adds more new keys to the map.
    ///		
    /// It adds more new keys to the map. It called by the observer notifier
    /// and it overrides the add() member function of the observer base.     
    virtual void add(const std::vector<Key>& keys) {
      int max = container.size() - 1;
      for (int i = 0; i < int(keys.size()); ++i) {
        int id = Parent::notifier()->id(keys[i]);
        if (id >= max) {
          max = id;
        }
      }
      container.resize(max + 1);
    }

    /// \brief Erase a key from the map.
    ///
    /// Erase a key from the map. It called by the observer notifier
    /// and it overrides the erase() member function of the observer base.     
    virtual void erase(const Key& key) {
      container[Parent::notifier()->id(key)] = Value();
    }

    /// \brief Erase more keys from the map.
    ///
    /// Erase more keys from the map. It called by the observer notifier
    /// and it overrides the erase() member function of the observer base.     
    virtual void erase(const std::vector<Key>& keys) {
      for (int i = 0; i < int(keys.size()); ++i) {
	container[Parent::notifier()->id(keys[i])] = Value();
      }
    }
    
    /// \brief Buildes the map.
    ///	
    /// It buildes the map. It called by the observer notifier
    /// and it overrides the build() member function of the observer base.
    virtual void build() { 
      int size = Parent::notifier()->maxId() + 1;
      container.reserve(size);
      container.resize(size);
    }

    /// \brief Clear the map.
    ///
    /// It erase all items from the map. It called by the observer notifier
    /// and it overrides the clear() member function of the observer base.     
    virtual void clear() { 
      container.clear();
    }
    
  private:
		
    Container container;

  };

}

#endif

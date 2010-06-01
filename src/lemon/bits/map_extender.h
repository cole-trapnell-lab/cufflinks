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

#ifndef LEMON_BITS_MAP_EXTENDER_H
#define LEMON_BITS_MAP_EXTENDER_H

#include <iterator>

#include <lemon/bits/traits.h>

#include <lemon/concept_check.h>
#include <lemon/concepts/maps.h>

///\file
///\brief Extenders for iterable maps.

namespace lemon {

  /// \ingroup graphbits
  /// 
  /// \brief Extender for maps
  template <typename _Map>
  class MapExtender : public _Map {
  public:

    typedef _Map Parent;
    typedef MapExtender Map;


    typedef typename Parent::Graph Graph;
    typedef typename Parent::Key Item;

    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    class MapIt;
    class ConstMapIt;

    friend class MapIt;
    friend class ConstMapIt;

  public:

    MapExtender(const Graph& graph) 
      : Parent(graph) {}

    MapExtender(const Graph& graph, const Value& value) 
      : Parent(graph, value) {}

    MapExtender& operator=(const MapExtender& cmap) {
      return operator=<MapExtender>(cmap);
    }

    template <typename CMap>
    MapExtender& operator=(const CMap& cmap) {
      Parent::operator=(cmap);
      return *this;
    } 

    class MapIt : public Item {
    public:
      
      typedef Item Parent;
      typedef typename Map::Value Value;
      
      MapIt() {}

      MapIt(Invalid i) : Parent(i) { }

      explicit MapIt(Map& _map) : map(_map) {
        map.notifier()->first(*this);
      }

      MapIt(const Map& _map, const Item& item) 
	: Parent(item), map(_map) {}

      MapIt& operator++() { 
	map.notifier()->next(*this);
	return *this; 
      }
      
      typename MapTraits<Map>::ConstReturnValue operator*() const {
	return map[*this];
      }

      typename MapTraits<Map>::ReturnValue operator*() {
	return map[*this];
      }
      
      void set(const Value& value) {
	map.set(*this, value);
      }
      
    protected:
      Map& map;
      
    };

    class ConstMapIt : public Item {
    public:

      typedef Item Parent;

      typedef typename Map::Value Value;
      
      ConstMapIt() {}

      ConstMapIt(Invalid i) : Parent(i) { }

      explicit ConstMapIt(Map& _map) : map(_map) {
        map.notifier()->first(*this);
      }

      ConstMapIt(const Map& _map, const Item& item) 
	: Parent(item), map(_map) {}

      ConstMapIt& operator++() { 
	map.notifier()->next(*this);
	return *this; 
      }

      typename MapTraits<Map>::ConstReturnValue operator*() const {
	return map[*this];
      }

    protected:
      const Map& map;
    };

    class ItemIt : public Item {
    public:
      
      typedef Item Parent;
      
      ItemIt() {}

      ItemIt(Invalid i) : Parent(i) { }

      explicit ItemIt(Map& _map) : map(_map) {
        map.notifier()->first(*this);
      }

      ItemIt(const Map& _map, const Item& item) 
	: Parent(item), map(_map) {}

      ItemIt& operator++() { 
	map.notifier()->next(*this);
	return *this; 
      }

    protected:
      const Map& map;
      
    };
  };

  /// \ingroup graphbits
  /// 
  /// \brief Extender for maps which use a subset of the items.
  template <typename _Graph, typename _Map>
  class SubMapExtender : public _Map {
  public:

    typedef _Map Parent;
    typedef SubMapExtender Map;

    typedef _Graph Graph;

    typedef typename Parent::Key Item;

    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    class MapIt;
    class ConstMapIt;

    friend class MapIt;
    friend class ConstMapIt;

  public:

    SubMapExtender(const Graph& _graph) 
      : Parent(_graph), graph(_graph) {}

    SubMapExtender(const Graph& _graph, const Value& _value) 
      : Parent(_graph, _value), graph(_graph) {}

    SubMapExtender& operator=(const SubMapExtender& cmap) {
      return operator=<MapExtender>(cmap);
    }

    template <typename CMap>
    SubMapExtender& operator=(const CMap& cmap) {
      checkConcept<concepts::ReadMap<Key, Value>, CMap>();
      Item it;
      for (graph.first(it); it != INVALID; graph.next(it)) {
        Parent::set(it, cmap[it]);
      }
      return *this;
    } 

    class MapIt : public Item {
    public:
      
      typedef Item Parent;
      typedef typename Map::Value Value;
      
      MapIt() {}

      MapIt(Invalid i) : Parent(i) { }

      explicit MapIt(Map& _map) : map(_map) {
        map.graph.first(*this);
      }

      MapIt(const Map& _map, const Item& item) 
	: Parent(item), map(_map) {}

      MapIt& operator++() { 
	map.graph.next(*this);
	return *this; 
      }
      
      typename MapTraits<Map>::ConstReturnValue operator*() const {
	return map[*this];
      }

      typename MapTraits<Map>::ReturnValue operator*() {
	return map[*this];
      }
      
      void set(const Value& value) {
	map.set(*this, value);
      }
      
    protected:
      Map& map;
      
    };

    class ConstMapIt : public Item {
    public:

      typedef Item Parent;

      typedef typename Map::Value Value;
      
      ConstMapIt() {}

      ConstMapIt(Invalid i) : Parent(i) { }

      explicit ConstMapIt(Map& _map) : map(_map) {
        map.graph.first(*this);
      }

      ConstMapIt(const Map& _map, const Item& item) 
	: Parent(item), map(_map) {}

      ConstMapIt& operator++() { 
	map.graph.next(*this);
	return *this; 
      }

      typename MapTraits<Map>::ConstReturnValue operator*() const {
	return map[*this];
      }

    protected:
      const Map& map;
    };

    class ItemIt : public Item {
    public:
      
      typedef Item Parent;
      
      ItemIt() {}

      ItemIt(Invalid i) : Parent(i) { }

      explicit ItemIt(Map& _map) : map(_map) {
        map.graph.first(*this);
      }

      ItemIt(const Map& _map, const Item& item) 
	: Parent(item), map(_map) {}

      ItemIt& operator++() { 
	map.graph.next(*this);
	return *this; 
      }

    protected:
      const Map& map;
      
    };
    
  private:

    const Graph& graph;
    
  };

}

#endif

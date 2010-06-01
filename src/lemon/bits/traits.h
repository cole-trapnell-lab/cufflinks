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

#ifndef LEMON_BITS_TRAITS_H
#define LEMON_BITS_TRAITS_H

#include <lemon/bits/utility.h>

///\file
///\brief Traits for graphs and maps
///

namespace lemon {
  template <typename _Graph, typename _Item>
  class ItemSetTraits {};
  

  template <typename Graph, typename Enable = void>
  struct NodeNotifierIndicator {
    typedef InvalidType Type;
  };
  template <typename Graph>
  struct NodeNotifierIndicator<
    Graph, 
    typename enable_if<typename Graph::NodeNotifier::Notifier, void>::type
  > { 
    typedef typename Graph::NodeNotifier Type;
  };

  template <typename _Graph>
  class ItemSetTraits<_Graph, typename _Graph::Node> {
  public:
    
    typedef _Graph Graph;

    typedef typename Graph::Node Item;
    typedef typename Graph::NodeIt ItemIt;

    typedef typename NodeNotifierIndicator<Graph>::Type ItemNotifier;

    template <typename _Value>
    class Map : public Graph::template NodeMap<_Value> {
    public:
      typedef typename Graph::template NodeMap<_Value> Parent; 
      typedef typename Graph::template NodeMap<_Value> Type; 
      typedef typename Parent::Value Value;

      Map(const Graph& _graph) : Parent(_graph) {}
      Map(const Graph& _graph, const Value& _value) 
	: Parent(_graph, _value) {}

     };

  };

  template <typename Graph, typename Enable = void>
  struct EdgeNotifierIndicator {
    typedef InvalidType Type;
  };
  template <typename Graph>
  struct EdgeNotifierIndicator<
    Graph, 
    typename enable_if<typename Graph::EdgeNotifier::Notifier, void>::type
  > { 
    typedef typename Graph::EdgeNotifier Type;
  };

  template <typename _Graph>
  class ItemSetTraits<_Graph, typename _Graph::Edge> {
  public:
    
    typedef _Graph Graph;

    typedef typename Graph::Edge Item;
    typedef typename Graph::EdgeIt ItemIt;

    typedef typename EdgeNotifierIndicator<Graph>::Type ItemNotifier;

    template <typename _Value>
    class Map : public Graph::template EdgeMap<_Value> {
    public:
      typedef typename Graph::template EdgeMap<_Value> Parent; 
      typedef typename Graph::template EdgeMap<_Value> Type; 
      typedef typename Parent::Value Value;

      Map(const Graph& _graph) : Parent(_graph) {}
      Map(const Graph& _graph, const Value& _value) 
	: Parent(_graph, _value) {}
    };

  };

  template <typename Graph, typename Enable = void>
  struct UEdgeNotifierIndicator {
    typedef InvalidType Type;
  };
  template <typename Graph>
  struct UEdgeNotifierIndicator<
    Graph, 
    typename enable_if<typename Graph::UEdgeNotifier::Notifier, void>::type
  > { 
    typedef typename Graph::UEdgeNotifier Type;
  };

  template <typename _Graph>
  class ItemSetTraits<_Graph, typename _Graph::UEdge> {
  public:
    
    typedef _Graph Graph;

    typedef typename Graph::UEdge Item;
    typedef typename Graph::UEdgeIt ItemIt;

    typedef typename UEdgeNotifierIndicator<Graph>::Type ItemNotifier;

    template <typename _Value>
    class Map : public Graph::template UEdgeMap<_Value> {
    public:
      typedef typename Graph::template UEdgeMap<_Value> Parent; 
      typedef typename Graph::template UEdgeMap<_Value> Type; 
      typedef typename Parent::Value Value;

      Map(const Graph& _graph) : Parent(_graph) {}
      Map(const Graph& _graph, const Value& _value) 
	: Parent(_graph, _value) {}
    };

  };

  template <typename Graph, typename Enable = void>
  struct ANodeNotifierIndicator {
    typedef InvalidType Type;
  };
  template <typename Graph>
  struct ANodeNotifierIndicator<
    Graph, 
    typename enable_if<typename Graph::ANodeNotifier::Notifier, void>::type
  > { 
    typedef typename Graph::ANodeNotifier Type;
  };

  template <typename _Graph>
  class ItemSetTraits<_Graph, typename _Graph::ANode> {
  public:
    
    typedef _Graph Graph;

    typedef typename Graph::ANode Item;
    typedef typename Graph::ANodeIt ItemIt;

    typedef typename ANodeNotifierIndicator<Graph>::Type ItemNotifier;

    template <typename _Value>
    class Map : public Graph::template ANodeMap<_Value> {
    public:
      typedef typename Graph::template ANodeMap<_Value> Parent; 
      typedef typename Graph::template ANodeMap<_Value> Type; 
      typedef typename Parent::Value Value;

      Map(const Graph& _graph) : Parent(_graph) {}
      Map(const Graph& _graph, const Value& _value) 
	: Parent(_graph, _value) {}
    };

  };

  template <typename Graph, typename Enable = void>
  struct BNodeNotifierIndicator {
    typedef InvalidType Type;
  };
  template <typename Graph>
  struct BNodeNotifierIndicator<
    Graph, 
    typename enable_if<typename Graph::BNodeNotifier::Notifier, void>::type
  > { 
    typedef typename Graph::BNodeNotifier Type;
  };

  template <typename _Graph>
  class ItemSetTraits<_Graph, typename _Graph::BNode> {
  public:
    
    typedef _Graph Graph;

    typedef typename Graph::BNode Item;
    typedef typename Graph::BNodeIt ItemIt;

    typedef typename BNodeNotifierIndicator<Graph>::Type ItemNotifier;

    template <typename _Value>
    class Map : public Graph::template BNodeMap<_Value> {
    public:
      typedef typename Graph::template BNodeMap<_Value> Parent; 
      typedef typename Graph::template BNodeMap<_Value> Type; 
      typedef typename Parent::Value Value;

      Map(const Graph& _graph) : Parent(_graph) {}
      Map(const Graph& _graph, const Value& _value) 
	: Parent(_graph, _value) {}
    };

  };


  template <typename Map, typename Enable = void>
  struct MapTraits {
    typedef False ReferenceMapTag;

    typedef typename Map::Key Key;
    typedef typename Map::Value Value;

    typedef const Value ConstReturnValue;
    typedef const Value ReturnValue;
  };

  template <typename Map>
  struct MapTraits<
    Map, typename enable_if<typename Map::ReferenceMapTag, void>::type > 
  {
    typedef True ReferenceMapTag;
    
    typedef typename Map::Key Key;
    typedef typename Map::Value Value;

    typedef typename Map::ConstReference ConstReturnValue;
    typedef typename Map::Reference ReturnValue;

    typedef typename Map::ConstReference ConstReference; 
    typedef typename Map::Reference Reference;
 };

  template <typename MatrixMap, typename Enable = void>
  struct MatrixMapTraits {
    typedef False ReferenceMapTag;

    typedef typename MatrixMap::FirstKey FirstKey;
    typedef typename MatrixMap::SecondKey SecondKey;
    typedef typename MatrixMap::Value Value;

    typedef const Value ConstReturnValue;
    typedef const Value ReturnValue;
  };

  template <typename MatrixMap>
  struct MatrixMapTraits<
    MatrixMap, typename enable_if<typename MatrixMap::ReferenceMapTag, 
                                  void>::type > 
  {
    typedef True ReferenceMapTag;
    
    typedef typename MatrixMap::FirstKey FirstKey;
    typedef typename MatrixMap::SecondKey SecondKey;
    typedef typename MatrixMap::Value Value;

    typedef typename MatrixMap::ConstReference ConstReturnValue;
    typedef typename MatrixMap::Reference ReturnValue;

    typedef typename MatrixMap::ConstReference ConstReference; 
    typedef typename MatrixMap::Reference Reference;
 };

  // Indicators for the tags

  template <typename Graph, typename Enable = void>
  struct NodeNumTagIndicator {
    static const bool value = false;
  };

  template <typename Graph>
  struct NodeNumTagIndicator<
    Graph, 
    typename enable_if<typename Graph::NodeNumTag, void>::type
  > {
    static const bool value = true;
  };

  template <typename Graph, typename Enable = void>
  struct EdgeNumTagIndicator {
    static const bool value = false;
  };

  template <typename Graph>
  struct EdgeNumTagIndicator<
    Graph, 
    typename enable_if<typename Graph::EdgeNumTag, void>::type
  > {
    static const bool value = true;
  };

  template <typename Graph, typename Enable = void>
  struct FindEdgeTagIndicator {
    static const bool value = false;
  };

  template <typename Graph>
  struct FindEdgeTagIndicator<
    Graph, 
    typename enable_if<typename Graph::FindEdgeTag, void>::type
  > {
    static const bool value = true;
  };

  template <typename Graph, typename Enable = void>
  struct UndirectedTagIndicator {
    static const bool value = false;
  };

  template <typename Graph>
  struct UndirectedTagIndicator<
    Graph, 
    typename enable_if<typename Graph::UndirectedTag, void>::type
  > {
    static const bool value = true;
  };

  template <typename Graph, typename Enable = void>
  struct BuildTagIndicator {
    static const bool value = false;
  };

  template <typename Graph>
  struct BuildTagIndicator<
    Graph, 
    typename enable_if<typename Graph::BuildTag, void>::type
  > {
    static const bool value = true;
  };

}

#endif

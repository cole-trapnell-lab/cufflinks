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

#ifndef LEMON_BITS_GRAPH_ADAPTOR_EXTENDER_H
#define LEMON_BITS_GRAPH_ADAPTOR_EXTENDER_H

#include <lemon/bits/invalid.h>
#include <lemon/error.h>

#include <lemon/bits/default_map.h>


///\ingroup graphbits
///\file
///\brief Extenders for the graph adaptor types
namespace lemon {

  /// \ingroup graphbits
  ///
  /// \brief Extender for the GraphAdaptors
  template <typename _Graph>
  class GraphAdaptorExtender : public _Graph {
  public:

    typedef _Graph Parent;
    typedef _Graph Graph;
    typedef GraphAdaptorExtender Adaptor;

    // Base extensions

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;

    int maxId(Node) const {
      return Parent::maxNodeId();
    }

    int maxId(Edge) const {
      return Parent::maxEdgeId();
    }

    Node fromId(int id, Node) const {
      return Parent::nodeFromId(id);
    }

    Edge fromId(int id, Edge) const {
      return Parent::edgeFromId(id);
    }

    Node oppositeNode(const Node &n, const Edge &e) const {
      if (n == Parent::source(e))
	return Parent::target(e);
      else if(n==Parent::target(e))
	return Parent::source(e);
      else
	return INVALID;
    }

    class NodeIt : public Node { 
      const Adaptor* graph;
    public:

      NodeIt() {}

      NodeIt(Invalid i) : Node(i) { }

      explicit NodeIt(const Adaptor& _graph) : graph(&_graph) {
	_graph.first(static_cast<Node&>(*this));
      }

      NodeIt(const Adaptor& _graph, const Node& node) 
	: Node(node), graph(&_graph) {}

      NodeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };


    class EdgeIt : public Edge { 
      const Adaptor* graph;
    public:

      EdgeIt() { }

      EdgeIt(Invalid i) : Edge(i) { }

      explicit EdgeIt(const Adaptor& _graph) : graph(&_graph) {
	_graph.first(static_cast<Edge&>(*this));
      }

      EdgeIt(const Adaptor& _graph, const Edge& e) : 
	Edge(e), graph(&_graph) { }

      EdgeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };


    class OutEdgeIt : public Edge { 
      const Adaptor* graph;
    public:

      OutEdgeIt() { }

      OutEdgeIt(Invalid i) : Edge(i) { }

      OutEdgeIt(const Adaptor& _graph, const Node& node) 
	: graph(&_graph) {
	_graph.firstOut(*this, node);
      }

      OutEdgeIt(const Adaptor& _graph, const Edge& edge) 
	: Edge(edge), graph(&_graph) {}

      OutEdgeIt& operator++() { 
	graph->nextOut(*this);
	return *this; 
      }

    };


    class InEdgeIt : public Edge { 
      const Adaptor* graph;
    public:

      InEdgeIt() { }

      InEdgeIt(Invalid i) : Edge(i) { }

      InEdgeIt(const Adaptor& _graph, const Node& node) 
	: graph(&_graph) {
	_graph.firstIn(*this, node);
      }

      InEdgeIt(const Adaptor& _graph, const Edge& edge) : 
	Edge(edge), graph(&_graph) {}

      InEdgeIt& operator++() { 
	graph->nextIn(*this);
	return *this; 
      }

    };

    /// \brief Base node of the iterator
    ///
    /// Returns the base node (ie. the source in this case) of the iterator
    Node baseNode(const OutEdgeIt &e) const {
      return Parent::source(e);
    }
    /// \brief Running node of the iterator
    ///
    /// Returns the running node (ie. the target in this case) of the
    /// iterator
    Node runningNode(const OutEdgeIt &e) const {
      return Parent::target(e);
    }

    /// \brief Base node of the iterator
    ///
    /// Returns the base node (ie. the target in this case) of the iterator
    Node baseNode(const InEdgeIt &e) const {
      return Parent::target(e);
    }
    /// \brief Running node of the iterator
    ///
    /// Returns the running node (ie. the source in this case) of the
    /// iterator
    Node runningNode(const InEdgeIt &e) const {
      return Parent::source(e);
    }

  };


  /// \ingroup graphbits
  ///
  /// \brief Extender for the UGraphAdaptors
  template <typename _UGraph> 
  class UGraphAdaptorExtender : public _UGraph {
  public:
    
    typedef _UGraph Parent;
    typedef _UGraph UGraph;
    typedef UGraphAdaptorExtender Adaptor;

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;
    typedef typename Parent::UEdge UEdge;

    // UGraph extension    

    int maxId(Node) const {
      return Parent::maxNodeId();
    }

    int maxId(Edge) const {
      return Parent::maxEdgeId();
    }

    int maxId(UEdge) const {
      return Parent::maxUEdgeId();
    }

    Node fromId(int id, Node) const {
      return Parent::nodeFromId(id);
    }

    Edge fromId(int id, Edge) const {
      return Parent::edgeFromId(id);
    }

    UEdge fromId(int id, UEdge) const {
      return Parent::uEdgeFromId(id);
    }

    Node oppositeNode(const Node &n, const UEdge &e) const {
      if( n == Parent::source(e))
	return Parent::target(e);
      else if( n == Parent::target(e))
	return Parent::source(e);
      else
	return INVALID;
    }

    Edge oppositeEdge(const Edge &e) const {
      return Parent::direct(e, !Parent::direction(e));
    }

    using Parent::direct;
    Edge direct(const UEdge &ue, const Node &s) const {
      return Parent::direct(ue, Parent::source(ue) == s);
    }


    class NodeIt : public Node { 
      const Adaptor* graph;
    public:

      NodeIt() {}

      NodeIt(Invalid i) : Node(i) { }

      explicit NodeIt(const Adaptor& _graph) : graph(&_graph) {
	_graph.first(static_cast<Node&>(*this));
      }

      NodeIt(const Adaptor& _graph, const Node& node) 
	: Node(node), graph(&_graph) {}

      NodeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };


    class EdgeIt : public Edge { 
      const Adaptor* graph;
    public:

      EdgeIt() { }

      EdgeIt(Invalid i) : Edge(i) { }

      explicit EdgeIt(const Adaptor& _graph) : graph(&_graph) {
	_graph.first(static_cast<Edge&>(*this));
      }

      EdgeIt(const Adaptor& _graph, const Edge& e) : 
	Edge(e), graph(&_graph) { }

      EdgeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };


    class OutEdgeIt : public Edge { 
      const Adaptor* graph;
    public:

      OutEdgeIt() { }

      OutEdgeIt(Invalid i) : Edge(i) { }

      OutEdgeIt(const Adaptor& _graph, const Node& node) 
	: graph(&_graph) {
	_graph.firstOut(*this, node);
      }

      OutEdgeIt(const Adaptor& _graph, const Edge& edge) 
	: Edge(edge), graph(&_graph) {}

      OutEdgeIt& operator++() { 
	graph->nextOut(*this);
	return *this; 
      }

    };


    class InEdgeIt : public Edge { 
      const Adaptor* graph;
    public:

      InEdgeIt() { }

      InEdgeIt(Invalid i) : Edge(i) { }

      InEdgeIt(const Adaptor& _graph, const Node& node) 
	: graph(&_graph) {
	_graph.firstIn(*this, node);
      }

      InEdgeIt(const Adaptor& _graph, const Edge& edge) : 
	Edge(edge), graph(&_graph) {}

      InEdgeIt& operator++() { 
	graph->nextIn(*this);
	return *this; 
      }

    };

    class UEdgeIt : public Parent::UEdge { 
      const Adaptor* graph;
    public:

      UEdgeIt() { }

      UEdgeIt(Invalid i) : UEdge(i) { }

      explicit UEdgeIt(const Adaptor& _graph) : graph(&_graph) {
	_graph.first(static_cast<UEdge&>(*this));
      }

      UEdgeIt(const Adaptor& _graph, const UEdge& e) : 
	UEdge(e), graph(&_graph) { }

      UEdgeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };

    class IncEdgeIt : public Parent::UEdge { 
      friend class UGraphAdaptorExtender;
      const Adaptor* graph;
      bool direction;
    public:

      IncEdgeIt() { }

      IncEdgeIt(Invalid i) : UEdge(i), direction(false) { }

      IncEdgeIt(const Adaptor& _graph, const Node &n) : graph(&_graph) {
	_graph.firstInc(static_cast<UEdge&>(*this), direction, n);
      }

      IncEdgeIt(const Adaptor& _graph, const UEdge &ue, const Node &n)
	: graph(&_graph), UEdge(ue) {
	direction = (_graph.source(ue) == n);
      }

      IncEdgeIt& operator++() {
	graph->nextInc(*this, direction);
	return *this; 
      }
    };

    /// \brief Base node of the iterator
    ///
    /// Returns the base node (ie. the source in this case) of the iterator
    Node baseNode(const OutEdgeIt &e) const {
      return Parent::source(static_cast<const Edge&>(e));
    }
    /// \brief Running node of the iterator
    ///
    /// Returns the running node (ie. the target in this case) of the
    /// iterator
    Node runningNode(const OutEdgeIt &e) const {
      return Parent::target(static_cast<const Edge&>(e));
    }

    /// \brief Base node of the iterator
    ///
    /// Returns the base node (ie. the target in this case) of the iterator
    Node baseNode(const InEdgeIt &e) const {
      return Parent::target(static_cast<const Edge&>(e));
    }
    /// \brief Running node of the iterator
    ///
    /// Returns the running node (ie. the source in this case) of the
    /// iterator
    Node runningNode(const InEdgeIt &e) const {
      return Parent::source(static_cast<const Edge&>(e));
    }

    /// Base node of the iterator
    ///
    /// Returns the base node of the iterator
    Node baseNode(const IncEdgeIt &e) const {
      return e.direction ? source(e) : target(e);
    }
    /// Running node of the iterator
    ///
    /// Returns the running node of the iterator
    Node runningNode(const IncEdgeIt &e) const {
      return e.direction ? target(e) : source(e);
    }

  };

  /// \ingroup graphbits
  ///
  /// \brief Extender for the BpUGraphAdaptors
  template <typename Base>
  class BpUGraphAdaptorExtender : public Base {
  public:
    typedef Base Parent;
    typedef BpUGraphAdaptorExtender Graph;

    typedef typename Parent::Node Node;
    typedef typename Parent::BNode BNode;
    typedef typename Parent::ANode ANode;
    typedef typename Parent::Edge Edge;
    typedef typename Parent::UEdge UEdge;


    int maxId(Node) const {
      return Parent::maxNodeId();
    }
    int maxId(BNode) const {
      return Parent::maxBNodeId();
    }
    int maxId(ANode) const {
      return Parent::maxANodeId();
    }
    int maxId(Edge) const {
      return Parent::maxEdgeId();
    }
    int maxId(UEdge) const {
      return Parent::maxUEdgeId();
    }


    Node fromId(int id, Node) const {
      return Parent::nodeFromId(id);
    }
    ANode fromId(int id, ANode) const {
      return Parent::nodeFromANodeId(id);
    }
    BNode fromId(int id, BNode) const {
      return Parent::nodeFromBNodeId(id);
    }
    Edge fromId(int id, Edge) const {
      return Parent::edgeFromId(id);
    }
    UEdge fromId(int id, UEdge) const {
      return Parent::uEdgeFromId(id);
    }  
  
    class NodeIt : public Node { 
      const Graph* graph;
    public:
    
      NodeIt() { }
    
      NodeIt(Invalid i) : Node(INVALID) { }
    
      explicit NodeIt(const Graph& _graph) : graph(&_graph) {
	graph->first(static_cast<Node&>(*this));
      }

      NodeIt(const Graph& _graph, const Node& node) 
	: Node(node), graph(&_graph) { }
    
      NodeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };

    class ANodeIt : public Node { 
      friend class BpUGraphAdaptorExtender;
      const Graph* graph;
    public:
    
      ANodeIt() { }
    
      ANodeIt(Invalid i) : Node(INVALID) { }
    
      explicit ANodeIt(const Graph& _graph) : graph(&_graph) {
	graph->firstANode(static_cast<Node&>(*this));
      }

      ANodeIt(const Graph& _graph, const Node& node) 
	: Node(node), graph(&_graph) {}
    
      ANodeIt& operator++() { 
	graph->nextANode(*this);
	return *this; 
      }
    };

    class BNodeIt : public Node { 
      friend class BpUGraphAdaptorExtender;
      const Graph* graph;
    public:
    
      BNodeIt() { }
    
      BNodeIt(Invalid i) : Node(INVALID) { }
    
      explicit BNodeIt(const Graph& _graph) : graph(&_graph) {
	graph->firstBNode(static_cast<Node&>(*this));
      }

      BNodeIt(const Graph& _graph, const Node& node) 
	: Node(node), graph(&_graph) {}
    
      BNodeIt& operator++() { 
	graph->nextBNode(*this);
	return *this; 
      }
    };

    class EdgeIt : public Edge { 
      friend class BpUGraphAdaptorExtender;
      const Graph* graph;
    public:
    
      EdgeIt() { }
    
      EdgeIt(Invalid i) : Edge(INVALID) { }
    
      explicit EdgeIt(const Graph& _graph) : graph(&_graph) {
	graph->first(static_cast<Edge&>(*this));
      }

      EdgeIt(const Graph& _graph, const Edge& edge) 
	: Edge(edge), graph(&_graph) { }
    
      EdgeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }

    };

    class UEdgeIt : public UEdge { 
      friend class BpUGraphAdaptorExtender;
      const Graph* graph;
    public:
    
      UEdgeIt() { }
    
      UEdgeIt(Invalid i) : UEdge(INVALID) { }
    
      explicit UEdgeIt(const Graph& _graph) : graph(&_graph) {
	graph->first(static_cast<UEdge&>(*this));
      }

      UEdgeIt(const Graph& _graph, const UEdge& edge) 
	: UEdge(edge), graph(&_graph) { }
    
      UEdgeIt& operator++() { 
	graph->next(*this);
	return *this; 
      }
    };

    class OutEdgeIt : public Edge { 
      friend class BpUGraphAdaptorExtender;
      const Graph* graph;
    public:
    
      OutEdgeIt() { }
    
      OutEdgeIt(Invalid i) : Edge(i) { }
    
      OutEdgeIt(const Graph& _graph, const Node& node) 
	: graph(&_graph) {
	graph->firstOut(*this, node);
      }
    
      OutEdgeIt(const Graph& _graph, const Edge& edge) 
	: Edge(edge), graph(&_graph) {}
    
      OutEdgeIt& operator++() { 
	graph->nextOut(*this);
	return *this; 
      }

    };


    class InEdgeIt : public Edge { 
      friend class BpUGraphAdaptorExtender;
      const Graph* graph;
    public:
    
      InEdgeIt() { }
    
      InEdgeIt(Invalid i) : Edge(i) { }
    
      InEdgeIt(const Graph& _graph, const Node& node) 
	: graph(&_graph) {
	graph->firstIn(*this, node);
      }
    
      InEdgeIt(const Graph& _graph, const Edge& edge) : 
	Edge(edge), graph(&_graph) {}
    
      InEdgeIt& operator++() { 
	graph->nextIn(*this);
	return *this; 
      }

    };
  
    /// \brief Base node of the iterator
    ///
    /// Returns the base node (ie. the source in this case) of the iterator
    Node baseNode(const OutEdgeIt &e) const {
      return Parent::source(static_cast<const Edge&>(e));
    }
    /// \brief Running node of the iterator
    ///
    /// Returns the running node (ie. the target in this case) of the
    /// iterator
    Node runningNode(const OutEdgeIt &e) const {
      return Parent::target(static_cast<const Edge&>(e));
    }
  
    /// \brief Base node of the iterator
    ///
    /// Returns the base node (ie. the target in this case) of the iterator
    Node baseNode(const InEdgeIt &e) const {
      return Parent::target(static_cast<const Edge&>(e));
    }
    /// \brief Running node of the iterator
    ///
    /// Returns the running node (ie. the source in this case) of the
    /// iterator
    Node runningNode(const InEdgeIt &e) const {
      return Parent::source(static_cast<const Edge&>(e));
    }
  
    class IncEdgeIt : public Parent::UEdge { 
      friend class BpUGraphAdaptorExtender;
      const Graph* graph;
      bool direction;
    public:
    
      IncEdgeIt() { }
    
      IncEdgeIt(Invalid i) : UEdge(i), direction(true) { }
    
      IncEdgeIt(const Graph& _graph, const Node &n) : graph(&_graph) {
	graph->firstInc(*this, direction, n);
      }

      IncEdgeIt(const Graph& _graph, const UEdge &ue, const Node &n)
	: graph(&_graph), UEdge(ue) {
	direction = (graph->source(ue) == n);
      }

      IncEdgeIt& operator++() {
	graph->nextInc(*this, direction);
	return *this; 
      }
    };
  

    /// Base node of the iterator
    ///
    /// Returns the base node of the iterator
    Node baseNode(const IncEdgeIt &e) const {
      return e.direction ? source(e) : target(e);
    }

    /// Running node of the iterator
    ///
    /// Returns the running node of the iterator
    Node runningNode(const IncEdgeIt &e) const {
      return e.direction ? target(e) : source(e);
    }

    Node oppositeNode(const Node &n, const UEdge &e) const {
      if( n == Parent::source(e))
	return Parent::target(e);
      else if( n == Parent::target(e))
	return Parent::source(e);
      else
	return INVALID;
    }

    Edge oppositeEdge(const Edge &e) const {
      return Parent::direct(e, !Parent::direction(e));
    }

    using Parent::direct;
    Edge direct(const UEdge &ue, const Node &s) const {
      return Parent::direct(ue, Parent::source(ue) == s);
    }

  };


}


#endif

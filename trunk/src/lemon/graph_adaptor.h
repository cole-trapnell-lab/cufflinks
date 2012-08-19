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

#ifndef LEMON_GRAPH_ADAPTOR_H
#define LEMON_GRAPH_ADAPTOR_H

///\ingroup graph_adaptors
///\file
///\brief Several graph adaptors.
///
///This file contains several useful graph adaptor functions.
///
///\author Marton Makai and Balazs Dezso

#include <lemon/bits/invalid.h>
#include <lemon/bits/variant.h>
#include <lemon/maps.h>

#include <lemon/bits/base_extender.h>
#include <lemon/bits/graph_adaptor_extender.h>
#include <lemon/bits/graph_extender.h>
#include <lemon/tolerance.h>

#include <algorithm>

namespace lemon {

  ///\brief Base type for the Graph Adaptors
  ///
  ///Base type for the Graph Adaptors
  ///
  ///This is the base type for most of LEMON graph adaptors. 
  ///This class implements a trivial graph adaptor i.e. it only wraps the 
  ///functions and types of the graph. The purpose of this class is to 
  ///make easier implementing graph adaptors. E.g. if an adaptor is 
  ///considered which differs from the wrapped graph only in some of its 
  ///functions or types, then it can be derived from GraphAdaptor,
  ///and only the 
  ///differences should be implemented.
  ///
  ///author Marton Makai 
  template<typename _Graph>
  class GraphAdaptorBase {
  public:
    typedef _Graph Graph;
    typedef GraphAdaptorBase Adaptor;
    typedef Graph ParentGraph;

  protected:
    Graph* graph;
    GraphAdaptorBase() : graph(0) { }
    void setGraph(Graph& _graph) { graph=&_graph; }

  public:
    GraphAdaptorBase(Graph& _graph) : graph(&_graph) { }

    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
   
    void first(Node& i) const { graph->first(i); }
    void first(Edge& i) const { graph->first(i); }
    void firstIn(Edge& i, const Node& n) const { graph->firstIn(i, n); }
    void firstOut(Edge& i, const Node& n ) const { graph->firstOut(i, n); }

    void next(Node& i) const { graph->next(i); }
    void next(Edge& i) const { graph->next(i); }
    void nextIn(Edge& i) const { graph->nextIn(i); }
    void nextOut(Edge& i) const { graph->nextOut(i); }

    Node source(const Edge& e) const { return graph->source(e); }
    Node target(const Edge& e) const { return graph->target(e); }

    typedef NodeNumTagIndicator<Graph> NodeNumTag;
    int nodeNum() const { return graph->nodeNum(); }
    
    typedef EdgeNumTagIndicator<Graph> EdgeNumTag;
    int edgeNum() const { return graph->edgeNum(); }

    typedef FindEdgeTagIndicator<Graph> FindEdgeTag;
    Edge findEdge(const Node& u, const Node& v, 
		  const Edge& prev = INVALID) {
      return graph->findEdge(u, v, prev);
    }
  
    Node addNode() const { 
      return Node(graph->addNode()); 
    }

    Edge addEdge(const Node& u, const Node& v) const { 
      return Edge(graph->addEdge(u, v)); 
    }

    void erase(const Node& i) const { graph->erase(i); }
    void erase(const Edge& i) const { graph->erase(i); }
  
    void clear() const { graph->clear(); }
    
    int id(const Node& v) const { return graph->id(v); }
    int id(const Edge& e) const { return graph->id(e); }

    Node fromNodeId(int ix) const {
      return graph->fromNodeId(ix);
    }

    Edge fromEdgeId(int ix) const {
      return graph->fromEdgeId(ix);
    }

    int maxNodeId() const {
      return graph->maxNodeId();
    }

    int maxEdgeId() const {
      return graph->maxEdgeId();
    }

    typedef typename ItemSetTraits<Graph, Node>::ItemNotifier NodeNotifier;

    NodeNotifier& notifier(Node) const {
      return graph->notifier(Node());
    } 

    typedef typename ItemSetTraits<Graph, Edge>::ItemNotifier EdgeNotifier;

    EdgeNotifier& notifier(Edge) const {
      return graph->notifier(Edge());
    } 
    
    template <typename _Value>
    class NodeMap : public Graph::template NodeMap<_Value> {
    public:

      typedef typename Graph::template NodeMap<_Value> Parent;

      explicit NodeMap(const Adaptor& ga) 
	: Parent(*ga.graph) {}

      NodeMap(const Adaptor& ga, const _Value& value)
	: Parent(*ga.graph, value) { }

      NodeMap& operator=(const NodeMap& cmap) {
        return operator=<NodeMap>(cmap);
      }

      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }
      
    };

    template <typename _Value>
    class EdgeMap : public Graph::template EdgeMap<_Value> {
    public:
      
      typedef typename Graph::template EdgeMap<_Value> Parent;
      
      explicit EdgeMap(const Adaptor& ga) 
	: Parent(*ga.graph) {}

      EdgeMap(const Adaptor& ga, const _Value& value)
	: Parent(*ga.graph, value) {}

      EdgeMap& operator=(const EdgeMap& cmap) {
        return operator=<EdgeMap>(cmap);
      }

      template <typename CMap>
      EdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }

    };

  };

  ///\ingroup graph_adaptors
  ///
  ///\brief Trivial Graph Adaptor
  /// 
  /// This class is an adaptor which does not change the adapted graph.
  /// It can be used only to test the graph adaptors.
  template <typename _Graph>
  class GraphAdaptor :
    public GraphAdaptorExtender<GraphAdaptorBase<_Graph> > { 
  public:
    typedef _Graph Graph;
    typedef GraphAdaptorExtender<GraphAdaptorBase<_Graph> > Parent;
  protected:
    GraphAdaptor() : Parent() { }

  public:
    explicit GraphAdaptor(Graph& _graph) { setGraph(_graph); }
  };

  /// \brief Just gives back a graph adaptor
  ///
  /// Just gives back a graph adaptor which 
  /// should be provide original graph
  template<typename Graph>
  GraphAdaptor<const Graph>
  graphAdaptor(const Graph& graph) {
    return GraphAdaptor<const Graph>(graph);
  }


  template <typename _Graph>
  class RevGraphAdaptorBase : public GraphAdaptorBase<_Graph> {
  public:
    typedef _Graph Graph;
    typedef GraphAdaptorBase<_Graph> Parent;
  protected:
    RevGraphAdaptorBase() : Parent() { }
  public:
    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;

    void firstIn(Edge& i, const Node& n) const { Parent::firstOut(i, n); }
    void firstOut(Edge& i, const Node& n ) const { Parent::firstIn(i, n); }

    void nextIn(Edge& i) const { Parent::nextOut(i); }
    void nextOut(Edge& i) const { Parent::nextIn(i); }

    Node source(const Edge& e) const { return Parent::target(e); }
    Node target(const Edge& e) const { return Parent::source(e); }

    typedef FindEdgeTagIndicator<Graph> FindEdgeTag;
    Edge findEdge(const Node& u, const Node& v, 
		  const Edge& prev = INVALID) {
      return Parent::findEdge(v, u, prev);
    }

  };
    

  ///\ingroup graph_adaptors
  ///
  ///\brief A graph adaptor which reverses the orientation of the edges.
  ///
  /// If \c g is defined as
  ///\code
  /// ListGraph g;
  ///\endcode
  /// then
  ///\code
  /// RevGraphAdaptor<ListGraph> ga(g);
  ///\endcode
  /// implements the graph obtained from \c g by 
  /// reversing the orientation of its edges.
  ///
  /// A good example of using RevGraphAdaptor is to decide that the
  /// directed graph is wheter strongly connected or not. If from one
  /// node each node is reachable and from each node is reachable this
  /// node then and just then the graph is strongly connected. Instead of
  /// this condition we use a little bit different. From one node each node
  /// ahould be reachable in the graph and in the reversed graph. Now this
  /// condition can be checked with the Dfs algorithm class and the
  /// RevGraphAdaptor algorithm class.
  ///
  /// And look at the code:
  ///
  ///\code
  /// bool stronglyConnected(const Graph& graph) {
  ///   if (NodeIt(graph) == INVALID) return true;
  ///   Dfs<Graph> dfs(graph);
  ///   dfs.run(NodeIt(graph));
  ///   for (NodeIt it(graph); it != INVALID; ++it) {
  ///     if (!dfs.reached(it)) {
  ///       return false;
  ///     }
  ///   }
  ///   typedef RevGraphAdaptor<const Graph> RGraph;
  ///   RGraph rgraph(graph);
  ///   DfsVisit<RGraph> rdfs(rgraph);
  ///   rdfs.run(NodeIt(graph));
  ///   for (NodeIt it(graph); it != INVALID; ++it) {
  ///     if (!rdfs.reached(it)) {
  ///       return false;
  ///     }
  ///   }
  ///   return true;
  /// }
  ///\endcode
  template<typename _Graph>
  class RevGraphAdaptor : 
    public GraphAdaptorExtender<RevGraphAdaptorBase<_Graph> > {
  public:
    typedef _Graph Graph;
    typedef GraphAdaptorExtender<
      RevGraphAdaptorBase<_Graph> > Parent;
  protected:
    RevGraphAdaptor() { }
  public:
    explicit RevGraphAdaptor(_Graph& _graph) { setGraph(_graph); }
  };

  /// \brief Just gives back a reverse graph adaptor
  ///
  /// Just gives back a reverse graph adaptor
  template<typename Graph>
  RevGraphAdaptor<const Graph>
  revGraphAdaptor(const Graph& graph) {
    return RevGraphAdaptor<const Graph>(graph);
  }

  template <typename _Graph, typename NodeFilterMap, 
	    typename EdgeFilterMap, bool checked = true>
  class SubGraphAdaptorBase : public GraphAdaptorBase<_Graph> {
  public:
    typedef _Graph Graph;
    typedef SubGraphAdaptorBase Adaptor;
    typedef GraphAdaptorBase<_Graph> Parent;
  protected:
    NodeFilterMap* node_filter_map;
    EdgeFilterMap* edge_filter_map;
    SubGraphAdaptorBase() : Parent(), 
			    node_filter_map(0), edge_filter_map(0) { }

    void setNodeFilterMap(NodeFilterMap& _node_filter_map) {
      node_filter_map=&_node_filter_map;
    }
    void setEdgeFilterMap(EdgeFilterMap& _edge_filter_map) {
      edge_filter_map=&_edge_filter_map;
    }

  public:

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;

    void first(Node& i) const { 
      Parent::first(i); 
      while (i!=INVALID && !(*node_filter_map)[i]) Parent::next(i); 
    }

    void first(Edge& i) const { 
      Parent::first(i); 
      while (i!=INVALID && (!(*edge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::source(i)]
	     || !(*node_filter_map)[Parent::target(i)])) Parent::next(i); 
    }

    void firstIn(Edge& i, const Node& n) const { 
      Parent::firstIn(i, n); 
      while (i!=INVALID && (!(*edge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::source(i)])) Parent::nextIn(i); 
    }

    void firstOut(Edge& i, const Node& n) const { 
      Parent::firstOut(i, n); 
      while (i!=INVALID && (!(*edge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::target(i)])) Parent::nextOut(i); 
    }

    void next(Node& i) const { 
      Parent::next(i); 
      while (i!=INVALID && !(*node_filter_map)[i]) Parent::next(i); 
    }

    void next(Edge& i) const { 
      Parent::next(i); 
      while (i!=INVALID && (!(*edge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::source(i)]
	     || !(*node_filter_map)[Parent::target(i)])) Parent::next(i); 
    }

    void nextIn(Edge& i) const { 
      Parent::nextIn(i); 
      while (i!=INVALID && (!(*edge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::source(i)])) Parent::nextIn(i); 
    }

    void nextOut(Edge& i) const { 
      Parent::nextOut(i); 
      while (i!=INVALID && (!(*edge_filter_map)[i] 
	     || !(*node_filter_map)[Parent::target(i)])) Parent::nextOut(i); 
    }

    ///\e

    /// This function hides \c n in the graph, i.e. the iteration 
    /// jumps over it. This is done by simply setting the value of \c n  
    /// to be false in the corresponding node-map.
    void hide(const Node& n) const { node_filter_map->set(n, false); }

    ///\e

    /// This function hides \c e in the graph, i.e. the iteration 
    /// jumps over it. This is done by simply setting the value of \c e  
    /// to be false in the corresponding edge-map.
    void hide(const Edge& e) const { edge_filter_map->set(e, false); }

    ///\e

    /// The value of \c n is set to be true in the node-map which stores 
    /// hide information. If \c n was hidden previuosly, then it is shown 
    /// again
     void unHide(const Node& n) const { node_filter_map->set(n, true); }

    ///\e

    /// The value of \c e is set to be true in the edge-map which stores 
    /// hide information. If \c e was hidden previuosly, then it is shown 
    /// again
    void unHide(const Edge& e) const { edge_filter_map->set(e, true); }

    /// Returns true if \c n is hidden.
    
    ///\e
    ///
    bool hidden(const Node& n) const { return !(*node_filter_map)[n]; }

    /// Returns true if \c n is hidden.
    
    ///\e
    ///
    bool hidden(const Edge& e) const { return !(*edge_filter_map)[e]; }

    typedef False NodeNumTag;
    typedef False EdgeNumTag;

    typedef FindEdgeTagIndicator<Graph> FindEdgeTag;
    Edge findEdge(const Node& source, const Node& target, 
		  const Edge& prev = INVALID) {
      if (!(*node_filter_map)[source] || !(*node_filter_map)[target]) {
        return INVALID;
      }
      Edge edge = Parent::findEdge(source, target, prev);
      while (edge != INVALID && !(*edge_filter_map)[edge]) {
        edge = Parent::findEdge(source, target, edge);
      }
      return edge;
    }

    template <typename _Value>
    class NodeMap 
      : public SubMapExtender<Adaptor, 
                              typename Parent::template NodeMap<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      //typedef SubMapExtender<Adaptor, typename Parent::
      //                       template NodeMap<_Value> > Parent;
    
      NodeMap(const Graph& g) 
	: Parent(g) {}
      NodeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
      NodeMap& operator=(const NodeMap& cmap) {
	return operator=<NodeMap>(cmap);
      }
    
      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };

    template <typename _Value>
    class EdgeMap 
      : public SubMapExtender<Adaptor, 
                              typename Parent::template EdgeMap<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      //typedef SubMapExtender<Adaptor, typename Parent::
      //                       template EdgeMap<_Value> > Parent;
    
      EdgeMap(const Graph& g) 
	: Parent(g) {}
      EdgeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
      EdgeMap& operator=(const EdgeMap& cmap) {
	return operator=<EdgeMap>(cmap);
      }
    
      template <typename CMap>
      EdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };

  };

  template <typename _Graph, typename NodeFilterMap, typename EdgeFilterMap>
  class SubGraphAdaptorBase<_Graph, NodeFilterMap, EdgeFilterMap, false> 
    : public GraphAdaptorBase<_Graph> {
  public:
    typedef _Graph Graph;
    typedef SubGraphAdaptorBase Adaptor;
    typedef GraphAdaptorBase<_Graph> Parent;
  protected:
    NodeFilterMap* node_filter_map;
    EdgeFilterMap* edge_filter_map;
    SubGraphAdaptorBase() : Parent(), 
			    node_filter_map(0), edge_filter_map(0) { }

    void setNodeFilterMap(NodeFilterMap& _node_filter_map) {
      node_filter_map=&_node_filter_map;
    }
    void setEdgeFilterMap(EdgeFilterMap& _edge_filter_map) {
      edge_filter_map=&_edge_filter_map;
    }

  public:

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;

    void first(Node& i) const { 
      Parent::first(i); 
      while (i!=INVALID && !(*node_filter_map)[i]) Parent::next(i); 
    }

    void first(Edge& i) const { 
      Parent::first(i); 
      while (i!=INVALID && !(*edge_filter_map)[i]) Parent::next(i); 
    }

    void firstIn(Edge& i, const Node& n) const { 
      Parent::firstIn(i, n); 
      while (i!=INVALID && !(*edge_filter_map)[i]) Parent::nextIn(i); 
    }

    void firstOut(Edge& i, const Node& n) const { 
      Parent::firstOut(i, n); 
      while (i!=INVALID && !(*edge_filter_map)[i]) Parent::nextOut(i); 
    }

    void next(Node& i) const { 
      Parent::next(i); 
      while (i!=INVALID && !(*node_filter_map)[i]) Parent::next(i); 
    }
    void next(Edge& i) const { 
      Parent::next(i); 
      while (i!=INVALID && !(*edge_filter_map)[i]) Parent::next(i); 
    }
    void nextIn(Edge& i) const { 
      Parent::nextIn(i); 
      while (i!=INVALID && !(*edge_filter_map)[i]) Parent::nextIn(i); 
    }

    void nextOut(Edge& i) const { 
      Parent::nextOut(i); 
      while (i!=INVALID && !(*edge_filter_map)[i]) Parent::nextOut(i); 
    }

    ///\e

    /// This function hides \c n in the graph, i.e. the iteration 
    /// jumps over it. This is done by simply setting the value of \c n  
    /// to be false in the corresponding node-map.
    void hide(const Node& n) const { node_filter_map->set(n, false); }

    ///\e

    /// This function hides \c e in the graph, i.e. the iteration 
    /// jumps over it. This is done by simply setting the value of \c e  
    /// to be false in the corresponding edge-map.
    void hide(const Edge& e) const { edge_filter_map->set(e, false); }

    ///\e

    /// The value of \c n is set to be true in the node-map which stores 
    /// hide information. If \c n was hidden previuosly, then it is shown 
    /// again
     void unHide(const Node& n) const { node_filter_map->set(n, true); }

    ///\e

    /// The value of \c e is set to be true in the edge-map which stores 
    /// hide information. If \c e was hidden previuosly, then it is shown 
    /// again
    void unHide(const Edge& e) const { edge_filter_map->set(e, true); }

    /// Returns true if \c n is hidden.
    
    ///\e
    ///
    bool hidden(const Node& n) const { return !(*node_filter_map)[n]; }

    /// Returns true if \c n is hidden.
    
    ///\e
    ///
    bool hidden(const Edge& e) const { return !(*edge_filter_map)[e]; }

    typedef False NodeNumTag;
    typedef False EdgeNumTag;

    typedef FindEdgeTagIndicator<Graph> FindEdgeTag;
    Edge findEdge(const Node& source, const Node& target, 
		  const Edge& prev = INVALID) {
      if (!(*node_filter_map)[source] || !(*node_filter_map)[target]) {
        return INVALID;
      }
      Edge edge = Parent::findEdge(source, target, prev);
      while (edge != INVALID && !(*edge_filter_map)[edge]) {
        edge = Parent::findEdge(source, target, edge);
      }
      return edge;
    }

    template <typename _Value>
    class NodeMap 
      : public SubMapExtender<Adaptor, 
                              typename Parent::template NodeMap<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      //typedef SubMapExtender<Adaptor, typename Parent::
      //                       template NodeMap<_Value> > Parent;
    
      NodeMap(const Graph& g) 
	: Parent(g) {}
      NodeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
      NodeMap& operator=(const NodeMap& cmap) {
	return operator=<NodeMap>(cmap);
      }
    
      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };

    template <typename _Value>
    class EdgeMap 
      : public SubMapExtender<Adaptor, 
                              typename Parent::template EdgeMap<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      //typedef SubMapExtender<Adaptor, typename Parent::
      //                       template EdgeMap<_Value> > Parent;
    
      EdgeMap(const Graph& g) 
	: Parent(g) {}
      EdgeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
      EdgeMap& operator=(const EdgeMap& cmap) {
	return operator=<EdgeMap>(cmap);
      }
    
      template <typename CMap>
      EdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };

  };

  /// \ingroup graph_adaptors
  ///
  /// \brief A graph adaptor for hiding nodes and edges from a graph.
  /// 
  /// SubGraphAdaptor shows the graph with filtered node-set and 
  /// edge-set. If the \c checked parameter is true then it filters the edgeset
  /// to do not get invalid edges without source or target.
  /// Let \f$ G=(V, A) \f$ be a directed graph
  /// and suppose that the graph instance \c g of type ListGraph
  /// implements \f$ G \f$.
  /// Let moreover \f$ b_V \f$ and \f$ b_A \f$ be bool-valued functions resp.
  /// on the node-set and edge-set.
  /// SubGraphAdaptor<...>::NodeIt iterates 
  /// on the node-set \f$ \{v\in V : b_V(v)=true\} \f$ and 
  /// SubGraphAdaptor<...>::EdgeIt iterates 
  /// on the edge-set \f$ \{e\in A : b_A(e)=true\} \f$. Similarly, 
  /// SubGraphAdaptor<...>::OutEdgeIt and
  /// SubGraphAdaptor<...>::InEdgeIt iterates 
  /// only on edges leaving and entering a specific node which have true value.
  /// 
  /// If the \c checked template parameter is false then we have to note that 
  /// the node-iterator cares only the filter on the node-set, and the 
  /// edge-iterator cares only the filter on the edge-set.
  /// This way the edge-map
  /// should filter all edges which's source or target is filtered by the 
  /// node-filter.
  ///\code
  /// typedef ListGraph Graph;
  /// Graph g;
  /// typedef Graph::Node Node;
  /// typedef Graph::Edge Edge;
  /// Node u=g.addNode(); //node of id 0
  /// Node v=g.addNode(); //node of id 1
  /// Node e=g.addEdge(u, v); //edge of id 0
  /// Node f=g.addEdge(v, u); //edge of id 1
  /// Graph::NodeMap<bool> nm(g, true);
  /// nm.set(u, false);
  /// Graph::EdgeMap<bool> em(g, true);
  /// em.set(e, false);
  /// typedef SubGraphAdaptor<Graph, Graph::NodeMap<bool>, Graph::EdgeMap<bool> > SubGA;
  /// SubGA ga(g, nm, em);
  /// for (SubGA::NodeIt n(ga); n!=INVALID; ++n) std::cout << g.id(n) << std::endl;
  /// std::cout << ":-)" << std::endl;
  /// for (SubGA::EdgeIt e(ga); e!=INVALID; ++e) std::cout << g.id(e) << std::endl;
  ///\endcode
  /// The output of the above code is the following.
  ///\code
  /// 1
  /// :-)
  /// 1
  ///\endcode
  /// Note that \c n is of type \c SubGA::NodeIt, but it can be converted to
  /// \c Graph::Node that is why \c g.id(n) can be applied.
  /// 
  /// For other examples see also the documentation of NodeSubGraphAdaptor and 
  /// EdgeSubGraphAdaptor.
  /// 
  /// \author Marton Makai

  template<typename _Graph, typename NodeFilterMap, 
	   typename EdgeFilterMap, bool checked = true>
  class SubGraphAdaptor : 
    public GraphAdaptorExtender<
    SubGraphAdaptorBase<_Graph, NodeFilterMap, EdgeFilterMap, checked> > {
  public:
    typedef _Graph Graph;
    typedef GraphAdaptorExtender< SubGraphAdaptorBase<_Graph, NodeFilterMap, 
                                                      EdgeFilterMap, checked> >
    Parent;

  protected:
    SubGraphAdaptor() { }
  public:

    SubGraphAdaptor(_Graph& _graph, NodeFilterMap& _node_filter_map, 
		    EdgeFilterMap& _edge_filter_map) { 
      setGraph(_graph);
      setNodeFilterMap(_node_filter_map);
      setEdgeFilterMap(_edge_filter_map);
    }

  };

  /// \brief Just gives back a sub graph adaptor
  ///
  /// Just gives back a sub graph adaptor
  template<typename Graph, typename NodeFilterMap, typename EdgeFilterMap>
  SubGraphAdaptor<const Graph, NodeFilterMap, EdgeFilterMap>
  subGraphAdaptor(const Graph& graph, 
                   NodeFilterMap& nfm, EdgeFilterMap& efm) {
    return SubGraphAdaptor<const Graph, NodeFilterMap, EdgeFilterMap>
      (graph, nfm, efm);
  }

  template<typename Graph, typename NodeFilterMap, typename EdgeFilterMap>
  SubGraphAdaptor<const Graph, const NodeFilterMap, EdgeFilterMap>
  subGraphAdaptor(const Graph& graph, 
                   NodeFilterMap& nfm, EdgeFilterMap& efm) {
    return SubGraphAdaptor<const Graph, const NodeFilterMap, EdgeFilterMap>
      (graph, nfm, efm);
  }

  template<typename Graph, typename NodeFilterMap, typename EdgeFilterMap>
  SubGraphAdaptor<const Graph, NodeFilterMap, const EdgeFilterMap>
  subGraphAdaptor(const Graph& graph, 
                   NodeFilterMap& nfm, EdgeFilterMap& efm) {
    return SubGraphAdaptor<const Graph, NodeFilterMap, const EdgeFilterMap>
      (graph, nfm, efm);
  }

  template<typename Graph, typename NodeFilterMap, typename EdgeFilterMap>
  SubGraphAdaptor<const Graph, const NodeFilterMap, const EdgeFilterMap>
  subGraphAdaptor(const Graph& graph, 
                   NodeFilterMap& nfm, EdgeFilterMap& efm) {
    return SubGraphAdaptor<const Graph, const NodeFilterMap, 
      const EdgeFilterMap>(graph, nfm, efm);
  }



  ///\ingroup graph_adaptors
  ///
  ///\brief An adaptor for hiding nodes from a graph.
  ///
  ///An adaptor for hiding nodes from a graph.
  ///This adaptor specializes SubGraphAdaptor in the way that only
  ///the node-set 
  ///can be filtered. In usual case the checked parameter is true, we get the
  ///induced subgraph. But if the checked parameter is false then we can
  ///filter only isolated nodes.
  ///\author Marton Makai
  template<typename Graph, typename NodeFilterMap, bool checked = true>
  class NodeSubGraphAdaptor : 
    public SubGraphAdaptor<Graph, NodeFilterMap, 
			   ConstMap<typename Graph::Edge,bool>, checked> {
  public:

    typedef SubGraphAdaptor<Graph, NodeFilterMap, 
			    ConstMap<typename Graph::Edge,bool>, checked > 
    Parent;

  protected:
    ConstMap<typename Graph::Edge, bool> const_true_map;

    NodeSubGraphAdaptor() : const_true_map(true) {
      Parent::setEdgeFilterMap(const_true_map);
    }

  public:

    NodeSubGraphAdaptor(Graph& _graph, NodeFilterMap& _node_filter_map) : 
      Parent(), const_true_map(true) { 
      Parent::setGraph(_graph);
      Parent::setNodeFilterMap(_node_filter_map);
      Parent::setEdgeFilterMap(const_true_map);
    }

  };


  /// \brief Just gives back a node sub graph adaptor
  ///
  /// Just gives back a node sub graph adaptor
  template<typename Graph, typename NodeFilterMap>
  NodeSubGraphAdaptor<const Graph, NodeFilterMap>
  nodeSubGraphAdaptor(const Graph& graph, NodeFilterMap& nfm) {
    return NodeSubGraphAdaptor<const Graph, NodeFilterMap>(graph, nfm);
  }

  template<typename Graph, typename NodeFilterMap>
  NodeSubGraphAdaptor<const Graph, const NodeFilterMap>
  nodeSubGraphAdaptor(const Graph& graph, const NodeFilterMap& nfm) {
    return NodeSubGraphAdaptor<const Graph, const NodeFilterMap>(graph, nfm);
  }

  ///\ingroup graph_adaptors
  ///
  ///\brief An adaptor for hiding edges from a graph.
  ///
  ///An adaptor for hiding edges from a graph.
  ///This adaptor specializes SubGraphAdaptor in the way that
  ///only the edge-set 
  ///can be filtered. The usefulness of this adaptor is demonstrated in the 
  ///problem of searching a maximum number of edge-disjoint shortest paths 
  ///between 
  ///two nodes \c s and \c t. Shortest here means being shortest w.r.t. 
  ///non-negative edge-lengths. Note that 
  ///the comprehension of the presented solution 
  ///need's some elementary knowledge from combinatorial optimization. 
  ///
  ///If a single shortest path is to be 
  ///searched between \c s and \c t, then this can be done easily by 
  ///applying the Dijkstra algorithm. What happens, if a maximum number of 
  ///edge-disjoint shortest paths is to be computed. It can be proved that an 
  ///edge can be in a shortest path if and only
  ///if it is tight with respect to 
  ///the potential function computed by Dijkstra.
  ///Moreover, any path containing 
  ///only such edges is a shortest one.
  ///Thus we have to compute a maximum number 
  ///of edge-disjoint paths between \c s and \c t in
  ///the graph which has edge-set 
  ///all the tight edges. The computation will be demonstrated
  ///on the following 
  ///graph, which is read from the dimacs file \c sub_graph_adaptor_demo.dim. 
  ///The full source code is available in \ref sub_graph_adaptor_demo.cc. 
  ///If you are interested in more demo programs, you can use 
  ///\ref dim_to_dot.cc to generate .dot files from dimacs files. 
  ///The .dot file of the following figure was generated by  
  ///the demo program \ref dim_to_dot.cc.
  ///
  ///\dot
  ///digraph lemon_dot_example {
  ///node [ shape=ellipse, fontname=Helvetica, fontsize=10 ];
  ///n0 [ label="0 (s)" ];
  ///n1 [ label="1" ];
  ///n2 [ label="2" ];
  ///n3 [ label="3" ];
  ///n4 [ label="4" ];
  ///n5 [ label="5" ];
  ///n6 [ label="6 (t)" ];
  ///edge [ shape=ellipse, fontname=Helvetica, fontsize=10 ];
  ///n5 ->  n6 [ label="9, length:4" ];
  ///n4 ->  n6 [ label="8, length:2" ];
  ///n3 ->  n5 [ label="7, length:1" ];
  ///n2 ->  n5 [ label="6, length:3" ];
  ///n2 ->  n6 [ label="5, length:5" ];
  ///n2 ->  n4 [ label="4, length:2" ];
  ///n1 ->  n4 [ label="3, length:3" ];
  ///n0 ->  n3 [ label="2, length:1" ];
  ///n0 ->  n2 [ label="1, length:2" ];
  ///n0 ->  n1 [ label="0, length:3" ];
  ///}
  ///\enddot
  ///
  ///\code
  ///Graph g;
  ///Node s, t;
  ///LengthMap length(g);
  ///
  ///readDimacs(std::cin, g, length, s, t);
  ///
  ///cout << "edges with lengths (of form id, source--length->target): " << endl;
  ///for(EdgeIt e(g); e!=INVALID; ++e) 
  ///  cout << g.id(e) << ", " << g.id(g.source(e)) << "--" 
  ///       << length[e] << "->" << g.id(g.target(e)) << endl;
  ///
  ///cout << "s: " << g.id(s) << " t: " << g.id(t) << endl;
  ///\endcode
  ///Next, the potential function is computed with Dijkstra.
  ///\code
  ///typedef Dijkstra<Graph, LengthMap> Dijkstra;
  ///Dijkstra dijkstra(g, length);
  ///dijkstra.run(s);
  ///\endcode
  ///Next, we consrtruct a map which filters the edge-set to the tight edges.
  ///\code
  ///typedef TightEdgeFilterMap<Graph, const Dijkstra::DistMap, LengthMap> 
  ///  TightEdgeFilter;
  ///TightEdgeFilter tight_edge_filter(g, dijkstra.distMap(), length);
  ///
  ///typedef EdgeSubGraphAdaptor<Graph, TightEdgeFilter> SubGA;
  ///SubGA ga(g, tight_edge_filter);
  ///\endcode
  ///Then, the maximum nimber of edge-disjoint \c s-\c t paths are computed 
  ///with a max flow algorithm Preflow.
  ///\code
  ///ConstMap<Edge, int> const_1_map(1);
  ///Graph::EdgeMap<int> flow(g, 0);
  ///
  ///Preflow<SubGA, ConstMap<Edge, int>, Graph::EdgeMap<int> > 
  ///  preflow(ga, const_1_map, s, t);
  ///preflow.run();
  ///\endcode
  ///Last, the output is:
  ///\code  
  ///cout << "maximum number of edge-disjoint shortest path: " 
  ///     << preflow.flowValue() << endl;
  ///cout << "edges of the maximum number of edge-disjoint shortest s-t paths: " 
  ///     << endl;
  ///for(EdgeIt e(g); e!=INVALID; ++e) 
  ///  if (preflow.flow(e))
  ///    cout << " " << g.id(g.source(e)) << "--"
  ///         << length[e] << "->" << g.id(g.target(e)) << endl;
  ///\endcode
  ///The program has the following (expected :-)) output:
  ///\code
  ///edges with lengths (of form id, source--length->target):
  /// 9, 5--4->6
  /// 8, 4--2->6
  /// 7, 3--1->5
  /// 6, 2--3->5
  /// 5, 2--5->6
  /// 4, 2--2->4
  /// 3, 1--3->4
  /// 2, 0--1->3
  /// 1, 0--2->2
  /// 0, 0--3->1
  ///s: 0 t: 6
  ///maximum number of edge-disjoint shortest path: 2
  ///edges of the maximum number of edge-disjoint shortest s-t paths:
  /// 9, 5--4->6
  /// 8, 4--2->6
  /// 7, 3--1->5
  /// 4, 2--2->4
  /// 2, 0--1->3
  /// 1, 0--2->2
  ///\endcode
  ///
  ///\author Marton Makai
  template<typename Graph, typename EdgeFilterMap>
  class EdgeSubGraphAdaptor : 
    public SubGraphAdaptor<Graph, ConstMap<typename Graph::Node,bool>, 
			   EdgeFilterMap, false> {
  public:
    typedef SubGraphAdaptor<Graph, ConstMap<typename Graph::Node,bool>, 
			    EdgeFilterMap, false> Parent;
  protected:
    ConstMap<typename Graph::Node, bool> const_true_map;

    EdgeSubGraphAdaptor() : const_true_map(true) {
      Parent::setNodeFilterMap(const_true_map);
    }

  public:

    EdgeSubGraphAdaptor(Graph& _graph, EdgeFilterMap& _edge_filter_map) : 
      Parent(), const_true_map(true) { 
      Parent::setGraph(_graph);
      Parent::setNodeFilterMap(const_true_map);
      Parent::setEdgeFilterMap(_edge_filter_map);
    }

  };

  /// \brief Just gives back an edge sub graph adaptor
  ///
  /// Just gives back an edge sub graph adaptor
  template<typename Graph, typename EdgeFilterMap>
  EdgeSubGraphAdaptor<const Graph, EdgeFilterMap>
  edgeSubGraphAdaptor(const Graph& graph, EdgeFilterMap& efm) {
    return EdgeSubGraphAdaptor<const Graph, EdgeFilterMap>(graph, efm);
  }

  template<typename Graph, typename EdgeFilterMap>
  EdgeSubGraphAdaptor<const Graph, const EdgeFilterMap>
  edgeSubGraphAdaptor(const Graph& graph, const EdgeFilterMap& efm) {
    return EdgeSubGraphAdaptor<const Graph, const EdgeFilterMap>(graph, efm);
  }

  template <typename _Graph>
  class UndirGraphAdaptorBase : 
    public UndirGraphExtender<GraphAdaptorBase<_Graph> > {
  public:
    typedef _Graph Graph;
    typedef UndirGraphAdaptorBase Adaptor;
    typedef UndirGraphExtender<GraphAdaptorBase<_Graph> > Parent;

  protected:

    UndirGraphAdaptorBase() : Parent() {}

  public:

    typedef typename Parent::UEdge UEdge;
    typedef typename Parent::Edge Edge;

  private:
    
    template <typename _Value>
    class EdgeMapBase {
    private:
      
      typedef typename _Graph::template EdgeMap<_Value> MapImpl;
      
    public:

      typedef typename MapTraits<MapImpl>::ReferenceMapTag ReferenceMapTag;

      typedef _Value Value;
      typedef Edge Key;
      
      EdgeMapBase(const Adaptor& adaptor) :
	forward_map(*adaptor.graph), backward_map(*adaptor.graph) {}

      EdgeMapBase(const Adaptor& adaptor, const Value& v) 
        : forward_map(*adaptor.graph, v), backward_map(*adaptor.graph, v) {}
      
      void set(const Edge& e, const Value& a) { 
	if (Parent::direction(e)) {
	  forward_map.set(e, a); 
        } else { 
	  backward_map.set(e, a);
        } 
      }

      typename MapTraits<MapImpl>::ConstReturnValue operator[](Edge e) const { 
	if (Parent::direction(e)) {
	  return forward_map[e]; 
	} else { 
	  return backward_map[e]; 
        }
      }

      typename MapTraits<MapImpl>::ReturnValue operator[](Edge e) { 
	if (Parent::direction(e)) {
	  return forward_map[e]; 
	} else { 
	  return backward_map[e]; 
        }
      }

    protected:

      MapImpl forward_map, backward_map; 

    };

  public:

    template <typename _Value>
    class EdgeMap 
      : public SubMapExtender<Adaptor, EdgeMapBase<_Value> > 
    {
    public:
      typedef Adaptor Graph;
      typedef SubMapExtender<Adaptor, EdgeMapBase<_Value> > Parent;
    
      EdgeMap(const Graph& g) 
	: Parent(g) {}
      EdgeMap(const Graph& g, const _Value& v) 
	: Parent(g, v) {}
    
      EdgeMap& operator=(const EdgeMap& cmap) {
	return operator=<EdgeMap>(cmap);
      }
    
      template <typename CMap>
      EdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
	return *this;
      }
    };
        
    template <typename _Value>
    class UEdgeMap : public Graph::template EdgeMap<_Value> {
    public:
      
      typedef typename Graph::template EdgeMap<_Value> Parent;
      
      explicit UEdgeMap(const Adaptor& ga) 
	: Parent(*ga.graph) {}

      UEdgeMap(const Adaptor& ga, const _Value& value)
	: Parent(*ga.graph, value) {}

      UEdgeMap& operator=(const UEdgeMap& cmap) {
        return operator=<UEdgeMap>(cmap);
      }

      template <typename CMap>
      UEdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }

    };
      
  };

  template <typename _Graph, typename Enable = void>
  class AlterableUndirGraphAdaptor 
    : public UGraphAdaptorExtender<UndirGraphAdaptorBase<_Graph> > {
  public:
    typedef UGraphAdaptorExtender<UndirGraphAdaptorBase<_Graph> > Parent;
    
  protected:

    AlterableUndirGraphAdaptor() : Parent() {}

  public:

    typedef typename Parent::EdgeNotifier UEdgeNotifier;
    typedef InvalidType EdgeNotifier;

  };

  template <typename _Graph>
  class AlterableUndirGraphAdaptor<
    _Graph, 
    typename enable_if<typename _Graph::EdgeNotifier::Notifier>::type > 
    : public UGraphAdaptorExtender<UndirGraphAdaptorBase<_Graph> > {
  public:

    typedef UGraphAdaptorExtender<UndirGraphAdaptorBase<_Graph> > Parent;
    typedef _Graph Graph;
    typedef typename _Graph::Edge GraphEdge;
    
  protected:

    AlterableUndirGraphAdaptor() 
      : Parent(), edge_notifier(*this), edge_notifier_proxy(*this) {}

    void setGraph(_Graph& g) {
      Parent::setGraph(g);
      edge_notifier_proxy.setNotifier(g.notifier(GraphEdge()));
    }

  public:

    ~AlterableUndirGraphAdaptor() {
      edge_notifier.clear();
    }

    typedef typename Parent::UEdge UEdge;
    typedef typename Parent::Edge Edge;

    typedef typename Parent::EdgeNotifier UEdgeNotifier;

    using Parent::notifier;

    typedef AlterationNotifier<AlterableUndirGraphAdaptor, 
                               Edge> EdgeNotifier;
    EdgeNotifier& notifier(Edge) const { return edge_notifier; }

  protected:

    class NotifierProxy : public Graph::EdgeNotifier::ObserverBase {
    public:

      typedef typename Graph::EdgeNotifier::ObserverBase Parent;
      typedef AlterableUndirGraphAdaptor AdaptorBase;
      
      NotifierProxy(const AdaptorBase& _adaptor)
        : Parent(), adaptor(&_adaptor) {
      }

      virtual ~NotifierProxy() {
        if (Parent::attached()) {
          Parent::detach();
        }
      }

      void setNotifier(typename Graph::EdgeNotifier& nf) {
        Parent::attach(nf);
      }

      
    protected:

      virtual void add(const GraphEdge& ge) {
        std::vector<Edge> edges;
        edges.push_back(AdaptorBase::Parent::direct(ge, true));
        edges.push_back(AdaptorBase::Parent::direct(ge, false));
        adaptor->notifier(Edge()).add(edges);
      }
      virtual void add(const std::vector<GraphEdge>& ge) {
        std::vector<Edge> edges;
        for (int i = 0; i < int(ge.size()); ++i) { 
          edges.push_back(AdaptorBase::Parent::direct(ge[i], true));
          edges.push_back(AdaptorBase::Parent::direct(ge[i], false));
        }
        adaptor->notifier(Edge()).add(edges);
      }
      virtual void erase(const GraphEdge& ge) {
        std::vector<Edge> edges;
        edges.push_back(AdaptorBase::Parent::direct(ge, true));
        edges.push_back(AdaptorBase::Parent::direct(ge, false));
        adaptor->notifier(Edge()).erase(edges);
      }
      virtual void erase(const std::vector<GraphEdge>& ge) {
        std::vector<Edge> edges;
        for (int i = 0; i < int(ge.size()); ++i) { 
          edges.push_back(AdaptorBase::Parent::direct(ge[i], true));
          edges.push_back(AdaptorBase::Parent::direct(ge[i], false));
        }
        adaptor->notifier(Edge()).erase(edges);
      }
      virtual void build() {
        adaptor->notifier(Edge()).build();
      }
      virtual void clear() {
        adaptor->notifier(Edge()).clear();
      }

      const AdaptorBase* adaptor;
    };


    mutable EdgeNotifier edge_notifier;
    NotifierProxy edge_notifier_proxy;

  };


  ///\ingroup graph_adaptors
  ///
  /// \brief An undirected graph is made from a directed graph by an adaptor
  ///
  /// This adaptor makes an undirected graph from a directed
  /// graph. All edge of the underlying will be showed in the adaptor
  /// as an undirected edge. Let's see an informal example about using
  /// this adaptor:
  ///
  /// There is a network of the streets of a town. Of course there are
  /// some one-way street in the town hence the network is a directed
  /// one. There is a crazy driver who go oppositely in the one-way
  /// street without moral sense. Of course he can pass this streets
  /// slower than the regular way, in fact his speed is half of the
  /// normal speed. How long should he drive to get from a source
  /// point to the target? Let see the example code which calculate it:
  ///
  ///\code
  /// typedef UndirGraphAdaptor<Graph> UGraph;
  /// UGraph ugraph(graph);
  ///
  /// typedef SimpleMap<LengthMap> FLengthMap;
  /// FLengthMap flength(length);
  ///
  /// typedef ScaleMap<LengthMap> RLengthMap;
  /// RLengthMap rlength(length, 2.0);
  ///
  /// typedef UGraph::CombinedEdgeMap<FLengthMap, RLengthMap > ULengthMap;
  /// ULengthMap ulength(flength, rlength);
  /// 
  /// Dijkstra<UGraph, ULengthMap> dijkstra(ugraph, ulength);
  /// std::cout << "Driving time : " << dijkstra.run(src, trg) << std::endl;
  ///\endcode
  ///
  /// The combined edge map makes the length map for the undirected
  /// graph. It is created from a forward and reverse map. The forward
  /// map is created from the original length map with a SimpleMap
  /// adaptor which just makes a read-write map from the reference map
  /// i.e. it forgets that it can be return reference to values. The
  /// reverse map is just the scaled original map with the ScaleMap
  /// adaptor. The combination solves that passing the reverse way
  /// takes double time than the original. To get the driving time we
  /// run the dijkstra algorithm on the undirected graph.
  ///
  /// \author Marton Makai and Balazs Dezso
  template<typename _Graph>
  class UndirGraphAdaptor : public AlterableUndirGraphAdaptor<_Graph> {
  public:
    typedef _Graph Graph;
    typedef AlterableUndirGraphAdaptor<_Graph> Parent;
  protected:
    UndirGraphAdaptor() { }
  public:

    /// \brief Constructor
    ///
    /// Constructor
    UndirGraphAdaptor(_Graph& _graph) { 
      setGraph(_graph);
    }

    /// \brief EdgeMap combined from two original EdgeMap
    ///
    /// This class adapts two original graph EdgeMap to
    /// get an edge map on the adaptor.
    template <typename _ForwardMap, typename _BackwardMap>
    class CombinedEdgeMap {
    public:
      
      typedef _ForwardMap ForwardMap;
      typedef _BackwardMap BackwardMap;

      typedef typename MapTraits<ForwardMap>::ReferenceMapTag ReferenceMapTag;

      typedef typename ForwardMap::Value Value;
      typedef typename Parent::Edge Key;

      /// \brief Constructor      
      ///
      /// Constructor      
      CombinedEdgeMap() : forward_map(0), backward_map(0) {}

      /// \brief Constructor      
      ///
      /// Constructor      
      CombinedEdgeMap(ForwardMap& _forward_map, BackwardMap& _backward_map) 
        : forward_map(&_forward_map), backward_map(&_backward_map) {}
      

      /// \brief Sets the value associated with a key.
      ///
      /// Sets the value associated with a key.
      void set(const Key& e, const Value& a) { 
	if (Parent::direction(e)) {
	  forward_map->set(e, a); 
        } else { 
	  backward_map->set(e, a);
        } 
      }

      /// \brief Returns the value associated with a key.
      ///
      /// Returns the value associated with a key.
      typename MapTraits<ForwardMap>::ConstReturnValue 
      operator[](const Key& e) const { 
	if (Parent::direction(e)) {
	  return (*forward_map)[e]; 
	} else { 
	  return (*backward_map)[e]; 
        }
      }

      /// \brief Returns the value associated with a key.
      ///
      /// Returns the value associated with a key.
      typename MapTraits<ForwardMap>::ReturnValue 
      operator[](const Key& e) { 
	if (Parent::direction(e)) {
	  return (*forward_map)[e]; 
	} else { 
	  return (*backward_map)[e]; 
        }
      }

      /// \brief Sets the forward map
      ///
      /// Sets the forward map
      void setForwardMap(ForwardMap& _forward_map) {
        forward_map = &_forward_map;
      }

      /// \brief Sets the backward map
      ///
      /// Sets the backward map
      void setBackwardMap(BackwardMap& _backward_map) {
        backward_map = &_backward_map;
      }

    protected:

      ForwardMap* forward_map;
      BackwardMap* backward_map; 

    };

  };

  /// \brief Just gives back an undir graph adaptor
  ///
  /// Just gives back an undir graph adaptor
  template<typename Graph>
  UndirGraphAdaptor<const Graph>
  undirGraphAdaptor(const Graph& graph) {
    return UndirGraphAdaptor<const Graph>(graph);
  }

  template<typename Graph, typename Number,  
           typename CapacityMap, typename FlowMap, 
           typename Tol = Tolerance<Number> >
  class ResForwardFilter {
    const CapacityMap* capacity;
    const FlowMap* flow;
    Tol tolerance;
  public:
    typedef typename Graph::Edge Key;
    typedef bool Value;

    ResForwardFilter(const CapacityMap& _capacity, const FlowMap& _flow,
                     const Tol& _tolerance = Tol()) 
      : capacity(&_capacity), flow(&_flow), tolerance(_tolerance) { }

    ResForwardFilter(const Tol& _tolerance) 
      : capacity(0), flow(0), tolerance(_tolerance)  { }

    void setCapacity(const CapacityMap& _capacity) { capacity = &_capacity; }
    void setFlow(const FlowMap& _flow) { flow = &_flow; }

    bool operator[](const typename Graph::Edge& e) const {
      return tolerance.positive((*capacity)[e] - (*flow)[e]);
    }
  };

  template<typename Graph, typename Number,
	   typename CapacityMap, typename FlowMap,
           typename Tol = Tolerance<Number> >
  class ResBackwardFilter {
    const CapacityMap* capacity;
    const FlowMap* flow;
    Tol tolerance;
  public:
    typedef typename Graph::Edge Key;
    typedef bool Value;

    ResBackwardFilter(const CapacityMap& _capacity, const FlowMap& _flow,
                      const Tol& _tolerance = Tol())
      : capacity(&_capacity), flow(&_flow), tolerance(_tolerance) { }
    ResBackwardFilter(const Tol& _tolerance = Tol())
      : capacity(0), flow(0), tolerance(_tolerance) { }
    void setCapacity(const CapacityMap& _capacity) { capacity = &_capacity; }
    void setFlow(const FlowMap& _flow) { flow = &_flow; }
    bool operator[](const typename Graph::Edge& e) const {
      return tolerance.positive((*flow)[e]);
    }
  };

  
  ///\ingroup graph_adaptors
  ///
  ///\brief An adaptor for composing the residual
  ///graph for directed flow and circulation problems.
  ///
  ///An adaptor for composing the residual graph for directed flow and
  ///circulation problems.  Let \f$ G=(V, A) \f$ be a directed graph
  ///and let \f$ F \f$ be a number type. Let moreover \f$ f,c:A\to F \f$,
  ///be functions on the edge-set.
  ///
  ///In the appications of ResGraphAdaptor, \f$ f \f$ usually stands
  ///for a flow and \f$ c \f$ for a capacity function.  Suppose that a
  ///graph instange \c g of type \c ListGraph implements \f$ G \f$.
  ///
  ///\code 
  ///  ListGraph g;
  ///\endcode 
  ///
  ///Then ResGraphAdaptor implements the graph structure with node-set
  /// \f$ V \f$ and edge-set \f$ A_{forward}\cup A_{backward} \f$,
  ///where \f$ A_{forward}=\{uv : uv\in A, f(uv)<c(uv)\} \f$ and 
  /// \f$ A_{backward}=\{vu : uv\in A, f(uv)>0\} \f$, i.e. the so called
  ///residual graph.  When we take the union 
  /// \f$ A_{forward}\cup A_{backward} \f$, multilicities are counted, i.e. 
  ///if an edge is in both \f$ A_{forward} \f$ and \f$ A_{backward} \f$, 
  ///then in the adaptor it appears twice. The following code shows how 
  ///such an instance can be constructed.
  ///
  ///\code 
  ///  typedef ListGraph Graph; 
  ///  Graph::EdgeMap<int> f(g);
  ///  Graph::EdgeMap<int> c(g); 
  ///  ResGraphAdaptor<Graph, int, Graph::EdgeMap<int>, Graph::EdgeMap<int> > ga(g); 
  ///\endcode
  ///\author Marton Makai
  ///
  template<typename Graph, typename Number, 
	   typename CapacityMap, typename FlowMap,
           typename Tol = Tolerance<Number> >
  class ResGraphAdaptor : 
    public EdgeSubGraphAdaptor< 
    UndirGraphAdaptor<const Graph>, 
    typename UndirGraphAdaptor<const Graph>::template CombinedEdgeMap<
    ResForwardFilter<const Graph, Number, CapacityMap, FlowMap>,  
    ResBackwardFilter<const Graph, Number, CapacityMap, FlowMap> > > {
  public:

    typedef UndirGraphAdaptor<const Graph> UGraph;

    typedef ResForwardFilter<const Graph, Number, CapacityMap, FlowMap> 
    ForwardFilter;

    typedef ResBackwardFilter<const Graph, Number, CapacityMap, FlowMap> 
    BackwardFilter;

    typedef typename UGraph::
    template CombinedEdgeMap<ForwardFilter, BackwardFilter>
    EdgeFilter;

    typedef EdgeSubGraphAdaptor<UGraph, EdgeFilter> Parent;

  protected:

    const CapacityMap* capacity;
    FlowMap* flow;

    UGraph ugraph;
    ForwardFilter forward_filter;
    BackwardFilter backward_filter;
    EdgeFilter edge_filter;

    void setCapacityMap(const CapacityMap& _capacity) {
      capacity=&_capacity;
      forward_filter.setCapacity(_capacity);
      backward_filter.setCapacity(_capacity);
    }

    void setFlowMap(FlowMap& _flow) {
      flow=&_flow;
      forward_filter.setFlow(_flow);
      backward_filter.setFlow(_flow);
    }

  public:

    /// \brief Constructor of the residual graph.
    ///
    /// Constructor of the residual graph. The parameters are the graph type,
    /// the flow map, the capacity map and a tolerance object.
    ResGraphAdaptor(const Graph& _graph, const CapacityMap& _capacity, 
                    FlowMap& _flow, const Tol& _tolerance = Tol()) 
      : Parent(), capacity(&_capacity), flow(&_flow), ugraph(_graph),
        forward_filter(_capacity, _flow, _tolerance), 
        backward_filter(_capacity, _flow, _tolerance),
        edge_filter(forward_filter, backward_filter)
    {
      Parent::setGraph(ugraph);
      Parent::setEdgeFilterMap(edge_filter);
    }

    typedef typename Parent::Edge Edge;

    /// \brief Gives back the residual capacity of the edge.
    ///
    /// Gives back the residual capacity of the edge.
    Number rescap(const Edge& edge) const {
      if (UGraph::direction(edge)) {
        return (*capacity)[edge]-(*flow)[edge]; 
      } else {
        return (*flow)[edge];
      }
    } 

    /// \brief Augment on the given edge in the residual graph.
    ///
    /// Augment on the given edge in the residual graph. It increase
    /// or decrease the flow on the original edge depend on the direction
    /// of the residual edge.
    void augment(const Edge& e, Number a) const {
      if (UGraph::direction(e)) {
        flow->set(e, (*flow)[e] + a);
      } else {  
        flow->set(e, (*flow)[e] - a);
      }
    }

    /// \brief Returns the direction of the edge.
    ///
    /// Returns true when the edge is same oriented as the original edge.
    static bool forward(const Edge& e) {
      return UGraph::direction(e);
    }

    /// \brief Returns the direction of the edge.
    ///
    /// Returns true when the edge is opposite oriented as the original edge.
    static bool backward(const Edge& e) {
      return !UGraph::direction(e);
    }

    /// \brief Gives back the forward oriented residual edge.
    ///
    /// Gives back the forward oriented residual edge.
    static Edge forward(const typename Graph::Edge& e) {
      return UGraph::direct(e, true);
    }

    /// \brief Gives back the backward oriented residual edge.
    ///
    /// Gives back the backward oriented residual edge.
    static Edge backward(const typename Graph::Edge& e) {
      return UGraph::direct(e, false);
    }

    /// \brief Residual capacity map.
    ///
    /// In generic residual graphs the residual capacity can be obtained 
    /// as a map. 
    class ResCap {
    protected:
      const ResGraphAdaptor* res_graph;
    public:
      typedef Number Value;
      typedef Edge Key;
      ResCap(const ResGraphAdaptor& _res_graph) 
        : res_graph(&_res_graph) {}
      
      Number operator[](const Edge& e) const {
        return res_graph->rescap(e);
      }
      
    };

  };



  template <typename _Graph, typename FirstOutEdgesMap>
  class ErasingFirstGraphAdaptorBase : public GraphAdaptorBase<_Graph> {
  public:
    typedef _Graph Graph;
    typedef GraphAdaptorBase<_Graph> Parent;
  protected:
    FirstOutEdgesMap* first_out_edges;
    ErasingFirstGraphAdaptorBase() : Parent(), 
				     first_out_edges(0) { }

    void setFirstOutEdgesMap(FirstOutEdgesMap& _first_out_edges) {
      first_out_edges=&_first_out_edges;
    }

  public:

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;

    void firstOut(Edge& i, const Node& n) const { 
      i=(*first_out_edges)[n];
    }

    void erase(const Edge& e) const {
      Node n=source(e);
      Edge f=e;
      Parent::nextOut(f);
      first_out_edges->set(n, f);
    }    
  };


  ///\ingroup graph_adaptors
  ///
  ///\brief For blocking flows.
  ///
  ///This graph adaptor is used for on-the-fly 
  ///Dinits blocking flow computations.
  ///For each node, an out-edge is stored which is used when the 
  ///\code
  ///OutEdgeIt& first(OutEdgeIt&, const Node&)
  ///\endcode
  ///is called. 
  ///
  ///\author Marton Makai
  ///
  template <typename _Graph, typename FirstOutEdgesMap>
  class ErasingFirstGraphAdaptor : 
    public GraphAdaptorExtender<
    ErasingFirstGraphAdaptorBase<_Graph, FirstOutEdgesMap> > {
  public:
    typedef _Graph Graph;
    typedef GraphAdaptorExtender<
      ErasingFirstGraphAdaptorBase<_Graph, FirstOutEdgesMap> > Parent;
    ErasingFirstGraphAdaptor(Graph& _graph, 
			     FirstOutEdgesMap& _first_out_edges) { 
      setGraph(_graph);
      setFirstOutEdgesMap(_first_out_edges);
    } 

  };

  /// \brief Base class for split graph adaptor
  ///
  /// Base class of split graph adaptor. In most case you do not need to
  /// use it directly but the documented member functions of this class can 
  /// be used with the SplitGraphAdaptor class.
  /// \sa SplitGraphAdaptor
  template <typename _Graph>
  class SplitGraphAdaptorBase 
    : public GraphAdaptorBase<const _Graph> {
  public:

    typedef _Graph Graph;

    typedef GraphAdaptorBase<const _Graph> Parent;

    typedef typename Graph::Node GraphNode;
    typedef typename Graph::Edge GraphEdge;

    class Node;
    class Edge;

    template <typename T> class NodeMap;
    template <typename T> class EdgeMap;
    

    class Node : public GraphNode {
      friend class SplitGraphAdaptorBase;
      template <typename T> friend class NodeMap;
    private:

      bool in_node;
      Node(GraphNode _node, bool _in_node)
	: GraphNode(_node), in_node(_in_node) {}
      
    public:

      Node() {}
      Node(Invalid) : GraphNode(INVALID), in_node(true) {}

      bool operator==(const Node& node) const {
	return GraphNode::operator==(node) && in_node == node.in_node;
      }
      
      bool operator!=(const Node& node) const {
	return !(*this == node);
      }
      
      bool operator<(const Node& node) const {
	return GraphNode::operator<(node) || 
	  (GraphNode::operator==(node) && in_node < node.in_node);
      }
    };

    class Edge {
      friend class SplitGraphAdaptorBase;
      template <typename T> friend class EdgeMap;
    private:
      typedef BiVariant<GraphEdge, GraphNode> EdgeImpl;

      explicit Edge(const GraphEdge& edge) : item(edge) {}
      explicit Edge(const GraphNode& node) : item(node) {}
      
      EdgeImpl item;

    public:
      Edge() {}
      Edge(Invalid) : item(GraphEdge(INVALID)) {}

      bool operator==(const Edge& edge) const {
        if (item.firstState()) {
          if (edge.item.firstState()) {
            return item.first() == edge.item.first();
          }
        } else {
          if (edge.item.secondState()) {
            return item.second() == edge.item.second();
          }
        }
        return false;
      }
      
      bool operator!=(const Edge& edge) const {
	return !(*this == edge);
      }
      
      bool operator<(const Edge& edge) const {
        if (item.firstState()) {
          if (edge.item.firstState()) {
            return item.first() < edge.item.first();
          }
          return false;
        } else {
          if (edge.item.secondState()) {
            return item.second() < edge.item.second();
          }
          return true;
        }
      }

      operator GraphEdge() const { return item.first(); }
      operator GraphNode() const { return item.second(); }

    };

    void first(Node& n) const {
      Parent::first(n);
      n.in_node = true;
    }

    void next(Node& n) const {
      if (n.in_node) {
	n.in_node = false;
      } else {
	n.in_node = true;
	Parent::next(n);
      }
    }

    void first(Edge& e) const {
      e.item.setSecond();
      Parent::first(e.item.second());
      if (e.item.second() == INVALID) {
        e.item.setFirst();
	Parent::first(e.item.first());
      }
    }

    void next(Edge& e) const {
      if (e.item.secondState()) {
	Parent::next(e.item.second());
        if (e.item.second() == INVALID) {
          e.item.setFirst();
          Parent::first(e.item.first());
        }
      } else {
	Parent::next(e.item.first());
      }      
    }

    void firstOut(Edge& e, const Node& n) const {
      if (n.in_node) {
        e.item.setSecond(n);
      } else {
        e.item.setFirst();
	Parent::firstOut(e.item.first(), n);
      }
    }

    void nextOut(Edge& e) const {
      if (!e.item.firstState()) {
	e.item.setFirst(INVALID);
      } else {
	Parent::nextOut(e.item.first());
      }      
    }

    void firstIn(Edge& e, const Node& n) const {
      if (!n.in_node) {
        e.item.setSecond(n);        
      } else {
        e.item.setFirst();
	Parent::firstIn(e.item.first(), n);
      }
    }

    void nextIn(Edge& e) const {
      if (!e.item.firstState()) {
	e.item.setFirst(INVALID);
      } else {
	Parent::nextIn(e.item.first());
      }
    }

    Node source(const Edge& e) const {
      if (e.item.firstState()) {
	return Node(Parent::source(e.item.first()), false);
      } else {
	return Node(e.item.second(), true);
      }
    }

    Node target(const Edge& e) const {
      if (e.item.firstState()) {
	return Node(Parent::target(e.item.first()), true);
      } else {
	return Node(e.item.second(), false);
      }
    }

    int id(const Node& n) const {
      return (Parent::id(n) << 1) | (n.in_node ? 0 : 1);
    }
    Node nodeFromId(int ix) const {
      return Node(Parent::nodeFromId(ix >> 1), (ix & 1) == 0);
    }
    int maxNodeId() const {
      return 2 * Parent::maxNodeId() + 1;
    }

    int id(const Edge& e) const {
      if (e.item.firstState()) {
        return Parent::id(e.item.first()) << 1;
      } else {
        return (Parent::id(e.item.second()) << 1) | 1;
      }
    }
    Edge edgeFromId(int ix) const {
      if ((ix & 1) == 0) {
        return Edge(Parent::edgeFromId(ix >> 1));
      } else {
        return Edge(Parent::nodeFromId(ix >> 1));
      }
    }
    int maxEdgeId() const {
      return std::max(Parent::maxNodeId() << 1, 
                      (Parent::maxEdgeId() << 1) | 1);
    }

    /// \brief Returns true when the node is in-node.
    ///
    /// Returns true when the node is in-node.
    static bool inNode(const Node& n) {
      return n.in_node;
    }

    /// \brief Returns true when the node is out-node.
    ///
    /// Returns true when the node is out-node.
    static bool outNode(const Node& n) {
      return !n.in_node;
    }

    /// \brief Returns true when the edge is edge in the original graph.
    ///
    /// Returns true when the edge is edge in the original graph.
    static bool origEdge(const Edge& e) {
      return e.item.firstState();
    }

    /// \brief Returns true when the edge binds an in-node and an out-node.
    ///
    /// Returns true when the edge binds an in-node and an out-node.
    static bool bindEdge(const Edge& e) {
      return e.item.secondState();
    }

    /// \brief Gives back the in-node created from the \c node.
    ///
    /// Gives back the in-node created from the \c node.
    static Node inNode(const GraphNode& n) {
      return Node(n, true);
    }

    /// \brief Gives back the out-node created from the \c node.
    ///
    /// Gives back the out-node created from the \c node.
    static Node outNode(const GraphNode& n) {
      return Node(n, false);
    }

    /// \brief Gives back the edge binds the two part of the node.
    /// 
    /// Gives back the edge binds the two part of the node.
    static Edge edge(const GraphNode& n) {
      return Edge(n);
    }

    /// \brief Gives back the edge of the original edge.
    /// 
    /// Gives back the edge of the original edge.
    static Edge edge(const GraphEdge& e) {
      return Edge(e);
    }

    typedef True NodeNumTag;

    int nodeNum() const {
      return  2 * countNodes(*Parent::graph);
    }

    typedef True EdgeNumTag;
    
    int edgeNum() const {
      return countEdges(*Parent::graph) + countNodes(*Parent::graph);
    }

    typedef True FindEdgeTag;

    Edge findEdge(const Node& u, const Node& v, 
		  const Edge& prev = INVALID) const {
      if (inNode(u)) {
        if (outNode(v)) {
          if (static_cast<const GraphNode&>(u) == 
              static_cast<const GraphNode&>(v) && prev == INVALID) {
            return Edge(u);
          }
        }
      } else {
        if (inNode(v)) {
          return Edge(findEdge(*Parent::graph, u, v, prev));
        }
      }
      return INVALID;
    }

    
    template <typename T>
    class NodeMap : public MapBase<Node, T> {
      typedef typename Parent::template NodeMap<T> NodeImpl;
    public:
      NodeMap(const SplitGraphAdaptorBase& _graph) 
	: inNodeMap(_graph), outNodeMap(_graph) {}
      NodeMap(const SplitGraphAdaptorBase& _graph, const T& t) 
	: inNodeMap(_graph, t), outNodeMap(_graph, t) {}
      NodeMap& operator=(const NodeMap& cmap) {
        return operator=<NodeMap>(cmap);
      }
      template <typename CMap>
      NodeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }

      void set(const Node& key, const T& val) {
	if (SplitGraphAdaptorBase::inNode(key)) { inNodeMap.set(key, val); }
	else {outNodeMap.set(key, val); }
      }
      
      typename MapTraits<NodeImpl>::ReturnValue 
      operator[](const Node& key) {
	if (SplitGraphAdaptorBase::inNode(key)) { return inNodeMap[key]; }
	else { return outNodeMap[key]; }
      }

      typename MapTraits<NodeImpl>::ConstReturnValue
      operator[](const Node& key) const {
	if (SplitGraphAdaptorBase::inNode(key)) { return inNodeMap[key]; }
	else { return outNodeMap[key]; }
      }

    private:
      NodeImpl inNodeMap, outNodeMap;
    };

    template <typename T>
    class EdgeMap : public MapBase<Edge, T> {
      typedef typename Parent::template EdgeMap<T> EdgeMapImpl;
      typedef typename Parent::template NodeMap<T> NodeMapImpl;
    public:

      EdgeMap(const SplitGraphAdaptorBase& _graph) 
	: edge_map(_graph), node_map(_graph) {}
      EdgeMap(const SplitGraphAdaptorBase& _graph, const T& t) 
	: edge_map(_graph, t), node_map(_graph, t) {}
      EdgeMap& operator=(const EdgeMap& cmap) {
        return operator=<EdgeMap>(cmap);
      }
      template <typename CMap>
      EdgeMap& operator=(const CMap& cmap) {
        Parent::operator=(cmap);
        return *this;
      }
      
      void set(const Edge& key, const T& val) {
	if (SplitGraphAdaptorBase::origEdge(key)) { 
          edge_map.set(key.item.first(), val); 
        } else {
          node_map.set(key.item.second(), val); 
        }
      }
      
      typename MapTraits<EdgeMapImpl>::ReturnValue
      operator[](const Edge& key) {
	if (SplitGraphAdaptorBase::origEdge(key)) { 
          return edge_map[key.item.first()];
        } else {
          return node_map[key.item.second()];
        }
      }

      typename MapTraits<EdgeMapImpl>::ConstReturnValue
      operator[](const Edge& key) const {
	if (SplitGraphAdaptorBase::origEdge(key)) { 
          return edge_map[key.item.first()];
        } else {
          return node_map[key.item.second()];
        }
      }

    private:
      typename Parent::template EdgeMap<T> edge_map;
      typename Parent::template NodeMap<T> node_map;
    };


  };

  template <typename _Graph, typename NodeEnable = void, 
            typename EdgeEnable = void>
  class AlterableSplitGraphAdaptor 
    : public GraphAdaptorExtender<SplitGraphAdaptorBase<_Graph> > {
  public:

    typedef GraphAdaptorExtender<SplitGraphAdaptorBase<_Graph> > Parent;
    typedef _Graph Graph;

    typedef typename Graph::Node GraphNode;
    typedef typename Graph::Node GraphEdge;

  protected:

    AlterableSplitGraphAdaptor() : Parent() {}

  public:
    
    typedef InvalidType NodeNotifier;
    typedef InvalidType EdgeNotifier;

  };

  template <typename _Graph, typename EdgeEnable>
  class AlterableSplitGraphAdaptor<
    _Graph,
    typename enable_if<typename _Graph::NodeNotifier::Notifier>::type,
    EdgeEnable> 
      : public GraphAdaptorExtender<SplitGraphAdaptorBase<_Graph> > {
  public:

    typedef GraphAdaptorExtender<SplitGraphAdaptorBase<_Graph> > Parent;
    typedef _Graph Graph;

    typedef typename Graph::Node GraphNode;
    typedef typename Graph::Edge GraphEdge;

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;
 
  protected:

    AlterableSplitGraphAdaptor() 
      : Parent(), node_notifier(*this), node_notifier_proxy(*this) {}

    void setGraph(_Graph& graph) {
      Parent::setGraph(graph);
      node_notifier_proxy.setNotifier(graph.notifier(GraphNode()));
    }

  public:

    ~AlterableSplitGraphAdaptor() {
      node_notifier.clear();
    }

    typedef AlterationNotifier<AlterableSplitGraphAdaptor, Node> NodeNotifier;
    typedef InvalidType EdgeNotifier;

    NodeNotifier& notifier(Node) const { return node_notifier; }

  protected:

    class NodeNotifierProxy : public Graph::NodeNotifier::ObserverBase {
    public:

      typedef typename Graph::NodeNotifier::ObserverBase Parent;
      typedef AlterableSplitGraphAdaptor AdaptorBase;
      
      NodeNotifierProxy(const AdaptorBase& _adaptor)
        : Parent(), adaptor(&_adaptor) {
      }

      virtual ~NodeNotifierProxy() {
        if (Parent::attached()) {
          Parent::detach();
        }
      }

      void setNotifier(typename Graph::NodeNotifier& graph_notifier) {
        Parent::attach(graph_notifier);
      }

      
    protected:

      virtual void add(const GraphNode& gn) {
        std::vector<Node> nodes;
        nodes.push_back(AdaptorBase::Parent::inNode(gn));
        nodes.push_back(AdaptorBase::Parent::outNode(gn));
        adaptor->notifier(Node()).add(nodes);
      }

      virtual void add(const std::vector<GraphNode>& gn) {
        std::vector<Node> nodes;
        for (int i = 0; i < int(gn.size()); ++i) {
          nodes.push_back(AdaptorBase::Parent::inNode(gn[i]));
          nodes.push_back(AdaptorBase::Parent::outNode(gn[i]));
        }
        adaptor->notifier(Node()).add(nodes);
      }

      virtual void erase(const GraphNode& gn) {
        std::vector<Node> nodes;
        nodes.push_back(AdaptorBase::Parent::inNode(gn));
        nodes.push_back(AdaptorBase::Parent::outNode(gn));
        adaptor->notifier(Node()).erase(nodes);
      }

      virtual void erase(const std::vector<GraphNode>& gn) {
        std::vector<Node> nodes;
        for (int i = 0; i < int(gn.size()); ++i) {
          nodes.push_back(AdaptorBase::Parent::inNode(gn[i]));
          nodes.push_back(AdaptorBase::Parent::outNode(gn[i]));
        }
        adaptor->notifier(Node()).erase(nodes);
      }
      virtual void build() {
        adaptor->notifier(Node()).build();
      }
      virtual void clear() {
        adaptor->notifier(Node()).clear();
      }

      const AdaptorBase* adaptor;
    };


    mutable NodeNotifier node_notifier;

    NodeNotifierProxy node_notifier_proxy;

  };

  template <typename _Graph>
  class AlterableSplitGraphAdaptor<
    _Graph,
    typename enable_if<typename _Graph::NodeNotifier::Notifier>::type,
    typename enable_if<typename _Graph::EdgeNotifier::Notifier>::type> 
      : public GraphAdaptorExtender<SplitGraphAdaptorBase<_Graph> > {
  public:

    typedef GraphAdaptorExtender<SplitGraphAdaptorBase<_Graph> > Parent;
    typedef _Graph Graph;

    typedef typename Graph::Node GraphNode;
    typedef typename Graph::Edge GraphEdge;

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;
 
  protected:
    
    AlterableSplitGraphAdaptor() 
      : Parent(), node_notifier(*this), edge_notifier(*this), 
        node_notifier_proxy(*this), edge_notifier_proxy(*this) {}
    
    void setGraph(_Graph& g) {
      Parent::setGraph(g);
      node_notifier_proxy.setNotifier(g.notifier(GraphNode()));
      edge_notifier_proxy.setNotifier(g.notifier(GraphEdge()));
    }

  public:

    ~AlterableSplitGraphAdaptor() {
      node_notifier.clear();
      edge_notifier.clear();
    }

    typedef AlterationNotifier<AlterableSplitGraphAdaptor, Node> NodeNotifier;
    typedef AlterationNotifier<AlterableSplitGraphAdaptor, Edge> EdgeNotifier;

    NodeNotifier& notifier(Node) const { return node_notifier; }
    EdgeNotifier& notifier(Edge) const { return edge_notifier; }

  protected:

    class NodeNotifierProxy : public Graph::NodeNotifier::ObserverBase {
    public:
      
      typedef typename Graph::NodeNotifier::ObserverBase Parent;
      typedef AlterableSplitGraphAdaptor AdaptorBase;
      
      NodeNotifierProxy(const AdaptorBase& _adaptor)
        : Parent(), adaptor(&_adaptor) {
      }

      virtual ~NodeNotifierProxy() {
        if (Parent::attached()) {
          Parent::detach();
        }
      }

      void setNotifier(typename Graph::NodeNotifier& graph_notifier) {
        Parent::attach(graph_notifier);
      }

      
    protected:

      virtual void add(const GraphNode& gn) {
        std::vector<Node> nodes;
        nodes.push_back(AdaptorBase::Parent::inNode(gn));
        nodes.push_back(AdaptorBase::Parent::outNode(gn));
        adaptor->notifier(Node()).add(nodes);
        adaptor->notifier(Edge()).add(AdaptorBase::Parent::edge(gn));
      }
      virtual void add(const std::vector<GraphNode>& gn) {
        std::vector<Node> nodes;
        std::vector<Edge> edges;
        for (int i = 0; i < int(gn.size()); ++i) {
          edges.push_back(AdaptorBase::Parent::edge(gn[i]));
          nodes.push_back(AdaptorBase::Parent::inNode(gn[i]));
          nodes.push_back(AdaptorBase::Parent::outNode(gn[i]));
        }
        adaptor->notifier(Node()).add(nodes);
        adaptor->notifier(Edge()).add(edges);
      }
      virtual void erase(const GraphNode& gn) {
        adaptor->notifier(Edge()).erase(AdaptorBase::Parent::edge(gn));
        std::vector<Node> nodes;
        nodes.push_back(AdaptorBase::Parent::inNode(gn));
        nodes.push_back(AdaptorBase::Parent::outNode(gn));
        adaptor->notifier(Node()).erase(nodes);
      }
      virtual void erase(const std::vector<GraphNode>& gn) {
        std::vector<Node> nodes;
        std::vector<Edge> edges;
        for (int i = 0; i < int(gn.size()); ++i) {
          edges.push_back(AdaptorBase::Parent::edge(gn[i]));
          nodes.push_back(AdaptorBase::Parent::inNode(gn[i]));
          nodes.push_back(AdaptorBase::Parent::outNode(gn[i]));
        }
        adaptor->notifier(Edge()).erase(edges);
        adaptor->notifier(Node()).erase(nodes);
      }
      virtual void build() {
        std::vector<Edge> edges;
        const typename Parent::Notifier* nf = Parent::notifier();
        GraphNode it;
        for (nf->first(it); it != INVALID; nf->next(it)) {
          edges.push_back(AdaptorBase::Parent::edge(it));
        }
        adaptor->notifier(Node()).build();
        adaptor->notifier(Edge()).add(edges);        
      }
      virtual void clear() {
        std::vector<Edge> edges;
        const typename Parent::Notifier* nf = Parent::notifier();
        GraphNode it;
        for (nf->first(it); it != INVALID; nf->next(it)) {
          edges.push_back(AdaptorBase::Parent::edge(it));
        }
        adaptor->notifier(Edge()).erase(edges);        
        adaptor->notifier(Node()).clear();
      }

      const AdaptorBase* adaptor;
    };

    class EdgeNotifierProxy : public Graph::EdgeNotifier::ObserverBase {
    public:

      typedef typename Graph::EdgeNotifier::ObserverBase Parent;
      typedef AlterableSplitGraphAdaptor AdaptorBase;
      
      EdgeNotifierProxy(const AdaptorBase& _adaptor)
        : Parent(), adaptor(&_adaptor) {
      }

      virtual ~EdgeNotifierProxy() {
        if (Parent::attached()) {
          Parent::detach();
        }
      }

      void setNotifier(typename Graph::EdgeNotifier& graph_notifier) {
        Parent::attach(graph_notifier);
      }

      
    protected:

      virtual void add(const GraphEdge& ge) {
        adaptor->notifier(Edge()).add(AdaptorBase::edge(ge));
      }
      virtual void add(const std::vector<GraphEdge>& ge) {
        std::vector<Edge> edges;
        for (int i = 0; i < int(ge.size()); ++i) {
          edges.push_back(AdaptorBase::edge(ge[i]));
        }
        adaptor->notifier(Edge()).add(edges);
      }
      virtual void erase(const GraphEdge& ge) {
        adaptor->notifier(Edge()).erase(AdaptorBase::edge(ge));
      }
      virtual void erase(const std::vector<GraphEdge>& ge) {
        std::vector<Edge> edges;
        for (int i = 0; i < int(ge.size()); ++i) {
          edges.push_back(AdaptorBase::edge(ge[i]));
        }
        adaptor->notifier(Edge()).erase(edges);
      }
      virtual void build() {
        std::vector<Edge> edges;
        const typename Parent::Notifier* nf = Parent::notifier();
        GraphEdge it;
        for (nf->first(it); it != INVALID; nf->next(it)) {
          edges.push_back(AdaptorBase::Parent::edge(it));
        }
        adaptor->notifier(Edge()).add(edges);
      }
      virtual void clear() {
        std::vector<Edge> edges;
        const typename Parent::Notifier* nf = Parent::notifier();
        GraphEdge it;
        for (nf->first(it); it != INVALID; nf->next(it)) {
          edges.push_back(AdaptorBase::Parent::edge(it));
        }
        adaptor->notifier(Edge()).erase(edges);
      }

      const AdaptorBase* adaptor;
    };


    mutable NodeNotifier node_notifier;
    mutable EdgeNotifier edge_notifier;

    NodeNotifierProxy node_notifier_proxy;
    EdgeNotifierProxy edge_notifier_proxy;

  };

  /// \ingroup graph_adaptors
  ///
  /// \brief Split graph adaptor class
  /// 
  /// This is an graph adaptor which splits all node into an in-node
  /// and an out-node. Formaly, the adaptor replaces each \f$ u \f$
  /// node in the graph with two node, \f$ u_{in} \f$ node and 
  /// \f$ u_{out} \f$ node. If there is an \f$ (v, u) \f$ edge in the 
  /// original graph the new target of the edge will be \f$ u_{in} \f$ and
  /// similarly the source of the original \f$ (u, v) \f$ edge will be
  /// \f$ u_{out} \f$.  The adaptor will add for each node in the 
  /// original graph an additional edge which will connect 
  /// \f$ (u_{in}, u_{out}) \f$.
  ///
  /// The aim of this class is to run algorithm with node costs if the 
  /// algorithm can use directly just edge costs. In this case we should use
  /// a \c SplitGraphAdaptor and set the node cost of the graph to the
  /// bind edge in the adapted graph.
  /// 
  /// By example a maximum flow algoritm can compute how many edge
  /// disjoint paths are in the graph. But we would like to know how
  /// many node disjoint paths are in the graph. First we have to
  /// adapt the graph with the \c SplitGraphAdaptor. Then run the flow
  /// algorithm on the adapted graph. The bottleneck of the flow will
  /// be the bind edges which bounds the flow with the count of the
  /// node disjoint paths.
  ///
  ///\code
  ///
  /// typedef SplitGraphAdaptor<SmartGraph> SGraph;
  ///
  /// SGraph sgraph(graph);
  ///
  /// typedef ConstMap<SGraph::Edge, int> SCapacity;
  /// SCapacity scapacity(1);
  ///
  /// SGraph::EdgeMap<int> sflow(sgraph);
  ///
  /// Preflow<SGraph, SCapacity> 
  ///   spreflow(sgraph, scapacity, 
  ///            SGraph::outNode(source), SGraph::inNode(target));
  ///                                            
  /// spreflow.run();
  ///
  ///\endcode
  ///
  /// The result of the mamixum flow on the original graph
  /// shows the next figure:
  ///
  /// \image html edge_disjoint.png
  /// \image latex edge_disjoint.eps "Edge disjoint paths" width=\textwidth
  /// 
  /// And the maximum flow on the adapted graph:
  ///
  /// \image html node_disjoint.png
  /// \image latex node_disjoint.eps "Node disjoint paths" width=\textwidth
  ///
  /// The second solution contains just 3 disjoint paths while the first 4.
  /// The full code can be found in the \ref disjoint_paths_demo.cc demo file.
  ///
  /// This graph adaptor is fully conform to the 
  /// \ref concepts::Graph "Graph" concept and
  /// contains some additional member functions and types. The 
  /// documentation of some member functions may be found just in the
  /// SplitGraphAdaptorBase class.
  ///
  /// \sa SplitGraphAdaptorBase
  template <typename _Graph>
  class SplitGraphAdaptor : public AlterableSplitGraphAdaptor<_Graph> {
  public:
    typedef AlterableSplitGraphAdaptor<_Graph> Parent;

    typedef typename Parent::Node Node;
    typedef typename Parent::Edge Edge;

    /// \brief Constructor of the adaptor.
    ///
    /// Constructor of the adaptor.
    SplitGraphAdaptor(_Graph& g) {
      Parent::setGraph(g);
    }

    /// \brief NodeMap combined from two original NodeMap
    ///
    /// This class adapt two of the original graph NodeMap to
    /// get a node map on the adapted graph.
    template <typename InNodeMap, typename OutNodeMap>
    class CombinedNodeMap {
    public:

      typedef Node Key;
      typedef typename InNodeMap::Value Value;

      /// \brief Constructor
      ///
      /// Constructor.
      CombinedNodeMap(InNodeMap& _inNodeMap, OutNodeMap& _outNodeMap) 
	: inNodeMap(_inNodeMap), outNodeMap(_outNodeMap) {}

      /// \brief The subscript operator.
      ///
      /// The subscript operator.
      Value& operator[](const Key& key) {
	if (Parent::inNode(key)) {
	  return inNodeMap[key];
	} else {
	  return outNodeMap[key];
	}
      }

      /// \brief The const subscript operator.
      ///
      /// The const subscript operator.
      Value operator[](const Key& key) const {
	if (Parent::inNode(key)) {
	  return inNodeMap[key];
	} else {
	  return outNodeMap[key];
	}
      }

      /// \brief The setter function of the map.
      /// 
      /// The setter function of the map.
      void set(const Key& key, const Value& value) {
	if (Parent::inNode(key)) {
	  inNodeMap.set(key, value);
	} else {
	  outNodeMap.set(key, value);
	}
      }
      
    private:
      
      InNodeMap& inNodeMap;
      OutNodeMap& outNodeMap;
      
    };


    /// \brief Just gives back a combined node map.
    /// 
    /// Just gives back a combined node map.
    template <typename InNodeMap, typename OutNodeMap>
    static CombinedNodeMap<InNodeMap, OutNodeMap> 
    combinedNodeMap(InNodeMap& in_map, OutNodeMap& out_map) {
      return CombinedNodeMap<InNodeMap, OutNodeMap>(in_map, out_map);
    }

    template <typename InNodeMap, typename OutNodeMap>
    static CombinedNodeMap<const InNodeMap, OutNodeMap> 
    combinedNodeMap(const InNodeMap& in_map, OutNodeMap& out_map) {
      return CombinedNodeMap<const InNodeMap, OutNodeMap>(in_map, out_map);
    }

    template <typename InNodeMap, typename OutNodeMap>
    static CombinedNodeMap<InNodeMap, const OutNodeMap> 
    combinedNodeMap(InNodeMap& in_map, const OutNodeMap& out_map) {
      return CombinedNodeMap<InNodeMap, const OutNodeMap>(in_map, out_map);
    }

    template <typename InNodeMap, typename OutNodeMap>
    static CombinedNodeMap<const InNodeMap, const OutNodeMap> 
    combinedNodeMap(const InNodeMap& in_map, const OutNodeMap& out_map) {
      return CombinedNodeMap<const InNodeMap, 
        const OutNodeMap>(in_map, out_map);
    }

    /// \brief EdgeMap combined from an original EdgeMap and NodeMap
    ///
    /// This class adapt an original graph EdgeMap and NodeMap to
    /// get an edge map on the adapted graph.
    template <typename GraphEdgeMap, typename GraphNodeMap>
    class CombinedEdgeMap {
    public:
      
      typedef Edge Key;
      typedef typename GraphEdgeMap::Value Value;
      
      /// \brief Constructor
      ///
      /// Constructor.
      CombinedEdgeMap(GraphEdgeMap& _edge_map, GraphNodeMap& _node_map) 
	: edge_map(_edge_map), node_map(_node_map) {}

      /// \brief The subscript operator.
      ///
      /// The subscript operator.
      void set(const Edge& edge, const Value& val) {
	if (Parent::origEdge(edge)) {
	  edge_map.set(edge, val);
	} else {
	  node_map.set(edge, val);
	}
      }

      /// \brief The const subscript operator.
      ///
      /// The const subscript operator.
      Value operator[](const Key& edge) const {
	if (Parent::origEdge(edge)) {
	  return edge_map[edge];
	} else {
	  return node_map[edge];
	}
      }      

      /// \brief The const subscript operator.
      ///
      /// The const subscript operator.
      Value& operator[](const Key& edge) {
	if (Parent::origEdge(edge)) {
	  return edge_map[edge];
	} else {
	  return node_map[edge];
	}
      }      
      
    private:
      GraphEdgeMap& edge_map;
      GraphNodeMap& node_map;
    };
                    
    /// \brief Just gives back a combined edge map.
    /// 
    /// Just gives back a combined edge map.
    template <typename GraphEdgeMap, typename GraphNodeMap>
    static CombinedEdgeMap<GraphEdgeMap, GraphNodeMap> 
    combinedEdgeMap(GraphEdgeMap& edge_map, GraphNodeMap& node_map) {
      return CombinedEdgeMap<GraphEdgeMap, GraphNodeMap>(edge_map, node_map);
    }

    template <typename GraphEdgeMap, typename GraphNodeMap>
    static CombinedEdgeMap<const GraphEdgeMap, GraphNodeMap> 
    combinedEdgeMap(const GraphEdgeMap& edge_map, GraphNodeMap& node_map) {
      return CombinedEdgeMap<const GraphEdgeMap, 
        GraphNodeMap>(edge_map, node_map);
    }

    template <typename GraphEdgeMap, typename GraphNodeMap>
    static CombinedEdgeMap<GraphEdgeMap, const GraphNodeMap> 
    combinedEdgeMap(GraphEdgeMap& edge_map, const GraphNodeMap& node_map) {
      return CombinedEdgeMap<GraphEdgeMap, 
        const GraphNodeMap>(edge_map, node_map);
    }

    template <typename GraphEdgeMap, typename GraphNodeMap>
    static CombinedEdgeMap<const GraphEdgeMap, const GraphNodeMap> 
    combinedEdgeMap(const GraphEdgeMap& edge_map, 
                    const GraphNodeMap& node_map) {
      return CombinedEdgeMap<const GraphEdgeMap, 
        const GraphNodeMap>(edge_map, node_map);
    }

  };

  /// \brief Just gives back a split graph adaptor
  ///
  /// Just gives back a split graph adaptor
  template<typename Graph>
  SplitGraphAdaptor<Graph>
  splitGraphAdaptor(const Graph& graph) {
    return SplitGraphAdaptor<Graph>(graph);
  }


} //namespace lemon

#endif //LEMON_GRAPH_ADAPTOR_H


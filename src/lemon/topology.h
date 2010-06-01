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

#ifndef LEMON_TOPOLOGY_H
#define LEMON_TOPOLOGY_H

#include <lemon/dfs.h>
#include <lemon/bfs.h>
#include <lemon/graph_utils.h>
#include <lemon/graph_adaptor.h>
#include <lemon/maps.h>

#include <lemon/concepts/graph.h>
#include <lemon/concepts/ugraph.h>
#include <lemon/concept_check.h>

#include <lemon/bin_heap.h>
#include <lemon/bucket_heap.h>

#include <stack>
#include <functional>

/// \ingroup graph_prop
/// \file
/// \brief Topology related algorithms
///
/// Topology related algorithms

namespace lemon {

  /// \ingroup graph_prop
  ///
  /// \brief Check that the given undirected graph is connected.
  ///
  /// Check that the given undirected graph connected.
  /// \param graph The undirected graph.
  /// \return %True when there is path between any two nodes in the graph.
  /// \note By definition, the empty graph is connected.
  template <typename UGraph>
  bool connected(const UGraph& graph) {
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::NodeIt NodeIt;
    if (NodeIt(graph) == INVALID) return true;
    Dfs<UGraph> dfs(graph);
    dfs.run(NodeIt(graph));
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	return false;
      }
    }
    return true;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Count the number of connected components of an undirected graph
  ///
  /// Count the number of connected components of an undirected graph
  ///
  /// \param graph The graph. It should be undirected.
  /// \return The number of components
  /// \note By definition, the empty graph consists
  /// of zero connected components.
  template <typename UGraph>
  int countConnectedComponents(const UGraph &graph) {
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::Node Node;
    typedef typename UGraph::Edge Edge;

    typedef NullMap<Node, Edge> PredMap;
    typedef NullMap<Node, int> DistMap;

    int compNum = 0;
    typename Bfs<UGraph>::
      template DefPredMap<PredMap>::
      template DefDistMap<DistMap>::
      Create bfs(graph);

    PredMap predMap;
    bfs.predMap(predMap);

    DistMap distMap;
    bfs.distMap(distMap);

    bfs.init();
    for(typename UGraph::NodeIt n(graph); n != INVALID; ++n) {
      if (!bfs.reached(n)) {
	bfs.addSource(n);
	bfs.start();
	++compNum;
      }
    }
    return compNum;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Find the connected components of an undirected graph
  ///
  /// Find the connected components of an undirected graph.
  ///
  /// \image html connected_components.png
  /// \image latex connected_components.eps "Connected components" width=\textwidth
  ///
  /// \param graph The graph. It should be undirected.
  /// \retval compMap A writable node map. The values will be set from 0 to
  /// the number of the connected components minus one. Each values of the map
  /// will be set exactly once, the values of a certain component will be
  /// set continuously.
  /// \return The number of components
  ///
  template <class UGraph, class NodeMap>
  int connectedComponents(const UGraph &graph, NodeMap &compMap) {
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::Node Node;
    typedef typename UGraph::Edge Edge;
    checkConcept<concepts::WriteMap<Node, int>, NodeMap>();

    typedef NullMap<Node, Edge> PredMap;
    typedef NullMap<Node, int> DistMap;

    int compNum = 0;
    typename Bfs<UGraph>::
      template DefPredMap<PredMap>::
      template DefDistMap<DistMap>::
      Create bfs(graph);

    PredMap predMap;
    bfs.predMap(predMap);

    DistMap distMap;
    bfs.distMap(distMap);
    
    bfs.init();
    for(typename UGraph::NodeIt n(graph); n != INVALID; ++n) {
      if(!bfs.reached(n)) {
	bfs.addSource(n);
	while (!bfs.emptyQueue()) {
	  compMap.set(bfs.nextNode(), compNum);
	  bfs.processNextNode();
	}
	++compNum;
      }
    }
    return compNum;
  }

  namespace _topology_bits {

    template <typename Graph, typename Iterator >
    struct LeaveOrderVisitor : public DfsVisitor<Graph> {
    public:
      typedef typename Graph::Node Node;
      LeaveOrderVisitor(Iterator it) : _it(it) {}

      void leave(const Node& node) {
	*(_it++) = node;
      }

    private:
      Iterator _it;
    };

    template <typename Graph, typename Map>
    struct FillMapVisitor : public DfsVisitor<Graph> {
    public:
      typedef typename Graph::Node Node;
      typedef typename Map::Value Value;

      FillMapVisitor(Map& map, Value& value) 
	: _map(map), _value(value) {}

      void reach(const Node& node) {
	_map.set(node, _value);
      }
    private:
      Map& _map;
      Value& _value;
    };

    template <typename Graph, typename EdgeMap>
    struct StronglyConnectedCutEdgesVisitor : public DfsVisitor<Graph> {
    public:
      typedef typename Graph::Node Node;
      typedef typename Graph::Edge Edge;

      StronglyConnectedCutEdgesVisitor(const Graph& graph, EdgeMap& cutMap, 
				       int& cutNum) 
	: _graph(graph), _cutMap(cutMap), _cutNum(cutNum), 
	  _compMap(graph), _num(0) {
      }

      void stop(const Node&) {
	++_num;
      }

      void reach(const Node& node) {
	_compMap.set(node, _num);
      }

      void examine(const Edge& edge) {
 	if (_compMap[_graph.source(edge)] != _compMap[_graph.target(edge)]) {
 	  _cutMap.set(edge, true);
 	  ++_cutNum;
 	}
      }
    private:
      const Graph& _graph;
      EdgeMap& _cutMap;
      int& _cutNum;

      typename Graph::template NodeMap<int> _compMap;
      int _num;
    };

  }


  /// \ingroup graph_prop
  ///
  /// \brief Check that the given directed graph is strongly connected.
  ///
  /// Check that the given directed graph is strongly connected. The
  /// graph is strongly connected when any two nodes of the graph are
  /// connected with directed paths in both direction.
  /// \return %False when the graph is not strongly connected.
  /// \see connected
  ///
  /// \note By definition, the empty graph is strongly connected.
  template <typename Graph>
  bool stronglyConnected(const Graph& graph) {
    checkConcept<concepts::Graph, Graph>();

    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;

    if (NodeIt(graph) == INVALID) return true;

    using namespace _topology_bits;

    typedef DfsVisitor<Graph> Visitor;
    Visitor visitor;

    DfsVisit<Graph, Visitor> dfs(graph, visitor);
    dfs.init();
    dfs.addSource(NodeIt(graph));
    dfs.start();

    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	return false;
      }
    }

    typedef RevGraphAdaptor<const Graph> RGraph;
    RGraph rgraph(graph);

    typedef DfsVisitor<Graph> RVisitor;
    RVisitor rvisitor;

    DfsVisit<RGraph, RVisitor> rdfs(rgraph, rvisitor);
    rdfs.init();    
    rdfs.addSource(NodeIt(graph));
    rdfs.start();

    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!rdfs.reached(it)) {
	return false;
      }
    }

    return true;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Count the strongly connected components of a directed graph
  ///
  /// Count the strongly connected components of a directed graph.
  /// The strongly connected components are the classes of an
  /// equivalence relation on the nodes of the graph. Two nodes are in
  /// the same class if they are connected with directed paths in both
  /// direction. 
  ///
  /// \param graph The graph.
  /// \return The number of components
  /// \note By definition, the empty graph has zero
  /// strongly connected components.
  template <typename Graph>
  int countStronglyConnectedComponents(const Graph& graph) {
    checkConcept<concepts::Graph, Graph>();

    using namespace _topology_bits;

    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::EdgeIt EdgeIt;
    
    typedef std::vector<Node> Container;
    typedef typename Container::iterator Iterator;

    Container nodes(countNodes(graph));
    typedef LeaveOrderVisitor<Graph, Iterator> Visitor;
    Visitor visitor(nodes.begin());
      
    DfsVisit<Graph, Visitor> dfs(graph, visitor);
    dfs.init();
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	dfs.start();
      }
    }

    typedef typename Container::reverse_iterator RIterator;
    typedef RevGraphAdaptor<const Graph> RGraph;

    RGraph rgraph(graph);

    typedef DfsVisitor<Graph> RVisitor;
    RVisitor rvisitor;

    DfsVisit<RGraph, RVisitor> rdfs(rgraph, rvisitor);

    int compNum = 0;

    rdfs.init();
    for (RIterator it = nodes.rbegin(); it != nodes.rend(); ++it) {
      if (!rdfs.reached(*it)) {
	rdfs.addSource(*it);
	rdfs.start();
	++compNum;
      }
    }
    return compNum;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Find the strongly connected components of a directed graph
  ///
  /// Find the strongly connected components of a directed graph.  The
  /// strongly connected components are the classes of an equivalence
  /// relation on the nodes of the graph. Two nodes are in
  /// relationship when there are directed paths between them in both
  /// direction. In addition, the numbering of components will satisfy
  /// that there is no edge going from a higher numbered component to
  /// a lower.
  ///
  /// \image html strongly_connected_components.png
  /// \image latex strongly_connected_components.eps "Strongly connected components" width=\textwidth
  ///
  /// \param graph The graph.
  /// \retval compMap A writable node map. The values will be set from 0 to
  /// the number of the strongly connected components minus one. Each value 
  /// of the map will be set exactly once, the values of a certain component 
  /// will be set continuously.
  /// \return The number of components
  ///
  template <typename Graph, typename NodeMap>
  int stronglyConnectedComponents(const Graph& graph, NodeMap& compMap) {
    checkConcept<concepts::Graph, Graph>();
    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;
    checkConcept<concepts::WriteMap<Node, int>, NodeMap>();

    using namespace _topology_bits;
    
    typedef std::vector<Node> Container;
    typedef typename Container::iterator Iterator;

    Container nodes(countNodes(graph));
    typedef LeaveOrderVisitor<Graph, Iterator> Visitor;
    Visitor visitor(nodes.begin());
      
    DfsVisit<Graph, Visitor> dfs(graph, visitor);
    dfs.init();
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	dfs.start();
      }
    }

    typedef typename Container::reverse_iterator RIterator;
    typedef RevGraphAdaptor<const Graph> RGraph;

    RGraph rgraph(graph);

    int compNum = 0;

    typedef FillMapVisitor<RGraph, NodeMap> RVisitor;
    RVisitor rvisitor(compMap, compNum);

    DfsVisit<RGraph, RVisitor> rdfs(rgraph, rvisitor);

    rdfs.init();
    for (RIterator it = nodes.rbegin(); it != nodes.rend(); ++it) {
      if (!rdfs.reached(*it)) {
	rdfs.addSource(*it);
	rdfs.start();
	++compNum;
      }
    }
    return compNum;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Find the cut edges of the strongly connected components.
  ///
  /// Find the cut edges of the strongly connected components.
  /// The strongly connected components are the classes of an equivalence
  /// relation on the nodes of the graph. Two nodes are in relationship
  /// when there are directed paths between them in both direction.
  /// The strongly connected components are separated by the cut edges.
  ///
  /// \param graph The graph.
  /// \retval cutMap A writable node map. The values will be set true when the
  /// edge is a cut edge.
  ///
  /// \return The number of cut edges
  template <typename Graph, typename EdgeMap>
  int stronglyConnectedCutEdges(const Graph& graph, EdgeMap& cutMap) {
    checkConcept<concepts::Graph, Graph>();
    typedef typename Graph::Node Node;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::NodeIt NodeIt;
    checkConcept<concepts::WriteMap<Edge, bool>, EdgeMap>();

    using namespace _topology_bits;
    
    typedef std::vector<Node> Container;
    typedef typename Container::iterator Iterator;

    Container nodes(countNodes(graph));
    typedef LeaveOrderVisitor<Graph, Iterator> Visitor;
    Visitor visitor(nodes.begin());
      
    DfsVisit<Graph, Visitor> dfs(graph, visitor);
    dfs.init();
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	dfs.start();
      }
    }

    typedef typename Container::reverse_iterator RIterator;
    typedef RevGraphAdaptor<const Graph> RGraph;

    RGraph rgraph(graph);

    int cutNum = 0;

    typedef StronglyConnectedCutEdgesVisitor<RGraph, EdgeMap> RVisitor;
    RVisitor rvisitor(rgraph, cutMap, cutNum);

    DfsVisit<RGraph, RVisitor> rdfs(rgraph, rvisitor);

    rdfs.init();
    for (RIterator it = nodes.rbegin(); it != nodes.rend(); ++it) {
      if (!rdfs.reached(*it)) {
	rdfs.addSource(*it);
	rdfs.start();
      }
    }
    return cutNum;
  }

  namespace _topology_bits {
    
    template <typename Graph>
    class CountBiNodeConnectedComponentsVisitor : public DfsVisitor<Graph> {
    public:
      typedef typename Graph::Node Node;
      typedef typename Graph::Edge Edge;
      typedef typename Graph::UEdge UEdge;

      CountBiNodeConnectedComponentsVisitor(const Graph& graph, int &compNum) 
	: _graph(graph), _compNum(compNum), 
	  _numMap(graph), _retMap(graph), _predMap(graph), _num(0) {}

      void start(const Node& node) {
	_predMap.set(node, INVALID);
      }
      
      void reach(const Node& node) {
	_numMap.set(node, _num);
	_retMap.set(node, _num);
	++_num;
      }

      void discover(const Edge& edge) {
	_predMap.set(_graph.target(edge), _graph.source(edge));
      }

      void examine(const Edge& edge) {
	if (_graph.source(edge) == _graph.target(edge) && 
	    _graph.direction(edge)) {
	  ++_compNum;
	  return;
	}
	if (_predMap[_graph.source(edge)] == _graph.target(edge)) {
	  return;
	}
	if (_retMap[_graph.source(edge)] > _numMap[_graph.target(edge)]) {
	  _retMap.set(_graph.source(edge), _numMap[_graph.target(edge)]);
	}
      }

      void backtrack(const Edge& edge) {
	if (_retMap[_graph.source(edge)] > _retMap[_graph.target(edge)]) {
	  _retMap.set(_graph.source(edge), _retMap[_graph.target(edge)]);
	}  
	if (_numMap[_graph.source(edge)] <= _retMap[_graph.target(edge)]) {
	  ++_compNum;
	}
      }
      
    private:
      const Graph& _graph;
      int& _compNum; 

      typename Graph::template NodeMap<int> _numMap;
      typename Graph::template NodeMap<int> _retMap;
      typename Graph::template NodeMap<Node> _predMap;
      int _num;
    };

    template <typename Graph, typename EdgeMap>
    class BiNodeConnectedComponentsVisitor : public DfsVisitor<Graph> {
    public:
      typedef typename Graph::Node Node;
      typedef typename Graph::Edge Edge;
      typedef typename Graph::UEdge UEdge;

      BiNodeConnectedComponentsVisitor(const Graph& graph, 
				       EdgeMap& compMap, int &compNum) 
	: _graph(graph), _compMap(compMap), _compNum(compNum), 
	  _numMap(graph), _retMap(graph), _predMap(graph), _num(0) {}

      void start(const Node& node) {
	_predMap.set(node, INVALID);
      }
      
      void reach(const Node& node) {
	_numMap.set(node, _num);
	_retMap.set(node, _num);
	++_num;
      }

      void discover(const Edge& edge) {
	Node target = _graph.target(edge);
	_predMap.set(target, edge);
	_edgeStack.push(edge);
      }

      void examine(const Edge& edge) {
	Node source = _graph.source(edge);
	Node target = _graph.target(edge);
	if (source == target && _graph.direction(edge)) {
	  _compMap.set(edge, _compNum);
	  ++_compNum;
	  return;
	}
	if (_numMap[target] < _numMap[source]) {
	  if (_predMap[source] != _graph.oppositeEdge(edge)) {
	    _edgeStack.push(edge);
	  }
	}
	if (_predMap[source] != INVALID && 
	    target == _graph.source(_predMap[source])) {
	  return;
	}
	if (_retMap[source] > _numMap[target]) {
	  _retMap.set(source, _numMap[target]);
	}
      }

      void backtrack(const Edge& edge) {
	Node source = _graph.source(edge);
	Node target = _graph.target(edge);
	if (_retMap[source] > _retMap[target]) {
	  _retMap.set(source, _retMap[target]);
	}  
	if (_numMap[source] <= _retMap[target]) {
	  while (_edgeStack.top() != edge) {
	    _compMap.set(_edgeStack.top(), _compNum);
	    _edgeStack.pop();
	  }
	  _compMap.set(edge, _compNum);
	  _edgeStack.pop();
	  ++_compNum;
	}
      }
      
    private:
      const Graph& _graph;
      EdgeMap& _compMap;
      int& _compNum; 

      typename Graph::template NodeMap<int> _numMap;
      typename Graph::template NodeMap<int> _retMap;
      typename Graph::template NodeMap<Edge> _predMap;
      std::stack<UEdge> _edgeStack;
      int _num;
    };


    template <typename Graph, typename NodeMap>
    class BiNodeConnectedCutNodesVisitor : public DfsVisitor<Graph> {
    public:
      typedef typename Graph::Node Node;
      typedef typename Graph::Edge Edge;
      typedef typename Graph::UEdge UEdge;

      BiNodeConnectedCutNodesVisitor(const Graph& graph, NodeMap& cutMap,
				     int& cutNum) 
	: _graph(graph), _cutMap(cutMap), _cutNum(cutNum),
	  _numMap(graph), _retMap(graph), _predMap(graph), _num(0) {}

      void start(const Node& node) {
	_predMap.set(node, INVALID);
	rootCut = false;
      }
      
      void reach(const Node& node) {
	_numMap.set(node, _num);
	_retMap.set(node, _num);
	++_num;
      }

      void discover(const Edge& edge) {
	_predMap.set(_graph.target(edge), _graph.source(edge));
      }

      void examine(const Edge& edge) {
	if (_graph.source(edge) == _graph.target(edge) && 
	    _graph.direction(edge)) {
	  if (!_cutMap[_graph.source(edge)]) {
	    _cutMap.set(_graph.source(edge), true);
	    ++_cutNum;
	  }
	  return;
	}
	if (_predMap[_graph.source(edge)] == _graph.target(edge)) return;
	if (_retMap[_graph.source(edge)] > _numMap[_graph.target(edge)]) {
	  _retMap.set(_graph.source(edge), _numMap[_graph.target(edge)]);
	}
      }

      void backtrack(const Edge& edge) {
	if (_retMap[_graph.source(edge)] > _retMap[_graph.target(edge)]) {
	  _retMap.set(_graph.source(edge), _retMap[_graph.target(edge)]);
	}  
	if (_numMap[_graph.source(edge)] <= _retMap[_graph.target(edge)]) {
	  if (_predMap[_graph.source(edge)] != INVALID) {
	    if (!_cutMap[_graph.source(edge)]) {
	      _cutMap.set(_graph.source(edge), true);
	      ++_cutNum;
	    }
	  } else if (rootCut) {
	    if (!_cutMap[_graph.source(edge)]) {
	      _cutMap.set(_graph.source(edge), true);
	      ++_cutNum;
	    }
	  } else {
	    rootCut = true;
	  }
	}
      }
      
    private:
      const Graph& _graph;
      NodeMap& _cutMap;
      int& _cutNum; 

      typename Graph::template NodeMap<int> _numMap;
      typename Graph::template NodeMap<int> _retMap;
      typename Graph::template NodeMap<Node> _predMap;
      std::stack<UEdge> _edgeStack;
      int _num;
      bool rootCut;
    };

  }

  template <typename UGraph>
  int countBiNodeConnectedComponents(const UGraph& graph);

  /// \ingroup graph_prop
  ///
  /// \brief Checks the graph is bi-node-connected.
  ///
  /// This function checks that the undirected graph is bi-node-connected  
  /// graph. The graph is bi-node-connected if any two undirected edge is 
  /// on same circle.
  ///
  /// \param graph The graph.
  /// \return %True when the graph bi-node-connected.
  template <typename UGraph>
  bool biNodeConnected(const UGraph& graph) {
    return countBiNodeConnectedComponents(graph) == 1;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Count the biconnected components.
  ///
  /// This function finds the bi-node-connected components in an undirected 
  /// graph. The biconnected components are the classes of an equivalence 
  /// relation on the undirected edges. Two undirected edge is in relationship
  /// when they are on same circle.
  ///
  /// \param graph The graph.
  /// \return The number of components.
  template <typename UGraph>
  int countBiNodeConnectedComponents(const UGraph& graph) {
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::NodeIt NodeIt;

    using namespace _topology_bits;

    typedef CountBiNodeConnectedComponentsVisitor<UGraph> Visitor;

    int compNum = 0;
    Visitor visitor(graph, compNum);

    DfsVisit<UGraph, Visitor> dfs(graph, visitor);
    dfs.init();
    
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	dfs.start();
      }
    }
    return compNum;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Find the bi-node-connected components.
  ///
  /// This function finds the bi-node-connected components in an undirected 
  /// graph. The bi-node-connected components are the classes of an equivalence
  /// relation on the undirected edges. Two undirected edge are in relationship
  /// when they are on same circle.
  ///
  /// \image html node_biconnected_components.png
  /// \image latex node_biconnected_components.eps "bi-node-connected components" width=\textwidth
  ///
  /// \param graph The graph.
  /// \retval compMap A writable uedge map. The values will be set from 0
  /// to the number of the biconnected components minus one. Each values 
  /// of the map will be set exactly once, the values of a certain component 
  /// will be set continuously.
  /// \return The number of components.
  ///
  template <typename UGraph, typename UEdgeMap>
  int biNodeConnectedComponents(const UGraph& graph, 
				UEdgeMap& compMap) {
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::NodeIt NodeIt;
    typedef typename UGraph::UEdge UEdge;
    checkConcept<concepts::WriteMap<UEdge, int>, UEdgeMap>();

    using namespace _topology_bits;

    typedef BiNodeConnectedComponentsVisitor<UGraph, UEdgeMap> Visitor;
    
    int compNum = 0;
    Visitor visitor(graph, compMap, compNum);

    DfsVisit<UGraph, Visitor> dfs(graph, visitor);
    dfs.init();
    
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	dfs.start();
      }
    }
    return compNum;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Find the bi-node-connected cut nodes.
  ///
  /// This function finds the bi-node-connected cut nodes in an undirected 
  /// graph. The bi-node-connected components are the classes of an equivalence
  /// relation on the undirected edges. Two undirected edges are in 
  /// relationship when they are on same circle. The biconnected components 
  /// are separted by nodes which are the cut nodes of the components.
  ///
  /// \param graph The graph.
  /// \retval cutMap A writable edge map. The values will be set true when
  /// the node separate two or more components.
  /// \return The number of the cut nodes.
  template <typename UGraph, typename NodeMap>
  int biNodeConnectedCutNodes(const UGraph& graph, NodeMap& cutMap) {
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::Node Node;
    typedef typename UGraph::NodeIt NodeIt;
    checkConcept<concepts::WriteMap<Node, bool>, NodeMap>();

    using namespace _topology_bits;

    typedef BiNodeConnectedCutNodesVisitor<UGraph, NodeMap> Visitor;
    
    int cutNum = 0;
    Visitor visitor(graph, cutMap, cutNum);

    DfsVisit<UGraph, Visitor> dfs(graph, visitor);
    dfs.init();
    
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	dfs.start();
      }
    }
    return cutNum;
  }

  namespace _topology_bits {
    
    template <typename Graph>
    class CountBiEdgeConnectedComponentsVisitor : public DfsVisitor<Graph> {
    public:
      typedef typename Graph::Node Node;
      typedef typename Graph::Edge Edge;
      typedef typename Graph::UEdge UEdge;

      CountBiEdgeConnectedComponentsVisitor(const Graph& graph, int &compNum) 
	: _graph(graph), _compNum(compNum), 
	  _numMap(graph), _retMap(graph), _predMap(graph), _num(0) {}

      void start(const Node& node) {
	_predMap.set(node, INVALID);
      }
      
      void reach(const Node& node) {
	_numMap.set(node, _num);
	_retMap.set(node, _num);
	++_num;
      }
      
      void leave(const Node& node) {
	if (_numMap[node] <= _retMap[node]) {
	  ++_compNum;
	}	
      }

      void discover(const Edge& edge) {
	_predMap.set(_graph.target(edge), edge);
      }

      void examine(const Edge& edge) {
	if (_predMap[_graph.source(edge)] == _graph.oppositeEdge(edge)) {
	  return;
	}
	if (_retMap[_graph.source(edge)] > _retMap[_graph.target(edge)]) {
	  _retMap.set(_graph.source(edge), _retMap[_graph.target(edge)]);
	}
      }

      void backtrack(const Edge& edge) {
	if (_retMap[_graph.source(edge)] > _retMap[_graph.target(edge)]) {
	  _retMap.set(_graph.source(edge), _retMap[_graph.target(edge)]);
	}  
      }
      
    private:
      const Graph& _graph;
      int& _compNum; 

      typename Graph::template NodeMap<int> _numMap;
      typename Graph::template NodeMap<int> _retMap;
      typename Graph::template NodeMap<Edge> _predMap;
      int _num;
    };

    template <typename Graph, typename NodeMap>
    class BiEdgeConnectedComponentsVisitor : public DfsVisitor<Graph> {
    public:
      typedef typename Graph::Node Node;
      typedef typename Graph::Edge Edge;
      typedef typename Graph::UEdge UEdge;

      BiEdgeConnectedComponentsVisitor(const Graph& graph, 
				       NodeMap& compMap, int &compNum) 
	: _graph(graph), _compMap(compMap), _compNum(compNum), 
	  _numMap(graph), _retMap(graph), _predMap(graph), _num(0) {}

      void start(const Node& node) {
	_predMap.set(node, INVALID);
      }
      
      void reach(const Node& node) {
	_numMap.set(node, _num);
	_retMap.set(node, _num);
	_nodeStack.push(node);
	++_num;
      }
      
      void leave(const Node& node) {
	if (_numMap[node] <= _retMap[node]) {
	  while (_nodeStack.top() != node) {
	    _compMap.set(_nodeStack.top(), _compNum);
	    _nodeStack.pop();
	  }
	  _compMap.set(node, _compNum);
	  _nodeStack.pop();
	  ++_compNum;
	}	
      }

      void discover(const Edge& edge) {
	_predMap.set(_graph.target(edge), edge);
      }

      void examine(const Edge& edge) {
	if (_predMap[_graph.source(edge)] == _graph.oppositeEdge(edge)) {
	  return;
	}
	if (_retMap[_graph.source(edge)] > _retMap[_graph.target(edge)]) {
	  _retMap.set(_graph.source(edge), _retMap[_graph.target(edge)]);
	}
      }

      void backtrack(const Edge& edge) {
	if (_retMap[_graph.source(edge)] > _retMap[_graph.target(edge)]) {
	  _retMap.set(_graph.source(edge), _retMap[_graph.target(edge)]);
	}  
      }
      
    private:
      const Graph& _graph;
      NodeMap& _compMap;
      int& _compNum; 

      typename Graph::template NodeMap<int> _numMap;
      typename Graph::template NodeMap<int> _retMap;
      typename Graph::template NodeMap<Edge> _predMap;
      std::stack<Node> _nodeStack;
      int _num;
    };


    template <typename Graph, typename EdgeMap>
    class BiEdgeConnectedCutEdgesVisitor : public DfsVisitor<Graph> {
    public:
      typedef typename Graph::Node Node;
      typedef typename Graph::Edge Edge;
      typedef typename Graph::UEdge UEdge;

      BiEdgeConnectedCutEdgesVisitor(const Graph& graph, 
				     EdgeMap& cutMap, int &cutNum) 
	: _graph(graph), _cutMap(cutMap), _cutNum(cutNum), 
	  _numMap(graph), _retMap(graph), _predMap(graph), _num(0) {}

      void start(const Node& node) {
	_predMap[node] = INVALID;
      }
      
      void reach(const Node& node) {
	_numMap.set(node, _num);
	_retMap.set(node, _num);
	++_num;
      }
      
      void leave(const Node& node) {
	if (_numMap[node] <= _retMap[node]) {
	  if (_predMap[node] != INVALID) {
	    _cutMap.set(_predMap[node], true);
	    ++_cutNum;
	  }
	}	
      }

      void discover(const Edge& edge) {
	_predMap.set(_graph.target(edge), edge);
      }

      void examine(const Edge& edge) {
	if (_predMap[_graph.source(edge)] == _graph.oppositeEdge(edge)) {
	  return;
	}
	if (_retMap[_graph.source(edge)] > _retMap[_graph.target(edge)]) {
	  _retMap.set(_graph.source(edge), _retMap[_graph.target(edge)]);
	}
      }

      void backtrack(const Edge& edge) {
	if (_retMap[_graph.source(edge)] > _retMap[_graph.target(edge)]) {
	  _retMap.set(_graph.source(edge), _retMap[_graph.target(edge)]);
	}  
      }
      
    private:
      const Graph& _graph;
      EdgeMap& _cutMap;
      int& _cutNum; 

      typename Graph::template NodeMap<int> _numMap;
      typename Graph::template NodeMap<int> _retMap;
      typename Graph::template NodeMap<Edge> _predMap;
      int _num;
    };
  }

  template <typename UGraph>
  int countBiEdgeConnectedComponents(const UGraph& graph);

  /// \ingroup graph_prop
  ///
  /// \brief Checks that the graph is bi-edge-connected.
  ///
  /// This function checks that the graph is bi-edge-connected. The undirected
  /// graph is bi-edge-connected when any two nodes are connected with two
  /// edge-disjoint paths.
  ///
  /// \param graph The undirected graph.
  /// \return The number of components.
  template <typename UGraph>
  bool biEdgeConnected(const UGraph& graph) { 
    return countBiEdgeConnectedComponents(graph) == 1;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Count the bi-edge-connected components.
  ///
  /// This function count the bi-edge-connected components in an undirected 
  /// graph. The bi-edge-connected components are the classes of an equivalence
  /// relation on the nodes. Two nodes are in relationship when they are  
  /// connected with at least two edge-disjoint paths.
  ///
  /// \param graph The undirected graph.
  /// \return The number of components.
  template <typename UGraph>
  int countBiEdgeConnectedComponents(const UGraph& graph) { 
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::NodeIt NodeIt;

    using namespace _topology_bits;

    typedef CountBiEdgeConnectedComponentsVisitor<UGraph> Visitor;
    
    int compNum = 0;
    Visitor visitor(graph, compNum);

    DfsVisit<UGraph, Visitor> dfs(graph, visitor);
    dfs.init();
    
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	dfs.start();
      }
    }
    return compNum;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Find the bi-edge-connected components.
  ///
  /// This function finds the bi-edge-connected components in an undirected 
  /// graph. The bi-edge-connected components are the classes of an equivalence
  /// relation on the nodes. Two nodes are in relationship when they are  
  /// connected at least two edge-disjoint paths.
  ///
  /// \image html edge_biconnected_components.png
  /// \image latex edge_biconnected_components.eps "bi-edge-connected components" width=\textwidth
  ///
  /// \param graph The graph.
  /// \retval compMap A writable node map. The values will be set from 0 to
  /// the number of the biconnected components minus one. Each values 
  /// of the map will be set exactly once, the values of a certain component 
  /// will be set continuously.
  /// \return The number of components.
  ///
  template <typename UGraph, typename NodeMap>
  int biEdgeConnectedComponents(const UGraph& graph, NodeMap& compMap) { 
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::NodeIt NodeIt;
    typedef typename UGraph::Node Node;
    checkConcept<concepts::WriteMap<Node, int>, NodeMap>();

    using namespace _topology_bits;

    typedef BiEdgeConnectedComponentsVisitor<UGraph, NodeMap> Visitor;
    
    int compNum = 0;
    Visitor visitor(graph, compMap, compNum);

    DfsVisit<UGraph, Visitor> dfs(graph, visitor);
    dfs.init();
    
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	dfs.start();
      }
    }
    return compNum;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Find the bi-edge-connected cut edges.
  ///
  /// This function finds the bi-edge-connected components in an undirected 
  /// graph. The bi-edge-connected components are the classes of an equivalence
  /// relation on the nodes. Two nodes are in relationship when they are 
  /// connected with at least two edge-disjoint paths. The bi-edge-connected 
  /// components are separted by edges which are the cut edges of the 
  /// components.
  ///
  /// \param graph The graph.
  /// \retval cutMap A writable node map. The values will be set true when the
  /// edge is a cut edge.
  /// \return The number of cut edges.
  template <typename UGraph, typename UEdgeMap>
  int biEdgeConnectedCutEdges(const UGraph& graph, UEdgeMap& cutMap) { 
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::NodeIt NodeIt;
    typedef typename UGraph::UEdge UEdge;
    checkConcept<concepts::WriteMap<UEdge, bool>, UEdgeMap>();

    using namespace _topology_bits;

    typedef BiEdgeConnectedCutEdgesVisitor<UGraph, UEdgeMap> Visitor;
    
    int cutNum = 0;
    Visitor visitor(graph, cutMap, cutNum);

    DfsVisit<UGraph, Visitor> dfs(graph, visitor);
    dfs.init();
    
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	dfs.start();
      }
    }
    return cutNum;
  }


  namespace _topology_bits {
    
    template <typename Graph, typename IntNodeMap>
    class TopologicalSortVisitor : public DfsVisitor<Graph> {
    public:
      typedef typename Graph::Node Node;
      typedef typename Graph::Edge edge;

      TopologicalSortVisitor(IntNodeMap& order, int num) 
	: _order(order), _num(num) {}
      
      void leave(const Node& node) {
	_order.set(node, --_num);
      }

    private:
      IntNodeMap& _order;
      int _num;
    };
    
  }

  /// \ingroup graph_prop
  ///
  /// \brief Sort the nodes of a DAG into topolgical order.
  ///
  /// Sort the nodes of a DAG into topolgical order.
  ///
  /// \param graph The graph. It should be directed and acyclic.
  /// \retval order A writable node map. The values will be set from 0 to
  /// the number of the nodes in the graph minus one. Each values of the map
  /// will be set exactly once, the values  will be set descending order.
  ///
  /// \see checkedTopologicalSort
  /// \see dag
  template <typename Graph, typename NodeMap>
  void topologicalSort(const Graph& graph, NodeMap& order) {
    using namespace _topology_bits;

    checkConcept<concepts::Graph, Graph>();
    checkConcept<concepts::WriteMap<typename Graph::Node, int>, NodeMap>();

    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::Edge Edge;

    TopologicalSortVisitor<Graph, NodeMap> 
      visitor(order, countNodes(graph)); 

    DfsVisit<Graph, TopologicalSortVisitor<Graph, NodeMap> >
      dfs(graph, visitor);

    dfs.init();
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	dfs.start();
      }
    }    
  }

  /// \ingroup graph_prop
  ///
  /// \brief Sort the nodes of a DAG into topolgical order.
  ///
  /// Sort the nodes of a DAG into topolgical order. It also checks
  /// that the given graph is DAG.
  ///
  /// \param graph The graph. It should be directed and acyclic.
  /// \retval order A readable - writable node map. The values will be set 
  /// from 0 to the number of the nodes in the graph minus one. Each values 
  /// of the map will be set exactly once, the values will be set descending 
  /// order.
  /// \return %False when the graph is not DAG.
  ///
  /// \see topologicalSort
  /// \see dag
  template <typename Graph, typename NodeMap>
  bool checkedTopologicalSort(const Graph& graph, NodeMap& order) {
    using namespace _topology_bits;

    checkConcept<concepts::Graph, Graph>();
    checkConcept<concepts::ReadWriteMap<typename Graph::Node, int>, NodeMap>();

    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::Edge Edge;

    order = constMap<Node, int, -1>();

    TopologicalSortVisitor<Graph, NodeMap> 
      visitor(order, countNodes(graph)); 

    DfsVisit<Graph, TopologicalSortVisitor<Graph, NodeMap> >
      dfs(graph, visitor);

    dfs.init();
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	while (!dfs.emptyQueue()) {
 	  Edge edge = dfs.nextEdge();
 	  Node target = graph.target(edge);
 	  if (dfs.reached(target) && order[target] == -1) {
 	    return false;
 	  }
 	  dfs.processNextEdge();
 	}
      }
    }
    return true;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Check that the given directed graph is a DAG.
  ///
  /// Check that the given directed graph is a DAG. The DAG is
  /// an Directed Acyclic Graph.
  /// \return %False when the graph is not DAG.
  /// \see acyclic
  template <typename Graph>
  bool dag(const Graph& graph) {

    checkConcept<concepts::Graph, Graph>();

    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::Edge Edge;

    typedef typename Graph::template NodeMap<bool> ProcessedMap;

    typename Dfs<Graph>::template DefProcessedMap<ProcessedMap>::
      Create dfs(graph);

    ProcessedMap processed(graph);
    dfs.processedMap(processed);

    dfs.init();
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	while (!dfs.emptyQueue()) {
	  Edge edge = dfs.nextEdge();
	  Node target = graph.target(edge);
	  if (dfs.reached(target) && !processed[target]) {
	    return false;
	  }
	  dfs.processNextEdge();
	}
      }
    }    
    return true;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Check that the given undirected graph is acyclic.
  ///
  /// Check that the given undirected graph acyclic.
  /// \param graph The undirected graph.
  /// \return %True when there is no circle in the graph.
  /// \see dag
  template <typename UGraph>
  bool acyclic(const UGraph& graph) {
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::Node Node;
    typedef typename UGraph::NodeIt NodeIt;
    typedef typename UGraph::Edge Edge;
    Dfs<UGraph> dfs(graph);
    dfs.init();
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	dfs.addSource(it);
	while (!dfs.emptyQueue()) {
	  Edge edge = dfs.nextEdge();
	  Node source = graph.source(edge);
	  Node target = graph.target(edge);
	  if (dfs.reached(target) && 
	      dfs.predEdge(source) != graph.oppositeEdge(edge)) {
	    return false;
	  }
	  dfs.processNextEdge();
	}
      }
    }
    return true;
  }

  /// \ingroup graph_prop
  ///
  /// \brief Check that the given undirected graph is tree.
  ///
  /// Check that the given undirected graph is tree.
  /// \param graph The undirected graph.
  /// \return %True when the graph is acyclic and connected.
  template <typename UGraph>
  bool tree(const UGraph& graph) {
    checkConcept<concepts::UGraph, UGraph>();
    typedef typename UGraph::Node Node;
    typedef typename UGraph::NodeIt NodeIt;
    typedef typename UGraph::Edge Edge;
    Dfs<UGraph> dfs(graph);
    dfs.init();
    dfs.addSource(NodeIt(graph));
    while (!dfs.emptyQueue()) {
      Edge edge = dfs.nextEdge();
      Node source = graph.source(edge);
      Node target = graph.target(edge);
      if (dfs.reached(target) && 
	  dfs.predEdge(source) != graph.oppositeEdge(edge)) {
	return false;
      }
      dfs.processNextEdge();
    }
    for (NodeIt it(graph); it != INVALID; ++it) {
      if (!dfs.reached(it)) {
	return false;
      }
    }    
    return true;
  }

  namespace _topology_bits {

    template <typename Graph>
    class BipartiteVisitor : public BfsVisitor<Graph> {
    public:
      typedef typename Graph::Edge Edge;
      typedef typename Graph::Node Node;

      BipartiteVisitor(const Graph& graph, bool& bipartite) 
        : _graph(graph), _part(graph), _bipartite(bipartite) {}
      
      void start(const Node& node) {
        _part[node] = true;
      }
      void discover(const Edge& edge) {
        _part.set(_graph.target(edge), !_part[_graph.source(edge)]);
      }
      void examine(const Edge& edge) {
        _bipartite = _bipartite && 
          _part[_graph.target(edge)] != _part[_graph.source(edge)];
      }

    private:

      const Graph& _graph;
      typename Graph::template NodeMap<bool> _part;
      bool& _bipartite;
    };

    template <typename Graph, typename PartMap>
    class BipartitePartitionsVisitor : public BfsVisitor<Graph> {
    public:
      typedef typename Graph::Edge Edge;
      typedef typename Graph::Node Node;

      BipartitePartitionsVisitor(const Graph& graph, 
                                 PartMap& part, bool& bipartite) 
        : _graph(graph), _part(part), _bipartite(bipartite) {}
      
      void start(const Node& node) {
        _part.set(node, true);
      }
      void discover(const Edge& edge) {
        _part.set(_graph.target(edge), !_part[_graph.source(edge)]);
      }
      void examine(const Edge& edge) {
        _bipartite = _bipartite && 
          _part[_graph.target(edge)] != _part[_graph.source(edge)];
      }

    private:

      const Graph& _graph;
      PartMap& _part;
      bool& _bipartite;
    };
  }

  /// \ingroup graph_prop
  ///
  /// \brief Check if the given undirected graph is bipartite or not
  ///
  /// The function checks if the given undirected \c graph graph is bipartite 
  /// or not. The \ref Bfs algorithm is used to calculate the result.
  /// \param graph The undirected graph.
  /// \return %True if \c graph is bipartite, %false otherwise.
  /// \sa bipartitePartitions
  ///
  /// \author Balazs Attila Mihaly  
  template<typename UGraph>
  inline bool bipartite(const UGraph &graph){
    using namespace _topology_bits;

    checkConcept<concepts::UGraph, UGraph>();
    
    typedef typename UGraph::NodeIt NodeIt;
    typedef typename UGraph::EdgeIt EdgeIt;
    
    bool bipartite = true;

    BipartiteVisitor<UGraph> 
      visitor(graph, bipartite);
    BfsVisit<UGraph, BipartiteVisitor<UGraph> > 
      bfs(graph, visitor);
    bfs.init();
    for(NodeIt it(graph); it != INVALID; ++it) {
      if(!bfs.reached(it)){
	bfs.addSource(it);
        while (!bfs.emptyQueue()) {
          bfs.processNextNode();
          if (!bipartite) return false;
        }
      }
    }
    return true;
  }
  
  /// \ingroup graph_prop
  ///
  /// \brief Check if the given undirected graph is bipartite or not
  ///
  /// The function checks if the given undirected graph is bipartite 
  /// or not. The  \ref  Bfs  algorithm  is   used  to  calculate the result. 
  /// During the execution, the \c partMap will be set as the two 
  /// partitions of the graph.
  /// \param graph The undirected graph.
  /// \retval partMap A writable bool map of nodes. It will be set as the
  /// two partitions of the graph. 
  /// \return %True if \c graph is bipartite, %false otherwise.
  ///
  /// \author Balazs Attila Mihaly  
  ///
  /// \image html bipartite_partitions.png
  /// \image latex bipartite_partitions.eps "Bipartite partititions" width=\textwidth
  template<typename UGraph, typename NodeMap>
  inline bool bipartitePartitions(const UGraph &graph, NodeMap &partMap){
    using namespace _topology_bits;

    checkConcept<concepts::UGraph, UGraph>();
    
    typedef typename UGraph::Node Node;
    typedef typename UGraph::NodeIt NodeIt;
    typedef typename UGraph::EdgeIt EdgeIt;

    bool bipartite = true;

    BipartitePartitionsVisitor<UGraph, NodeMap> 
      visitor(graph, partMap, bipartite);
    BfsVisit<UGraph, BipartitePartitionsVisitor<UGraph, NodeMap> > 
      bfs(graph, visitor);
    bfs.init();
    for(NodeIt it(graph); it != INVALID; ++it) {
      if(!bfs.reached(it)){
	bfs.addSource(it);
        while (!bfs.emptyQueue()) {
          bfs.processNextNode();
          if (!bipartite) return false;
        }
      }
    }
    return true;
  }

  /// \brief Returns true when there is not loop edge in the graph.
  ///
  /// Returns true when there is not loop edge in the graph.
  template <typename Graph>
  bool loopFree(const Graph& graph) {
    for (typename Graph::EdgeIt it(graph); it != INVALID; ++it) {
      if (graph.source(it) == graph.target(it)) return false;
    }
    return true;
  }

  /// \brief Returns true when there is not parallel edges in the graph.
  ///
  /// Returns true when there is not parallel edges in the graph.
  template <typename Graph>
  bool parallelFree(const Graph& graph) {
    typename Graph::template NodeMap<bool> reached(graph, false);
    for (typename Graph::NodeIt n(graph); n != INVALID; ++n) {
      for (typename Graph::OutEdgeIt e(graph, n); e != INVALID; ++e) {
        if (reached[graph.target(e)]) return false;
        reached.set(graph.target(e), true);
      }
      for (typename Graph::OutEdgeIt e(graph, n); e != INVALID; ++e) {
        reached.set(graph.target(e), false);
      }
    }
    return true;
  }

  /// \brief Returns true when there is not loop edge and parallel
  /// edges in the graph.
  ///
  /// Returns true when there is not loop edge and parallel edges in
  /// the graph.
  template <typename Graph>
  bool simpleGraph(const Graph& graph) {
    typename Graph::template NodeMap<bool> reached(graph, false);
    for (typename Graph::NodeIt n(graph); n != INVALID; ++n) {
      reached.set(n, true);
      for (typename Graph::OutEdgeIt e(graph, n); e != INVALID; ++e) {
        if (reached[graph.target(e)]) return false;
        reached.set(graph.target(e), true);
      }
      for (typename Graph::OutEdgeIt e(graph, n); e != INVALID; ++e) {
        reached.set(graph.target(e), false);
      }
      reached.set(n, false);
    }
    return true;
  }
   
} //namespace lemon

#endif //LEMON_TOPOLOGY_H

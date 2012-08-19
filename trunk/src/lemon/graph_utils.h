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

#ifndef LEMON_GRAPH_UTILS_H
#define LEMON_GRAPH_UTILS_H

#include <iterator>
#include <vector>
#include <map>
#include <lemon/math.h>
#include <algorithm>

#include <lemon/bits/invalid.h>
#include <lemon/bits/utility.h>
#include <lemon/maps.h>
#include <lemon/bits/traits.h>

#include <lemon/bits/alteration_notifier.h>
#include <lemon/bits/default_map.h>

///\ingroup gutils
///\file
///\brief Graph utilities.

namespace lemon {

  /// \addtogroup gutils
  /// @{

  ///Creates convenience typedefs for the graph types and iterators

  ///This \c \#define creates convenience typedefs for the following types
  ///of \c Graph: \c Node,  \c NodeIt, \c Edge, \c EdgeIt, \c InEdgeIt,
  ///\c OutEdgeIt
  ///\note If \c G it a template parameter, it should be used in this way.
  ///\code
  ///  GRAPH_TYPEDEFS(typename G);
  ///\endcode
  ///
  ///\warning There are no typedefs for the graph maps because of the lack of
  ///template typedefs in C++.
#define GRAPH_TYPEDEFS(Graph)				\
  typedef Graph::     Node      Node;			\
    typedef Graph::   NodeIt    NodeIt;			\
    typedef Graph::   Edge      Edge;			\
    typedef Graph::   EdgeIt    EdgeIt;			\
    typedef Graph:: InEdgeIt  InEdgeIt;			\
    typedef Graph::OutEdgeIt OutEdgeIt

  ///Creates convenience typedefs for the undirected graph types and iterators

  ///This \c \#define creates the same convenience typedefs as defined by
  ///\ref GRAPH_TYPEDEFS(Graph) and three more, namely it creates
  ///\c UEdge, \c UEdgeIt, \c IncEdgeIt,
  ///
  ///\note If \c G it a template parameter, it should be used in this way.
  ///\code
  ///  UGRAPH_TYPEDEFS(typename G);
  ///\endcode
  ///
  ///\warning There are no typedefs for the graph maps because of the lack of
  ///template typedefs in C++.
#define UGRAPH_TYPEDEFS(Graph)				\
  GRAPH_TYPEDEFS(Graph);				\
    typedef Graph:: UEdge   UEdge;			\
    typedef Graph:: UEdgeIt UEdgeIt;			\
    typedef Graph:: IncEdgeIt   IncEdgeIt

  ///\brief Creates convenience typedefs for the bipartite undirected graph 
  ///types and iterators

  ///This \c \#define creates the same convenience typedefs as defined by
  ///\ref UGRAPH_TYPEDEFS(Graph) and two more, namely it creates
  ///\c ANodeIt, \c BNodeIt, 
  ///
  ///\note If \c G it a template parameter, it should be used in this way.
  ///\code
  ///  BPUGRAPH_TYPEDEFS(typename G);
  ///\endcode
  ///
  ///\warning There are no typedefs for the graph maps because of the lack of
  ///template typedefs in C++.
#define BPUGRAPH_TYPEDEFS(Graph)            \
  UGRAPH_TYPEDEFS(Graph);		    \
    typedef Graph::ANode ANode;             \
    typedef Graph::BNode BNode;             \
    typedef Graph::ANodeIt ANodeIt;	    \
    typedef Graph::BNodeIt BNodeIt

  /// \brief Function to count the items in the graph.
  ///
  /// This function counts the items (nodes, edges etc) in the graph.
  /// The complexity of the function is O(n) because
  /// it iterates on all of the items.

  template <typename Graph, typename Item>
  inline int countItems(const Graph& g) {
    typedef typename ItemSetTraits<Graph, Item>::ItemIt ItemIt;
    int num = 0;
    for (ItemIt it(g); it != INVALID; ++it) {
      ++num;
    }
    return num;
  }

  // Node counting:

  namespace _graph_utils_bits {
    
    template <typename Graph, typename Enable = void>
    struct CountNodesSelector {
      static int count(const Graph &g) {
        return countItems<Graph, typename Graph::Node>(g);
      }
    };

    template <typename Graph>
    struct CountNodesSelector<
      Graph, typename 
      enable_if<typename Graph::NodeNumTag, void>::type> 
    {
      static int count(const Graph &g) {
        return g.nodeNum();
      }
    };    
  }

  /// \brief Function to count the nodes in the graph.
  ///
  /// This function counts the nodes in the graph.
  /// The complexity of the function is O(n) but for some
  /// graph structures it is specialized to run in O(1).
  ///
  /// If the graph contains a \e nodeNum() member function and a 
  /// \e NodeNumTag tag then this function calls directly the member
  /// function to query the cardinality of the node set.
  template <typename Graph>
  inline int countNodes(const Graph& g) {
    return _graph_utils_bits::CountNodesSelector<Graph>::count(g);
  }

  namespace _graph_utils_bits {
    
    template <typename Graph, typename Enable = void>
    struct CountANodesSelector {
      static int count(const Graph &g) {
        return countItems<Graph, typename Graph::ANode>(g);
      }
    };

    template <typename Graph>
    struct CountANodesSelector<
      Graph, typename 
      enable_if<typename Graph::NodeNumTag, void>::type> 
    {
      static int count(const Graph &g) {
        return g.aNodeNum();
      }
    };    
  }

  /// \brief Function to count the anodes in the graph.
  ///
  /// This function counts the anodes in the graph.
  /// The complexity of the function is O(an) but for some
  /// graph structures it is specialized to run in O(1).
  ///
  /// If the graph contains an \e aNodeNum() member function and a 
  /// \e NodeNumTag tag then this function calls directly the member
  /// function to query the cardinality of the A-node set.
  template <typename Graph>
  inline int countANodes(const Graph& g) {
    return _graph_utils_bits::CountANodesSelector<Graph>::count(g);
  }

  namespace _graph_utils_bits {
    
    template <typename Graph, typename Enable = void>
    struct CountBNodesSelector {
      static int count(const Graph &g) {
        return countItems<Graph, typename Graph::BNode>(g);
      }
    };

    template <typename Graph>
    struct CountBNodesSelector<
      Graph, typename 
      enable_if<typename Graph::NodeNumTag, void>::type> 
    {
      static int count(const Graph &g) {
        return g.bNodeNum();
      }
    };    
  }

  /// \brief Function to count the bnodes in the graph.
  ///
  /// This function counts the bnodes in the graph.
  /// The complexity of the function is O(bn) but for some
  /// graph structures it is specialized to run in O(1).
  ///
  /// If the graph contains a \e bNodeNum() member function and a 
  /// \e NodeNumTag tag then this function calls directly the member
  /// function to query the cardinality of the B-node set.
  template <typename Graph>
  inline int countBNodes(const Graph& g) {
    return _graph_utils_bits::CountBNodesSelector<Graph>::count(g);
  }


  // Edge counting:

  namespace _graph_utils_bits {
    
    template <typename Graph, typename Enable = void>
    struct CountEdgesSelector {
      static int count(const Graph &g) {
        return countItems<Graph, typename Graph::Edge>(g);
      }
    };

    template <typename Graph>
    struct CountEdgesSelector<
      Graph, 
      typename enable_if<typename Graph::EdgeNumTag, void>::type> 
    {
      static int count(const Graph &g) {
        return g.edgeNum();
      }
    };    
  }

  /// \brief Function to count the edges in the graph.
  ///
  /// This function counts the edges in the graph.
  /// The complexity of the function is O(e) but for some
  /// graph structures it is specialized to run in O(1).
  ///
  /// If the graph contains a \e edgeNum() member function and a 
  /// \e EdgeNumTag tag then this function calls directly the member
  /// function to query the cardinality of the edge set.
  template <typename Graph>
  inline int countEdges(const Graph& g) {
    return _graph_utils_bits::CountEdgesSelector<Graph>::count(g);
  }

  // Undirected edge counting:
  namespace _graph_utils_bits {
    
    template <typename Graph, typename Enable = void>
    struct CountUEdgesSelector {
      static int count(const Graph &g) {
        return countItems<Graph, typename Graph::UEdge>(g);
      }
    };

    template <typename Graph>
    struct CountUEdgesSelector<
      Graph, 
      typename enable_if<typename Graph::EdgeNumTag, void>::type> 
    {
      static int count(const Graph &g) {
        return g.uEdgeNum();
      }
    };    
  }

  /// \brief Function to count the undirected edges in the graph.
  ///
  /// This function counts the undirected edges in the graph.
  /// The complexity of the function is O(e) but for some
  /// graph structures it is specialized to run in O(1).
  ///
  /// If the graph contains a \e uEdgeNum() member function and a 
  /// \e EdgeNumTag tag then this function calls directly the member
  /// function to query the cardinality of the undirected edge set.
  template <typename Graph>
  inline int countUEdges(const Graph& g) {
    return _graph_utils_bits::CountUEdgesSelector<Graph>::count(g);

  }


  template <typename Graph, typename DegIt>
  inline int countNodeDegree(const Graph& _g, const typename Graph::Node& _n) {
    int num = 0;
    for (DegIt it(_g, _n); it != INVALID; ++it) {
      ++num;
    }
    return num;
  }

  /// \brief Function to count the number of the out-edges from node \c n.
  ///
  /// This function counts the number of the out-edges from node \c n
  /// in the graph.  
  template <typename Graph>
  inline int countOutEdges(const Graph& _g,  const typename Graph::Node& _n) {
    return countNodeDegree<Graph, typename Graph::OutEdgeIt>(_g, _n);
  }

  /// \brief Function to count the number of the in-edges to node \c n.
  ///
  /// This function counts the number of the in-edges to node \c n
  /// in the graph.  
  template <typename Graph>
  inline int countInEdges(const Graph& _g,  const typename Graph::Node& _n) {
    return countNodeDegree<Graph, typename Graph::InEdgeIt>(_g, _n);
  }

  /// \brief Function to count the number of the inc-edges to node \c n.
  ///
  /// This function counts the number of the inc-edges to node \c n
  /// in the graph.  
  template <typename Graph>
  inline int countIncEdges(const Graph& _g,  const typename Graph::Node& _n) {
    return countNodeDegree<Graph, typename Graph::IncEdgeIt>(_g, _n);
  }

  namespace _graph_utils_bits {
    
    template <typename Graph, typename Enable = void>
    struct FindEdgeSelector {
      typedef typename Graph::Node Node;
      typedef typename Graph::Edge Edge;
      static Edge find(const Graph &g, Node u, Node v, Edge e) {
        if (e == INVALID) {
          g.firstOut(e, u);
        } else {
          g.nextOut(e);
        }
        while (e != INVALID && g.target(e) != v) {
          g.nextOut(e);
        }
        return e;
      }
    };

    template <typename Graph>
    struct FindEdgeSelector<
      Graph, 
      typename enable_if<typename Graph::FindEdgeTag, void>::type> 
    {
      typedef typename Graph::Node Node;
      typedef typename Graph::Edge Edge;
      static Edge find(const Graph &g, Node u, Node v, Edge prev) {
        return g.findEdge(u, v, prev);
      }
    };    
  }

  /// \brief Finds an edge between two nodes of a graph.
  ///
  /// Finds an edge from node \c u to node \c v in graph \c g.
  ///
  /// If \c prev is \ref INVALID (this is the default value), then
  /// it finds the first edge from \c u to \c v. Otherwise it looks for
  /// the next edge from \c u to \c v after \c prev.
  /// \return The found edge or \ref INVALID if there is no such an edge.
  ///
  /// Thus you can iterate through each edge from \c u to \c v as it follows.
  ///\code
  /// for(Edge e=findEdge(g,u,v);e!=INVALID;e=findEdge(g,u,v,e)) {
  ///   ...
  /// }
  ///\endcode
  ///
  ///\sa EdgeLookUp
  ///\sa AllEdgeLookUp
  ///\sa DynEdgeLookUp
  ///\sa ConEdgeIt
  template <typename Graph>
  inline typename Graph::Edge 
  findEdge(const Graph &g, typename Graph::Node u, typename Graph::Node v,
           typename Graph::Edge prev = INVALID) {
    return _graph_utils_bits::FindEdgeSelector<Graph>::find(g, u, v, prev);
  }

  /// \brief Iterator for iterating on edges connected the same nodes.
  ///
  /// Iterator for iterating on edges connected the same nodes. It is 
  /// higher level interface for the findEdge() function. You can
  /// use it the following way:
  ///\code
  /// for (ConEdgeIt<Graph> it(g, src, trg); it != INVALID; ++it) {
  ///   ...
  /// }
  ///\endcode
  /// 
  ///\sa findEdge()
  ///\sa EdgeLookUp
  ///\sa AllEdgeLookUp
  ///\sa DynEdgeLookUp
  ///
  /// \author Balazs Dezso 
  template <typename _Graph>
  class ConEdgeIt : public _Graph::Edge {
  public:

    typedef _Graph Graph;
    typedef typename Graph::Edge Parent;

    typedef typename Graph::Edge Edge;
    typedef typename Graph::Node Node;

    /// \brief Constructor.
    ///
    /// Construct a new ConEdgeIt iterating on the edges which
    /// connects the \c u and \c v node.
    ConEdgeIt(const Graph& g, Node u, Node v) : graph(g) {
      Parent::operator=(findEdge(graph, u, v));
    }

    /// \brief Constructor.
    ///
    /// Construct a new ConEdgeIt which continues the iterating from 
    /// the \c e edge.
    ConEdgeIt(const Graph& g, Edge e) : Parent(e), graph(g) {}
    
    /// \brief Increment operator.
    ///
    /// It increments the iterator and gives back the next edge.
    ConEdgeIt& operator++() {
      Parent::operator=(findEdge(graph, graph.source(*this), 
				 graph.target(*this), *this));
      return *this;
    }
  private:
    const Graph& graph;
  };

  namespace _graph_utils_bits {
    
    template <typename Graph, typename Enable = void>
    struct FindUEdgeSelector {
      typedef typename Graph::Node Node;
      typedef typename Graph::UEdge UEdge;
      static UEdge find(const Graph &g, Node u, Node v, UEdge e) {
        bool b;
        if (u != v) {
          if (e == INVALID) {
            g.firstInc(e, b, u);
          } else {
            b = g.source(e) == u;
            g.nextInc(e, b);
          }
          while (e != INVALID && (b ? g.target(e) : g.source(e)) != v) {
            g.nextInc(e, b);
          }
        } else {
          if (e == INVALID) {
            g.firstInc(e, b, u);
          } else {
            b = true;
            g.nextInc(e, b);
          }
          while (e != INVALID && (!b || g.target(e) != v)) {
            g.nextInc(e, b);
          }
        }
        return e;
      }
    };

    template <typename Graph>
    struct FindUEdgeSelector<
      Graph, 
      typename enable_if<typename Graph::FindEdgeTag, void>::type> 
    {
      typedef typename Graph::Node Node;
      typedef typename Graph::UEdge UEdge;
      static UEdge find(const Graph &g, Node u, Node v, UEdge prev) {
        return g.findUEdge(u, v, prev);
      }
    };    
  }

  /// \brief Finds an uedge between two nodes of a graph.
  ///
  /// Finds an uedge from node \c u to node \c v in graph \c g.
  /// If the node \c u and node \c v is equal then each loop edge
  /// will be enumerated.
  ///
  /// If \c prev is \ref INVALID (this is the default value), then
  /// it finds the first edge from \c u to \c v. Otherwise it looks for
  /// the next edge from \c u to \c v after \c prev.
  /// \return The found edge or \ref INVALID if there is no such an edge.
  ///
  /// Thus you can iterate through each edge from \c u to \c v as it follows.
  ///\code
  /// for(UEdge e = findUEdge(g,u,v); e != INVALID; 
  ///     e = findUEdge(g,u,v,e)) {
  ///   ...
  /// }
  ///\endcode
  ///
  ///\sa ConEdgeIt

  template <typename Graph>
  inline typename Graph::UEdge 
  findUEdge(const Graph &g, typename Graph::Node u, typename Graph::Node v,
            typename Graph::UEdge p = INVALID) {
    return _graph_utils_bits::FindUEdgeSelector<Graph>::find(g, u, v, p);
  }

  /// \brief Iterator for iterating on uedges connected the same nodes.
  ///
  /// Iterator for iterating on uedges connected the same nodes. It is 
  /// higher level interface for the findUEdge() function. You can
  /// use it the following way:
  ///\code
  /// for (ConUEdgeIt<Graph> it(g, src, trg); it != INVALID; ++it) {
  ///   ...
  /// }
  ///\endcode
  ///
  ///\sa findUEdge()
  ///
  /// \author Balazs Dezso 
  template <typename _Graph>
  class ConUEdgeIt : public _Graph::UEdge {
  public:

    typedef _Graph Graph;
    typedef typename Graph::UEdge Parent;

    typedef typename Graph::UEdge UEdge;
    typedef typename Graph::Node Node;

    /// \brief Constructor.
    ///
    /// Construct a new ConUEdgeIt iterating on the edges which
    /// connects the \c u and \c v node.
    ConUEdgeIt(const Graph& g, Node u, Node v) : graph(g) {
      Parent::operator=(findUEdge(graph, u, v));
    }

    /// \brief Constructor.
    ///
    /// Construct a new ConUEdgeIt which continues the iterating from 
    /// the \c e edge.
    ConUEdgeIt(const Graph& g, UEdge e) : Parent(e), graph(g) {}
    
    /// \brief Increment operator.
    ///
    /// It increments the iterator and gives back the next edge.
    ConUEdgeIt& operator++() {
      Parent::operator=(findUEdge(graph, graph.source(*this), 
				      graph.target(*this), *this));
      return *this;
    }
  private:
    const Graph& graph;
  };

  /// \brief Copy a map.
  ///
  /// This function copies the \c from map to the \c to map. It uses the
  /// given iterator to iterate on the data structure and it uses the \c ref
  /// mapping to convert the from's keys to the to's keys.
  template <typename To, typename From, 
	    typename ItemIt, typename Ref>	    
  void copyMap(To& to, const From& from, 
	       ItemIt it, const Ref& ref) {
    for (; it != INVALID; ++it) {
      to[ref[it]] = from[it];
    }
  }

  /// \brief Copy the from map to the to map.
  ///
  /// Copy the \c from map to the \c to map. It uses the given iterator
  /// to iterate on the data structure.
  template <typename To, typename From, typename ItemIt>	    
  void copyMap(To& to, const From& from, ItemIt it) {
    for (; it != INVALID; ++it) {
      to[it] = from[it];
    }
  }

  namespace _graph_utils_bits {

    template <typename Graph, typename Item, typename RefMap>
    class MapCopyBase {
    public:
      virtual void copy(const Graph& from, const RefMap& refMap) = 0;
      
      virtual ~MapCopyBase() {}
    };

    template <typename Graph, typename Item, typename RefMap, 
              typename ToMap, typename FromMap>
    class MapCopy : public MapCopyBase<Graph, Item, RefMap> {
    public:

      MapCopy(ToMap& tmap, const FromMap& map) 
        : _tmap(tmap), _map(map) {}
      
      virtual void copy(const Graph& graph, const RefMap& refMap) {
        typedef typename ItemSetTraits<Graph, Item>::ItemIt ItemIt;
        for (ItemIt it(graph); it != INVALID; ++it) {
          _tmap.set(refMap[it], _map[it]);
        }
      }

    private:
      ToMap& _tmap;
      const FromMap& _map;
    };

    template <typename Graph, typename Item, typename RefMap, typename It>
    class ItemCopy : public MapCopyBase<Graph, Item, RefMap> {
    public:

      ItemCopy(It& it, const Item& item) : _it(it), _item(item) {}
      
      virtual void copy(const Graph&, const RefMap& refMap) {
        _it = refMap[_item];
      }

    private:
      It& _it;
      Item _item;
    };

    template <typename Graph, typename Item, typename RefMap, typename Ref>
    class RefCopy : public MapCopyBase<Graph, Item, RefMap> {
    public:

      RefCopy(Ref& map) : _map(map) {}
      
      virtual void copy(const Graph& graph, const RefMap& refMap) {
        typedef typename ItemSetTraits<Graph, Item>::ItemIt ItemIt;
        for (ItemIt it(graph); it != INVALID; ++it) {
          _map.set(it, refMap[it]);
        }
      }

    private:
      Ref& _map;
    };

    template <typename Graph, typename Item, typename RefMap, 
              typename CrossRef>
    class CrossRefCopy : public MapCopyBase<Graph, Item, RefMap> {
    public:

      CrossRefCopy(CrossRef& cmap) : _cmap(cmap) {}
      
      virtual void copy(const Graph& graph, const RefMap& refMap) {
        typedef typename ItemSetTraits<Graph, Item>::ItemIt ItemIt;
        for (ItemIt it(graph); it != INVALID; ++it) {
          _cmap.set(refMap[it], it);
        }
      }

    private:
      CrossRef& _cmap;
    };

    template <typename Graph, typename Enable = void>
    struct GraphCopySelector {
      template <typename From, typename NodeRefMap, typename EdgeRefMap>
      static void copy(Graph &to, const From& from,
                       NodeRefMap& nodeRefMap, EdgeRefMap& edgeRefMap) {
        for (typename From::NodeIt it(from); it != INVALID; ++it) {
          nodeRefMap[it] = to.addNode();
        }
        for (typename From::EdgeIt it(from); it != INVALID; ++it) {
          edgeRefMap[it] = to.addEdge(nodeRefMap[from.source(it)], 
                                          nodeRefMap[from.target(it)]);
        }
      }
    };

    template <typename Graph>
    struct GraphCopySelector<
      Graph, 
      typename enable_if<typename Graph::BuildTag, void>::type> 
    {
      template <typename From, typename NodeRefMap, typename EdgeRefMap>
      static void copy(Graph &to, const From& from,
                       NodeRefMap& nodeRefMap, EdgeRefMap& edgeRefMap) {
        to.build(from, nodeRefMap, edgeRefMap);
      }
    };

    template <typename UGraph, typename Enable = void>
    struct UGraphCopySelector {
      template <typename From, typename NodeRefMap, typename UEdgeRefMap>
      static void copy(UGraph &to, const From& from,
                       NodeRefMap& nodeRefMap, UEdgeRefMap& uEdgeRefMap) {
        for (typename From::NodeIt it(from); it != INVALID; ++it) {
          nodeRefMap[it] = to.addNode();
        }
        for (typename From::UEdgeIt it(from); it != INVALID; ++it) {
          uEdgeRefMap[it] = to.addEdge(nodeRefMap[from.source(it)], 
				       nodeRefMap[from.target(it)]);
        }
      }
    };

    template <typename UGraph>
    struct UGraphCopySelector<
      UGraph, 
      typename enable_if<typename UGraph::BuildTag, void>::type> 
    {
      template <typename From, typename NodeRefMap, typename UEdgeRefMap>
      static void copy(UGraph &to, const From& from,
                       NodeRefMap& nodeRefMap, UEdgeRefMap& uEdgeRefMap) {
        to.build(from, nodeRefMap, uEdgeRefMap);
      }
    };

    template <typename BpUGraph, typename Enable = void>
    struct BpUGraphCopySelector {
      template <typename From, typename ANodeRefMap, 
                typename BNodeRefMap, typename UEdgeRefMap>
      static void copy(BpUGraph &to, const From& from,
                       ANodeRefMap& aNodeRefMap, BNodeRefMap& bNodeRefMap,
                       UEdgeRefMap& uEdgeRefMap) {
        for (typename From::ANodeIt it(from); it != INVALID; ++it) {
          aNodeRefMap[it] = to.addANode();
        }
        for (typename From::BNodeIt it(from); it != INVALID; ++it) {
          bNodeRefMap[it] = to.addBNode();
        }
        for (typename From::UEdgeIt it(from); it != INVALID; ++it) {
          uEdgeRefMap[it] = to.addEdge(aNodeRefMap[from.aNode(it)], 
                                           bNodeRefMap[from.bNode(it)]);
        }
      }
    };

    template <typename BpUGraph>
    struct BpUGraphCopySelector<
      BpUGraph, 
      typename enable_if<typename BpUGraph::BuildTag, void>::type> 
    {
      template <typename From, typename ANodeRefMap, 
                typename BNodeRefMap, typename UEdgeRefMap>
      static void copy(BpUGraph &to, const From& from,
                       ANodeRefMap& aNodeRefMap, BNodeRefMap& bNodeRefMap,
                       UEdgeRefMap& uEdgeRefMap) {
        to.build(from, aNodeRefMap, bNodeRefMap, uEdgeRefMap);
      }
    };
    

  }

  /// \brief Class to copy a graph.
  ///
  /// Class to copy a graph to another graph (duplicate a graph). The
  /// simplest way of using it is through the \c copyGraph() function.
  template <typename To, typename From>
  class GraphCopy {
  private:

    typedef typename From::Node Node;
    typedef typename From::NodeIt NodeIt;
    typedef typename From::Edge Edge;
    typedef typename From::EdgeIt EdgeIt;

    typedef typename To::Node TNode;
    typedef typename To::Edge TEdge;

    typedef typename From::template NodeMap<TNode> NodeRefMap;
    typedef typename From::template EdgeMap<TEdge> EdgeRefMap;
    
    
  public: 


    /// \brief Constructor for the GraphCopy.
    ///
    /// It copies the content of the \c _from graph into the
    /// \c _to graph.
    GraphCopy(To& _to, const From& _from) 
      : from(_from), to(_to) {}

    /// \brief Destructor of the GraphCopy
    ///
    /// Destructor of the GraphCopy
    ~GraphCopy() {
      for (int i = 0; i < int(nodeMapCopies.size()); ++i) {
        delete nodeMapCopies[i];
      }
      for (int i = 0; i < int(edgeMapCopies.size()); ++i) {
        delete edgeMapCopies[i];
      }

    }

    /// \brief Copies the node references into the given map.
    ///
    /// Copies the node references into the given map.
    template <typename NodeRef>
    GraphCopy& nodeRef(NodeRef& map) {
      nodeMapCopies.push_back(new _graph_utils_bits::RefCopy<From, Node, 
                              NodeRefMap, NodeRef>(map));
      return *this;
    }

    /// \brief Copies the node cross references into the given map.
    ///
    ///  Copies the node cross references (reverse references) into
    ///  the given map.
    template <typename NodeCrossRef>
    GraphCopy& nodeCrossRef(NodeCrossRef& map) {
      nodeMapCopies.push_back(new _graph_utils_bits::CrossRefCopy<From, Node,
                              NodeRefMap, NodeCrossRef>(map));
      return *this;
    }

    /// \brief Make copy of the given map.
    ///
    /// Makes copy of the given map for the newly created graph. 
    /// The new map's key type is the to graph's node type,
    /// and the copied map's key type is the from graph's node
    /// type.  
    template <typename ToMap, typename FromMap>
    GraphCopy& nodeMap(ToMap& tmap, const FromMap& map) {
      nodeMapCopies.push_back(new _graph_utils_bits::MapCopy<From, Node, 
                              NodeRefMap, ToMap, FromMap>(tmap, map));
      return *this;
    }

    /// \brief Make a copy of the given node.
    ///
    /// Make a copy of the given node.
    GraphCopy& node(TNode& tnode, const Node& snode) {
      nodeMapCopies.push_back(new _graph_utils_bits::ItemCopy<From, Node, 
                              NodeRefMap, TNode>(tnode, snode));
      return *this;
    }

    /// \brief Copies the edge references into the given map.
    ///
    /// Copies the edge references into the given map.
    template <typename EdgeRef>
    GraphCopy& edgeRef(EdgeRef& map) {
      edgeMapCopies.push_back(new _graph_utils_bits::RefCopy<From, Edge, 
                              EdgeRefMap, EdgeRef>(map));
      return *this;
    }

    /// \brief Copies the edge cross references into the given map.
    ///
    ///  Copies the edge cross references (reverse references) into
    ///  the given map.
    template <typename EdgeCrossRef>
    GraphCopy& edgeCrossRef(EdgeCrossRef& map) {
      edgeMapCopies.push_back(new _graph_utils_bits::CrossRefCopy<From, Edge,
                              EdgeRefMap, EdgeCrossRef>(map));
      return *this;
    }

    /// \brief Make copy of the given map.
    ///
    /// Makes copy of the given map for the newly created graph. 
    /// The new map's key type is the to graph's edge type,
    /// and the copied map's key type is the from graph's edge
    /// type.  
    template <typename ToMap, typename FromMap>
    GraphCopy& edgeMap(ToMap& tmap, const FromMap& map) {
      edgeMapCopies.push_back(new _graph_utils_bits::MapCopy<From, Edge, 
                              EdgeRefMap, ToMap, FromMap>(tmap, map));
      return *this;
    }

    /// \brief Make a copy of the given edge.
    ///
    /// Make a copy of the given edge.
    GraphCopy& edge(TEdge& tedge, const Edge& sedge) {
      edgeMapCopies.push_back(new _graph_utils_bits::ItemCopy<From, Edge, 
                              EdgeRefMap, TEdge>(tedge, sedge));
      return *this;
    }

    /// \brief Executes the copies.
    ///
    /// Executes the copies.
    void run() {
      NodeRefMap nodeRefMap(from);
      EdgeRefMap edgeRefMap(from);
      _graph_utils_bits::GraphCopySelector<To>::
        copy(to, from, nodeRefMap, edgeRefMap);
      for (int i = 0; i < int(nodeMapCopies.size()); ++i) {
        nodeMapCopies[i]->copy(from, nodeRefMap);
      }
      for (int i = 0; i < int(edgeMapCopies.size()); ++i) {
        edgeMapCopies[i]->copy(from, edgeRefMap);
      }      
    }

  protected:


    const From& from;
    To& to;

    std::vector<_graph_utils_bits::MapCopyBase<From, Node, NodeRefMap>* > 
    nodeMapCopies;

    std::vector<_graph_utils_bits::MapCopyBase<From, Edge, EdgeRefMap>* > 
    edgeMapCopies;

  };

  /// \brief Copy a graph to another graph.
  ///
  /// Copy a graph to another graph.
  /// The usage of the function:
  /// 
  ///\code
  /// copyGraph(trg, src).nodeRef(nr).edgeCrossRef(ecr).run();
  ///\endcode
  /// 
  /// After the copy the \c nr map will contain the mapping from the
  /// nodes of the \c from graph to the nodes of the \c to graph and
  /// \c ecr will contain the mapping from the edges of the \c to graph
  /// to the edges of the \c from graph.
  ///
  /// \see GraphCopy 
  template <typename To, typename From>
  GraphCopy<To, From> copyGraph(To& to, const From& from) {
    return GraphCopy<To, From>(to, from);
  }

  /// \brief Class to copy an undirected graph.
  ///
  /// Class to copy an undirected graph to another graph (duplicate a graph).
  /// The simplest way of using it is through the \c copyUGraph() function.
  template <typename To, typename From>
  class UGraphCopy {
  private:

    typedef typename From::Node Node;
    typedef typename From::NodeIt NodeIt;
    typedef typename From::Edge Edge;
    typedef typename From::EdgeIt EdgeIt;
    typedef typename From::UEdge UEdge;
    typedef typename From::UEdgeIt UEdgeIt;

    typedef typename To::Node TNode;
    typedef typename To::Edge TEdge;
    typedef typename To::UEdge TUEdge;

    typedef typename From::template NodeMap<TNode> NodeRefMap;
    typedef typename From::template UEdgeMap<TUEdge> UEdgeRefMap;

    struct EdgeRefMap {
      EdgeRefMap(const To& _to, const From& _from,
                 const UEdgeRefMap& _uedge_ref, const NodeRefMap& _node_ref) 
        : to(_to), from(_from), 
          uedge_ref(_uedge_ref), node_ref(_node_ref) {}

      typedef typename From::Edge Key;
      typedef typename To::Edge Value;

      Value operator[](const Key& key) const {
        bool forward = 
          (from.direction(key) == 
           (node_ref[from.source(static_cast<const UEdge&>(key))] == 
            to.source(uedge_ref[static_cast<const UEdge&>(key)])));
	return to.direct(uedge_ref[key], forward); 
      }
      
      const To& to;
      const From& from;
      const UEdgeRefMap& uedge_ref;
      const NodeRefMap& node_ref;
    };

    
  public: 


    /// \brief Constructor for the GraphCopy.
    ///
    /// It copies the content of the \c _from graph into the
    /// \c _to graph.
    UGraphCopy(To& _to, const From& _from) 
      : from(_from), to(_to) {}

    /// \brief Destructor of the GraphCopy
    ///
    /// Destructor of the GraphCopy
    ~UGraphCopy() {
      for (int i = 0; i < int(nodeMapCopies.size()); ++i) {
        delete nodeMapCopies[i];
      }
      for (int i = 0; i < int(edgeMapCopies.size()); ++i) {
        delete edgeMapCopies[i];
      }
      for (int i = 0; i < int(uEdgeMapCopies.size()); ++i) {
        delete uEdgeMapCopies[i];
      }

    }

    /// \brief Copies the node references into the given map.
    ///
    /// Copies the node references into the given map.
    template <typename NodeRef>
    UGraphCopy& nodeRef(NodeRef& map) {
      nodeMapCopies.push_back(new _graph_utils_bits::RefCopy<From, Node, 
                              NodeRefMap, NodeRef>(map));
      return *this;
    }

    /// \brief Copies the node cross references into the given map.
    ///
    ///  Copies the node cross references (reverse references) into
    ///  the given map.
    template <typename NodeCrossRef>
    UGraphCopy& nodeCrossRef(NodeCrossRef& map) {
      nodeMapCopies.push_back(new _graph_utils_bits::CrossRefCopy<From, Node,
                              NodeRefMap, NodeCrossRef>(map));
      return *this;
    }

    /// \brief Make copy of the given map.
    ///
    /// Makes copy of the given map for the newly created graph. 
    /// The new map's key type is the to graph's node type,
    /// and the copied map's key type is the from graph's node
    /// type.  
    template <typename ToMap, typename FromMap>
    UGraphCopy& nodeMap(ToMap& tmap, const FromMap& map) {
      nodeMapCopies.push_back(new _graph_utils_bits::MapCopy<From, Node, 
                              NodeRefMap, ToMap, FromMap>(tmap, map));
      return *this;
    }

    /// \brief Make a copy of the given node.
    ///
    /// Make a copy of the given node.
    UGraphCopy& node(TNode& tnode, const Node& snode) {
      nodeMapCopies.push_back(new _graph_utils_bits::ItemCopy<From, Node, 
                              NodeRefMap, TNode>(tnode, snode));
      return *this;
    }

    /// \brief Copies the edge references into the given map.
    ///
    /// Copies the edge references into the given map.
    template <typename EdgeRef>
    UGraphCopy& edgeRef(EdgeRef& map) {
      edgeMapCopies.push_back(new _graph_utils_bits::RefCopy<From, Edge, 
                              EdgeRefMap, EdgeRef>(map));
      return *this;
    }

    /// \brief Copies the edge cross references into the given map.
    ///
    ///  Copies the edge cross references (reverse references) into
    ///  the given map.
    template <typename EdgeCrossRef>
    UGraphCopy& edgeCrossRef(EdgeCrossRef& map) {
      edgeMapCopies.push_back(new _graph_utils_bits::CrossRefCopy<From, Edge,
                              EdgeRefMap, EdgeCrossRef>(map));
      return *this;
    }

    /// \brief Make copy of the given map.
    ///
    /// Makes copy of the given map for the newly created graph. 
    /// The new map's key type is the to graph's edge type,
    /// and the copied map's key type is the from graph's edge
    /// type.  
    template <typename ToMap, typename FromMap>
    UGraphCopy& edgeMap(ToMap& tmap, const FromMap& map) {
      edgeMapCopies.push_back(new _graph_utils_bits::MapCopy<From, Edge, 
                              EdgeRefMap, ToMap, FromMap>(tmap, map));
      return *this;
    }

    /// \brief Make a copy of the given edge.
    ///
    /// Make a copy of the given edge.
    UGraphCopy& edge(TEdge& tedge, const Edge& sedge) {
      edgeMapCopies.push_back(new _graph_utils_bits::ItemCopy<From, Edge, 
                              EdgeRefMap, TEdge>(tedge, sedge));
      return *this;
    }

    /// \brief Copies the undirected edge references into the given map.
    ///
    /// Copies the undirected edge references into the given map.
    template <typename UEdgeRef>
    UGraphCopy& uEdgeRef(UEdgeRef& map) {
      uEdgeMapCopies.push_back(new _graph_utils_bits::RefCopy<From, UEdge, 
                               UEdgeRefMap, UEdgeRef>(map));
      return *this;
    }

    /// \brief Copies the undirected edge cross references into the given map.
    ///
    /// Copies the undirected edge cross references (reverse
    /// references) into the given map.
    template <typename UEdgeCrossRef>
    UGraphCopy& uEdgeCrossRef(UEdgeCrossRef& map) {
      uEdgeMapCopies.push_back(new _graph_utils_bits::CrossRefCopy<From, 
                               UEdge, UEdgeRefMap, UEdgeCrossRef>(map));
      return *this;
    }

    /// \brief Make copy of the given map.
    ///
    /// Makes copy of the given map for the newly created graph. 
    /// The new map's key type is the to graph's undirected edge type,
    /// and the copied map's key type is the from graph's undirected edge
    /// type.  
    template <typename ToMap, typename FromMap>
    UGraphCopy& uEdgeMap(ToMap& tmap, const FromMap& map) {
      uEdgeMapCopies.push_back(new _graph_utils_bits::MapCopy<From, UEdge, 
                               UEdgeRefMap, ToMap, FromMap>(tmap, map));
      return *this;
    }

    /// \brief Make a copy of the given undirected edge.
    ///
    /// Make a copy of the given undirected edge.
    UGraphCopy& uEdge(TUEdge& tuedge, const UEdge& suedge) {
      uEdgeMapCopies.push_back(new _graph_utils_bits::ItemCopy<From, UEdge, 
                               UEdgeRefMap, TUEdge>(tuedge, suedge));
      return *this;
    }

    /// \brief Executes the copies.
    ///
    /// Executes the copies.
    void run() {
      NodeRefMap nodeRefMap(from);
      UEdgeRefMap uEdgeRefMap(from);
      EdgeRefMap edgeRefMap(to, from, uEdgeRefMap, nodeRefMap);
      _graph_utils_bits::UGraphCopySelector<To>::
        copy(to, from, nodeRefMap, uEdgeRefMap);
      for (int i = 0; i < int(nodeMapCopies.size()); ++i) {
        nodeMapCopies[i]->copy(from, nodeRefMap);
      }
      for (int i = 0; i < int(uEdgeMapCopies.size()); ++i) {
        uEdgeMapCopies[i]->copy(from, uEdgeRefMap);
      }
      for (int i = 0; i < int(edgeMapCopies.size()); ++i) {
        edgeMapCopies[i]->copy(from, edgeRefMap);
      }
    }

  private:
    
    const From& from;
    To& to;

    std::vector<_graph_utils_bits::MapCopyBase<From, Node, NodeRefMap>* > 
    nodeMapCopies;

    std::vector<_graph_utils_bits::MapCopyBase<From, Edge, EdgeRefMap>* > 
    edgeMapCopies;

    std::vector<_graph_utils_bits::MapCopyBase<From, UEdge, UEdgeRefMap>* > 
    uEdgeMapCopies;

  };

  /// \brief Copy an undirected graph to another graph.
  ///
  /// Copy an undirected graph to another graph.
  /// The usage of the function:
  /// 
  ///\code
  /// copyUGraph(trg, src).nodeRef(nr).edgeCrossRef(ecr).run();
  ///\endcode
  /// 
  /// After the copy the \c nr map will contain the mapping from the
  /// nodes of the \c from graph to the nodes of the \c to graph and
  /// \c ecr will contain the mapping from the edges of the \c to graph
  /// to the edges of the \c from graph.
  ///
  /// \see UGraphCopy 
  template <typename To, typename From>
  UGraphCopy<To, From> 
  copyUGraph(To& to, const From& from) {
    return UGraphCopy<To, From>(to, from);
  }

  /// \brief Class to copy a bipartite undirected graph.
  ///
  /// Class to copy a bipartite undirected graph to another graph
  /// (duplicate a graph).  The simplest way of using it is through
  /// the \c copyBpUGraph() function.
  template <typename To, typename From>
  class BpUGraphCopy {
  private:

    typedef typename From::Node Node;
    typedef typename From::ANode ANode;
    typedef typename From::BNode BNode;
    typedef typename From::NodeIt NodeIt;
    typedef typename From::Edge Edge;
    typedef typename From::EdgeIt EdgeIt;
    typedef typename From::UEdge UEdge;
    typedef typename From::UEdgeIt UEdgeIt;

    typedef typename To::Node TNode;
    typedef typename To::Edge TEdge;
    typedef typename To::UEdge TUEdge;

    typedef typename From::template ANodeMap<TNode> ANodeRefMap;
    typedef typename From::template BNodeMap<TNode> BNodeRefMap;
    typedef typename From::template UEdgeMap<TUEdge> UEdgeRefMap;

    struct NodeRefMap {
      NodeRefMap(const From& _from, const ANodeRefMap& _anode_ref,
                 const BNodeRefMap& _bnode_ref)
        : from(_from), anode_ref(_anode_ref), bnode_ref(_bnode_ref) {}

      typedef typename From::Node Key;
      typedef typename To::Node Value;

      Value operator[](const Key& key) const {
	return from.aNode(key) ? anode_ref[key] : bnode_ref[key]; 
      }
      
      const From& from;
      const ANodeRefMap& anode_ref;
      const BNodeRefMap& bnode_ref;
    };

    struct EdgeRefMap {
      EdgeRefMap(const To& _to, const From& _from,
                 const UEdgeRefMap& _uedge_ref, const NodeRefMap& _node_ref) 
        : to(_to), from(_from), 
          uedge_ref(_uedge_ref), node_ref(_node_ref) {}

      typedef typename From::Edge Key;
      typedef typename To::Edge Value;

      Value operator[](const Key& key) const {
        bool forward = 
          (from.direction(key) == 
           (node_ref[from.source(static_cast<const UEdge&>(key))] == 
            to.source(uedge_ref[static_cast<const UEdge&>(key)])));
	return to.direct(uedge_ref[key], forward); 
      }
      
      const To& to;
      const From& from;
      const UEdgeRefMap& uedge_ref;
      const NodeRefMap& node_ref;
    };
    
  public: 


    /// \brief Constructor for the GraphCopy.
    ///
    /// It copies the content of the \c _from graph into the
    /// \c _to graph.
    BpUGraphCopy(To& _to, const From& _from) 
      : from(_from), to(_to) {}

    /// \brief Destructor of the GraphCopy
    ///
    /// Destructor of the GraphCopy
    ~BpUGraphCopy() {
      for (int i = 0; i < int(aNodeMapCopies.size()); ++i) {
        delete aNodeMapCopies[i];
      }
      for (int i = 0; i < int(bNodeMapCopies.size()); ++i) {
        delete bNodeMapCopies[i];
      }
      for (int i = 0; i < int(nodeMapCopies.size()); ++i) {
        delete nodeMapCopies[i];
      }
      for (int i = 0; i < int(edgeMapCopies.size()); ++i) {
        delete edgeMapCopies[i];
      }
      for (int i = 0; i < int(uEdgeMapCopies.size()); ++i) {
        delete uEdgeMapCopies[i];
      }

    }

    /// \brief Copies the A-node references into the given map.
    ///
    /// Copies the A-node references into the given map.
    template <typename ANodeRef>
    BpUGraphCopy& aNodeRef(ANodeRef& map) {
      aNodeMapCopies.push_back(new _graph_utils_bits::RefCopy<From, ANode, 
                               ANodeRefMap, ANodeRef>(map));
      return *this;
    }

    /// \brief Copies the A-node cross references into the given map.
    ///
    /// Copies the A-node cross references (reverse references) into
    /// the given map.
    template <typename ANodeCrossRef>
    BpUGraphCopy& aNodeCrossRef(ANodeCrossRef& map) {
      aNodeMapCopies.push_back(new _graph_utils_bits::CrossRefCopy<From, 
                               ANode, ANodeRefMap, ANodeCrossRef>(map));
      return *this;
    }

    /// \brief Make copy of the given A-node map.
    ///
    /// Makes copy of the given map for the newly created graph. 
    /// The new map's key type is the to graph's node type,
    /// and the copied map's key type is the from graph's node
    /// type.  
    template <typename ToMap, typename FromMap>
    BpUGraphCopy& aNodeMap(ToMap& tmap, const FromMap& map) {
      aNodeMapCopies.push_back(new _graph_utils_bits::MapCopy<From, ANode, 
                               ANodeRefMap, ToMap, FromMap>(tmap, map));
      return *this;
    }

    /// \brief Copies the B-node references into the given map.
    ///
    /// Copies the B-node references into the given map.
    template <typename BNodeRef>
    BpUGraphCopy& bNodeRef(BNodeRef& map) {
      bNodeMapCopies.push_back(new _graph_utils_bits::RefCopy<From, BNode, 
                               BNodeRefMap, BNodeRef>(map));
      return *this;
    }

    /// \brief Copies the B-node cross references into the given map.
    ///
    ///  Copies the B-node cross references (reverse references) into
    ///  the given map.
    template <typename BNodeCrossRef>
    BpUGraphCopy& bNodeCrossRef(BNodeCrossRef& map) {
      bNodeMapCopies.push_back(new _graph_utils_bits::CrossRefCopy<From, 
                              BNode, BNodeRefMap, BNodeCrossRef>(map));
      return *this;
    }

    /// \brief Make copy of the given B-node map.
    ///
    /// Makes copy of the given map for the newly created graph. 
    /// The new map's key type is the to graph's node type,
    /// and the copied map's key type is the from graph's node
    /// type.  
    template <typename ToMap, typename FromMap>
    BpUGraphCopy& bNodeMap(ToMap& tmap, const FromMap& map) {
      bNodeMapCopies.push_back(new _graph_utils_bits::MapCopy<From, BNode, 
                               BNodeRefMap, ToMap, FromMap>(tmap, map));
      return *this;
    }
    /// \brief Copies the node references into the given map.
    ///
    /// Copies the node references into the given map.
    template <typename NodeRef>
    BpUGraphCopy& nodeRef(NodeRef& map) {
      nodeMapCopies.push_back(new _graph_utils_bits::RefCopy<From, Node, 
                              NodeRefMap, NodeRef>(map));
      return *this;
    }

    /// \brief Copies the node cross references into the given map.
    ///
    ///  Copies the node cross references (reverse references) into
    ///  the given map.
    template <typename NodeCrossRef>
    BpUGraphCopy& nodeCrossRef(NodeCrossRef& map) {
      nodeMapCopies.push_back(new _graph_utils_bits::CrossRefCopy<From, Node,
                              NodeRefMap, NodeCrossRef>(map));
      return *this;
    }

    /// \brief Make copy of the given map.
    ///
    /// Makes copy of the given map for the newly created graph. 
    /// The new map's key type is the to graph's node type,
    /// and the copied map's key type is the from graph's node
    /// type.  
    template <typename ToMap, typename FromMap>
    BpUGraphCopy& nodeMap(ToMap& tmap, const FromMap& map) {
      nodeMapCopies.push_back(new _graph_utils_bits::MapCopy<From, Node, 
                              NodeRefMap, ToMap, FromMap>(tmap, map));
      return *this;
    }

    /// \brief Make a copy of the given node.
    ///
    /// Make a copy of the given node.
    BpUGraphCopy& node(TNode& tnode, const Node& snode) {
      nodeMapCopies.push_back(new _graph_utils_bits::ItemCopy<From, Node, 
                              NodeRefMap, TNode>(tnode, snode));
      return *this;
    }

    /// \brief Copies the edge references into the given map.
    ///
    /// Copies the edge references into the given map.
    template <typename EdgeRef>
    BpUGraphCopy& edgeRef(EdgeRef& map) {
      edgeMapCopies.push_back(new _graph_utils_bits::RefCopy<From, Edge, 
                              EdgeRefMap, EdgeRef>(map));
      return *this;
    }

    /// \brief Copies the edge cross references into the given map.
    ///
    ///  Copies the edge cross references (reverse references) into
    ///  the given map.
    template <typename EdgeCrossRef>
    BpUGraphCopy& edgeCrossRef(EdgeCrossRef& map) {
      edgeMapCopies.push_back(new _graph_utils_bits::CrossRefCopy<From, Edge,
                              EdgeRefMap, EdgeCrossRef>(map));
      return *this;
    }

    /// \brief Make copy of the given map.
    ///
    /// Makes copy of the given map for the newly created graph. 
    /// The new map's key type is the to graph's edge type,
    /// and the copied map's key type is the from graph's edge
    /// type.  
    template <typename ToMap, typename FromMap>
    BpUGraphCopy& edgeMap(ToMap& tmap, const FromMap& map) {
      edgeMapCopies.push_back(new _graph_utils_bits::MapCopy<From, Edge, 
                              EdgeRefMap, ToMap, FromMap>(tmap, map));
      return *this;
    }

    /// \brief Make a copy of the given edge.
    ///
    /// Make a copy of the given edge.
    BpUGraphCopy& edge(TEdge& tedge, const Edge& sedge) {
      edgeMapCopies.push_back(new _graph_utils_bits::ItemCopy<From, Edge, 
                              EdgeRefMap, TEdge>(tedge, sedge));
      return *this;
    }

    /// \brief Copies the undirected edge references into the given map.
    ///
    /// Copies the undirected edge references into the given map.
    template <typename UEdgeRef>
    BpUGraphCopy& uEdgeRef(UEdgeRef& map) {
      uEdgeMapCopies.push_back(new _graph_utils_bits::RefCopy<From, UEdge, 
                               UEdgeRefMap, UEdgeRef>(map));
      return *this;
    }

    /// \brief Copies the undirected edge cross references into the given map.
    ///
    /// Copies the undirected edge cross references (reverse
    /// references) into the given map.
    template <typename UEdgeCrossRef>
    BpUGraphCopy& uEdgeCrossRef(UEdgeCrossRef& map) {
      uEdgeMapCopies.push_back(new _graph_utils_bits::CrossRefCopy<From, 
                               UEdge, UEdgeRefMap, UEdgeCrossRef>(map));
      return *this;
    }

    /// \brief Make copy of the given map.
    ///
    /// Makes copy of the given map for the newly created graph. 
    /// The new map's key type is the to graph's undirected edge type,
    /// and the copied map's key type is the from graph's undirected edge
    /// type.  
    template <typename ToMap, typename FromMap>
    BpUGraphCopy& uEdgeMap(ToMap& tmap, const FromMap& map) {
      uEdgeMapCopies.push_back(new _graph_utils_bits::MapCopy<From, UEdge, 
                               UEdgeRefMap, ToMap, FromMap>(tmap, map));
      return *this;
    }

    /// \brief Make a copy of the given undirected edge.
    ///
    /// Make a copy of the given undirected edge.
    BpUGraphCopy& uEdge(TUEdge& tuedge, const UEdge& suedge) {
      uEdgeMapCopies.push_back(new _graph_utils_bits::ItemCopy<From, UEdge, 
                               UEdgeRefMap, TUEdge>(tuedge, suedge));
      return *this;
    }

    /// \brief Executes the copies.
    ///
    /// Executes the copies.
    void run() {
      ANodeRefMap aNodeRefMap(from);
      BNodeRefMap bNodeRefMap(from);
      NodeRefMap nodeRefMap(from, aNodeRefMap, bNodeRefMap);
      UEdgeRefMap uEdgeRefMap(from);
      EdgeRefMap edgeRefMap(to, from, uEdgeRefMap, nodeRefMap);
      _graph_utils_bits::BpUGraphCopySelector<To>::
        copy(to, from, aNodeRefMap, bNodeRefMap, uEdgeRefMap);
      for (int i = 0; i < int(aNodeMapCopies.size()); ++i) {
        aNodeMapCopies[i]->copy(from, aNodeRefMap);
      }
      for (int i = 0; i < int(bNodeMapCopies.size()); ++i) {
        bNodeMapCopies[i]->copy(from, bNodeRefMap);
      }
      for (int i = 0; i < int(nodeMapCopies.size()); ++i) {
        nodeMapCopies[i]->copy(from, nodeRefMap);
      }
      for (int i = 0; i < int(uEdgeMapCopies.size()); ++i) {
        uEdgeMapCopies[i]->copy(from, uEdgeRefMap);
      }
      for (int i = 0; i < int(edgeMapCopies.size()); ++i) {
        edgeMapCopies[i]->copy(from, edgeRefMap);
      }
    }

  private:
    
    const From& from;
    To& to;

    std::vector<_graph_utils_bits::MapCopyBase<From, ANode, ANodeRefMap>* > 
    aNodeMapCopies;

    std::vector<_graph_utils_bits::MapCopyBase<From, BNode, BNodeRefMap>* > 
    bNodeMapCopies;

    std::vector<_graph_utils_bits::MapCopyBase<From, Node, NodeRefMap>* > 
    nodeMapCopies;

    std::vector<_graph_utils_bits::MapCopyBase<From, Edge, EdgeRefMap>* > 
    edgeMapCopies;

    std::vector<_graph_utils_bits::MapCopyBase<From, UEdge, UEdgeRefMap>* > 
    uEdgeMapCopies;

  };

  /// \brief Copy a bipartite undirected graph to another graph.
  ///
  /// Copy a bipartite undirected graph to another graph.
  /// The usage of the function:
  /// 
  ///\code
  /// copyBpUGraph(trg, src).aNodeRef(anr).edgeCrossRef(ecr).run();
  ///\endcode
  /// 
  /// After the copy the \c nr map will contain the mapping from the
  /// nodes of the \c from graph to the nodes of the \c to graph and
  /// \c ecr will contain the mapping from the edges of the \c to graph
  /// to the edges of the \c from graph.
  ///
  /// \see BpUGraphCopy
  template <typename To, typename From>
  BpUGraphCopy<To, From> 
  copyBpUGraph(To& to, const From& from) {
    return BpUGraphCopy<To, From>(to, from);
  }


  /// @}

  /// \addtogroup graph_maps
  /// @{

  /// Provides an immutable and unique id for each item in the graph.

  /// The IdMap class provides a unique and immutable id for each item of the
  /// same type (e.g. node) in the graph. This id is <ul><li>\b unique:
  /// different items (nodes) get different ids <li>\b immutable: the id of an
  /// item (node) does not change (even if you delete other nodes).  </ul>
  /// Through this map you get access (i.e. can read) the inner id values of
  /// the items stored in the graph. This map can be inverted with its member
  /// class \c InverseMap.
  ///
  template <typename _Graph, typename _Item>
  class IdMap {
  public:
    typedef _Graph Graph;
    typedef int Value;
    typedef _Item Item;
    typedef _Item Key;

    /// \brief Constructor.
    ///
    /// Constructor of the map.
    explicit IdMap(const Graph& _graph) : graph(&_graph) {}

    /// \brief Gives back the \e id of the item.
    ///
    /// Gives back the immutable and unique \e id of the item.
    int operator[](const Item& item) const { return graph->id(item);}

    /// \brief Gives back the item by its id.
    ///
    /// Gives back the item by its id.
    Item operator()(int id) { return graph->fromId(id, Item()); }

  private:
    const Graph* graph;

  public:

    /// \brief The class represents the inverse of its owner (IdMap).
    ///
    /// The class represents the inverse of its owner (IdMap).
    /// \see inverse()
    class InverseMap {
    public:

      /// \brief Constructor.
      ///
      /// Constructor for creating an id-to-item map.
      explicit InverseMap(const Graph& _graph) : graph(&_graph) {}

      /// \brief Constructor.
      ///
      /// Constructor for creating an id-to-item map.
      explicit InverseMap(const IdMap& idMap) : graph(idMap.graph) {}

      /// \brief Gives back the given item from its id.
      ///
      /// Gives back the given item from its id.
      /// 
      Item operator[](int id) const { return graph->fromId(id, Item());}

    private:
      const Graph* graph;
    };

    /// \brief Gives back the inverse of the map.
    ///
    /// Gives back the inverse of the IdMap.
    InverseMap inverse() const { return InverseMap(*graph);} 

  };

  
  /// \brief General invertable graph-map type.

  /// This type provides simple invertable graph-maps. 
  /// The InvertableMap wraps an arbitrary ReadWriteMap 
  /// and if a key is set to a new value then store it
  /// in the inverse map.
  ///
  /// The values of the map can be accessed
  /// with stl compatible forward iterator.
  ///
  /// \param _Graph The graph type.
  /// \param _Item The item type of the graph.
  /// \param _Value The value type of the map.
  ///
  /// \see IterableValueMap
  template <typename _Graph, typename _Item, typename _Value>
  class InvertableMap : protected DefaultMap<_Graph, _Item, _Value> {
  private:
    
    typedef DefaultMap<_Graph, _Item, _Value> Map;
    typedef _Graph Graph;

    typedef std::map<_Value, _Item> Container;
    Container invMap;    

  public:
 
    /// The key type of InvertableMap (Node, Edge, UEdge).
    typedef typename Map::Key Key;
    /// The value type of the InvertableMap.
    typedef typename Map::Value Value;



    /// \brief Constructor.
    ///
    /// Construct a new InvertableMap for the graph.
    ///
    explicit InvertableMap(const Graph& graph) : Map(graph) {} 

    /// \brief Forward iterator for values.
    ///
    /// This iterator is an stl compatible forward
    /// iterator on the values of the map. The values can
    /// be accessed in the [beginValue, endValue) range.
    ///
    class ValueIterator 
      : public std::iterator<std::forward_iterator_tag, Value> {
      friend class InvertableMap;
    private:
      ValueIterator(typename Container::const_iterator _it) 
        : it(_it) {}
    public:
      
      ValueIterator() {}

      ValueIterator& operator++() { ++it; return *this; }
      ValueIterator operator++(int) { 
        ValueIterator tmp(*this); 
        operator++();
        return tmp; 
      }

      const Value& operator*() const { return it->first; }
      const Value* operator->() const { return &(it->first); }

      bool operator==(ValueIterator jt) const { return it == jt.it; }
      bool operator!=(ValueIterator jt) const { return it != jt.it; }
      
    private:
      typename Container::const_iterator it;
    };

    /// \brief Returns an iterator to the first value.
    ///
    /// Returns an stl compatible iterator to the 
    /// first value of the map. The values of the
    /// map can be accessed in the [beginValue, endValue)
    /// range.
    ValueIterator beginValue() const {
      return ValueIterator(invMap.begin());
    }

    /// \brief Returns an iterator after the last value.
    ///
    /// Returns an stl compatible iterator after the 
    /// last value of the map. The values of the
    /// map can be accessed in the [beginValue, endValue)
    /// range.
    ValueIterator endValue() const {
      return ValueIterator(invMap.end());
    }
    
    /// \brief The setter function of the map.
    ///
    /// Sets the mapped value.
    void set(const Key& key, const Value& val) {
      Value oldval = Map::operator[](key);
      typename Container::iterator it = invMap.find(oldval);
      if (it != invMap.end() && it->second == key) {
	invMap.erase(it);
      }      
      invMap.insert(make_pair(val, key));
      Map::set(key, val);
    }

    /// \brief The getter function of the map.
    ///
    /// It gives back the value associated with the key.
    typename MapTraits<Map>::ConstReturnValue 
    operator[](const Key& key) const {
      return Map::operator[](key);
    }

    /// \brief Gives back the item by its value.
    ///
    /// Gives back the item by its value.
    Key operator()(const Value& key) const {
      typename Container::const_iterator it = invMap.find(key);
      return it != invMap.end() ? it->second : INVALID;
    }

  protected:

    /// \brief Erase the key from the map.
    ///
    /// Erase the key to the map. It is called by the
    /// \c AlterationNotifier.
    virtual void erase(const Key& key) {
      Value val = Map::operator[](key);
      typename Container::iterator it = invMap.find(val);
      if (it != invMap.end() && it->second == key) {
	invMap.erase(it);
      }
      Map::erase(key);
    }

    /// \brief Erase more keys from the map.
    ///
    /// Erase more keys from the map. It is called by the
    /// \c AlterationNotifier.
    virtual void erase(const std::vector<Key>& keys) {
      for (int i = 0; i < int(keys.size()); ++i) {
	Value val = Map::operator[](keys[i]);
	typename Container::iterator it = invMap.find(val);
	if (it != invMap.end() && it->second == keys[i]) {
	  invMap.erase(it);
	}
      }
      Map::erase(keys);
    }

    /// \brief Clear the keys from the map and inverse map.
    ///
    /// Clear the keys from the map and inverse map. It is called by the
    /// \c AlterationNotifier.
    virtual void clear() {
      invMap.clear();
      Map::clear();
    }

  public:

    /// \brief The inverse map type.
    ///
    /// The inverse of this map. The subscript operator of the map
    /// gives back always the item what was last assigned to the value. 
    class InverseMap {
    public:
      /// \brief Constructor of the InverseMap.
      ///
      /// Constructor of the InverseMap.
      explicit InverseMap(const InvertableMap& _inverted) 
        : inverted(_inverted) {}

      /// The value type of the InverseMap.
      typedef typename InvertableMap::Key Value;
      /// The key type of the InverseMap.
      typedef typename InvertableMap::Value Key; 

      /// \brief Subscript operator. 
      ///
      /// Subscript operator. It gives back always the item 
      /// what was last assigned to the value.
      Value operator[](const Key& key) const {
	return inverted(key);
      }
      
    private:
      const InvertableMap& inverted;
    };

    /// \brief It gives back the just readable inverse map.
    ///
    /// It gives back the just readable inverse map.
    InverseMap inverse() const {
      return InverseMap(*this);
    } 


    
  };

  /// \brief Provides a mutable, continuous and unique descriptor for each 
  /// item in the graph.
  ///
  /// The DescriptorMap class provides a unique and continuous (but mutable)
  /// descriptor (id) for each item of the same type (e.g. node) in the
  /// graph. This id is <ul><li>\b unique: different items (nodes) get
  /// different ids <li>\b continuous: the range of the ids is the set of
  /// integers between 0 and \c n-1, where \c n is the number of the items of
  /// this type (e.g. nodes) (so the id of a node can change if you delete an
  /// other node, i.e. this id is mutable).  </ul> This map can be inverted
  /// with its member class \c InverseMap.
  ///
  /// \param _Graph The graph class the \c DescriptorMap belongs to.
  /// \param _Item The Item is the Key of the Map. It may be Node, Edge or 
  /// UEdge.
  template <typename _Graph, typename _Item>
  class DescriptorMap : protected DefaultMap<_Graph, _Item, int> {

    typedef _Item Item;
    typedef DefaultMap<_Graph, _Item, int> Map;

  public:
    /// The graph class of DescriptorMap.
    typedef _Graph Graph;

    /// The key type of DescriptorMap (Node, Edge, UEdge).
    typedef typename Map::Key Key;
    /// The value type of DescriptorMap.
    typedef typename Map::Value Value;

    /// \brief Constructor.
    ///
    /// Constructor for descriptor map.
    explicit DescriptorMap(const Graph& _graph) : Map(_graph) {
      Item it;
      const typename Map::Notifier* nf = Map::notifier(); 
      for (nf->first(it); it != INVALID; nf->next(it)) {
	Map::set(it, invMap.size());
	invMap.push_back(it);	
      }      
    }

  protected:

    /// \brief Add a new key to the map.
    ///
    /// Add a new key to the map. It is called by the
    /// \c AlterationNotifier.
    virtual void add(const Item& item) {
      Map::add(item);
      Map::set(item, invMap.size());
      invMap.push_back(item);
    }

    /// \brief Add more new keys to the map.
    ///
    /// Add more new keys to the map. It is called by the
    /// \c AlterationNotifier.
    virtual void add(const std::vector<Item>& items) {
      Map::add(items);
      for (int i = 0; i < int(items.size()); ++i) {
	Map::set(items[i], invMap.size());
	invMap.push_back(items[i]);
      }
    }

    /// \brief Erase the key from the map.
    ///
    /// Erase the key from the map. It is called by the
    /// \c AlterationNotifier.
    virtual void erase(const Item& item) {
      Map::set(invMap.back(), Map::operator[](item));
      invMap[Map::operator[](item)] = invMap.back();
      invMap.pop_back();
      Map::erase(item);
    }

    /// \brief Erase more keys from the map.
    ///
    /// Erase more keys from the map. It is called by the
    /// \c AlterationNotifier.
    virtual void erase(const std::vector<Item>& items) {
      for (int i = 0; i < int(items.size()); ++i) {
	Map::set(invMap.back(), Map::operator[](items[i]));
	invMap[Map::operator[](items[i])] = invMap.back();
	invMap.pop_back();
      }
      Map::erase(items);
    }

    /// \brief Build the unique map.
    ///
    /// Build the unique map. It is called by the
    /// \c AlterationNotifier.
    virtual void build() {
      Map::build();
      Item it;
      const typename Map::Notifier* nf = Map::notifier(); 
      for (nf->first(it); it != INVALID; nf->next(it)) {
	Map::set(it, invMap.size());
	invMap.push_back(it);	
      }      
    }
    
    /// \brief Clear the keys from the map.
    ///
    /// Clear the keys from the map. It is called by the
    /// \c AlterationNotifier.
    virtual void clear() {
      invMap.clear();
      Map::clear();
    }

  public:

    /// \brief Returns the maximal value plus one.
    ///
    /// Returns the maximal value plus one in the map.
    unsigned int size() const {
      return invMap.size();
    }

    /// \brief Swaps the position of the two items in the map.
    ///
    /// Swaps the position of the two items in the map.
    void swap(const Item& p, const Item& q) {
      int pi = Map::operator[](p);
      int qi = Map::operator[](q);
      Map::set(p, qi);
      invMap[qi] = p;
      Map::set(q, pi);
      invMap[pi] = q;
    }

    /// \brief Gives back the \e descriptor of the item.
    ///
    /// Gives back the mutable and unique \e descriptor of the map.
    int operator[](const Item& item) const {
      return Map::operator[](item);
    }

    /// \brief Gives back the item by its descriptor.
    ///
    /// Gives back th item by its descriptor.
    Item operator()(int id) const {
      return invMap[id];
    }
    
  private:

    typedef std::vector<Item> Container;
    Container invMap;

  public:
    /// \brief The inverse map type of DescriptorMap.
    ///
    /// The inverse map type of DescriptorMap.
    class InverseMap {
    public:
      /// \brief Constructor of the InverseMap.
      ///
      /// Constructor of the InverseMap.
      explicit InverseMap(const DescriptorMap& _inverted) 
	: inverted(_inverted) {}


      /// The value type of the InverseMap.
      typedef typename DescriptorMap::Key Value;
      /// The key type of the InverseMap.
      typedef typename DescriptorMap::Value Key; 

      /// \brief Subscript operator. 
      ///
      /// Subscript operator. It gives back the item 
      /// that the descriptor belongs to currently.
      Value operator[](const Key& key) const {
	return inverted(key);
      }

      /// \brief Size of the map.
      ///
      /// Returns the size of the map.
      unsigned int size() const {
	return inverted.size();
      }
      
    private:
      const DescriptorMap& inverted;
    };

    /// \brief Gives back the inverse of the map.
    ///
    /// Gives back the inverse of the map.
    const InverseMap inverse() const {
      return InverseMap(*this);
    }
  };

  /// \brief Returns the source of the given edge.
  ///
  /// The SourceMap gives back the source Node of the given edge. 
  /// \see TargetMap
  /// \author Balazs Dezso
  template <typename Graph>
  class SourceMap {
  public:

    typedef typename Graph::Node Value;
    typedef typename Graph::Edge Key;

    /// \brief Constructor
    ///
    /// Constructor
    /// \param _graph The graph that the map belongs to.
    explicit SourceMap(const Graph& _graph) : graph(_graph) {}

    /// \brief The subscript operator.
    ///
    /// The subscript operator.
    /// \param edge The edge 
    /// \return The source of the edge 
    Value operator[](const Key& edge) const {
      return graph.source(edge);
    }

  private:
    const Graph& graph;
  };

  /// \brief Returns a \ref SourceMap class.
  ///
  /// This function just returns an \ref SourceMap class.
  /// \relates SourceMap
  template <typename Graph>
  inline SourceMap<Graph> sourceMap(const Graph& graph) {
    return SourceMap<Graph>(graph);
  } 

  /// \brief Returns the target of the given edge.
  ///
  /// The TargetMap gives back the target Node of the given edge. 
  /// \see SourceMap
  /// \author Balazs Dezso
  template <typename Graph>
  class TargetMap {
  public:

    typedef typename Graph::Node Value;
    typedef typename Graph::Edge Key;

    /// \brief Constructor
    ///
    /// Constructor
    /// \param _graph The graph that the map belongs to.
    explicit TargetMap(const Graph& _graph) : graph(_graph) {}

    /// \brief The subscript operator.
    ///
    /// The subscript operator.
    /// \param e The edge 
    /// \return The target of the edge 
    Value operator[](const Key& e) const {
      return graph.target(e);
    }

  private:
    const Graph& graph;
  };

  /// \brief Returns a \ref TargetMap class.
  ///
  /// This function just returns a \ref TargetMap class.
  /// \relates TargetMap
  template <typename Graph>
  inline TargetMap<Graph> targetMap(const Graph& graph) {
    return TargetMap<Graph>(graph);
  }

  /// \brief Returns the "forward" directed edge view of an undirected edge.
  ///
  /// Returns the "forward" directed edge view of an undirected edge.
  /// \see BackwardMap
  /// \author Balazs Dezso
  template <typename Graph>
  class ForwardMap {
  public:

    typedef typename Graph::Edge Value;
    typedef typename Graph::UEdge Key;

    /// \brief Constructor
    ///
    /// Constructor
    /// \param _graph The graph that the map belongs to.
    explicit ForwardMap(const Graph& _graph) : graph(_graph) {}

    /// \brief The subscript operator.
    ///
    /// The subscript operator.
    /// \param key An undirected edge 
    /// \return The "forward" directed edge view of undirected edge 
    Value operator[](const Key& key) const {
      return graph.direct(key, true);
    }

  private:
    const Graph& graph;
  };

  /// \brief Returns a \ref ForwardMap class.
  ///
  /// This function just returns an \ref ForwardMap class.
  /// \relates ForwardMap
  template <typename Graph>
  inline ForwardMap<Graph> forwardMap(const Graph& graph) {
    return ForwardMap<Graph>(graph);
  }

  /// \brief Returns the "backward" directed edge view of an undirected edge.
  ///
  /// Returns the "backward" directed edge view of an undirected edge.
  /// \see ForwardMap
  /// \author Balazs Dezso
  template <typename Graph>
  class BackwardMap {
  public:

    typedef typename Graph::Edge Value;
    typedef typename Graph::UEdge Key;

    /// \brief Constructor
    ///
    /// Constructor
    /// \param _graph The graph that the map belongs to.
    explicit BackwardMap(const Graph& _graph) : graph(_graph) {}

    /// \brief The subscript operator.
    ///
    /// The subscript operator.
    /// \param key An undirected edge 
    /// \return The "backward" directed edge view of undirected edge 
    Value operator[](const Key& key) const {
      return graph.direct(key, false);
    }

  private:
    const Graph& graph;
  };

  /// \brief Returns a \ref BackwardMap class

  /// This function just returns a \ref BackwardMap class.
  /// \relates BackwardMap
  template <typename Graph>
  inline BackwardMap<Graph> backwardMap(const Graph& graph) {
    return BackwardMap<Graph>(graph);
  }

  /// \brief Potential difference map
  ///
  /// If there is an potential map on the nodes then we
  /// can get an edge map as we get the substraction of the
  /// values of the target and source.
  template <typename Graph, typename NodeMap>
  class PotentialDifferenceMap {
  public:
    typedef typename Graph::Edge Key;
    typedef typename NodeMap::Value Value;

    /// \brief Constructor
    ///
    /// Contructor of the map
    explicit PotentialDifferenceMap(const Graph& _graph, 
                                    const NodeMap& _potential) 
      : graph(_graph), potential(_potential) {}

    /// \brief Const subscription operator
    ///
    /// Const subscription operator
    Value operator[](const Key& edge) const {
      return potential[graph.target(edge)] - potential[graph.source(edge)];
    }

  private:
    const Graph& graph;
    const NodeMap& potential;
  };

  /// \brief Returns a PotentialDifferenceMap.
  ///
  /// This function just returns a PotentialDifferenceMap.
  /// \relates PotentialDifferenceMap
  template <typename Graph, typename NodeMap>
  PotentialDifferenceMap<Graph, NodeMap> 
  potentialDifferenceMap(const Graph& graph, const NodeMap& potential) {
    return PotentialDifferenceMap<Graph, NodeMap>(graph, potential);
  }

  /// \brief Map of the node in-degrees.
  ///
  /// This map returns the in-degree of a node. Once it is constructed,
  /// the degrees are stored in a standard NodeMap, so each query is done
  /// in constant time. On the other hand, the values are updated automatically
  /// whenever the graph changes.
  ///
  /// \warning Besides addNode() and addEdge(), a graph structure may provide
  /// alternative ways to modify the graph. The correct behavior of InDegMap
  /// is not guarantied if these additional features are used. For example
  /// the functions \ref ListGraph::changeSource() "changeSource()",
  /// \ref ListGraph::changeTarget() "changeTarget()" and
  /// \ref ListGraph::reverseEdge() "reverseEdge()"
  /// of \ref ListGraph will \e not update the degree values correctly.
  ///
  /// \sa OutDegMap

  template <typename _Graph>
  class InDegMap  
    : protected ItemSetTraits<_Graph, typename _Graph::Edge>
      ::ItemNotifier::ObserverBase {

  public:
    
    typedef _Graph Graph;
    typedef int Value;
    typedef typename Graph::Node Key;

    typedef typename ItemSetTraits<_Graph, typename _Graph::Edge>
    ::ItemNotifier::ObserverBase Parent;

  private:

    class AutoNodeMap : public DefaultMap<_Graph, Key, int> {
    public:

      typedef DefaultMap<_Graph, Key, int> Parent;
      typedef typename Parent::Graph Graph;

      AutoNodeMap(const Graph& graph) : Parent(graph, 0) {}
      
      virtual void add(const Key& key) {
	Parent::add(key);
	Parent::set(key, 0);
      }

      virtual void add(const std::vector<Key>& keys) {
	Parent::add(keys);
	for (int i = 0; i < int(keys.size()); ++i) {
	  Parent::set(keys[i], 0);
	}
      }

      virtual void build() {
	Parent::build();
	Key it;
	typename Parent::Notifier* nf = Parent::notifier();
	for (nf->first(it); it != INVALID; nf->next(it)) {
	  Parent::set(it, 0);
	}
      }
    };

  public:

    /// \brief Constructor.
    ///
    /// Constructor for creating in-degree map.
    explicit InDegMap(const Graph& _graph) : graph(_graph), deg(_graph) {
      Parent::attach(graph.notifier(typename _Graph::Edge()));
      
      for(typename _Graph::NodeIt it(graph); it != INVALID; ++it) {
	deg[it] = countInEdges(graph, it);
      }
    }
    
    /// Gives back the in-degree of a Node.
    int operator[](const Key& key) const {
      return deg[key];
    }

  protected:
    
    typedef typename Graph::Edge Edge;

    virtual void add(const Edge& edge) {
      ++deg[graph.target(edge)];
    }

    virtual void add(const std::vector<Edge>& edges) {
      for (int i = 0; i < int(edges.size()); ++i) {
        ++deg[graph.target(edges[i])];
      }
    }

    virtual void erase(const Edge& edge) {
      --deg[graph.target(edge)];
    }

    virtual void erase(const std::vector<Edge>& edges) {
      for (int i = 0; i < int(edges.size()); ++i) {
        --deg[graph.target(edges[i])];
      }
    }

    virtual void build() {
      for(typename _Graph::NodeIt it(graph); it != INVALID; ++it) {
	deg[it] = countInEdges(graph, it);
      }      
    }

    virtual void clear() {
      for(typename _Graph::NodeIt it(graph); it != INVALID; ++it) {
	deg[it] = 0;
      }
    }
  private:
    
    const _Graph& graph;
    AutoNodeMap deg;
  };

  /// \brief Map of the node out-degrees.
  ///
  /// This map returns the out-degree of a node. Once it is constructed,
  /// the degrees are stored in a standard NodeMap, so each query is done
  /// in constant time. On the other hand, the values are updated automatically
  /// whenever the graph changes.
  ///
  /// \warning Besides addNode() and addEdge(), a graph structure may provide
  /// alternative ways to modify the graph. The correct behavior of OutDegMap
  /// is not guarantied if these additional features are used. For example
  /// the functions \ref ListGraph::changeSource() "changeSource()",
  /// \ref ListGraph::changeTarget() "changeTarget()" and
  /// \ref ListGraph::reverseEdge() "reverseEdge()"
  /// of \ref ListGraph will \e not update the degree values correctly.
  ///
  /// \sa InDegMap

  template <typename _Graph>
  class OutDegMap  
    : protected ItemSetTraits<_Graph, typename _Graph::Edge>
      ::ItemNotifier::ObserverBase {

  public:

    typedef typename ItemSetTraits<_Graph, typename _Graph::Edge>
    ::ItemNotifier::ObserverBase Parent;
    
    typedef _Graph Graph;
    typedef int Value;
    typedef typename Graph::Node Key;

  private:

    class AutoNodeMap : public DefaultMap<_Graph, Key, int> {
    public:

      typedef DefaultMap<_Graph, Key, int> Parent;
      typedef typename Parent::Graph Graph;

      AutoNodeMap(const Graph& graph) : Parent(graph, 0) {}
      
      virtual void add(const Key& key) {
	Parent::add(key);
	Parent::set(key, 0);
      }
      virtual void add(const std::vector<Key>& keys) {
	Parent::add(keys);
	for (int i = 0; i < int(keys.size()); ++i) {
	  Parent::set(keys[i], 0);
	}
      }
      virtual void build() {
	Parent::build();
	Key it;
	typename Parent::Notifier* nf = Parent::notifier();
	for (nf->first(it); it != INVALID; nf->next(it)) {
	  Parent::set(it, 0);
	}
      }
    };

  public:

    /// \brief Constructor.
    ///
    /// Constructor for creating out-degree map.
    explicit OutDegMap(const Graph& _graph) : graph(_graph), deg(_graph) {
      Parent::attach(graph.notifier(typename _Graph::Edge()));
      
      for(typename _Graph::NodeIt it(graph); it != INVALID; ++it) {
	deg[it] = countOutEdges(graph, it);
      }
    }

    /// Gives back the out-degree of a Node.
    int operator[](const Key& key) const {
      return deg[key];
    }

  protected:
    
    typedef typename Graph::Edge Edge;

    virtual void add(const Edge& edge) {
      ++deg[graph.source(edge)];
    }

    virtual void add(const std::vector<Edge>& edges) {
      for (int i = 0; i < int(edges.size()); ++i) {
        ++deg[graph.source(edges[i])];
      }
    }

    virtual void erase(const Edge& edge) {
      --deg[graph.source(edge)];
    }

    virtual void erase(const std::vector<Edge>& edges) {
      for (int i = 0; i < int(edges.size()); ++i) {
        --deg[graph.source(edges[i])];
      }
    }

    virtual void build() {
      for(typename _Graph::NodeIt it(graph); it != INVALID; ++it) {
	deg[it] = countOutEdges(graph, it);
      }      
    }

    virtual void clear() {
      for(typename _Graph::NodeIt it(graph); it != INVALID; ++it) {
	deg[it] = 0;
      }
    }
  private:
    
    const _Graph& graph;
    AutoNodeMap deg;
  };


  ///Dynamic edge look up between given endpoints.
  
  ///\ingroup gutils
  ///Using this class, you can find an edge in a graph from a given
  ///source to a given target in amortized time <em>O(log d)</em>,
  ///where <em>d</em> is the out-degree of the source node.
  ///
  ///It is possible to find \e all parallel edges between two nodes with
  ///the \c findFirst() and \c findNext() members.
  ///
  ///See the \ref EdgeLookUp and \ref AllEdgeLookUp classes if your
  ///graph do not changed so frequently.
  ///
  ///This class uses a self-adjusting binary search tree, Sleator's
  ///and Tarjan's Splay tree for guarantee the logarithmic amortized
  ///time bound for edge lookups. This class also guarantees the
  ///optimal time bound in a constant factor for any distribution of
  ///queries.
  ///
  ///\param G The type of the underlying graph.  
  ///
  ///\sa EdgeLookUp  
  ///\sa AllEdgeLookUp  
  template<class G>
  class DynEdgeLookUp 
    : protected ItemSetTraits<G, typename G::Edge>::ItemNotifier::ObserverBase
  {
  public:
    typedef typename ItemSetTraits<G, typename G::Edge>
    ::ItemNotifier::ObserverBase Parent;

    GRAPH_TYPEDEFS(typename G);
    typedef G Graph;

  protected:

    class AutoNodeMap : public DefaultMap<G, Node, Edge> {
    public:

      typedef DefaultMap<G, Node, Edge> Parent;

      AutoNodeMap(const G& graph) : Parent(graph, INVALID) {}
      
      virtual void add(const Node& node) {
	Parent::add(node);
	Parent::set(node, INVALID);
      }

      virtual void add(const std::vector<Node>& nodes) {
	Parent::add(nodes);
	for (int i = 0; i < int(nodes.size()); ++i) {
	  Parent::set(nodes[i], INVALID);
	}
      }

      virtual void build() {
	Parent::build();
	Node it;
	typename Parent::Notifier* nf = Parent::notifier();
	for (nf->first(it); it != INVALID; nf->next(it)) {
	  Parent::set(it, INVALID);
	}
      }
    };

    const Graph &_g;
    AutoNodeMap _head;
    typename Graph::template EdgeMap<Edge> _parent;
    typename Graph::template EdgeMap<Edge> _left;
    typename Graph::template EdgeMap<Edge> _right;
    
    class EdgeLess {
      const Graph &g;
    public:
      EdgeLess(const Graph &_g) : g(_g) {}
      bool operator()(Edge a,Edge b) const 
      {
	return g.target(a)<g.target(b);
      }
    };
    
  public:
    
    ///Constructor

    ///Constructor.
    ///
    ///It builds up the search database.
    DynEdgeLookUp(const Graph &g) 
      : _g(g),_head(g),_parent(g),_left(g),_right(g) 
    { 
      Parent::attach(_g.notifier(typename Graph::Edge()));
      refresh(); 
    }
    
  protected:

    virtual void add(const Edge& edge) {
      insert(edge);
    }

    virtual void add(const std::vector<Edge>& edges) {
      for (int i = 0; i < int(edges.size()); ++i) {
	insert(edges[i]);
      }
    }

    virtual void erase(const Edge& edge) {
      remove(edge);
    }

    virtual void erase(const std::vector<Edge>& edges) {
      for (int i = 0; i < int(edges.size()); ++i) {
	remove(edges[i]);
      }     
    }

    virtual void build() {
      refresh();
    }

    virtual void clear() {
      for(NodeIt n(_g);n!=INVALID;++n) {
	_head.set(n, INVALID);
      }
    }

    void insert(Edge edge) {
      Node s = _g.source(edge);
      Node t = _g.target(edge);
      _left.set(edge, INVALID);
      _right.set(edge, INVALID);
      
      Edge e = _head[s];
      if (e == INVALID) {
	_head.set(s, edge);
	_parent.set(edge, INVALID);
	return;
      }
      while (true) {
	if (t < _g.target(e)) {
	  if (_left[e] == INVALID) {
	    _left.set(e, edge);
	    _parent.set(edge, e);
	    splay(edge);
	    return;
	  } else {
	    e = _left[e];
	  }
	} else {
	  if (_right[e] == INVALID) {
	    _right.set(e, edge);
	    _parent.set(edge, e);
	    splay(edge);
	    return;
	  } else {
	    e = _right[e];
	  }
	}
      }
    }

    void remove(Edge edge) {
      if (_left[edge] == INVALID) {
	if (_right[edge] != INVALID) {
	  _parent.set(_right[edge], _parent[edge]);
	}
	if (_parent[edge] != INVALID) {
	  if (_left[_parent[edge]] == edge) {
	    _left.set(_parent[edge], _right[edge]);
	  } else {
	    _right.set(_parent[edge], _right[edge]);
	  }
	} else {
	  _head.set(_g.source(edge), _right[edge]);
	}
      } else if (_right[edge] == INVALID) {
	_parent.set(_left[edge], _parent[edge]);
	if (_parent[edge] != INVALID) {
	  if (_left[_parent[edge]] == edge) {
	    _left.set(_parent[edge], _left[edge]);
	  } else {
	    _right.set(_parent[edge], _left[edge]);
	  }
	} else {
	  _head.set(_g.source(edge), _left[edge]);
	}
      } else {
	Edge e = _left[edge];
	if (_right[e] != INVALID) {
	  e = _right[e];	  
	  while (_right[e] != INVALID) {
	    e = _right[e];
	  }
	  Edge s = _parent[e];
	  _right.set(_parent[e], _left[e]);
	  if (_left[e] != INVALID) {
	    _parent.set(_left[e], _parent[e]);
	  }
	  
	  _left.set(e, _left[edge]);
	  _parent.set(_left[edge], e);
	  _right.set(e, _right[edge]);
	  _parent.set(_right[edge], e);

	  _parent.set(e, _parent[edge]);
	  if (_parent[edge] != INVALID) {
	    if (_left[_parent[edge]] == edge) {
	      _left.set(_parent[edge], e);
	    } else {
	      _right.set(_parent[edge], e);
	    }
	  }
	  splay(s);
	} else {
	  _right.set(e, _right[edge]);
	  _parent.set(_right[edge], e);

	  if (_parent[edge] != INVALID) {
	    if (_left[_parent[edge]] == edge) {
	      _left.set(_parent[edge], e);
	    } else {
	      _right.set(_parent[edge], e);
	    }
	  } else {
	    _head.set(_g.source(edge), e);
	  }
	}
      }
    }

    Edge refreshRec(std::vector<Edge> &v,int a,int b) 
    {
      int m=(a+b)/2;
      Edge me=v[m];
      if (a < m) {
	Edge left = refreshRec(v,a,m-1);
	_left.set(me, left);
	_parent.set(left, me);
      } else {
	_left.set(me, INVALID);
      }
      if (m < b) {
	Edge right = refreshRec(v,m+1,b);
	_right.set(me, right);
	_parent.set(right, me);
      } else {
	_right.set(me, INVALID);
      }
      return me;
    }

    void refresh() {
      for(NodeIt n(_g);n!=INVALID;++n) {
	std::vector<Edge> v;
	for(OutEdgeIt e(_g,n);e!=INVALID;++e) v.push_back(e);
	if(v.size()) {
	  std::sort(v.begin(),v.end(),EdgeLess(_g));
	  Edge head = refreshRec(v,0,v.size()-1);
	  _head.set(n, head);
	  _parent.set(head, INVALID);
	}
	else _head.set(n, INVALID);
      }
    }

    void zig(Edge v) {        
      Edge w = _parent[v];
      _parent.set(v, _parent[w]);
      _parent.set(w, v);
      _left.set(w, _right[v]);
      _right.set(v, w);
      if (_parent[v] != INVALID) {
	if (_right[_parent[v]] == w) {
	  _right.set(_parent[v], v);
	} else {
	  _left.set(_parent[v], v);
	}
      }
      if (_left[w] != INVALID){
	_parent.set(_left[w], w);
      }
    }

    void zag(Edge v) {        
      Edge w = _parent[v];
      _parent.set(v, _parent[w]);
      _parent.set(w, v);
      _right.set(w, _left[v]);
      _left.set(v, w);
      if (_parent[v] != INVALID){
	if (_left[_parent[v]] == w) {
	  _left.set(_parent[v], v);
	} else {
	  _right.set(_parent[v], v);
	}
      }
      if (_right[w] != INVALID){
	_parent.set(_right[w], w);
      }
    }

    void splay(Edge v) {
      while (_parent[v] != INVALID) {
	if (v == _left[_parent[v]]) {
	  if (_parent[_parent[v]] == INVALID) {
	    zig(v);
	  } else {
	    if (_parent[v] == _left[_parent[_parent[v]]]) {
	      zig(_parent[v]);
	      zig(v);
	    } else {
	      zig(v);
	      zag(v);
	    }
	  }
	} else {
	  if (_parent[_parent[v]] == INVALID) {
	    zag(v);
	  } else {
	    if (_parent[v] == _left[_parent[_parent[v]]]) {
	      zag(v);
	      zig(v);
	    } else {
	      zag(_parent[v]);
	      zag(v);
	    }
	  }
	}
      }
      _head[_g.source(v)] = v;
    }


  public:
    
    ///Find an edge between two nodes.
    
    ///Find an edge between two nodes in time <em>O(</em>log<em>d)</em>, where
    /// <em>d</em> is the number of outgoing edges of \c s.
    ///\param s The source node
    ///\param t The target node
    ///\return An edge from \c s to \c t if there exists,
    ///\ref INVALID otherwise.
    Edge operator()(Node s, Node t) const
    {
      Edge e = _head[s];
      while (true) {
	if (_g.target(e) == t) {
	  const_cast<DynEdgeLookUp&>(*this).splay(e);
	  return e;
	} else if (t < _g.target(e)) {
	  if (_left[e] == INVALID) {
	    const_cast<DynEdgeLookUp&>(*this).splay(e);
	    return INVALID;
	  } else {
	    e = _left[e];
	  }
	} else  {
	  if (_right[e] == INVALID) {
	    const_cast<DynEdgeLookUp&>(*this).splay(e);
	    return INVALID;
	  } else {
	    e = _right[e];
	  }
	}
      }
    }

    ///Find the first edge between two nodes.
    
    ///Find the first edge between two nodes in time
    /// <em>O(</em>log<em>d)</em>, where <em>d</em> is the number of
    /// outgoing edges of \c s.  
    ///\param s The source node 
    ///\param t The target node
    ///\return An edge from \c s to \c t if there exists, \ref INVALID
    /// otherwise.
    Edge findFirst(Node s, Node t) const
    {
      Edge e = _head[s];
      Edge r = INVALID;
      while (true) {
	if (_g.target(e) < t) {
	  if (_right[e] == INVALID) {
	    const_cast<DynEdgeLookUp&>(*this).splay(e);
	    return r;
	  } else {
	    e = _right[e];
	  }
	} else {
	  if (_g.target(e) == t) {
	    r = e;
	  }
	  if (_left[e] == INVALID) {
	    const_cast<DynEdgeLookUp&>(*this).splay(e);
	    return r;
	  } else {
	    e = _left[e];
	  }
	}
      }
    }

    ///Find the next edge between two nodes.
    
    ///Find the next edge between two nodes in time
    /// <em>O(</em>log<em>d)</em>, where <em>d</em> is the number of
    /// outgoing edges of \c s.  
    ///\param s The source node 
    ///\param t The target node
    ///\return An edge from \c s to \c t if there exists, \ref INVALID
    /// otherwise.

    ///\note If \c e is not the result of the previous \c findFirst()
    ///operation then the amorized time bound can not be guaranteed.
#ifdef DOXYGEN
    Edge findNext(Node s, Node t, Edge e) const
#else
    Edge findNext(Node, Node t, Edge e) const
#endif
    {
      if (_right[e] != INVALID) {
	e = _right[e];
	while (_left[e] != INVALID) {
	  e = _left[e];
	}
	const_cast<DynEdgeLookUp&>(*this).splay(e);
      } else {
	while (_parent[e] != INVALID && _right[_parent[e]] ==  e) {
	  e = _parent[e];
	}
	if (_parent[e] == INVALID) {
	  return INVALID;
	} else {
	  e = _parent[e];
	  const_cast<DynEdgeLookUp&>(*this).splay(e);
	}
      }
      if (_g.target(e) == t) return e;
      else return INVALID;    
    }

  };

  ///Fast edge look up between given endpoints.
  
  ///\ingroup gutils
  ///Using this class, you can find an edge in a graph from a given
  ///source to a given target in time <em>O(log d)</em>,
  ///where <em>d</em> is the out-degree of the source node.
  ///
  ///It is not possible to find \e all parallel edges between two nodes.
  ///Use \ref AllEdgeLookUp for this purpose.
  ///
  ///\warning This class is static, so you should refresh() (or at least
  ///refresh(Node)) this data structure
  ///whenever the graph changes. This is a time consuming (superlinearly
  ///proportional (<em>O(m</em>log<em>m)</em>) to the number of edges).
  ///
  ///\param G The type of the underlying graph.
  ///
  ///\sa DynEdgeLookUp
  ///\sa AllEdgeLookUp  
  template<class G>
  class EdgeLookUp 
  {
  public:
    GRAPH_TYPEDEFS(typename G);
    typedef G Graph;

  protected:
    const Graph &_g;
    typename Graph::template NodeMap<Edge> _head;
    typename Graph::template EdgeMap<Edge> _left;
    typename Graph::template EdgeMap<Edge> _right;
    
    class EdgeLess {
      const Graph &g;
    public:
      EdgeLess(const Graph &_g) : g(_g) {}
      bool operator()(Edge a,Edge b) const 
      {
	return g.target(a)<g.target(b);
      }
    };
    
  public:
    
    ///Constructor

    ///Constructor.
    ///
    ///It builds up the search database, which remains valid until the graph
    ///changes.
    EdgeLookUp(const Graph &g) :_g(g),_head(g),_left(g),_right(g) {refresh();}
    
  private:
    Edge refreshRec(std::vector<Edge> &v,int a,int b) 
    {
      int m=(a+b)/2;
      Edge me=v[m];
      _left[me] = a<m?refreshRec(v,a,m-1):INVALID;
      _right[me] = m<b?refreshRec(v,m+1,b):INVALID;
      return me;
    }
  public:
    ///Refresh the data structure at a node.

    ///Build up the search database of node \c n.
    ///
    ///It runs in time <em>O(d</em>log<em>d)</em>, where <em>d</em> is
    ///the number of the outgoing edges of \c n.
    void refresh(Node n) 
    {
      std::vector<Edge> v;
      for(OutEdgeIt e(_g,n);e!=INVALID;++e) v.push_back(e);
      if(v.size()) {
	std::sort(v.begin(),v.end(),EdgeLess(_g));
	_head[n]=refreshRec(v,0,v.size()-1);
      }
      else _head[n]=INVALID;
    }
    ///Refresh the full data structure.

    ///Build up the full search database. In fact, it simply calls
    ///\ref refresh(Node) "refresh(n)" for each node \c n.
    ///
    ///It runs in time <em>O(m</em>log<em>D)</em>, where <em>m</em> is
    ///the number of the edges of \c n and <em>D</em> is the maximum
    ///out-degree of the graph.

    void refresh() 
    {
      for(NodeIt n(_g);n!=INVALID;++n) refresh(n);
    }
    
    ///Find an edge between two nodes.
    
    ///Find an edge between two nodes in time <em>O(</em>log<em>d)</em>, where
    /// <em>d</em> is the number of outgoing edges of \c s.
    ///\param s The source node
    ///\param t The target node
    ///\return An edge from \c s to \c t if there exists,
    ///\ref INVALID otherwise.
    ///
    ///\warning If you change the graph, refresh() must be called before using
    ///this operator. If you change the outgoing edges of
    ///a single node \c n, then
    ///\ref refresh(Node) "refresh(n)" is enough.
    ///
    Edge operator()(Node s, Node t) const
    {
      Edge e;
      for(e=_head[s];
	  e!=INVALID&&_g.target(e)!=t;
	  e = t < _g.target(e)?_left[e]:_right[e]) ;
      return e;
    }

  };

  ///Fast look up of all edges between given endpoints.
  
  ///\ingroup gutils
  ///This class is the same as \ref EdgeLookUp, with the addition
  ///that it makes it possible to find all edges between given endpoints.
  ///
  ///\warning This class is static, so you should refresh() (or at least
  ///refresh(Node)) this data structure
  ///whenever the graph changes. This is a time consuming (superlinearly
  ///proportional (<em>O(m</em>log<em>m)</em>) to the number of edges).
  ///
  ///\param G The type of the underlying graph.
  ///
  ///\sa DynEdgeLookUp
  ///\sa EdgeLookUp  
  template<class G>
  class AllEdgeLookUp : public EdgeLookUp<G>
  {
    using EdgeLookUp<G>::_g;
    using EdgeLookUp<G>::_right;
    using EdgeLookUp<G>::_left;
    using EdgeLookUp<G>::_head;

    GRAPH_TYPEDEFS(typename G);
    typedef G Graph;
    
    typename Graph::template EdgeMap<Edge> _next;
    
    Edge refreshNext(Edge head,Edge next=INVALID)
    {
      if(head==INVALID) return next;
      else {
	next=refreshNext(_right[head],next);
// 	_next[head]=next;
	_next[head]=( next!=INVALID && _g.target(next)==_g.target(head))
	  ? next : INVALID;
	return refreshNext(_left[head],head);
      }
    }
    
    void refreshNext()
    {
      for(NodeIt n(_g);n!=INVALID;++n) refreshNext(_head[n]);
    }
    
  public:
    ///Constructor

    ///Constructor.
    ///
    ///It builds up the search database, which remains valid until the graph
    ///changes.
    AllEdgeLookUp(const Graph &g) : EdgeLookUp<G>(g), _next(g) {refreshNext();}

    ///Refresh the data structure at a node.

    ///Build up the search database of node \c n.
    ///
    ///It runs in time <em>O(d</em>log<em>d)</em>, where <em>d</em> is
    ///the number of the outgoing edges of \c n.
    
    void refresh(Node n) 
    {
      EdgeLookUp<G>::refresh(n);
      refreshNext(_head[n]);
    }
    
    ///Refresh the full data structure.

    ///Build up the full search database. In fact, it simply calls
    ///\ref refresh(Node) "refresh(n)" for each node \c n.
    ///
    ///It runs in time <em>O(m</em>log<em>D)</em>, where <em>m</em> is
    ///the number of the edges of \c n and <em>D</em> is the maximum
    ///out-degree of the graph.

    void refresh() 
    {
      for(NodeIt n(_g);n!=INVALID;++n) refresh(_head[n]);
    }
    
    ///Find an edge between two nodes.
    
    ///Find an edge between two nodes.
    ///\param s The source node
    ///\param t The target node
    ///\param prev The previous edge between \c s and \c t. It it is INVALID or
    ///not given, the operator finds the first appropriate edge.
    ///\return An edge from \c s to \c t after \c prev or
    ///\ref INVALID if there is no more.
    ///
    ///For example, you can count the number of edges from \c u to \c v in the
    ///following way.
    ///\code
    ///AllEdgeLookUp<ListGraph> ae(g);
    ///...
    ///int n=0;
    ///for(Edge e=ae(u,v);e!=INVALID;e=ae(u,v,e)) n++;
    ///\endcode
    ///
    ///Finding the first edge take <em>O(</em>log<em>d)</em> time, where
    /// <em>d</em> is the number of outgoing edges of \c s. Then, the
    ///consecutive edges are found in constant time.
    ///
    ///\warning If you change the graph, refresh() must be called before using
    ///this operator. If you change the outgoing edges of
    ///a single node \c n, then
    ///\ref refresh(Node) "refresh(n)" is enough.
    ///
#ifdef DOXYGEN
    Edge operator()(Node s, Node t, Edge prev=INVALID) const {}
#else
    using EdgeLookUp<G>::operator() ;
    Edge operator()(Node s, Node t, Edge prev) const
    {
      return prev==INVALID?(*this)(s,t):_next[prev];
    }
#endif
      
  };

  /// @}

} //END OF NAMESPACE LEMON

#endif

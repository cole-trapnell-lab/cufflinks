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

///\ingroup graph_concepts
///\file
///\brief The concept of graph components.


#ifndef LEMON_CONCEPT_GRAPH_COMPONENTS_H
#define LEMON_CONCEPT_GRAPH_COMPONENTS_H

#include <lemon/bits/invalid.h>
#include <lemon/concepts/maps.h>

#include <lemon/bits/alteration_notifier.h>

namespace lemon {
  namespace concepts {

    /// \brief Skeleton class for graph Node and Edge types
    ///
    /// This class describes the interface of Node and Edge (and UEdge
    /// in undirected graphs) subtypes of graph types.
    ///
    /// \note This class is a template class so that we can use it to
    /// create graph skeleton classes. The reason for this is than Node
    /// and Edge types should \em not derive from the same base class.
    /// For Node you should instantiate it with character 'n' and for Edge
    /// with 'e'.

#ifndef DOXYGEN
    template <char _selector = '0'>
#endif
    class GraphItem {
    public:
      /// \brief Default constructor.
      ///      
      /// \warning The default constructor is not required to set
      /// the item to some well-defined value. So you should consider it
      /// as uninitialized.
      GraphItem() {}
      /// \brief Copy constructor.
      ///
      /// Copy constructor.
      ///
      GraphItem(const GraphItem &) {}
      /// \brief Invalid constructor \& conversion.
      ///
      /// This constructor initializes the item to be invalid.
      /// \sa Invalid for more details.
      GraphItem(Invalid) {}
      /// \brief Assign operator for nodes.
      ///
      /// The nodes are assignable. 
      ///
      GraphItem& operator=(GraphItem const&) { return *this; }
      /// \brief Equality operator.
      ///
      /// Two iterators are equal if and only if they represents the
      /// same node in the graph or both are invalid.
      bool operator==(GraphItem) const { return false; }
      /// \brief Inequality operator.
      ///
      /// \sa operator==(const Node& n)
      ///
      bool operator!=(GraphItem) const { return false; }

      /// \brief Artificial ordering operator.
      ///
      /// To allow the use of graph descriptors as key type in std::map or
      /// similar associative container we require this.
      ///
      /// \note This operator only have to define some strict ordering of
      /// the items; this order has nothing to do with the iteration
      /// ordering of the items.
      bool operator<(GraphItem) const { return false; }

      template<typename _GraphItem>
      struct Constraints {
	void constraints() {
	  _GraphItem i1;
	  _GraphItem i2 = i1;
	  _GraphItem i3 = INVALID;
	  
	  i1 = i2 = i3;

	  bool b;
	  //	  b = (ia == ib) && (ia != ib) && (ia < ib);
	  b = (ia == ib) && (ia != ib);
	  b = (ia == INVALID) && (ib != INVALID);
          b = (ia < ib);
	}

	const _GraphItem &ia;
	const _GraphItem &ib;
      };
    };

    /// \brief An empty base graph class.
    ///  
    /// This class provides the minimal set of features needed for a graph
    /// structure. All graph concepts have to be conform to this base
    /// graph. It just provides types for nodes and edges and functions to
    /// get the source and the target of the edges.
    class BaseGraphComponent {
    public:

      typedef BaseGraphComponent Graph;
      
      /// \brief Node class of the graph.
      ///
      /// This class represents the Nodes of the graph. 
      ///
      typedef GraphItem<'n'> Node;

      /// \brief Edge class of the graph.
      ///
      /// This class represents the Edges of the graph. 
      ///
      typedef GraphItem<'e'> Edge;

      /// \brief Gives back the target node of an edge.
      ///
      /// Gives back the target node of an edge.
      ///
      Node target(const Edge&) const { return INVALID;}

      /// \brief Gives back the source node of an edge.
      ///
      /// Gives back the source node of an edge.
      ///
      Node source(const Edge&) const { return INVALID;}

      /// \brief Gives back the opposite node on the given edge.
      ///
      /// Gives back the opposite node on the given edge.
      Node oppositeNode(const Node&, const Edge&) const {
        return INVALID;
      }

      template <typename _Graph>
      struct Constraints {
	typedef typename _Graph::Node Node;
	typedef typename _Graph::Edge Edge;
      
	void constraints() {
	  checkConcept<GraphItem<'n'>, Node>();
	  checkConcept<GraphItem<'e'>, Edge>();
	  {
	    Node n;
	    Edge e(INVALID);
	    n = graph.source(e);
	    n = graph.target(e);
            n = graph.oppositeNode(n, e);
	  }      
	}
      
	const _Graph& graph;
      };
    };

    /// \brief An empty base undirected graph class.
    ///  
    /// This class provides the minimal set of features needed for an
    /// undirected graph structure. All undirected graph concepts have
    /// to be conform to this base graph. It just provides types for
    /// nodes, edges and undirected edges and functions to get the
    /// source and the target of the edges and undirected edges,
    /// conversion from edges to undirected edges and function to get
    /// both direction of the undirected edges.
    class BaseUGraphComponent : public BaseGraphComponent {
    public:
      typedef BaseGraphComponent::Node Node;
      typedef BaseGraphComponent::Edge Edge;
      /// \brief Undirected edge class of the graph.
      ///
      /// This class represents the undirected edges of the graph.
      /// The undirected graphs can be used as a directed graph which
      /// for each edge contains the opposite edge too so the graph is
      /// bidirected. The undirected edge represents two opposite
      /// directed edges.
      class UEdge : public GraphItem<'u'> {
      public:
        typedef GraphItem<'u'> Parent;
        /// \brief Default constructor.
        ///      
        /// \warning The default constructor is not required to set
        /// the item to some well-defined value. So you should consider it
        /// as uninitialized.
        UEdge() {}
        /// \brief Copy constructor.
        ///
        /// Copy constructor.
        ///
        UEdge(const UEdge &) : Parent() {}
        /// \brief Invalid constructor \& conversion.
        ///
        /// This constructor initializes the item to be invalid.
        /// \sa Invalid for more details.
        UEdge(Invalid) {}
        /// \brief Converter from edge to undirected edge.
        ///
        /// Besides the core graph item functionality each edge should
        /// be convertible to the represented undirected edge. 
        UEdge(const Edge&) {}
        /// \brief Assign edge to undirected edge.
        ///
        /// Besides the core graph item functionality each edge should
        /// be convertible to the represented undirected edge. 
        UEdge& operator=(const Edge&) { return *this; }
      };

      /// \brief Returns the direction of the edge.
      ///
      /// Returns the direction of the edge. Each edge represents an
      /// undirected edge with a direction. It gives back the
      /// direction.
      bool direction(const Edge&) const { return true; }

      /// \brief Returns the directed edge.
      ///
      /// Returns the directed edge from its direction and the
      /// represented undirected edge.
      Edge direct(const UEdge&, bool) const { return INVALID;} 

      /// \brief Returns the directed edge.
      ///
      /// Returns the directed edge from its source and the
      /// represented undirected edge.
      Edge direct(const UEdge&, const Node&) const { return INVALID;} 

      /// \brief Returns the opposite edge.
      ///
      /// Returns the opposite edge. It is the edge representing the
      /// same undirected edge and has opposite direction.
      Edge oppositeEdge(const Edge&) const { return INVALID;}

      /// \brief Gives back the target node of an undirected edge.
      ///
      /// Gives back the target node of an undirected edge. The name
      /// target is a little confusing because the undirected edge
      /// does not have target but it just means that one of the end 
      /// node.
      Node target(const UEdge&) const { return INVALID;}

      /// \brief Gives back the source node of an undirected edge.
      ///
      /// Gives back the source node of an undirected edge. The name
      /// source is a little confusing because the undirected edge
      /// does not have source but it just means that one of the end 
      /// node.
      Node source(const UEdge&) const { return INVALID;}
      
      template <typename _Graph>
      struct Constraints {
	typedef typename _Graph::Node Node;
	typedef typename _Graph::Edge Edge;
	typedef typename _Graph::UEdge UEdge;
      
	void constraints() {
          checkConcept<BaseGraphComponent, _Graph>();
	  checkConcept<GraphItem<'u'>, UEdge>();
	  {
	    Node n;
	    UEdge ue(INVALID);
            Edge e;
	    n = graph.source(ue);
	    n = graph.target(ue);
            e = graph.direct(ue, true);
            e = graph.direct(ue, n);
            e = graph.oppositeEdge(e);
            ue = e;
            bool d = graph.direction(e);
            ignore_unused_variable_warning(d);
	  }      
	}
      
	const _Graph& graph;
      };

    };

    /// \brief An empty base bipartite undirected graph class.
    ///  
    /// This class provides the minimal set of features needed for an
    /// bipartite undirected graph structure. All bipartite undirected
    /// graph concepts have to be conform to this base graph. It just
    /// provides types for nodes, A-nodes, B-nodes, edges and
    /// undirected edges and functions to get the source and the
    /// target of the edges and undirected edges, conversion from
    /// edges to undirected edges and function to get both direction
    /// of the undirected edges.
    class BaseBpUGraphComponent : public BaseUGraphComponent {
    public:
      typedef BaseUGraphComponent::Node Node;
      typedef BaseUGraphComponent::Edge Edge;
      typedef BaseUGraphComponent::UEdge UEdge;

      /// \brief Helper class for A-nodes.
      ///
      /// This class is just a helper class for A-nodes, it is not
      /// suggested to use it directly. It can be converted easily to
      /// node and vice versa. The usage of this class is limited
      /// to use just as template parameters for special map types. 
      class ANode : public Node {
      public:
        typedef Node Parent;

        /// \brief Default constructor.
        ///      
        /// \warning The default constructor is not required to set
        /// the item to some well-defined value. So you should consider it
        /// as uninitialized.
        ANode() {}
        /// \brief Copy constructor.
        ///
        /// Copy constructor.
        ///
        ANode(const ANode &) : Parent() {}
        /// \brief Invalid constructor \& conversion.
        ///
        /// This constructor initializes the item to be invalid.
        /// \sa Invalid for more details.
        ANode(Invalid) {}
        /// \brief Converter from node to A-node.
        ///
        /// Besides the core graph item functionality each node should
        /// be convertible to the represented A-node if it is it possible. 
        ANode(const Node&) {}
        /// \brief Assign node to A-node.
        ///
        /// Besides the core graph item functionality each node should
        /// be convertible to the represented A-node if it is it possible. 
        ANode& operator=(const Node&) { return *this; }
      };

      /// \brief Helper class for B-nodes.
      ///
      /// This class is just a helper class for B-nodes, it is not
      /// suggested to use it directly. It can be converted easily to
      /// node and vice versa. The usage of this class is limited
      /// to use just as template parameters for special map types. 
      class BNode : public Node {
      public:
        typedef Node Parent;

        /// \brief Default constructor.
        ///      
        /// \warning The default constructor is not required to set
        /// the item to some well-defined value. So you should consider it
        /// as uninitialized.
        BNode() {}
        /// \brief Copy constructor.
        ///
        /// Copy constructor.
        ///
        BNode(const BNode &) : Parent() {}
        /// \brief Invalid constructor \& conversion.
        ///
        /// This constructor initializes the item to be invalid.
        /// \sa Invalid for more details.
        BNode(Invalid) {}
        /// \brief Converter from node to B-node.
        ///
        /// Besides the core graph item functionality each node should
        /// be convertible to the represented B-node if it is it possible. 
        BNode(const Node&) {}
        /// \brief Assign node to B-node.
        ///
        /// Besides the core graph item functionality each node should
        /// be convertible to the represented B-node if it is it possible. 
        BNode& operator=(const Node&) { return *this; }
      };
      
      /// \brief Gives back %true when the node is A-node.
      ///
      /// Gives back %true when the node is A-node.
      bool aNode(const Node&) const { return false; }

      /// \brief Gives back %true when the node is B-node.
      ///
      /// Gives back %true when the node is B-node.
      bool bNode(const Node&) const { return false; }

      /// \brief Gives back the A-node of the undirected edge.
      ///
      /// Gives back the A-node of the undirected edge.
      Node aNode(const UEdge&) const { return INVALID; }

      /// \brief Gives back the B-node of the undirected edge.
      ///
      /// Gives back the B-node of the undirected edge.
      Node bNode(const UEdge&) const { return INVALID; }
      
      template <typename _Graph>
      struct Constraints {
	typedef typename _Graph::Node Node;
	typedef typename _Graph::ANode ANode;
	typedef typename _Graph::BNode BNode;
	typedef typename _Graph::Edge Edge;
	typedef typename _Graph::UEdge UEdge;
      
	void constraints() {
          checkConcept<BaseUGraphComponent, _Graph>();
	  checkConcept<GraphItem<'a'>, ANode>();
	  checkConcept<GraphItem<'b'>, BNode>();
	  {
	    Node n;
	    UEdge ue(INVALID);
            bool b;
	    n = graph.aNode(ue);
	    n = graph.bNode(ue);
            b = graph.aNode(n);
            b = graph.bNode(n);
            ANode an;
            an = n; n = an;
            BNode bn;
            bn = n; n = bn;            
            ignore_unused_variable_warning(b);
	  }      
	}
      
	const _Graph& graph;
      };

    };

    /// \brief An empty idable base graph class.
    ///  
    /// This class provides beside the core graph features
    /// core id functions for the graph structure.
    /// The most of the base graphs should be conform to this concept.
    /// The id's are unique and immutable.
    template <typename _Base = BaseGraphComponent>
    class IDableGraphComponent : public _Base {
    public:

      typedef _Base Base;
      typedef typename Base::Node Node;
      typedef typename Base::Edge Edge;

      /// \brief Gives back an unique integer id for the Node. 
      ///
      /// Gives back an unique integer id for the Node. 
      ///
      int id(const Node&) const { return -1;}

      /// \brief Gives back the node by the unique id.
      ///
      /// Gives back the node by the unique id.
      /// If the graph does not contain node with the given id
      /// then the result of the function is undetermined. 
      Node nodeFromId(int) const { return INVALID;}

      /// \brief Gives back an unique integer id for the Edge. 
      ///
      /// Gives back an unique integer id for the Edge. 
      ///
      int id(const Edge&) const { return -1;}

      /// \brief Gives back the edge by the unique id.
      ///
      /// Gives back the edge by the unique id.
      /// If the graph does not contain edge with the given id
      /// then the result of the function is undetermined. 
      Edge edgeFromId(int) const { return INVALID;}

      /// \brief Gives back an integer greater or equal to the maximum
      /// Node id.
      ///
      /// Gives back an integer greater or equal to the maximum Node
      /// id.
      int maxNodeId() const { return -1;}

      /// \brief Gives back an integer greater or equal to the maximum
      /// Edge id.
      ///
      /// Gives back an integer greater or equal to the maximum Edge
      /// id.
      int maxEdgeId() const { return -1;}

      template <typename _Graph>
      struct Constraints {

	void constraints() {
	  checkConcept<Base, _Graph >();
	  typename _Graph::Node node;
	  int nid = graph.id(node);
	  nid = graph.id(node);
	  node = graph.nodeFromId(nid);
	  typename _Graph::Edge edge;
	  int eid = graph.id(edge);
	  eid = graph.id(edge);
	  edge = graph.edgeFromId(eid);

	  nid = graph.maxNodeId();
	  ignore_unused_variable_warning(nid);
	  eid = graph.maxEdgeId();
	  ignore_unused_variable_warning(eid);
	}

	const _Graph& graph;
      };
    };

    /// \brief An empty idable base undirected graph class.
    ///  
    /// This class provides beside the core undirected graph features
    /// core id functions for the undirected graph structure.  The
    /// most of the base undirected graphs should be conform to this
    /// concept.  The id's are unique and immutable.
    template <typename _Base = BaseUGraphComponent>
    class IDableUGraphComponent : public IDableGraphComponent<_Base> {
    public:

      typedef _Base Base;
      typedef typename Base::UEdge UEdge;

      using IDableGraphComponent<_Base>::id;

      /// \brief Gives back an unique integer id for the UEdge. 
      ///
      /// Gives back an unique integer id for the UEdge. 
      ///
      int id(const UEdge&) const { return -1;}

      /// \brief Gives back the undirected edge by the unique id.
      ///
      /// Gives back the undirected edge by the unique id.  If the
      /// graph does not contain edge with the given id then the
      /// result of the function is undetermined.
      UEdge uEdgeFromId(int) const { return INVALID;}

      /// \brief Gives back an integer greater or equal to the maximum
      /// UEdge id.
      ///
      /// Gives back an integer greater or equal to the maximum UEdge
      /// id.
      int maxUEdgeId() const { return -1;}

      template <typename _Graph>
      struct Constraints {

	void constraints() {
	  checkConcept<Base, _Graph >();
	  checkConcept<IDableGraphComponent<Base>, _Graph >();
	  typename _Graph::UEdge uedge;
	  int ueid = graph.id(uedge);
	  ueid = graph.id(uedge);
	  uedge = graph.uEdgeFromId(ueid);
	  ueid = graph.maxUEdgeId();
	  ignore_unused_variable_warning(ueid);
	}

	const _Graph& graph;
      };
    };

    /// \brief An empty idable base bipartite undirected graph class.
    ///  
    /// This class provides beside the core bipartite undirected graph
    /// features core id functions for the bipartite undirected graph
    /// structure.  The most of the base undirected graphs should be
    /// conform to this concept.
    template <typename _Base = BaseBpUGraphComponent>
    class IDableBpUGraphComponent : public IDableUGraphComponent<_Base> {
    public:

      typedef _Base Base;
      typedef typename Base::Node Node;

      using IDableUGraphComponent<_Base>::id;

      /// \brief Gives back an unique integer id for the ANode. 
      ///
      /// Gives back an unique integer id for the ANode. 
      ///
      int aNodeId(const Node&) const { return -1;}

      /// \brief Gives back the undirected edge by the unique id.
      ///
      /// Gives back the undirected edge by the unique id.  If the
      /// graph does not contain edge with the given id then the
      /// result of the function is undetermined.
      Node nodeFromANodeId(int) const { return INVALID;}

      /// \brief Gives back an integer greater or equal to the maximum
      /// ANode id.
      ///
      /// Gives back an integer greater or equal to the maximum ANode
      /// id.
      int maxANodeId() const { return -1;}

      /// \brief Gives back an unique integer id for the BNode. 
      ///
      /// Gives back an unique integer id for the BNode. 
      ///
      int bNodeId(const Node&) const { return -1;}

      /// \brief Gives back the undirected edge by the unique id.
      ///
      /// Gives back the undirected edge by the unique id.  If the
      /// graph does not contain edge with the given id then the
      /// result of the function is undetermined.
      Node nodeFromBNodeId(int) const { return INVALID;}

      /// \brief Gives back an integer greater or equal to the maximum
      /// BNode id.
      ///
      /// Gives back an integer greater or equal to the maximum BNode
      /// id.
      int maxBNodeId() const { return -1;}

      template <typename _Graph>
      struct Constraints {

	void constraints() {
	  checkConcept<Base, _Graph >();
	  checkConcept<IDableGraphComponent<Base>, _Graph >();
	  typename _Graph::Node node(INVALID);
	  int id;
          id = graph.aNodeId(node);
          id = graph.bNodeId(node);
          node = graph.nodeFromANodeId(id);
          node = graph.nodeFromBNodeId(id);
          id = graph.maxANodeId();
          id = graph.maxBNodeId();
	}

	const _Graph& graph;
      };
    };

    /// \brief Skeleton class for graph NodeIt and EdgeIt
    ///
    /// Skeleton class for graph NodeIt and EdgeIt.
    ///
    template <typename _Graph, typename _Item>
    class GraphItemIt : public _Item {
    public:
      /// \brief Default constructor.
      ///
      /// @warning The default constructor sets the iterator
      /// to an undefined value.
      GraphItemIt() {}
      /// \brief Copy constructor.
      ///
      /// Copy constructor.
      ///
      GraphItemIt(const GraphItemIt& ) {}
      /// \brief Sets the iterator to the first item.
      ///
      /// Sets the iterator to the first item of \c the graph.
      ///
      explicit GraphItemIt(const _Graph&) {}
      /// \brief Invalid constructor \& conversion.
      ///
      /// This constructor initializes the item to be invalid.
      /// \sa Invalid for more details.
      GraphItemIt(Invalid) {}
      /// \brief Assign operator for items.
      ///
      /// The items are assignable. 
      ///
      GraphItemIt& operator=(const GraphItemIt&) { return *this; }      
      /// \brief Next item.
      /// 
      /// Assign the iterator to the next item.
      ///
      GraphItemIt& operator++() { return *this; }
      /// \brief Equality operator
      /// 
      /// Two iterators are equal if and only if they point to the
      /// same object or both are invalid.
      bool operator==(const GraphItemIt&) const { return true;}
      /// \brief Inequality operator
      ///	
      /// \sa operator==(Node n)
      ///
      bool operator!=(const GraphItemIt&) const { return true;}
      
      template<typename _GraphItemIt>
      struct Constraints {
	void constraints() {
	  _GraphItemIt it1(g);	
	  _GraphItemIt it2;

	  it2 = ++it1;
	  ++it2 = it1;
	  ++(++it1);

	  _Item bi = it1;
	  bi = it2;
	}
	_Graph& g;
      };
    };

    /// \brief Skeleton class for graph InEdgeIt and OutEdgeIt
    ///
    /// \note Because InEdgeIt and OutEdgeIt may not inherit from the same
    /// base class, the _selector is a additional template parameter. For 
    /// InEdgeIt you should instantiate it with character 'i' and for 
    /// OutEdgeIt with 'o'.
    template <typename _Graph,
	      typename _Item = typename _Graph::Edge,
              typename _Base = typename _Graph::Node, 
	      char _selector = '0'>
    class GraphIncIt : public _Item {
    public:
      /// \brief Default constructor.
      ///
      /// @warning The default constructor sets the iterator
      /// to an undefined value.
      GraphIncIt() {}
      /// \brief Copy constructor.
      ///
      /// Copy constructor.
      ///
      GraphIncIt(GraphIncIt const& gi) : _Item(gi) {}
      /// \brief Sets the iterator to the first edge incoming into or outgoing 
      /// from the node.
      ///
      /// Sets the iterator to the first edge incoming into or outgoing 
      /// from the node.
      ///
      explicit GraphIncIt(const _Graph&, const _Base&) {}
      /// \brief Invalid constructor \& conversion.
      ///
      /// This constructor initializes the item to be invalid.
      /// \sa Invalid for more details.
      GraphIncIt(Invalid) {}
      /// \brief Assign operator for iterators.
      ///
      /// The iterators are assignable. 
      ///
      GraphIncIt& operator=(GraphIncIt const&) { return *this; }      
      /// \brief Next item.
      ///
      /// Assign the iterator to the next item.
      ///
      GraphIncIt& operator++() { return *this; }

      /// \brief Equality operator
      ///
      /// Two iterators are equal if and only if they point to the
      /// same object or both are invalid.
      bool operator==(const GraphIncIt&) const { return true;}

      /// \brief Inequality operator
      ///
      /// \sa operator==(Node n)
      ///
      bool operator!=(const GraphIncIt&) const { return true;}

      template <typename _GraphIncIt>
      struct Constraints {
	void constraints() {
	  checkConcept<GraphItem<_selector>, _GraphIncIt>();
	  _GraphIncIt it1(graph, node);
	  _GraphIncIt it2;

	  it2 = ++it1;
	  ++it2 = it1;
	  ++(++it1);
	  _Item e = it1;
	  e = it2;

	}

	_Item edge;
	_Base node;
	_Graph graph;
	_GraphIncIt it;
      };
    };


    /// \brief An empty iterable graph class.
    ///
    /// This class provides beside the core graph features
    /// iterator based iterable interface for the graph structure.
    /// This concept is part of the Graph concept.
    template <typename _Base = BaseGraphComponent>
    class IterableGraphComponent : public _Base {

    public:
    
      typedef _Base Base;
      typedef typename Base::Node Node;
      typedef typename Base::Edge Edge;

      typedef IterableGraphComponent Graph;

      /// \name Base iteration
      /// 
      /// This interface provides functions for iteration on graph items
      ///
      /// @{  

      /// \brief Gives back the first node in the iterating order.
      ///      
      /// Gives back the first node in the iterating order.
      ///     
      void first(Node&) const {}

      /// \brief Gives back the next node in the iterating order.
      ///
      /// Gives back the next node in the iterating order.
      ///     
      void next(Node&) const {}

      /// \brief Gives back the first edge in the iterating order.
      ///
      /// Gives back the first edge in the iterating order.
      ///     
      void first(Edge&) const {}

      /// \brief Gives back the next edge in the iterating order.
      ///
      /// Gives back the next edge in the iterating order.
      ///     
      void next(Edge&) const {}


      /// \brief Gives back the first of the edges point to the given
      /// node.
      ///
      /// Gives back the first of the edges point to the given node.
      ///     
      void firstIn(Edge&, const Node&) const {}

      /// \brief Gives back the next of the edges points to the given
      /// node.
      ///
      /// Gives back the next of the edges points to the given node.
      ///
      void nextIn(Edge&) const {}

      /// \brief Gives back the first of the edges start from the
      /// given node.
      ///      
      /// Gives back the first of the edges start from the given node.
      ///     
      void firstOut(Edge&, const Node&) const {}

      /// \brief Gives back the next of the edges start from the given
      /// node.
      ///
      /// Gives back the next of the edges start from the given node.
      ///     
      void nextOut(Edge&) const {}

      /// @}

      /// \name Class based iteration
      /// 
      /// This interface provides functions for iteration on graph items
      ///
      /// @{

      /// \brief This iterator goes through each node.
      ///
      /// This iterator goes through each node.
      ///
      typedef GraphItemIt<Graph, Node> NodeIt;

      /// \brief This iterator goes through each node.
      ///
      /// This iterator goes through each node.
      ///
      typedef GraphItemIt<Graph, Edge> EdgeIt;

      /// \brief This iterator goes trough the incoming edges of a node.
      ///
      /// This iterator goes trough the \e inccoming edges of a certain node
      /// of a graph.
      typedef GraphIncIt<Graph, Edge, Node, 'i'> InEdgeIt;

      /// \brief This iterator goes trough the outgoing edges of a node.
      ///
      /// This iterator goes trough the \e outgoing edges of a certain node
      /// of a graph.
      typedef GraphIncIt<Graph, Edge, Node, 'o'> OutEdgeIt;

      /// \brief The base node of the iterator.
      ///
      /// Gives back the base node of the iterator.
      /// It is always the target of the pointed edge.
      Node baseNode(const InEdgeIt&) const { return INVALID; }

      /// \brief The running node of the iterator.
      ///
      /// Gives back the running node of the iterator.
      /// It is always the source of the pointed edge.
      Node runningNode(const InEdgeIt&) const { return INVALID; }

      /// \brief The base node of the iterator.
      ///
      /// Gives back the base node of the iterator.
      /// It is always the source of the pointed edge.
      Node baseNode(const OutEdgeIt&) const { return INVALID; }

      /// \brief The running node of the iterator.
      ///
      /// Gives back the running node of the iterator.
      /// It is always the target of the pointed edge.
      Node runningNode(const OutEdgeIt&) const { return INVALID; }

      /// @}

      template <typename _Graph> 
      struct Constraints {
	void constraints() {
	  checkConcept<Base, _Graph>();

          {
            typename _Graph::Node node(INVALID);      
            typename _Graph::Edge edge(INVALID);
            {
              graph.first(node);
              graph.next(node);
            }
            {
              graph.first(edge);
              graph.next(edge);
            }
            {
              graph.firstIn(edge, node);
              graph.nextIn(edge);
            }
            {
              graph.firstOut(edge, node);
              graph.nextOut(edge);
            }
          }           

          {
            checkConcept<GraphItemIt<_Graph, typename _Graph::Edge>,
              typename _Graph::EdgeIt >();
            checkConcept<GraphItemIt<_Graph, typename _Graph::Node>,
              typename _Graph::NodeIt >();
            checkConcept<GraphIncIt<_Graph, typename _Graph::Edge, 
              typename _Graph::Node, 'i'>, typename _Graph::InEdgeIt>();
            checkConcept<GraphIncIt<_Graph, typename _Graph::Edge, 
              typename _Graph::Node, 'o'>, typename _Graph::OutEdgeIt>();

            typename _Graph::Node n;
            typename _Graph::InEdgeIt ieit(INVALID);
            typename _Graph::OutEdgeIt oeit(INVALID);
            n = graph.baseNode(ieit);
            n = graph.runningNode(ieit);
            n = graph.baseNode(oeit);
            n = graph.runningNode(oeit);
            ignore_unused_variable_warning(n);
          }
        }
	
	const _Graph& graph;
	
      };
    };

    /// \brief An empty iterable undirected graph class.
    ///
    /// This class provides beside the core graph features iterator
    /// based iterable interface for the undirected graph structure.
    /// This concept is part of the UGraph concept.
    template <typename _Base = BaseUGraphComponent>
    class IterableUGraphComponent : public IterableGraphComponent<_Base> {
    public:

      typedef _Base Base;
      typedef typename Base::Node Node;
      typedef typename Base::Edge Edge;
      typedef typename Base::UEdge UEdge;

    
      typedef IterableUGraphComponent Graph;

      /// \name Base iteration
      /// 
      /// This interface provides functions for iteration on graph items
      /// @{  

      using IterableGraphComponent<_Base>::first;
      using IterableGraphComponent<_Base>::next;

      /// \brief Gives back the first undirected edge in the iterating
      /// order.
      ///
      /// Gives back the first undirected edge in the iterating order.
      ///     
      void first(UEdge&) const {}

      /// \brief Gives back the next undirected edge in the iterating
      /// order.
      ///
      /// Gives back the next undirected edge in the iterating order.
      ///     
      void next(UEdge&) const {}


      /// \brief Gives back the first of the undirected edges from the
      /// given node.
      ///
      /// Gives back the first of the undirected edges from the given
      /// node. The bool parameter gives back that direction which
      /// gives a good direction of the uedge so the source of the
      /// directed edge is the given node.
      void firstInc(UEdge&, bool&, const Node&) const {}

      /// \brief Gives back the next of the undirected edges from the
      /// given node.
      ///
      /// Gives back the next of the undirected edges from the given
      /// node. The bool parameter should be used as the \c firstInc()
      /// use it.
      void nextInc(UEdge&, bool&) const {}

      using IterableGraphComponent<_Base>::baseNode;
      using IterableGraphComponent<_Base>::runningNode;

      /// @}

      /// \name Class based iteration
      /// 
      /// This interface provides functions for iteration on graph items
      ///
      /// @{

      /// \brief This iterator goes through each node.
      ///
      /// This iterator goes through each node.
      typedef GraphItemIt<Graph, UEdge> UEdgeIt;
      /// \brief This iterator goes trough the incident edges of a
      /// node.
      ///
      /// This iterator goes trough the incident edges of a certain
      /// node of a graph.
      typedef GraphIncIt<Graph, UEdge, Node, 'u'> IncEdgeIt;
      /// \brief The base node of the iterator.
      ///
      /// Gives back the base node of the iterator.
      Node baseNode(const IncEdgeIt&) const { return INVALID; }

      /// \brief The running node of the iterator.
      ///
      /// Gives back the running node of the iterator.
      Node runningNode(const IncEdgeIt&) const { return INVALID; }

      /// @}

      template <typename _Graph> 
      struct Constraints {
	void constraints() {
	  checkConcept<IterableGraphComponent<Base>, _Graph>();

          {
            typename _Graph::Node node(INVALID);
            typename _Graph::UEdge uedge(INVALID);
            bool dir;
            {
              graph.first(uedge);
              graph.next(uedge);
            }
            {
              graph.firstInc(uedge, dir, node);
              graph.nextInc(uedge, dir);
            }
            
          }	
  
          {
            checkConcept<GraphItemIt<_Graph, typename _Graph::UEdge>,
              typename _Graph::UEdgeIt >();
            checkConcept<GraphIncIt<_Graph, typename _Graph::UEdge, 
              typename _Graph::Node, 'u'>, typename _Graph::IncEdgeIt>();
            
            typename _Graph::Node n;
            typename _Graph::IncEdgeIt ueit(INVALID);
            n = graph.baseNode(ueit);
            n = graph.runningNode(ueit);
          }
        }
	
	const _Graph& graph;
	
      };
    };

    /// \brief An empty iterable bipartite undirected graph class.
    ///
    /// This class provides beside the core graph features iterator
    /// based iterable interface for the bipartite undirected graph
    /// structure. This concept is part of the BpUGraph concept.
    template <typename _Base = BaseUGraphComponent>
    class IterableBpUGraphComponent : public IterableUGraphComponent<_Base> {
    public:

      typedef _Base Base;
      typedef typename Base::Node Node;
      typedef typename Base::UEdge UEdge;
    
      typedef IterableBpUGraphComponent Graph;

      /// \name Base iteration
      /// 
      /// This interface provides functions for iteration on graph items
      /// @{  

      using IterableUGraphComponent<_Base>::first;
      using IterableUGraphComponent<_Base>::next;

      /// \brief Gives back the first A-node in the iterating order.
      ///
      /// Gives back the first undirected A-node in the iterating
      /// order.
      ///     
      void firstANode(Node&) const {}

      /// \brief Gives back the next A-node in the iterating order.
      ///
      /// Gives back the next A-node in the iterating order.
      ///     
      void nextANode(Node&) const {}

      /// \brief Gives back the first B-node in the iterating order.
      ///
      /// Gives back the first undirected B-node in the iterating
      /// order.
      ///     
      void firstBNode(Node&) const {}

      /// \brief Gives back the next B-node in the iterating order.
      ///
      /// Gives back the next B-node in the iterating order.
      ///     
      void nextBNode(Node&) const {}


      /// \brief Gives back the first of the undirected edges start
      /// from the given A-node.
      ///      
      /// Gives back the first of the undirected edges start from the
      /// given A-node.
      void firstFromANode(UEdge&, const Node&) const {}

      /// \brief Gives back the next of the undirected edges start
      /// from the given A-node.
      ///      
      /// Gives back the next of the undirected edges start from the
      /// given A-node.
      void nextFromANode(UEdge&) const {}

      /// \brief Gives back the first of the undirected edges start
      /// from the given B-node.
      ///      
      /// Gives back the first of the undirected edges start from the
      /// given B-node.
      void firstFromBNode(UEdge&, const Node&) const {}

      /// \brief Gives back the next of the undirected edges start
      /// from the given B-node.
      ///      
      /// Gives back the next of the undirected edges start from the
      /// given B-node.
      void nextFromBNode(UEdge&) const {}


      /// @}

      /// \name Class based iteration
      /// 
      /// This interface provides functions for iteration on graph items
      ///
      /// @{

      /// \brief This iterator goes through each A-node.
      ///
      /// This iterator goes through each A-node.
      typedef GraphItemIt<Graph, Node> ANodeIt;

      /// \brief This iterator goes through each B-node.
      ///
      /// This iterator goes through each B-node.
      typedef GraphItemIt<Graph, Node> BNodeIt;

      /// @}

      template <typename _Graph> 
      struct Constraints {
	void constraints() {
	  checkConcept<IterableUGraphComponent<Base>, _Graph>();

          {
            typename _Graph::Node node(INVALID);
            typename _Graph::UEdge uedge(INVALID);
            graph.firstANode(node);
            graph.nextANode(node);
            graph.firstBNode(node);
            graph.nextBNode(node);

            graph.firstFromANode(uedge, node);
            graph.nextFromANode(uedge);
            graph.firstFromBNode(uedge, node);
            graph.nextFromBNode(uedge);
          }
          {
            checkConcept<GraphItemIt<_Graph, typename _Graph::Node>,
              typename _Graph::ANodeIt >();
            checkConcept<GraphItemIt<_Graph, typename _Graph::Node>,
              typename _Graph::BNodeIt >();
          }

	}
	
	const _Graph& graph;
	
      };
    };

    /// \brief An empty alteration notifier graph class.
    ///  
    /// This class provides beside the core graph features alteration
    /// notifier interface for the graph structure.  This implements
    /// an observer-notifier pattern for each graph item. More
    /// obsevers can be registered into the notifier and whenever an
    /// alteration occured in the graph all the observers will
    /// notified about it.
    template <typename _Base = BaseGraphComponent>
    class AlterableGraphComponent : public _Base {
    public:

      typedef _Base Base;
      typedef typename Base::Node Node;
      typedef typename Base::Edge Edge;


      /// The node observer registry.
      typedef AlterationNotifier<AlterableGraphComponent, Node> 
      NodeNotifier;
      /// The edge observer registry.
      typedef AlterationNotifier<AlterableGraphComponent, Edge> 
      EdgeNotifier;
      
      /// \brief Gives back the node alteration notifier.
      ///
      /// Gives back the node alteration notifier.
      NodeNotifier& notifier(Node) const {
	return NodeNotifier();
      }
      
      /// \brief Gives back the edge alteration notifier.
      ///
      /// Gives back the edge alteration notifier.
      EdgeNotifier& notifier(Edge) const {
	return EdgeNotifier();
      }

      template <typename _Graph> 
      struct Constraints {
	void constraints() {
	  checkConcept<Base, _Graph>();
          typename _Graph::NodeNotifier& nn 
            = graph.notifier(typename _Graph::Node());

          typename _Graph::EdgeNotifier& en 
            = graph.notifier(typename _Graph::Edge());
          
          ignore_unused_variable_warning(nn);
          ignore_unused_variable_warning(en);
	}
	
	const _Graph& graph;
	
      };
      
    };

    /// \brief An empty alteration notifier undirected graph class.
    ///  
    /// This class provides beside the core graph features alteration
    /// notifier interface for the graph structure.  This implements
    /// an observer-notifier pattern for each graph item. More
    /// obsevers can be registered into the notifier and whenever an
    /// alteration occured in the graph all the observers will
    /// notified about it.
    template <typename _Base = BaseUGraphComponent>
    class AlterableUGraphComponent : public AlterableGraphComponent<_Base> {
    public:

      typedef _Base Base;
      typedef typename Base::UEdge UEdge;


      /// The edge observer registry.
      typedef AlterationNotifier<AlterableUGraphComponent, UEdge> 
      UEdgeNotifier;
      
      /// \brief Gives back the edge alteration notifier.
      ///
      /// Gives back the edge alteration notifier.
      UEdgeNotifier& notifier(UEdge) const {
	return UEdgeNotifier();
      }

      template <typename _Graph> 
      struct Constraints {
	void constraints() {
	  checkConcept<AlterableGraphComponent<Base>, _Graph>();
          typename _Graph::UEdgeNotifier& uen 
            = graph.notifier(typename _Graph::UEdge());
          ignore_unused_variable_warning(uen);
	}
	
	const _Graph& graph;
	
      };
      
    };

    /// \brief An empty alteration notifier bipartite undirected graph
    /// class.
    ///  
    /// This class provides beside the core graph features alteration
    /// notifier interface for the graph structure.  This implements
    /// an observer-notifier pattern for each graph item. More
    /// obsevers can be registered into the notifier and whenever an
    /// alteration occured in the graph all the observers will
    /// notified about it.
    template <typename _Base = BaseUGraphComponent>
    class AlterableBpUGraphComponent : public AlterableUGraphComponent<_Base> {
    public:

      typedef _Base Base;
      typedef typename Base::ANode ANode;
      typedef typename Base::BNode BNode;


      /// The A-node observer registry.
      typedef AlterationNotifier<AlterableBpUGraphComponent, ANode> 
      ANodeNotifier;

      /// The B-node observer registry.
      typedef AlterationNotifier<AlterableBpUGraphComponent, BNode> 
      BNodeNotifier;
      
      /// \brief Gives back the A-node alteration notifier.
      ///
      /// Gives back the A-node alteration notifier.
      ANodeNotifier& notifier(ANode) const {
	return ANodeNotifier();
      }

      /// \brief Gives back the B-node alteration notifier.
      ///
      /// Gives back the B-node alteration notifier.
      BNodeNotifier& notifier(BNode) const {
	return BNodeNotifier();
      }

      template <typename _Graph> 
      struct Constraints {
	void constraints() {
          checkConcept<AlterableUGraphComponent<Base>, _Graph>();
          typename _Graph::ANodeNotifier& ann 
            = graph.notifier(typename _Graph::ANode());
          typename _Graph::BNodeNotifier& bnn 
            = graph.notifier(typename _Graph::BNode());
          ignore_unused_variable_warning(ann);
          ignore_unused_variable_warning(bnn);
	}
	
	const _Graph& graph;
	
      };
      
    };


    /// \brief Class describing the concept of graph maps
    /// 
    /// This class describes the common interface of the graph maps
    /// (NodeMap, EdgeMap), that is \ref maps-page "maps" which can be used to
    /// associate data to graph descriptors (nodes or edges).
    template <typename _Graph, typename _Item, typename _Value>
    class GraphMap : public ReadWriteMap<_Item, _Value> {
    public:

      typedef ReadWriteMap<_Item, _Value> Parent;

      /// The graph type of the map.
      typedef _Graph Graph;
      /// The key type of the map.
      typedef _Item Key;
      /// The value type of the map.
      typedef _Value Value;

      /// \brief Construct a new map.
      ///
      /// Construct a new map for the graph.
      explicit GraphMap(const Graph&) {}
      /// \brief Construct a new map with default value.
      ///
      /// Construct a new map for the graph and initalise the values.
      GraphMap(const Graph&, const Value&) {}
      /// \brief Copy constructor.
      ///
      /// Copy Constructor.
      GraphMap(const GraphMap&) : Parent() {}
      
      /// \brief Assign operator.
      ///
      /// Assign operator. It does not mofify the underlying graph,
      /// it just iterates on the current item set and set the  map
      /// with the value returned by the assigned map. 
      template <typename CMap>
      GraphMap& operator=(const CMap&) { 
        checkConcept<ReadMap<Key, Value>, CMap>();
        return *this;
      }

      template<typename _Map>
      struct Constraints {
	void constraints() {
	  checkConcept<ReadWriteMap<Key, Value>, _Map >();
	  // Construction with a graph parameter
	  _Map a(g);
	  // Constructor with a graph and a default value parameter
	  _Map a2(g,t);
	  // Copy constructor.
	  _Map b(c);
          
          ReadMap<Key, Value> cmap;
          b = cmap;

	  ignore_unused_variable_warning(a2);
	  ignore_unused_variable_warning(b);
	}

	const _Map &c;
	const Graph &g;
	const typename GraphMap::Value &t;
      };

    };

    /// \brief An empty mappable graph class.
    ///
    /// This class provides beside the core graph features
    /// map interface for the graph structure.
    /// This concept is part of the Graph concept.
    template <typename _Base = BaseGraphComponent>
    class MappableGraphComponent : public _Base  {
    public:

      typedef _Base Base;
      typedef typename Base::Node Node;
      typedef typename Base::Edge Edge;

      typedef MappableGraphComponent Graph;

      /// \brief ReadWrite map of the nodes.
      ///
      /// ReadWrite map of the nodes.
      ///
      template <typename _Value>
      class NodeMap : public GraphMap<Graph, Node, _Value> {
      public:
        typedef GraphMap<MappableGraphComponent, Node, _Value> Parent;

	/// \brief Construct a new map.
	///
	/// Construct a new map for the graph.
	explicit NodeMap(const MappableGraphComponent& graph) 
          : Parent(graph) {}

	/// \brief Construct a new map with default value.
	///
	/// Construct a new map for the graph and initalise the values.
	NodeMap(const MappableGraphComponent& graph, const _Value& value)
          : Parent(graph, value) {}

	/// \brief Copy constructor.
	///
	/// Copy Constructor.
	NodeMap(const NodeMap& nm) : Parent(nm) {}

	/// \brief Assign operator.
	///
	/// Assign operator.
        template <typename CMap>
        NodeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<Node, _Value>, CMap>();
          return *this;
        }

      };

      /// \brief ReadWrite map of the edges.
      ///
      /// ReadWrite map of the edges.
      ///
      template <typename _Value>
      class EdgeMap : public GraphMap<Graph, Edge, _Value> {
      public:
        typedef GraphMap<MappableGraphComponent, Edge, _Value> Parent;

	/// \brief Construct a new map.
	///
	/// Construct a new map for the graph.
	explicit EdgeMap(const MappableGraphComponent& graph) 
          : Parent(graph) {}

	/// \brief Construct a new map with default value.
	///
	/// Construct a new map for the graph and initalise the values.
	EdgeMap(const MappableGraphComponent& graph, const _Value& value)
          : Parent(graph, value) {}

	/// \brief Copy constructor.
	///
	/// Copy Constructor.
	EdgeMap(const EdgeMap& nm) : Parent(nm) {}

	/// \brief Assign operator.
	///
	/// Assign operator.
        template <typename CMap>
        EdgeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<Edge, _Value>, CMap>();
          return *this;
        }

      };


      template <typename _Graph>
      struct Constraints {

	struct Dummy {
	  int value;
	  Dummy() : value(0) {}
	  Dummy(int _v) : value(_v) {}
	};

	void constraints() {
	  checkConcept<Base, _Graph>();
	  { // int map test
	    typedef typename _Graph::template NodeMap<int> IntNodeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::Node, int>, 
	      IntNodeMap >();
	  } { // bool map test
	    typedef typename _Graph::template NodeMap<bool> BoolNodeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::Node, bool>,
	      BoolNodeMap >();
	  } { // Dummy map test
	    typedef typename _Graph::template NodeMap<Dummy> DummyNodeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::Node, Dummy>,
	      DummyNodeMap >();
	  } 

	  { // int map test
	    typedef typename _Graph::template EdgeMap<int> IntEdgeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::Edge, int>,
	      IntEdgeMap >();
	  } { // bool map test
	    typedef typename _Graph::template EdgeMap<bool> BoolEdgeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::Edge, bool>,
	      BoolEdgeMap >();
	  } { // Dummy map test
	    typedef typename _Graph::template EdgeMap<Dummy> DummyEdgeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::Edge, Dummy>, 
	      DummyEdgeMap >();
	  } 
	}

	_Graph& graph;
      };
    };

    /// \brief An empty mappable base bipartite undirected graph class.
    ///
    /// This class provides beside the core graph features
    /// map interface for the graph structure.
    /// This concept is part of the UGraph concept.
    template <typename _Base = BaseUGraphComponent>
    class MappableUGraphComponent : public MappableGraphComponent<_Base>  {
    public:

      typedef _Base Base;
      typedef typename Base::UEdge UEdge;

      typedef MappableUGraphComponent Graph;

      /// \brief ReadWrite map of the uedges.
      ///
      /// ReadWrite map of the uedges.
      ///
      template <typename _Value>
      class UEdgeMap : public GraphMap<Graph, UEdge, _Value> {  
      public:
        typedef GraphMap<MappableUGraphComponent, UEdge, _Value> Parent;

	/// \brief Construct a new map.
	///
	/// Construct a new map for the graph.
	explicit UEdgeMap(const MappableUGraphComponent& graph) 
          : Parent(graph) {}

	/// \brief Construct a new map with default value.
	///
	/// Construct a new map for the graph and initalise the values.
	UEdgeMap(const MappableUGraphComponent& graph, const _Value& value)
          : Parent(graph, value) {}

	/// \brief Copy constructor.
	///
	/// Copy Constructor.
	UEdgeMap(const UEdgeMap& nm) : Parent(nm) {}

	/// \brief Assign operator.
	///
	/// Assign operator.
        template <typename CMap>
        UEdgeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<UEdge, _Value>, CMap>();
          return *this;
        }

      };


      template <typename _Graph>
      struct Constraints {

	struct Dummy {
	  int value;
	  Dummy() : value(0) {}
	  Dummy(int _v) : value(_v) {}
	};

	void constraints() {
	  checkConcept<MappableGraphComponent<Base>, _Graph>();

	  { // int map test
	    typedef typename _Graph::template UEdgeMap<int> IntUEdgeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::UEdge, int>,
	      IntUEdgeMap >();
	  } { // bool map test
	    typedef typename _Graph::template UEdgeMap<bool> BoolUEdgeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::UEdge, bool>,
	      BoolUEdgeMap >();
	  } { // Dummy map test
	    typedef typename _Graph::template UEdgeMap<Dummy> DummyUEdgeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::UEdge, Dummy>, 
	      DummyUEdgeMap >();
	  } 
	}

	_Graph& graph;
      };
    };

    /// \brief An empty mappable base bipartite undirected graph
    /// class.
    ///
    /// This class provides beside the core graph features
    /// map interface for the graph structure.
    /// This concept is part of the BpUGraph concept.
    template <typename _Base = BaseBpUGraphComponent>
    class MappableBpUGraphComponent : public MappableUGraphComponent<_Base>  {
    public:

      typedef _Base Base;
      typedef typename Base::Node Node;

      typedef MappableBpUGraphComponent Graph;

      /// \brief ReadWrite map of the A-nodes.
      ///
      /// ReadWrite map of the A-nodes.
      ///
      template <typename _Value>
      class ANodeMap : public GraphMap<Graph, Node, _Value> {  
      public:
        typedef GraphMap<MappableBpUGraphComponent, Node, _Value> Parent;

	/// \brief Construct a new map.
	///
	/// Construct a new map for the graph.
	explicit ANodeMap(const MappableBpUGraphComponent& graph) 
          : Parent(graph) {}

	/// \brief Construct a new map with default value.
	///
	/// Construct a new map for the graph and initalise the values.
	ANodeMap(const MappableBpUGraphComponent& graph, const _Value& value)
          : Parent(graph, value) {}

	/// \brief Copy constructor.
	///
	/// Copy Constructor.
	ANodeMap(const ANodeMap& nm) : Parent(nm) {}

	/// \brief Assign operator.
	///
	/// Assign operator.
        template <typename CMap>
        ANodeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<Node, _Value>, CMap>();
          return *this;
        }

      };

      /// \brief ReadWrite map of the B-nodes.
      ///
      /// ReadWrite map of the A-nodes.
      ///
      template <typename _Value>
      class BNodeMap : public GraphMap<Graph, Node, _Value> {  
      public:
        typedef GraphMap<MappableBpUGraphComponent, Node, _Value> Parent;

	/// \brief Construct a new map.
	///
	/// Construct a new map for the graph.
	explicit BNodeMap(const MappableBpUGraphComponent& graph) 
          : Parent(graph) {}

	/// \brief Construct a new map with default value.
	///
	/// Construct a new map for the graph and initalise the values.
	BNodeMap(const MappableBpUGraphComponent& graph, const _Value& value)
          : Parent(graph, value) {}

	/// \brief Copy constructor.
	///
	/// Copy Constructor.
	BNodeMap(const BNodeMap& nm) : Parent(nm) {}

	/// \brief Assign operator.
	///
	/// Assign operator.
        template <typename CMap>
        BNodeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<Node, _Value>, CMap>();
          return *this;
        }

      };


      template <typename _Graph>
      struct Constraints {

	struct Dummy {
	  int value;
	  Dummy() : value(0) {}
	  Dummy(int _v) : value(_v) {}
	};

	void constraints() {
	  checkConcept<MappableUGraphComponent<Base>, _Graph>();

	  { // int map test
	    typedef typename _Graph::template ANodeMap<int> IntANodeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::ANode, int>,
	      IntANodeMap >();
	  } { // bool map test
	    typedef typename _Graph::template ANodeMap<bool> BoolANodeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::ANode, bool>,
	      BoolANodeMap >();
	  } { // Dummy map test
	    typedef typename _Graph::template ANodeMap<Dummy> DummyANodeMap;
	    checkConcept<GraphMap<_Graph, typename _Graph::ANode, Dummy>, 
	      DummyANodeMap >();
	  } 
	}

	_Graph& graph;
      };
    };


    /// \brief An empty extendable graph class.
    ///
    /// This class provides beside the core graph features graph
    /// extendable interface for the graph structure.  The main
    /// difference between the base and this interface is that the
    /// graph alterations should handled already on this level.
    template <typename _Base = BaseGraphComponent>
    class ExtendableGraphComponent : public _Base {
    public:
      typedef _Base Base;

      typedef typename _Base::Node Node;
      typedef typename _Base::Edge Edge;

      /// \brief Adds a new node to the graph.
      ///
      /// Adds a new node to the graph.
      ///
      Node addNode() {
	return INVALID;
      }
    
      /// \brief Adds a new edge connects the given two nodes.
      ///
      /// Adds a new edge connects the the given two nodes.
      Edge addEdge(const Node&, const Node&) {
	return INVALID;
      }

      template <typename _Graph>
      struct Constraints {
	void constraints() {
          checkConcept<Base, _Graph>();
	  typename _Graph::Node node_a, node_b;
	  node_a = graph.addNode();
	  node_b = graph.addNode();
	  typename _Graph::Edge edge;
	  edge = graph.addEdge(node_a, node_b);
	}

	_Graph& graph;
      };
    };

    /// \brief An empty extendable base undirected graph class.
    ///
    /// This class provides beside the core undirected graph features
    /// core undircted graph extend interface for the graph structure.
    /// The main difference between the base and this interface is
    /// that the graph alterations should handled already on this
    /// level.
    template <typename _Base = BaseUGraphComponent>
    class ExtendableUGraphComponent : public _Base {
    public:

      typedef _Base Base;
      typedef typename _Base::Node Node;
      typedef typename _Base::UEdge UEdge;

      /// \brief Adds a new node to the graph.
      ///
      /// Adds a new node to the graph.
      ///
      Node addNode() {
	return INVALID;
      }
    
      /// \brief Adds a new edge connects the given two nodes.
      ///
      /// Adds a new edge connects the the given two nodes.
      UEdge addEdge(const Node&, const Node&) {
	return INVALID;
      }

      template <typename _Graph>
      struct Constraints {
	void constraints() {
	  checkConcept<Base, _Graph>();
	  typename _Graph::Node node_a, node_b;
	  node_a = graph.addNode();
	  node_b = graph.addNode();
	  typename _Graph::UEdge uedge;
	  uedge = graph.addUEdge(node_a, node_b);
	}

	_Graph& graph;
      };
    };

    /// \brief An empty extendable base undirected graph class.
    ///
    /// This class provides beside the core bipartite undirected graph
    /// features core undircted graph extend interface for the graph
    /// structure.  The main difference between the base and this
    /// interface is that the graph alterations should handled already
    /// on this level.
    template <typename _Base = BaseBpUGraphComponent>
    class ExtendableBpUGraphComponent 
      : public ExtendableUGraphComponent<_Base> {

      typedef _Base Base;

      template <typename _Graph>
      struct Constraints {
	void constraints() {
          checkConcept<ExtendableUGraphComponent<Base>, _Graph>();
	}
      };
    };

    /// \brief An empty erasable graph class.
    ///  
    /// This class provides beside the core graph features core erase
    /// functions for the graph structure. The main difference between
    /// the base and this interface is that the graph alterations
    /// should handled already on this level.
    template <typename _Base = BaseGraphComponent>
    class ErasableGraphComponent : public _Base {
    public:

      typedef _Base Base;
      typedef typename Base::Node Node;
      typedef typename Base::Edge Edge;

      /// \brief Erase a node from the graph.
      ///
      /// Erase a node from the graph. This function should 
      /// erase all edges connecting to the node.
      void erase(const Node&) {}    

      /// \brief Erase an edge from the graph.
      ///
      /// Erase an edge from the graph.
      ///
      void erase(const Edge&) {}

      template <typename _Graph>
      struct Constraints {
	void constraints() {
          checkConcept<Base, _Graph>();
	  typename _Graph::Node node;
	  graph.erase(node);
	  typename _Graph::Edge edge;
	  graph.erase(edge);
	}

	_Graph& graph;
      };
    };

    /// \brief An empty erasable base undirected graph class.
    ///  
    /// This class provides beside the core undirected graph features
    /// core erase functions for the undirceted graph structure. The
    /// main difference between the base and this interface is that
    /// the graph alterations should handled already on this level.
    template <typename _Base = BaseUGraphComponent>
    class ErasableUGraphComponent : public _Base {
    public:

      typedef _Base Base;
      typedef typename Base::Node Node;
      typedef typename Base::UEdge UEdge;

      /// \brief Erase a node from the graph.
      ///
      /// Erase a node from the graph. This function should erase
      /// edges connecting to the node.
      void erase(const Node&) {}    

      /// \brief Erase an edge from the graph.
      ///
      /// Erase an edge from the graph.
      ///
      void erase(const UEdge&) {}

      template <typename _Graph>
      struct Constraints {
	void constraints() {
          checkConcept<Base, _Graph>();
	  typename _Graph::Node node;
	  graph.erase(node);
	  typename _Graph::Edge edge;
	  graph.erase(edge);
	}

	_Graph& graph;
      };
    };

    /// \brief An empty erasable base bipartite undirected graph class.
    ///  
    /// This class provides beside the core bipartite undirected graph
    /// features core erase functions for the undirceted graph
    /// structure. The main difference between the base and this
    /// interface is that the graph alterations should handled already
    /// on this level.
    template <typename _Base = BaseBpUGraphComponent>
    class ErasableBpUGraphComponent : public ErasableUGraphComponent<_Base> {
    public:

      typedef _Base Base;

      template <typename _Graph>
      struct Constraints {
	void constraints() {
          checkConcept<ErasableUGraphComponent<Base>, _Graph>();
	}
      };
    };

    /// \brief An empty clearable base graph class.
    ///
    /// This class provides beside the core graph features core clear
    /// functions for the graph structure. The main difference between
    /// the base and this interface is that the graph alterations
    /// should handled already on this level.
    template <typename _Base = BaseGraphComponent>
    class ClearableGraphComponent : public _Base {
    public:

      typedef _Base Base;

      /// \brief Erase all nodes and edges from the graph.
      ///
      /// Erase all nodes and edges from the graph.
      ///
      void clear() {}    

      template <typename _Graph>
      struct Constraints {
	void constraints() {
          checkConcept<Base, _Graph>();
	  graph.clear();
	}

	_Graph graph;
      };
    };

    /// \brief An empty clearable base undirected graph class.
    ///
    /// This class provides beside the core undirected graph features
    /// core clear functions for the undirected graph structure. The
    /// main difference between the base and this interface is that
    /// the graph alterations should handled already on this level.
    template <typename _Base = BaseUGraphComponent>
    class ClearableUGraphComponent : public ClearableGraphComponent<_Base> {
    public:

      typedef _Base Base;

      template <typename _Graph>
      struct Constraints {
	void constraints() {
          checkConcept<ClearableUGraphComponent<Base>, _Graph>();
	}

	_Graph graph;
      };
    };

    /// \brief An empty clearable base bipartite undirected graph
    /// class.
    ///
    /// This class provides beside the core bipartite undirected graph
    /// features core clear functions for the undirected graph
    /// structure. The main difference between the base and this
    /// interface is that the graph alterations should handled already
    /// on this level.
    template <typename _Base = BaseUGraphComponent>
    class ClearableBpUGraphComponent : public ClearableUGraphComponent<_Base> {
    public:

      typedef _Base Base;

      template <typename _Graph>
      struct Constraints {
	void constraints() {
          checkConcept<ClearableBpUGraphComponent<Base>, _Graph>();
	}

      };

    };

  }

}

#endif

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

#ifndef LEMON_CONCEPT_GRAPH_H
#define LEMON_CONCEPT_GRAPH_H

///\ingroup graph_concepts
///\file
///\brief The concept of Directed Graphs.

#include <lemon/bits/invalid.h>
#include <lemon/bits/utility.h>
#include <lemon/concepts/maps.h>
#include <lemon/concept_check.h>
#include <lemon/concepts/graph_components.h>

namespace lemon {
  namespace concepts {

    /// \ingroup graph_concepts
    ///
    /// \brief Class describing the concept of Directed Graphs.
    ///
    /// This class describes the \ref concept "concept" of the
    /// immutable directed graphs.
    ///
    /// Note that actual graph implementation like @ref ListGraph or
    /// @ref SmartGraph may have several additional functionality.
    ///
    /// \sa concept
    class Graph {
    private:
      ///Graphs are \e not copy constructible. Use GraphCopy() instead.
      
      ///Graphs are \e not copy constructible. Use GraphCopy() instead.
      ///
      Graph(const Graph &) {};
      ///\brief Assignment of \ref Graph "Graph"s to another ones are
      ///\e not allowed. Use GraphCopy() instead.
      
      ///Assignment of \ref Graph "Graph"s to another ones are
      ///\e not allowed.  Use GraphCopy() instead.

      void operator=(const Graph &) {}
    public:
      ///\e

      /// Defalult constructor.

      /// Defalult constructor.
      ///
      Graph() { }
      /// Class for identifying a node of the graph

      /// This class identifies a node of the graph. It also serves
      /// as a base class of the node iterators,
      /// thus they will convert to this type.
      class Node {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        Node() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        Node(const Node&) { }

        /// Invalid constructor \& conversion.

        /// This constructor initializes the iterator to be invalid.
        /// \sa Invalid for more details.
        Node(Invalid) { }
        /// Equality operator

        /// Two iterators are equal if and only if they point to the
        /// same object or both are invalid.
        bool operator==(Node) const { return true; }

        /// Inequality operator
        
        /// \sa operator==(Node n)
        ///
        bool operator!=(Node) const { return true; }

	/// Artificial ordering operator.
	
	/// To allow the use of graph descriptors as key type in std::map or
	/// similar associative container we require this.
	///
	/// \note This operator only have to define some strict ordering of
	/// the items; this order has nothing to do with the iteration
	/// ordering of the items.
	bool operator<(Node) const { return false; }

      };
    
      /// This iterator goes through each node.

      /// This iterator goes through each node.
      /// Its usage is quite simple, for example you can count the number
      /// of nodes in graph \c g of type \c Graph like this:
      ///\code
      /// int count=0;
      /// for (Graph::NodeIt n(g); n!=INVALID; ++n) ++count;
      ///\endcode
      class NodeIt : public Node {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        NodeIt() { }
        /// Copy constructor.
        
        /// Copy constructor.
        ///
        NodeIt(const NodeIt& n) : Node(n) { }
        /// Invalid constructor \& conversion.

        /// Initialize the iterator to be invalid.
        /// \sa Invalid for more details.
        NodeIt(Invalid) { }
        /// Sets the iterator to the first node.

        /// Sets the iterator to the first node of \c g.
        ///
        NodeIt(const Graph&) { }
        /// Node -> NodeIt conversion.

        /// Sets the iterator to the node of \c the graph pointed by 
	/// the trivial iterator.
        /// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        NodeIt(const Graph&, const Node&) { }
        /// Next node.

        /// Assign the iterator to the next node.
        ///
        NodeIt& operator++() { return *this; }
      };
    
    
      /// Class for identifying an edge of the graph

      /// This class identifies an edge of the graph. It also serves
      /// as a base class of the edge iterators,
      /// thus they will convert to this type.
      class Edge {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        Edge() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        Edge(const Edge&) { }
        /// Initialize the iterator to be invalid.

        /// Initialize the iterator to be invalid.
        ///
        Edge(Invalid) { }
        /// Equality operator

        /// Two iterators are equal if and only if they point to the
        /// same object or both are invalid.
        bool operator==(Edge) const { return true; }
        /// Inequality operator

        /// \sa operator==(Edge n)
        ///
        bool operator!=(Edge) const { return true; }

	/// Artificial ordering operator.
	
	/// To allow the use of graph descriptors as key type in std::map or
	/// similar associative container we require this.
	///
	/// \note This operator only have to define some strict ordering of
	/// the items; this order has nothing to do with the iteration
	/// ordering of the items.
	bool operator<(Edge) const { return false; }
      };
    
      /// This iterator goes trough the outgoing edges of a node.

      /// This iterator goes trough the \e outgoing edges of a certain node
      /// of a graph.
      /// Its usage is quite simple, for example you can count the number
      /// of outgoing edges of a node \c n
      /// in graph \c g of type \c Graph as follows.
      ///\code
      /// int count=0;
      /// for (Graph::OutEdgeIt e(g, n); e!=INVALID; ++e) ++count;
      ///\endcode
    
      class OutEdgeIt : public Edge {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        OutEdgeIt() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        OutEdgeIt(const OutEdgeIt& e) : Edge(e) { }
        /// Initialize the iterator to be invalid.

        /// Initialize the iterator to be invalid.
        ///
        OutEdgeIt(Invalid) { }
        /// This constructor sets the iterator to the first outgoing edge.
    
        /// This constructor sets the iterator to the first outgoing edge of
        /// the node.
        OutEdgeIt(const Graph&, const Node&) { }
        /// Edge -> OutEdgeIt conversion

        /// Sets the iterator to the value of the trivial iterator.
	/// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        OutEdgeIt(const Graph&, const Edge&) { }
        ///Next outgoing edge
        
        /// Assign the iterator to the next 
        /// outgoing edge of the corresponding node.
        OutEdgeIt& operator++() { return *this; }
      };

      /// This iterator goes trough the incoming edges of a node.

      /// This iterator goes trough the \e incoming edges of a certain node
      /// of a graph.
      /// Its usage is quite simple, for example you can count the number
      /// of outgoing edges of a node \c n
      /// in graph \c g of type \c Graph as follows.
      ///\code
      /// int count=0;
      /// for(Graph::InEdgeIt e(g, n); e!=INVALID; ++e) ++count;
      ///\endcode

      class InEdgeIt : public Edge {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        InEdgeIt() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        InEdgeIt(const InEdgeIt& e) : Edge(e) { }
        /// Initialize the iterator to be invalid.

        /// Initialize the iterator to be invalid.
        ///
        InEdgeIt(Invalid) { }
        /// This constructor sets the iterator to first incoming edge.
    
        /// This constructor set the iterator to the first incoming edge of
        /// the node.
        InEdgeIt(const Graph&, const Node&) { }
        /// Edge -> InEdgeIt conversion

        /// Sets the iterator to the value of the trivial iterator \c e.
        /// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        InEdgeIt(const Graph&, const Edge&) { }
        /// Next incoming edge

        /// Assign the iterator to the next inedge of the corresponding node.
        ///
        InEdgeIt& operator++() { return *this; }
      };
      /// This iterator goes through each edge.

      /// This iterator goes through each edge of a graph.
      /// Its usage is quite simple, for example you can count the number
      /// of edges in a graph \c g of type \c Graph as follows:
      ///\code
      /// int count=0;
      /// for(Graph::EdgeIt e(g); e!=INVALID; ++e) ++count;
      ///\endcode
      class EdgeIt : public Edge {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        EdgeIt() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        EdgeIt(const EdgeIt& e) : Edge(e) { }
        /// Initialize the iterator to be invalid.

        /// Initialize the iterator to be invalid.
        ///
        EdgeIt(Invalid) { }
        /// This constructor sets the iterator to the first edge.
    
        /// This constructor sets the iterator to the first edge of \c g.
        ///@param g the graph
        EdgeIt(const Graph& g) { ignore_unused_variable_warning(g); }
        /// Edge -> EdgeIt conversion

        /// Sets the iterator to the value of the trivial iterator \c e.
        /// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        EdgeIt(const Graph&, const Edge&) { } 
        ///Next edge
        
        /// Assign the iterator to the next edge.
        EdgeIt& operator++() { return *this; }
      };
      ///Gives back the target node of an edge.

      ///Gives back the target node of an edge.
      ///
      Node target(Edge) const { return INVALID; }
      ///Gives back the source node of an edge.

      ///Gives back the source node of an edge.
      ///
      Node source(Edge) const { return INVALID; }

      void first(Node&) const {}
      void next(Node&) const {}

      void first(Edge&) const {}
      void next(Edge&) const {}


      void firstIn(Edge&, const Node&) const {}
      void nextIn(Edge&) const {}

      void firstOut(Edge&, const Node&) const {}
      void nextOut(Edge&) const {}

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

      /// \brief The opposite node on the given edge.
      ///
      /// Gives back the opposite node on the given edge.
      Node oppositeNode(const Node&, const Edge&) const { return INVALID; }

      /// \brief Read write map of the nodes to type \c T.
      /// 
      /// ReadWrite map of the nodes to type \c T.
      /// \sa Reference
      template<class T> 
      class NodeMap : public ReadWriteMap< Node, T > {
      public:

        ///\e
        NodeMap(const Graph&) { }
        ///\e
        NodeMap(const Graph&, T) { }

        ///Copy constructor
        NodeMap(const NodeMap& nm) : ReadWriteMap< Node, T >(nm) { }
        ///Assignment operator
        template <typename CMap>
        NodeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<Node, T>, CMap>();
          return *this; 
        }
      };

      /// \brief Read write map of the edges to type \c T.
      ///
      /// Reference map of the edges to type \c T.
      /// \sa Reference
      template<class T> 
      class EdgeMap : public ReadWriteMap<Edge,T> {
      public:

        ///\e
        EdgeMap(const Graph&) { }
        ///\e
        EdgeMap(const Graph&, T) { }
        ///Copy constructor
        EdgeMap(const EdgeMap& em) : ReadWriteMap<Edge,T>(em) { }
        ///Assignment operator
        template <typename CMap>
        EdgeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<Edge, T>, CMap>();
          return *this; 
        }
      };

      template <typename RGraph>
      struct Constraints {
        void constraints() {
          checkConcept<IterableGraphComponent<>, Graph>();
          checkConcept<MappableGraphComponent<>, Graph>();
        }
      };

    };
    
  } //namespace concepts  
} //namespace lemon



#endif // LEMON_CONCEPT_GRAPH_H

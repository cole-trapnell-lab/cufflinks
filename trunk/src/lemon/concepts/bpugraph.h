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

/// \ingroup graph_concepts
/// \file
/// \brief The concept of Bipartite Undirected Graphs.

#ifndef LEMON_CONCEPT_BPUGRAPH_H
#define LEMON_CONCEPT_BPUGRAPH_H

#include <lemon/concepts/graph_components.h>

#include <lemon/concepts/graph.h>
#include <lemon/concepts/ugraph.h>

#include <lemon/bits/utility.h>

namespace lemon {
  namespace concepts {

    /// \ingroup graph_concepts
    ///
    /// \brief Class describing the concept of Bipartite Undirected Graphs.
    ///
    /// This class describes the common interface of all 
    /// Undirected Bipartite Graphs.
    ///
    /// As all concept describing classes it provides only interface
    /// without any sensible implementation. So any algorithm for
    /// bipartite undirected graph should compile with this class, but it 
    /// will not run properly, of course.
    ///
    /// In LEMON bipartite undirected graphs also fulfill the concept of 
    /// the undirected graphs (\ref lemon::concepts::UGraph "UGraph Concept"). 
    ///
    /// You can assume that all undirected bipartite graph can be handled
    /// as an undirected graph and consequently as a static graph.
    ///
    /// The bipartite graph stores two types of nodes which are named
    /// ANode and BNode. The graph type contains two types ANode and
    /// BNode which are inherited from Node type. Moreover they have
    /// constructor which converts Node to either ANode or BNode when
    /// it is possible. Therefor everywhere the Node type can be used
    /// instead of ANode and BNode. So the usage of the ANode and
    /// BNode is not suggested.
    ///
    /// The iteration on the partition can be done with the ANodeIt and 
    /// BNodeIt classes. The node map can be used to map values to the nodes
    /// and similarly we can use to map values for just the ANodes and
    /// BNodes the ANodeMap and BNodeMap template classes.

    class BpUGraph {
    public:
      /// \brief The undirected graph should be tagged by the
      /// UndirectedTag.
      ///
      /// The undirected graph should be tagged by the UndirectedTag. This
      /// tag helps the enable_if technics to make compile time 
      /// specializations for undirected graphs.  
      typedef True UndirectedTag;

      /// \brief The base type of node iterators, 
      /// or in other words, the trivial node iterator.
      ///
      /// This is the base type of each node iterator,
      /// thus each kind of node iterator converts to this.
      /// More precisely each kind of node iterator should be inherited 
      /// from the trivial node iterator. The Node class represents
      /// both of two types of nodes. 
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

      /// \brief Helper class for ANodes.
      ///
      /// This class is just a helper class for ANodes, it is not
      /// suggested to use it directly. It can be converted easily to
      /// node and vice versa. The usage of this class is limited
      /// to use just as template parameters for special map types. 
      class ANode : public Node {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        ANode() : Node() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        ANode(const ANode&) : Node() { }

        /// Construct the same node as ANode.

        /// Construct the same node as ANode. It may throws assertion
        /// when the given node is from the BNode set.
        ANode(const Node&) : Node() { }

        /// Assign node to A-node.

        /// Besides the core graph item functionality each node should
        /// be convertible to the represented A-node if it is it possible. 
        ANode& operator=(const Node&) { return *this; }

        /// Invalid constructor \& conversion.

        /// This constructor initializes the iterator to be invalid.
        /// \sa Invalid for more details.
        ANode(Invalid) { }
        /// Equality operator

        /// Two iterators are equal if and only if they point to the
        /// same object or both are invalid.
        bool operator==(ANode) const { return true; }

        /// Inequality operator
        
        /// \sa operator==(ANode n)
        ///
        bool operator!=(ANode) const { return true; }

	/// Artificial ordering operator.
	
	/// To allow the use of graph descriptors as key type in std::map or
	/// similar associative container we require this.
	///
	/// \note This operator only have to define some strict ordering of
	/// the items; this order has nothing to do with the iteration
	/// ordering of the items.
	bool operator<(ANode) const { return false; }

      };

      /// \brief Helper class for BNodes.
      ///
      /// This class is just a helper class for BNodes, it is not
      /// suggested to use it directly. It can be converted easily to
      /// node and vice versa. The usage of this class is limited
      /// to use just as template parameters for special map types. 
      class BNode : public Node {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        BNode() : Node() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        BNode(const BNode&) : Node() { }

        /// Construct the same node as BNode.

        /// Construct the same node as BNode. It may throws assertion
        /// when the given node is from the ANode set.
        BNode(const Node&) : Node() { }

        /// Assign node to B-node.

        /// Besides the core graph item functionality each node should
        /// be convertible to the represented B-node if it is it possible. 
        BNode& operator=(const Node&) { return *this; }

        /// Invalid constructor \& conversion.

        /// This constructor initializes the iterator to be invalid.
        /// \sa Invalid for more details.
        BNode(Invalid) { }
        /// Equality operator

        /// Two iterators are equal if and only if they point to the
        /// same object or both are invalid.
        bool operator==(BNode) const { return true; }

        /// Inequality operator
        
        /// \sa operator==(BNode n)
        ///
        bool operator!=(BNode) const { return true; }

	/// Artificial ordering operator.
	
	/// To allow the use of graph descriptors as key type in std::map or
	/// similar associative container we require this.
	///
	/// \note This operator only have to define some strict ordering of
	/// the items; this order has nothing to do with the iteration
	/// ordering of the items.
	bool operator<(BNode) const { return false; }

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
        NodeIt(const BpUGraph&) { }
        /// Node -> NodeIt conversion.

        /// Sets the iterator to the node of \c the graph pointed by 
	/// the trivial iterator.
        /// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        NodeIt(const BpUGraph&, const Node&) { }
        /// Next node.

        /// Assign the iterator to the next node.
        ///
        NodeIt& operator++() { return *this; }
      };

      /// This iterator goes through each ANode.

      /// This iterator goes through each ANode.
      /// Its usage is quite simple, for example you can count the number
      /// of nodes in graph \c g of type \c Graph like this:
      ///\code
      /// int count=0;
      /// for (Graph::ANodeIt n(g); n!=INVALID; ++n) ++count;
      ///\endcode
      class ANodeIt : public Node {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        ANodeIt() { }
        /// Copy constructor.
        
        /// Copy constructor.
        ///
        ANodeIt(const ANodeIt& n) : Node(n) { }
        /// Invalid constructor \& conversion.

        /// Initialize the iterator to be invalid.
        /// \sa Invalid for more details.
        ANodeIt(Invalid) { }
        /// Sets the iterator to the first node.

        /// Sets the iterator to the first node of \c g.
        ///
        ANodeIt(const BpUGraph&) { }
        /// Node -> ANodeIt conversion.

        /// Sets the iterator to the node of \c the graph pointed by 
	/// the trivial iterator.
        /// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        ANodeIt(const BpUGraph&, const Node&) { }
        /// Next node.

        /// Assign the iterator to the next node.
        ///
        ANodeIt& operator++() { return *this; }
      };

      /// This iterator goes through each BNode.

      /// This iterator goes through each BNode.
      /// Its usage is quite simple, for example you can count the number
      /// of nodes in graph \c g of type \c Graph like this:
      ///\code
      /// int count=0;
      /// for (Graph::BNodeIt n(g); n!=INVALID; ++n) ++count;
      ///\endcode
      class BNodeIt : public Node {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        BNodeIt() { }
        /// Copy constructor.
        
        /// Copy constructor.
        ///
        BNodeIt(const BNodeIt& n) : Node(n) { }
        /// Invalid constructor \& conversion.

        /// Initialize the iterator to be invalid.
        /// \sa Invalid for more details.
        BNodeIt(Invalid) { }
        /// Sets the iterator to the first node.

        /// Sets the iterator to the first node of \c g.
        ///
        BNodeIt(const BpUGraph&) { }
        /// Node -> BNodeIt conversion.

        /// Sets the iterator to the node of \c the graph pointed by 
	/// the trivial iterator.
        /// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        BNodeIt(const BpUGraph&, const Node&) { }
        /// Next node.

        /// Assign the iterator to the next node.
        ///
        BNodeIt& operator++() { return *this; }
      };
    
    
      /// The base type of the undirected edge iterators.

      /// The base type of the undirected edge iterators.
      ///
      class UEdge {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        UEdge() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        UEdge(const UEdge&) { }
        /// Initialize the iterator to be invalid.

        /// Initialize the iterator to be invalid.
        ///
        UEdge(Invalid) { }
        /// Equality operator

        /// Two iterators are equal if and only if they point to the
        /// same object or both are invalid.
        bool operator==(UEdge) const { return true; }
        /// Inequality operator

        /// \sa operator==(UEdge n)
        ///
        bool operator!=(UEdge) const { return true; }

	/// Artificial ordering operator.
	
	/// To allow the use of graph descriptors as key type in std::map or
	/// similar associative container we require this.
	///
	/// \note This operator only have to define some strict ordering of
	/// the items; this order has nothing to do with the iteration
	/// ordering of the items.
	bool operator<(UEdge) const { return false; }
      };

      /// This iterator goes through each undirected edge.

      /// This iterator goes through each undirected edge of a graph.
      /// Its usage is quite simple, for example you can count the number
      /// of undirected edges in a graph \c g of type \c Graph as follows:
      ///\code
      /// int count=0;
      /// for(Graph::UEdgeIt e(g); e!=INVALID; ++e) ++count;
      ///\endcode
      class UEdgeIt : public UEdge {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        UEdgeIt() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        UEdgeIt(const UEdgeIt& e) : UEdge(e) { }
        /// Initialize the iterator to be invalid.

        /// Initialize the iterator to be invalid.
        ///
        UEdgeIt(Invalid) { }
        /// This constructor sets the iterator to the first undirected edge.
    
        /// This constructor sets the iterator to the first undirected edge.
        UEdgeIt(const BpUGraph&) { }
        /// UEdge -> UEdgeIt conversion

        /// Sets the iterator to the value of the trivial iterator.
        /// This feature necessitates that each time we
        /// iterate the undirected edge-set, the iteration order is the 
	/// same.
        UEdgeIt(const BpUGraph&, const UEdge&) { } 
        /// Next undirected edge
        
        /// Assign the iterator to the next undirected edge.
        UEdgeIt& operator++() { return *this; }
      };

      /// \brief This iterator goes trough the incident undirected 
      /// edges of a node.
      ///
      /// This iterator goes trough the incident undirected edges
      /// of a certain node
      /// of a graph.
      /// Its usage is quite simple, for example you can compute the
      /// degree (i.e. count the number
      /// of incident edges of a node \c n
      /// in graph \c g of type \c Graph as follows.
      ///\code
      /// int count=0;
      /// for(Graph::IncEdgeIt e(g, n); e!=INVALID; ++e) ++count;
      ///\endcode
      class IncEdgeIt : public UEdge {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        IncEdgeIt() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        IncEdgeIt(const IncEdgeIt& e) : UEdge(e) { }
        /// Initialize the iterator to be invalid.

        /// Initialize the iterator to be invalid.
        ///
        IncEdgeIt(Invalid) { }
        /// This constructor sets the iterator to first incident edge.
    
        /// This constructor set the iterator to the first incident edge of
        /// the node.
        IncEdgeIt(const BpUGraph&, const Node&) { }
        /// UEdge -> IncEdgeIt conversion

        /// Sets the iterator to the value of the trivial iterator \c e.
        /// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        IncEdgeIt(const BpUGraph&, const UEdge&) { }
        /// Next incident edge

        /// Assign the iterator to the next incident edge
	/// of the corresponding node.
        IncEdgeIt& operator++() { return *this; }
      };

      /// The directed edge type.

      /// The directed edge type. It can be converted to the
      /// undirected edge.
      class Edge : public UEdge {
      public:
        /// Default constructor

        /// @warning The default constructor sets the iterator
        /// to an undefined value.
        Edge() { }
        /// Copy constructor.

        /// Copy constructor.
        ///
        Edge(const Edge& e) : UEdge(e) { }
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
      /// This iterator goes through each directed edge.

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
        EdgeIt(const BpUGraph &g) { ignore_unused_variable_warning(g); }
        /// Edge -> EdgeIt conversion

        /// Sets the iterator to the value of the trivial iterator \c e.
        /// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        EdgeIt(const BpUGraph&, const Edge&) { } 
        ///Next edge
        
        /// Assign the iterator to the next edge.
        EdgeIt& operator++() { return *this; }
      };
   
      /// This iterator goes trough the outgoing directed edges of a node.

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
        ///@param n the node
        ///@param g the graph
        OutEdgeIt(const BpUGraph& n, const Node& g) {
	  ignore_unused_variable_warning(n);
	  ignore_unused_variable_warning(g);
	}
        /// Edge -> OutEdgeIt conversion

        /// Sets the iterator to the value of the trivial iterator.
	/// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        OutEdgeIt(const BpUGraph&, const Edge&) { }
        ///Next outgoing edge
        
        /// Assign the iterator to the next 
        /// outgoing edge of the corresponding node.
        OutEdgeIt& operator++() { return *this; }
      };

      /// This iterator goes trough the incoming directed edges of a node.

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
        ///@param n the node
        ///@param g the graph
        InEdgeIt(const BpUGraph& g, const Node& n) { 
	  ignore_unused_variable_warning(n);
	  ignore_unused_variable_warning(g);
	}
        /// Edge -> InEdgeIt conversion

        /// Sets the iterator to the value of the trivial iterator \c e.
        /// This feature necessitates that each time we 
        /// iterate the edge-set, the iteration order is the same.
        InEdgeIt(const BpUGraph&, const Edge&) { }
        /// Next incoming edge

        /// Assign the iterator to the next inedge of the corresponding node.
        ///
        InEdgeIt& operator++() { return *this; }
      };

      /// \brief Read write map of the nodes to type \c T.
      /// 
      /// ReadWrite map of the nodes to type \c T.
      /// \sa Reference
      /// \todo Wrong documentation
      template<class T> 
      class NodeMap : public ReadWriteMap< Node, T >
      {
      public:

        ///\e
        NodeMap(const BpUGraph&) { }
        ///\e
        NodeMap(const BpUGraph&, T) { }

        ///Copy constructor
        NodeMap(const NodeMap& nm) : ReadWriteMap< Node, T >(nm) { }
        ///Assignment operator
        NodeMap& operator=(const NodeMap&) { return *this; }
        ///Assignment operator
        template <typename CMap>
        NodeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<Node, T>, CMap>();
          return *this; 
        }
      };

      /// \brief Read write map of the ANodes to type \c T.
      /// 
      /// ReadWrite map of the ANodes to type \c T.
      /// \sa Reference
      /// \todo Wrong documentation
      template<class T> 
      class ANodeMap : public ReadWriteMap< Node, T >
      {
      public:

        ///\e
        ANodeMap(const BpUGraph&) { }
        ///\e
        ANodeMap(const BpUGraph&, T) { }

        ///Copy constructor
        ANodeMap(const ANodeMap& nm) : ReadWriteMap< Node, T >(nm) { }
        ///Assignment operator
        ANodeMap& operator=(const ANodeMap&) { return *this; }
        ///Assignment operator
        template <typename CMap>
        ANodeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<Node, T>, CMap>();
          return *this; 
        }
      };

      /// \brief Read write map of the BNodes to type \c T.
      /// 
      /// ReadWrite map of the BNodes to type \c T.
      /// \sa Reference
      /// \todo Wrong documentation
      template<class T> 
      class BNodeMap : public ReadWriteMap< Node, T >
      {
      public:

        ///\e
        BNodeMap(const BpUGraph&) { }
        ///\e
        BNodeMap(const BpUGraph&, T) { }

        ///Copy constructor
        BNodeMap(const BNodeMap& nm) : ReadWriteMap< Node, T >(nm) { }
        ///Assignment operator
        BNodeMap& operator=(const BNodeMap&) { return *this; }
        ///Assignment operator
        template <typename CMap>
        BNodeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<Node, T>, CMap>();
          return *this; 
        }
      };

      /// \brief Read write map of the directed edges to type \c T.
      ///
      /// Reference map of the directed edges to type \c T.
      /// \sa Reference
      /// \todo Wrong documentation
      template<class T> 
      class EdgeMap : public ReadWriteMap<Edge,T>
      {
      public:

        ///\e
        EdgeMap(const BpUGraph&) { }
        ///\e
        EdgeMap(const BpUGraph&, T) { }
        ///Copy constructor
        EdgeMap(const EdgeMap& em) : ReadWriteMap<Edge,T>(em) { }
        ///Assignment operator
        EdgeMap& operator=(const EdgeMap&) { return *this; }
        ///Assignment operator
        template <typename CMap>
        EdgeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<Edge, T>, CMap>();
          return *this; 
        }
      };

      /// Read write map of the undirected edges to type \c T.

      /// Reference map of the edges to type \c T.
      /// \sa Reference
      /// \todo Wrong documentation
      template<class T> 
      class UEdgeMap : public ReadWriteMap<UEdge,T>
      {
      public:

        ///\e
        UEdgeMap(const BpUGraph&) { }
        ///\e
        UEdgeMap(const BpUGraph&, T) { }
        ///Copy constructor
        UEdgeMap(const UEdgeMap& em) : ReadWriteMap<UEdge,T>(em) {}
        ///Assignment operator
        UEdgeMap &operator=(const UEdgeMap&) { return *this; }
        ///Assignment operator
        template <typename CMap>
        UEdgeMap& operator=(const CMap&) { 
          checkConcept<ReadMap<UEdge, T>, CMap>();
          return *this; 
        }
      };

      /// \brief Direct the given undirected edge.
      ///
      /// Direct the given undirected edge. The returned edge source
      /// will be the given node.
      Edge direct(const UEdge&, const Node&) const {
	return INVALID;
      }

      /// \brief Direct the given undirected edge.
      ///
      /// Direct the given undirected edge. The returned edge
      /// represents the given undirected edge and the direction comes
      /// from the given bool.  The source of the undirected edge and
      /// the directed edge is the same when the given bool is true.
      Edge direct(const UEdge&, bool) const {
	return INVALID;
      }

      /// \brief Returns true when the given node is an ANode.
      ///
      /// Returns true when the given node is an ANode.
      bool aNode(Node) const { return true;}

      /// \brief Returns true when the given node is an BNode.
      ///
      /// Returns true when the given node is an BNode.
      bool bNode(Node) const { return true;}

      /// \brief Returns the edge's end node which is in the ANode set.
      ///
      /// Returns the edge's end node which is in the ANode set.
      Node aNode(UEdge) const { return INVALID;}

      /// \brief Returns the edge's end node which is in the BNode set.
      ///
      /// Returns the edge's end node which is in the BNode set.
      Node bNode(UEdge) const { return INVALID;}

      /// \brief Returns true if the edge has default orientation.
      ///
      /// Returns whether the given directed edge is same orientation as
      /// the corresponding undirected edge's default orientation.
      bool direction(Edge) const { return true; }

      /// \brief Returns the opposite directed edge.
      ///
      /// Returns the opposite directed edge.
      Edge oppositeEdge(Edge) const { return INVALID; }

      /// \brief Opposite node on an edge
      ///
      /// \return the opposite of the given Node on the given UEdge
      Node oppositeNode(Node, UEdge) const { return INVALID; }

      /// \brief First node of the undirected edge.
      ///
      /// \return the first node of the given UEdge.
      ///
      /// Naturally undirected edges don't have direction and thus
      /// don't have source and target node. But we use these two methods
      /// to query the two endnodes of the edge. The direction of the edge
      /// which arises this way is called the inherent direction of the
      /// undirected edge, and is used to define the "default" direction
      /// of the directed versions of the edges.
      /// \sa direction
      Node source(UEdge) const { return INVALID; }

      /// \brief Second node of the undirected edge.
      Node target(UEdge) const { return INVALID; }

      /// \brief Source node of the directed edge.
      Node source(Edge) const { return INVALID; }

      /// \brief Target node of the directed edge.
      Node target(Edge) const { return INVALID; }

      /// \brief Base node of the iterator
      ///
      /// Returns the base node (the source in this case) of the iterator
      Node baseNode(OutEdgeIt e) const {
	return source(e);
      }

      /// \brief Running node of the iterator
      ///
      /// Returns the running node (the target in this case) of the
      /// iterator
      Node runningNode(OutEdgeIt e) const {
	return target(e);
      }

      /// \brief Base node of the iterator
      ///
      /// Returns the base node (the target in this case) of the iterator
      Node baseNode(InEdgeIt e) const {
	return target(e);
      }
      /// \brief Running node of the iterator
      ///
      /// Returns the running node (the source in this case) of the
      /// iterator
      Node runningNode(InEdgeIt e) const {
	return source(e);
      }

      /// \brief Base node of the iterator
      ///
      /// Returns the base node of the iterator
      Node baseNode(IncEdgeIt) const {
	return INVALID;
      }
      
      /// \brief Running node of the iterator
      ///
      /// Returns the running node of the iterator
      Node runningNode(IncEdgeIt) const {
	return INVALID;
      }

      void first(Node&) const {}
      void next(Node&) const {}

      void first(Edge&) const {}
      void next(Edge&) const {}

      void first(UEdge&) const {}
      void next(UEdge&) const {}

      void firstANode(Node&) const {}
      void nextANode(Node&) const {}

      void firstBNode(Node&) const {}
      void nextBNode(Node&) const {}

      void firstIn(Edge&, const Node&) const {}
      void nextIn(Edge&) const {}

      void firstOut(Edge&, const Node&) const {}
      void nextOut(Edge&) const {}

      void firstInc(UEdge &, bool &, const Node &) const {}
      void nextInc(UEdge &, bool &) const {}

      void firstFromANode(UEdge&, const Node&) const {}
      void nextFromANode(UEdge&) const {}

      void firstFromBNode(UEdge&, const Node&) const {}
      void nextFromBNode(UEdge&) const {}

      template <typename Graph>
      struct Constraints {
	void constraints() {
	  checkConcept<IterableBpUGraphComponent<>, Graph>();
	  checkConcept<MappableBpUGraphComponent<>, Graph>();
	}
      };

    };

  }

}

#endif

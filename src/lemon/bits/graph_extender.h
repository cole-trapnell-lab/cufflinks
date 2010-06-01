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

#ifndef LEMON_BITS_GRAPH_EXTENDER_H
#define LEMON_BITS_GRAPH_EXTENDER_H

#include <lemon/bits/invalid.h>
#include <lemon/bits/utility.h>
#include <lemon/error.h>

#include <lemon/bits/map_extender.h>
#include <lemon/bits/default_map.h>

#include <lemon/concept_check.h>
#include <lemon/concepts/maps.h>

///\ingroup graphbits
///\file
///\brief Extenders for the graph types
namespace lemon {
	
	/// \ingroup graphbits
	///
	/// \brief Extender for the Graphs
	template <typename Base>
	class GraphExtender : public Base {
	public:
		
		typedef Base Parent;
		typedef GraphExtender Graph;
		
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
		
		// Alterable extension
		
		typedef AlterationNotifier<GraphExtender, Node> NodeNotifier;
		typedef AlterationNotifier<GraphExtender, Edge> EdgeNotifier;
		
		
	protected:
		
		mutable NodeNotifier node_notifier;
		mutable EdgeNotifier edge_notifier;
		
	public:
		
		NodeNotifier& notifier(Node) const {
			return node_notifier;
		}
		
		EdgeNotifier& notifier(Edge) const {
			return edge_notifier;
		}
		
		class NodeIt : public Node { 
			const Graph* graph;
		public:
			
			NodeIt() {}
			
			NodeIt(Invalid i) : Node(i) { }
			
			explicit NodeIt(const Graph& _graph) : graph(&_graph) {
				_graph.first(static_cast<Node&>(*this));
			}
			
			NodeIt(const Graph& _graph, const Node& node) 
			: Node(node), graph(&_graph) {}
			
			NodeIt& operator++() { 
				graph->next(*this);
				return *this; 
			}
			
		};
		
		
		class EdgeIt : public Edge { 
			const Graph* graph;
		public:
			
			EdgeIt() { }
			
			EdgeIt(Invalid i) : Edge(i) { }
			
			explicit EdgeIt(const Graph& _graph) : graph(&_graph) {
				_graph.first(static_cast<Edge&>(*this));
			}
			
			EdgeIt(const Graph& _graph, const Edge& e) : 
			Edge(e), graph(&_graph) { }
			
			EdgeIt& operator++() { 
				graph->next(*this);
				return *this; 
			}
			
		};
		
		
		class OutEdgeIt : public Edge { 
			const Graph* graph;
		public:
			
			OutEdgeIt() { }
			
			OutEdgeIt(Invalid i) : Edge(i) { }
			
			OutEdgeIt(const Graph& _graph, const Node& node) 
			: graph(&_graph) {
				_graph.firstOut(*this, node);
			}
			
			OutEdgeIt(const Graph& _graph, const Edge& edge) 
			: Edge(edge), graph(&_graph) {}
			
			OutEdgeIt& operator++() { 
				graph->nextOut(*this);
				return *this; 
			}
			
		};
		
		
		class InEdgeIt : public Edge { 
			const Graph* graph;
		public:
			
			InEdgeIt() { }
			
			InEdgeIt(Invalid i) : Edge(i) { }
			
			InEdgeIt(const Graph& _graph, const Node& node) 
			: graph(&_graph) {
				_graph.firstIn(*this, node);
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
		/// Returns the base node (i.e. the source in this case) of the iterator
		Node baseNode(const OutEdgeIt &e) const {
			return Parent::source(e);
		}
		/// \brief Running node of the iterator
		///
		/// Returns the running node (i.e. the target in this case) of the
		/// iterator
		Node runningNode(const OutEdgeIt &e) const {
			return Parent::target(e);
		}
		
		/// \brief Base node of the iterator
		///
		/// Returns the base node (i.e. the target in this case) of the iterator
		Node baseNode(const InEdgeIt &e) const {
			return Parent::target(e);
		}
		/// \brief Running node of the iterator
		///
		/// Returns the running node (i.e. the source in this case) of the
		/// iterator
		Node runningNode(const InEdgeIt &e) const {
			return Parent::source(e);
		}
		
		
		template <typename _Value>
		class NodeMap 
		: public MapExtender<DefaultMap<Graph, Node, _Value> > {
		public:
			typedef GraphExtender Graph;
			typedef MapExtender<DefaultMap<Graph, Node, _Value> > Parent;
			
			explicit NodeMap(const Graph& graph) 
			: Parent(graph) {}
			NodeMap(const Graph& graph, const _Value& value) 
			: Parent(graph, value) {}
			
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
		: public MapExtender<DefaultMap<Graph, Edge, _Value> > {
		public:
			typedef GraphExtender Graph;
			typedef MapExtender<DefaultMap<Graph, Edge, _Value> > Parent;
			
			explicit EdgeMap(const Graph& graph) 
			: Parent(graph) {}
			EdgeMap(const Graph& graph, const _Value& value) 
			: Parent(graph, value) {}
			
			EdgeMap& operator=(const EdgeMap& cmap) {
				return operator=<EdgeMap>(cmap);
			}
			
			template <typename CMap>
			EdgeMap& operator=(const CMap& cmap) {
				Parent::operator=(cmap);
				return *this;
			}
		};
		
		
		Node addNode() {
			Node node = Parent::addNode();
			notifier(Node()).add(node);
			return node;
		}
		
		Edge addEdge(const Node& from, const Node& to) {
			Edge edge = Parent::addEdge(from, to);
			notifier(Edge()).add(edge);
			return edge;
		}
		
		void clear() {
			notifier(Edge()).clear();
			notifier(Node()).clear();
			Parent::clear();
		}
		
		template <typename Graph, typename NodeRefMap, typename EdgeRefMap>
		void build(const Graph& graph, NodeRefMap& nodeRef, EdgeRefMap& edgeRef) {
			Parent::build(graph, nodeRef, edgeRef);
			notifier(Node()).build();
			notifier(Edge()).build();
		}
		
		void erase(const Node& node) {
			Edge edge;
			Parent::firstOut(edge, node);
			while (edge != INVALID ) {
				erase(edge);
				Parent::firstOut(edge, node);
			} 
			
			Parent::firstIn(edge, node);
			while (edge != INVALID ) {
				erase(edge);
				Parent::firstIn(edge, node);
			}
			
			notifier(Node()).erase(node);
			Parent::erase(node);
		}
		
		void erase(const Edge& edge) {
			notifier(Edge()).erase(edge);
			Parent::erase(edge);
		}
		
		GraphExtender() {
			node_notifier.setContainer(*this);
			edge_notifier.setContainer(*this);
		} 
		
		
		~GraphExtender() {
			edge_notifier.clear();
			node_notifier.clear();
		}
	};
	
	/// \ingroup graphbits
	///
	/// \brief Extender for the UGraphs
	template <typename Base> 
	class UGraphExtender : public Base {
	public:
		
		typedef Base Parent;
		typedef UGraphExtender Graph;
		
		typedef True UndirectedTag;
		
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
		
		// Alterable extension
		
		typedef AlterationNotifier<UGraphExtender, Node> NodeNotifier;
		typedef AlterationNotifier<UGraphExtender, Edge> EdgeNotifier;
		typedef AlterationNotifier<UGraphExtender, UEdge> UEdgeNotifier;
		
		
	protected:
		
		mutable NodeNotifier node_notifier;
		mutable EdgeNotifier edge_notifier;
		mutable UEdgeNotifier uedge_notifier;
		
	public:
		
		NodeNotifier& notifier(Node) const {
			return node_notifier;
		}
		
		EdgeNotifier& notifier(Edge) const {
			return edge_notifier;
		}
		
		UEdgeNotifier& notifier(UEdge) const {
			return uedge_notifier;
		}
		
		
		
		class NodeIt : public Node { 
			const Graph* graph;
		public:
			
			NodeIt() {}
			
			NodeIt(Invalid i) : Node(i) { }
			
			explicit NodeIt(const Graph& _graph) : graph(&_graph) {
				_graph.first(static_cast<Node&>(*this));
			}
			
			NodeIt(const Graph& _graph, const Node& node) 
			: Node(node), graph(&_graph) {}
			
			NodeIt& operator++() { 
				graph->next(*this);
				return *this; 
			}
			
		};
		
		
		class EdgeIt : public Edge { 
			const Graph* graph;
		public:
			
			EdgeIt() { }
			
			EdgeIt(Invalid i) : Edge(i) { }
			
			explicit EdgeIt(const Graph& _graph) : graph(&_graph) {
				_graph.first(static_cast<Edge&>(*this));
			}
			
			EdgeIt(const Graph& _graph, const Edge& e) : 
			Edge(e), graph(&_graph) { }
			
			EdgeIt& operator++() { 
				graph->next(*this);
				return *this; 
			}
			
		};
		
		
		class OutEdgeIt : public Edge { 
			const Graph* graph;
		public:
			
			OutEdgeIt() { }
			
			OutEdgeIt(Invalid i) : Edge(i) { }
			
			OutEdgeIt(const Graph& _graph, const Node& node) 
			: graph(&_graph) {
				_graph.firstOut(*this, node);
			}
			
			OutEdgeIt(const Graph& _graph, const Edge& edge) 
			: Edge(edge), graph(&_graph) {}
			
			OutEdgeIt& operator++() { 
				graph->nextOut(*this);
				return *this; 
			}
			
		};
		
		
		class InEdgeIt : public Edge { 
			const Graph* graph;
		public:
			
			InEdgeIt() { }
			
			InEdgeIt(Invalid i) : Edge(i) { }
			
			InEdgeIt(const Graph& _graph, const Node& node) 
			: graph(&_graph) {
				_graph.firstIn(*this, node);
			}
			
			InEdgeIt(const Graph& _graph, const Edge& edge) : 
			Edge(edge), graph(&_graph) {}
			
			InEdgeIt& operator++() { 
				graph->nextIn(*this);
				return *this; 
			}
			
		};
		
		
		class UEdgeIt : public Parent::UEdge { 
			const Graph* graph;
		public:
			
			UEdgeIt() { }
			
			UEdgeIt(Invalid i) : UEdge(i) { }
			
			explicit UEdgeIt(const Graph& _graph) : graph(&_graph) {
				_graph.first(static_cast<UEdge&>(*this));
			}
			
			UEdgeIt(const Graph& _graph, const UEdge& e) : 
			UEdge(e), graph(&_graph) { }
			
			UEdgeIt& operator++() { 
				graph->next(*this);
				return *this; 
			}
			
		};
		
		class IncEdgeIt : public Parent::UEdge {
			friend class UGraphExtender;
			const Graph* graph;
			bool direction;
		public:
			
			IncEdgeIt() { }
			
			IncEdgeIt(Invalid i) : UEdge(i), direction(false) { }
			
			IncEdgeIt(const Graph& _graph, const Node &n) : graph(&_graph) {
				_graph.firstInc(*this, direction, n);
			}
			
			IncEdgeIt(const Graph& _graph, const UEdge &ue, const Node &n)
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
		
		// Mappable extension
		
		template <typename _Value>
		class NodeMap 
		: public MapExtender<DefaultMap<Graph, Node, _Value> > {
		public:
			typedef UGraphExtender Graph;
			typedef MapExtender<DefaultMap<Graph, Node, _Value> > Parent;
			
			NodeMap(const Graph& graph) 
			: Parent(graph) {}
			NodeMap(const Graph& graph, const _Value& value) 
			: Parent(graph, value) {}
			
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
		: public MapExtender<DefaultMap<Graph, Edge, _Value> > {
		public:
			typedef UGraphExtender Graph;
			typedef MapExtender<DefaultMap<Graph, Edge, _Value> > Parent;
			
			EdgeMap(const Graph& graph) 
			: Parent(graph) {}
			EdgeMap(const Graph& graph, const _Value& value) 
			: Parent(graph, value) {}
			
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
		class UEdgeMap 
		: public MapExtender<DefaultMap<Graph, UEdge, _Value> > {
		public:
			typedef UGraphExtender Graph;
			typedef MapExtender<DefaultMap<Graph, UEdge, _Value> > Parent;
			
			UEdgeMap(const Graph& graph) 
			: Parent(graph) {}
			
			UEdgeMap(const Graph& graph, const _Value& value) 
			: Parent(graph, value) {}
			
			UEdgeMap& operator=(const UEdgeMap& cmap) {
				return operator=<UEdgeMap>(cmap);
			}
			
			template <typename CMap>
			UEdgeMap& operator=(const CMap& cmap) {
				Parent::operator=(cmap);
				return *this;
			}
			
		};
		
		// Alteration extension
		
		Node addNode() {
			Node node = Parent::addNode();
			notifier(Node()).add(node);
			return node;
		}
		
		UEdge addEdge(const Node& from, const Node& to) {
			UEdge uedge = Parent::addEdge(from, to);
			notifier(UEdge()).add(uedge);
			std::vector<Edge> ev;
			ev.push_back(Parent::direct(uedge, true));
			ev.push_back(Parent::direct(uedge, false));      
			notifier(Edge()).add(ev);
			return uedge;
		}
		
		void clear() {
			notifier(Edge()).clear();
			notifier(UEdge()).clear();
			notifier(Node()).clear();
			Parent::clear();
		}
		
		template <typename Graph, typename NodeRefMap, typename UEdgeRefMap>
		void build(const Graph& graph, NodeRefMap& nodeRef, 
				   UEdgeRefMap& uEdgeRef) {
			Parent::build(graph, nodeRef, uEdgeRef);
			notifier(Node()).build();
			notifier(UEdge()).build();
			notifier(Edge()).build();
		}
		
		void erase(const Node& node) {
			Edge edge;
			Parent::firstOut(edge, node);
			while (edge != INVALID ) {
				erase(edge);
				Parent::firstOut(edge, node);
			} 
			
			Parent::firstIn(edge, node);
			while (edge != INVALID ) {
				erase(edge);
				Parent::firstIn(edge, node);
			}
			
			notifier(Node()).erase(node);
			Parent::erase(node);
		}
		
		void erase(const UEdge& uedge) {
			std::vector<Edge> ev;
			ev.push_back(Parent::direct(uedge, true));
			ev.push_back(Parent::direct(uedge, false));      
			notifier(Edge()).erase(ev);
			notifier(UEdge()).erase(uedge);
			Parent::erase(uedge);
		}
		
		UGraphExtender() {
			node_notifier.setContainer(*this); 
			edge_notifier.setContainer(*this);
			uedge_notifier.setContainer(*this);
		} 
		
		~UGraphExtender() {
			uedge_notifier.clear();
			edge_notifier.clear();
			node_notifier.clear(); 
		} 
		
	};
	
	/// \ingroup graphbits
	///
	/// \brief Extender for the BpUGraphs
	template <typename Base>
	class BpUGraphExtender : public Base {
	public:
		
		typedef Base Parent;
		typedef BpUGraphExtender Graph;
		
		typedef True UndirectedTag;
		
		typedef typename Parent::Node Node;
		typedef typename Parent::ANode ANode;
		typedef typename Parent::BNode BNode;
		typedef typename Parent::Edge Edge;
		typedef typename Parent::UEdge UEdge;
		
		
		Node oppositeNode(const Node& node, const UEdge& edge) const {
			return Parent::aNode(edge) == node ? 
			Parent::bNode(edge) : Parent::aNode(edge);
		}
		
		using Parent::direct;
		Edge direct(const UEdge& edge, const Node& node) const {
			return Parent::direct(edge, node == Parent::source(edge));
		}
		
		Edge oppositeEdge(const Edge& edge) const {
			return direct(edge, !Parent::direction(edge));
		}
		
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
		
		typedef AlterationNotifier<BpUGraphExtender, ANode> ANodeNotifier;
		typedef AlterationNotifier<BpUGraphExtender, BNode> BNodeNotifier;
		typedef AlterationNotifier<BpUGraphExtender, Node> NodeNotifier;
		typedef AlterationNotifier<BpUGraphExtender, Edge> EdgeNotifier;
		typedef AlterationNotifier<BpUGraphExtender, UEdge> UEdgeNotifier;
		
	protected:
		
		mutable ANodeNotifier anode_notifier;
		mutable BNodeNotifier bnode_notifier;
		mutable NodeNotifier node_notifier;
		mutable EdgeNotifier edge_notifier;
		mutable UEdgeNotifier uedge_notifier;
		
	public:
		
		NodeNotifier& notifier(Node) const {
			return node_notifier;
		}
		
		ANodeNotifier& notifier(ANode) const {
			return anode_notifier;
		}
		
		BNodeNotifier& notifier(BNode) const {
			return bnode_notifier;
		}
		
		EdgeNotifier& notifier(Edge) const {
			return edge_notifier;
		}
		
		UEdgeNotifier& notifier(UEdge) const {
			return uedge_notifier;
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
			friend class BpUGraphExtender;
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
			friend class BpUGraphExtender;
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
			friend class BpUGraphExtender;
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
			friend class BpUGraphExtender;
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
			friend class BpUGraphExtender;
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
			friend class BpUGraphExtender;
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
			friend class BpUGraphExtender;
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
		
		template <typename _Value>
		class ANodeMap 
		: public MapExtender<DefaultMap<Graph, ANode, _Value> > {
		public:
			typedef BpUGraphExtender Graph;
			typedef MapExtender<DefaultMap<Graph, ANode, _Value> > Parent;
			
			ANodeMap(const Graph& graph) 
			: Parent(graph) {}
			ANodeMap(const Graph& graph, const _Value& value) 
			: Parent(graph, value) {}
			
			ANodeMap& operator=(const ANodeMap& cmap) {
				return operator=<ANodeMap>(cmap);
			}
			
			template <typename CMap>
			ANodeMap& operator=(const CMap& cmap) {
				Parent::operator=(cmap);
				return *this;
			}
			
		};
		
		template <typename _Value>
		class BNodeMap 
		: public MapExtender<DefaultMap<Graph, BNode, _Value> > {
		public:
			typedef BpUGraphExtender Graph;
			typedef MapExtender<DefaultMap<Graph, BNode, _Value> > Parent;
			
			BNodeMap(const Graph& graph) 
			: Parent(graph) {}
			BNodeMap(const Graph& graph, const _Value& value) 
			: Parent(graph, value) {}
			
			BNodeMap& operator=(const BNodeMap& cmap) {
				return operator=<BNodeMap>(cmap);
			}
			
			template <typename CMap>
			BNodeMap& operator=(const CMap& cmap) {
				Parent::operator=(cmap);
				return *this;
			}
			
		};
		
	public:
		
		template <typename _Value>
		class NodeMap {
		public:
			typedef BpUGraphExtender Graph;
			
			typedef Node Key;
			typedef _Value Value;
			
			/// The reference type of the map;
			typedef typename ANodeMap<_Value>::Reference Reference;
			/// The const reference type of the map;
			typedef typename ANodeMap<_Value>::ConstReference ConstReference;
			
			typedef True ReferenceMapTag;
			
			NodeMap(const Graph& _graph) 
			: graph(_graph), aNodeMap(_graph), bNodeMap(_graph) {}
			NodeMap(const Graph& _graph, const _Value& _value) 
			: graph(_graph), aNodeMap(_graph, _value), bNodeMap(_graph, _value) {}
			
			NodeMap& operator=(const NodeMap& cmap) {
				return operator=<NodeMap>(cmap);
			}
			
			template <typename CMap>
			NodeMap& operator=(const CMap& cmap) {
				checkConcept<concepts::ReadMap<Node, _Value>, CMap>();
				aNodeMap = cmap;
				bNodeMap = cmap;
				return *this;
			}
			
			ConstReference operator[](const Key& node) const {
				if (Parent::aNode(node)) {
					return aNodeMap[node];
				} else {
					return bNodeMap[node];
				}
			} 
			
			Reference operator[](const Key& node) {
				if (Parent::aNode(node)) {
					return aNodeMap[node];
				} else {
					return bNodeMap[node];
				}
			}
			
			void set(const Key& node, const Value& value) {
				if (Parent::aNode(node)) {
					aNodeMap.set(node, value);
				} else {
					bNodeMap.set(node, value);
				}
			}
			
			class MapIt : public NodeIt {
			public:
				
				typedef NodeIt Parent;
				
				explicit MapIt(NodeMap& _map) 
				: Parent(_map.graph), map(_map) {}
				
				typename MapTraits<NodeMap>::ConstReturnValue operator*() const {
					return map[*this];
				}
				
				typename MapTraits<NodeMap>::ReturnValue operator*() {
					return map[*this];
				}
				
				void set(const Value& value) {
					map.set(*this, value);
				}
				
			private:
				NodeMap& map;
			};
			
			class ConstMapIt : public NodeIt {
			public:
				
				typedef NodeIt Parent;
				
				explicit ConstMapIt(const NodeMap& _map) 
				: Parent(_map.graph), map(_map) {}
				
				typename MapTraits<NodeMap>::ConstReturnValue operator*() const {
					return map[*this];
				}
				
			private:
				const NodeMap& map;
			};
			
			class ItemIt : public NodeIt {
			public:
				
				typedef NodeIt Parent;
				
				explicit ItemIt(const NodeMap& _map)
				: Parent(_map.graph) {}
				
			};
			
		private:
			const Graph& graph;
			ANodeMap<_Value> aNodeMap;
			BNodeMap<_Value> bNodeMap;
		};
		
		
		template <typename _Value>
		class EdgeMap 
		: public MapExtender<DefaultMap<Graph, Edge, _Value> > {
		public:
			typedef BpUGraphExtender Graph;
			typedef MapExtender<DefaultMap<Graph, Edge, _Value> > Parent;
			
			EdgeMap(const Graph& graph) 
			: Parent(graph) {}
			EdgeMap(const Graph& graph, const _Value& value) 
			: Parent(graph, value) {}
			
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
		class UEdgeMap 
		: public MapExtender<DefaultMap<Graph, UEdge, _Value> > {
		public:
			typedef BpUGraphExtender Graph;
			typedef MapExtender<DefaultMap<Graph, UEdge, _Value> > Parent;
			
			UEdgeMap(const Graph& graph) 
			: Parent(graph) {}
			UEdgeMap(const Graph& graph, const _Value& value) 
			: Parent(graph, value) {}
			
			UEdgeMap& operator=(const UEdgeMap& cmap) {
				return operator=<UEdgeMap>(cmap);
			}
			
			template <typename CMap>
			UEdgeMap& operator=(const CMap& cmap) {
				Parent::operator=(cmap);
				return *this;
			}
		};
		
		
		Node addANode() {
			Node node = Parent::addANode();
			notifier(ANode()).add(node);
			notifier(Node()).add(node);
			return node;
		}
		
		Node addBNode() {
			Node node = Parent::addBNode();
			notifier(BNode()).add(node);
			notifier(Node()).add(node);
			return node;
		}
		
		UEdge addEdge(const Node& s, const Node& t) {
			UEdge uedge = Parent::addEdge(s, t);
			notifier(UEdge()).add(uedge);
			
			std::vector<Edge> ev;
			ev.push_back(Parent::direct(uedge, true));
			ev.push_back(Parent::direct(uedge, false));
			notifier(Edge()).add(ev);
			
			return uedge;
		}
		
		void clear() {
			notifier(Edge()).clear();
			notifier(UEdge()).clear();
			notifier(Node()).clear();
			notifier(BNode()).clear();
			notifier(ANode()).clear();
			Parent::clear();
		}
		
		template <typename Graph, typename ANodeRefMap, 
		typename BNodeRefMap, typename UEdgeRefMap>
		void build(const Graph& graph, ANodeRefMap& aNodeRef, 
				   BNodeRefMap& bNodeRef, UEdgeRefMap& uEdgeRef) {
			Parent::build(graph, aNodeRef, bNodeRef, uEdgeRef);
			notifier(ANode()).build();
			notifier(BNode()).build();
			notifier(Node()).build();
			notifier(UEdge()).build();
			notifier(Edge()).build();
		}
		
		void erase(const Node& node) {
			UEdge uedge;
			if (Parent::aNode(node)) {
				Parent::firstFromANode(uedge, node);
				while (uedge != INVALID) {
					erase(uedge);
					Parent::firstFromANode(uedge, node);
				}
				notifier(ANode()).erase(node);
			} else {
				Parent::firstFromBNode(uedge, node);
				while (uedge != INVALID) {
					erase(uedge);
					Parent::firstFromBNode(uedge, node);
				}
				notifier(BNode()).erase(node);
			}
			
			notifier(Node()).erase(node);
			Parent::erase(node);
		}
		
		void erase(const UEdge& uedge) {
			std::vector<Edge> ev;
			ev.push_back(Parent::direct(uedge, true));
			ev.push_back(Parent::direct(uedge, false));
			notifier(Edge()).erase(ev);
			notifier(UEdge()).erase(uedge);
			Parent::erase(uedge);
		}
		
		
		BpUGraphExtender() {
			anode_notifier.setContainer(*this); 
			bnode_notifier.setContainer(*this); 
			node_notifier.setContainer(*this); 
			edge_notifier.setContainer(*this); 
			uedge_notifier.setContainer(*this);
		} 
		
		~BpUGraphExtender() {
			uedge_notifier.clear();
			edge_notifier.clear(); 
			node_notifier.clear(); 
			anode_notifier.clear(); 
			bnode_notifier.clear(); 
		}
		
		Edge findEdge(Node u, Node v, Edge prev = INVALID) const {
			UEdge uedge = Parent::findUEdge(u, v, prev);
			if (uedge != INVALID) {
				return Parent::direct(uedge, Parent::aNode(u));
			} else {
				return INVALID;
			}
		}
		
	};
	
}

#endif

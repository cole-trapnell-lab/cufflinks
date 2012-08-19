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

#ifndef LEMON_SMART_GRAPH_H
#define LEMON_SMART_GRAPH_H

///\ingroup graphs
///\file
///\brief SmartGraph and SmartUGraph classes.

#include <vector>

#include <lemon/bits/invalid.h>

#include <lemon/bits/base_extender.h>
#include <lemon/bits/graph_extender.h>

#include <lemon/bits/utility.h>
#include <lemon/error.h>

#include <lemon/bits/graph_extender.h>

namespace lemon {
	
	class SmartGraph;
	///Base of SmartGraph
	
	///Base of SmartGraph
	///
	class SmartGraphBase {
	protected:
		
		struct NodeT 
		{
			int first_in, first_out;      
			NodeT() {}
		};
		struct EdgeT 
		{
			int target, source, next_in, next_out;      
			EdgeT() {}  
		};
		
		std::vector<NodeT> nodes;
		
		std::vector<EdgeT> edges;
		
		
	public:
		
		typedef SmartGraphBase Graph;
		
		class Node;
		class Edge;
		
		
	public:
		
		SmartGraphBase() : nodes(), edges() { }
		SmartGraphBase(const SmartGraphBase &_g) 
		: nodes(_g.nodes), edges(_g.edges) { }
		
		typedef True NodeNumTag;
		typedef True EdgeNumTag;
		
		int nodeNum() const { return nodes.size(); }
		int edgeNum() const { return edges.size(); }
		
		int maxNodeId() const { return nodes.size()-1; }
		int maxEdgeId() const { return edges.size()-1; }
		
		Node addNode() {
			int n = nodes.size();     
			nodes.push_back(NodeT());
			nodes[n].first_in = -1;
			nodes[n].first_out = -1;
			return Node(n);
		}
		
		Edge addEdge(Node u, Node v) {
			int n = edges.size(); 
			edges.push_back(EdgeT());
			edges[n].source = u.id; 
			edges[n].target = v.id;
			edges[n].next_out = nodes[u.id].first_out;
			edges[n].next_in = nodes[v.id].first_in;
			nodes[u.id].first_out = nodes[v.id].first_in = n;
			
			return Edge(n);
		}
		
		void clear() {
			edges.clear();
			nodes.clear();
		}
		
		Node source(Edge e) const { return Node(edges[e.id].source); }
		Node target(Edge e) const { return Node(edges[e.id].target); }
		
		static int id(Node v) { return v.id; }
		static int id(Edge e) { return e.id; }
		
		static Node nodeFromId(int id) { return Node(id);}
		static Edge edgeFromId(int id) { return Edge(id);}
		
		class Node {
			friend class SmartGraphBase;
			friend class SmartGraph;
			
		protected:
			int id;
			explicit Node(int _id) : id(_id) {}
		public:
			Node() {}
			Node (Invalid) : id(-1) {}
			bool operator==(const Node i) const {return id == i.id;}
			bool operator!=(const Node i) const {return id != i.id;}
			bool operator<(const Node i) const {return id < i.id;}
		};
		
		
		class Edge {
			friend class SmartGraphBase;
			friend class SmartGraph;
			
		protected:
			int id;
			explicit Edge(int _id) : id(_id) {}
		public:
			Edge() { }
			Edge (Invalid) : id(-1) {}
			bool operator==(const Edge i) const {return id == i.id;}
			bool operator!=(const Edge i) const {return id != i.id;}
			bool operator<(const Edge i) const {return id < i.id;}
		};
		
		void first(Node& node) const {
			node.id = nodes.size() - 1;
		}
		
		static void next(Node& node) {
			--node.id;
		}
		
		void first(Edge& edge) const {
			edge.id = edges.size() - 1;
		}
		
		static void next(Edge& edge) {
			--edge.id;
		}
		
		void firstOut(Edge& edge, const Node& node) const {
			edge.id = nodes[node.id].first_out;
		}
		
		void nextOut(Edge& edge) const {
			edge.id = edges[edge.id].next_out;
		}
		
		void firstIn(Edge& edge, const Node& node) const {
			edge.id = nodes[node.id].first_in;
		}
		
		void nextIn(Edge& edge) const {
			edge.id = edges[edge.id].next_in;
		}
		
	};
	
	typedef GraphExtender<SmartGraphBase> ExtendedSmartGraphBase;
	
	///\ingroup graphs
	///
	///\brief A smart graph class.
	///
	///This is a simple and fast graph implementation.
	///It is also quite memory efficient, but at the price
	///that <b> it does support only limited (only stack-like)
	///node and edge deletions</b>.
	///It conforms to 
	///the \ref concepts::Graph "Graph concept" with an
	///important extra feature that
	///its maps are real \ref concepts::ReferenceMap "reference map"s.
	///
	///\sa concepts::Graph.
	///
	///\author Alpar Juttner
	class SmartGraph : public ExtendedSmartGraphBase {
	public:
		
		typedef ExtendedSmartGraphBase Parent;
		
	private:
		
		///SmartGraph is \e not copy constructible. Use GraphCopy() instead.
		
		///SmartGraph is \e not copy constructible. Use GraphCopy() instead.
		///
		SmartGraph(const SmartGraph &) : ExtendedSmartGraphBase() {};
		///\brief Assignment of SmartGraph to another one is \e not allowed.
		///Use GraphCopy() instead.
		
		///Assignment of SmartGraph to another one is \e not allowed.
		///Use GraphCopy() instead.
		void operator=(const SmartGraph &) {}
		
	public:
		
		/// Constructor
		
		/// Constructor.
		///
		SmartGraph() {};
		
		///Add a new node to the graph.
		
		/// \return the new node.
		///
		Node addNode() { return Parent::addNode(); }
		
		///Add a new edge to the graph.
		
		///Add a new edge to the graph with source node \c s
		///and target node \c t.
		///\return the new edge.
		Edge addEdge(const Node& s, const Node& t) { 
			return Parent::addEdge(s, t); 
		}
		
		/// \brief Using this it is possible to avoid the superfluous memory
		/// allocation.
		
		/// Using this it is possible to avoid the superfluous memory
		/// allocation: if you know that the graph you want to build will
		/// be very large (e.g. it will contain millions of nodes and/or edges)
		/// then it is worth reserving space for this amount before starting
		/// to build the graph.
		/// \sa reserveEdge
		void reserveNode(int n) { nodes.reserve(n); };
		
		/// \brief Using this it is possible to avoid the superfluous memory
		/// allocation.
		
		/// Using this it is possible to avoid the superfluous memory
		/// allocation: if you know that the graph you want to build will
		/// be very large (e.g. it will contain millions of nodes and/or edges)
		/// then it is worth reserving space for this amount before starting
		/// to build the graph.
		/// \sa reserveNode
		void reserveEdge(int m) { edges.reserve(m); };
		
		///Clear the graph.
		
		///Erase all the nodes and edges from the graph.
		///
		void clear() {
			Parent::clear();
		}
		
		///Split a node.
		
		///This function splits a node. First a new node is added to the graph,
		///then the source of each outgoing edge of \c n is moved to this new node.
		///If \c connect is \c true (this is the default value), then a new edge
		///from \c n to the newly created node is also added.
		///\return The newly created node.
		///
		///\note The <tt>Edge</tt>s
		///referencing a moved edge remain
		///valid. However <tt>InEdge</tt>'s and <tt>OutEdge</tt>'s
		///may be invalidated.
		///\warning This functionality cannot be used together with the Snapshot
		///feature.
		///\todo It could be implemented in a bit faster way.
		Node split(Node n, bool connect = true)
		{
			Node b = addNode();
			nodes[b.id].first_out=nodes[n.id].first_out;
			nodes[n.id].first_out=-1;
			for(int i=nodes[b.id].first_out;i!=-1;i++) edges[i].source=b.id;
			if(connect) addEdge(n,b);
			return b;
		}
		
	public:
		
		class Snapshot;
		
	protected:
		
		void restoreSnapshot(const Snapshot &s)
		{
			while(s.edge_num<edges.size()) {
				Edge edge = edgeFromId(edges.size()-1);
				Parent::notifier(Edge()).erase(edge);
				nodes[edges.back().source].first_out=edges.back().next_out;
				nodes[edges.back().target].first_in=edges.back().next_in;
				edges.pop_back();
			}
			while(s.node_num<nodes.size()) {
				Node node = nodeFromId(nodes.size()-1);
				Parent::notifier(Node()).erase(node);
				nodes.pop_back();
			}
		}    
		
	public:
		
		///Class to make a snapshot of the graph and to restrore to it later.
		
		///Class to make a snapshot of the graph and to restrore to it later.
		///
		///The newly added nodes and edges can be removed using the
		///restore() function.
		///\note After you restore a state, you cannot restore
		///a later state, in other word you cannot add again the edges deleted
		///by restore() using another one Snapshot instance.
		///
		///\warning If you do not use correctly the snapshot that can cause
		///either broken program, invalid state of the graph, valid but
		///not the restored graph or no change. Because the runtime performance
		///the validity of the snapshot is not stored.
		class Snapshot 
		{
			SmartGraph *g;
		protected:
			friend class SmartGraph;
			unsigned int node_num;
			unsigned int edge_num;
		public:
			///Default constructor.
			
			///Default constructor.
			///To actually make a snapshot you must call save().
			///
			Snapshot() : g(0) {}
			///Constructor that immediately makes a snapshot
			
			///This constructor immediately makes a snapshot of the graph.
			///\param _g The graph we make a snapshot of.
			Snapshot(SmartGraph &_g) :g(&_g) {
				node_num=g->nodes.size();
				edge_num=g->edges.size();
			}
			
			///Make a snapshot.
			
			///Make a snapshot of the graph.
			///
			///This function can be called more than once. In case of a repeated
			///call, the previous snapshot gets lost.
			///\param _g The graph we make the snapshot of.
			void save(SmartGraph &_g) 
			{
				g=&_g;
				node_num=g->nodes.size();
				edge_num=g->edges.size();
			}
			
			///Undo the changes until a snapshot.
			
			///Undo the changes until a snapshot created by save().
			///
			///\note After you restored a state, you cannot restore
			///a later state, in other word you cannot add again the edges deleted
			///by restore().
			void restore()
			{
				g->restoreSnapshot(*this);
			}
		};
	};
	
	
	class SmartUGraphBase {
		
	protected:
		
		struct NodeT {
			int first_out;
		};
		
		struct EdgeT {
			int target;
			int next_out;
		};
		
		std::vector<NodeT> nodes;
		std::vector<EdgeT> edges;
		
		int first_free_edge;
		
	public:
		
		typedef SmartUGraphBase Graph;
		
		class Node;
		class Edge;
		class UEdge;
		
		class Node {
			friend class SmartUGraphBase;
		protected:
			
			int id;
			explicit Node(int pid) { id = pid;}
			
		public:
			Node() {}
			Node (Invalid) { id = -1; }
			bool operator==(const Node& node) const {return id == node.id;}
			bool operator!=(const Node& node) const {return id != node.id;}
			bool operator<(const Node& node) const {return id < node.id;}
		};
		
		class UEdge {
			friend class SmartUGraphBase;
		protected:
			
			int id;
			explicit UEdge(int pid) { id = pid;}
			
		public:
			UEdge() {}
			UEdge (Invalid) { id = -1; }
			bool operator==(const UEdge& edge) const {return id == edge.id;}
			bool operator!=(const UEdge& edge) const {return id != edge.id;}
			bool operator<(const UEdge& edge) const {return id < edge.id;}
		};
		
		class Edge {
			friend class SmartUGraphBase;
		protected:
			
			int id;
			explicit Edge(int pid) { id = pid;}
			
		public:
			operator UEdge() const { return uEdgeFromId(id / 2); }
			
			Edge() {}
			Edge (Invalid) { id = -1; }
			bool operator==(const Edge& edge) const {return id == edge.id;}
			bool operator!=(const Edge& edge) const {return id != edge.id;}
			bool operator<(const Edge& edge) const {return id < edge.id;}
		};
		
		
		
		SmartUGraphBase()
		: nodes(), edges() {}
		
		
		int maxNodeId() const { return nodes.size()-1; } 
		int maxUEdgeId() const { return edges.size() / 2 - 1; }
		int maxEdgeId() const { return edges.size()-1; }
		
		Node source(Edge e) const { return Node(edges[e.id ^ 1].target); }
		Node target(Edge e) const { return Node(edges[e.id].target); }
		
		Node source(UEdge e) const { return Node(edges[2 * e.id].target); }
		Node target(UEdge e) const { return Node(edges[2 * e.id + 1].target); }
		
		static bool direction(Edge e) {
			return (e.id & 1) == 1;
		}
		
		static Edge direct(UEdge e, bool d) {
			return Edge(e.id * 2 + (d ? 1 : 0));
		}
		
		void first(Node& node) const { 
			node.id = nodes.size() - 1;
		}
		
		void next(Node& node) const {
			--node.id;
		}
		
		void first(Edge& edge) const { 
			edge.id = edges.size() - 1;
		}
		
		void next(Edge& edge) const {
			--edge.id;
		}
		
		void first(UEdge& edge) const { 
			edge.id = edges.size() / 2 - 1;
		}
		
		void next(UEdge& edge) const {
			--edge.id;
		}
		
		void firstOut(Edge &edge, const Node& v) const {
			edge.id = nodes[v.id].first_out;
		}
		void nextOut(Edge &edge) const {
			edge.id = edges[edge.id].next_out;
		}
		
		void firstIn(Edge &edge, const Node& v) const {
			edge.id = ((nodes[v.id].first_out) ^ 1);
			if (edge.id == -2) edge.id = -1;
		}
		void nextIn(Edge &edge) const {
			edge.id = ((edges[edge.id ^ 1].next_out) ^ 1);
			if (edge.id == -2) edge.id = -1;
		}
		
		void firstInc(UEdge &edge, bool& d, const Node& v) const {
			int de = nodes[v.id].first_out;
			if (de != -1) {
				edge.id = de / 2;
				d = ((de & 1) == 1);
			} else {
				edge.id = -1;
				d = true;
			}
		}
		void nextInc(UEdge &edge, bool& d) const {
			int de = (edges[(edge.id * 2) | (d ? 1 : 0)].next_out);
			if (de != -1) {
				edge.id = de / 2;
				d = ((de & 1) == 1);
			} else {
				edge.id = -1;
				d = true;      
			}
		}
		
		static int id(Node v) { return v.id; }
		static int id(Edge e) { return e.id; }
		static int id(UEdge e) { return e.id; }
		
		static Node nodeFromId(int id) { return Node(id);}
		static Edge edgeFromId(int id) { return Edge(id);}
		static UEdge uEdgeFromId(int id) { return UEdge(id);}
		
		Node addNode() {     
			int n = nodes.size();
			nodes.push_back(NodeT());
			nodes[n].first_out = -1;
			
			return Node(n);
		}
		
		UEdge addEdge(Node u, Node v) {
			int n = edges.size();
			edges.push_back(EdgeT());
			edges.push_back(EdgeT());
			
			edges[n].target = u.id;
			edges[n | 1].target = v.id;
			
			edges[n].next_out = nodes[v.id].first_out;
			nodes[v.id].first_out = n;
			
			edges[n | 1].next_out = nodes[u.id].first_out;	
			nodes[u.id].first_out = (n | 1);
			
			return UEdge(n / 2);
		}
		
		void clear() {
			edges.clear();
			nodes.clear();
		}
		
	};
	
	typedef UGraphExtender<SmartUGraphBase> ExtendedSmartUGraphBase;
	
	/// \ingroup graphs
	///
	/// \brief A smart undirected graph class.
	///
	/// This is a simple and fast undirected graph implementation.
	/// It is also quite memory efficient, but at the price
	/// that <b> it does support only limited (only stack-like)
	/// node and edge deletions</b>.
	/// Except from this it conforms to 
	/// the \ref concepts::UGraph "UGraph concept".
	///
	///It also has an
	///important extra feature that
	///its maps are real \ref concepts::ReferenceMap "reference map"s.
	///
	/// \sa concepts::UGraph.
	///
	class SmartUGraph : public ExtendedSmartUGraphBase {
	private:
		
		///SmartUGraph is \e not copy constructible. Use UGraphCopy() instead.
		
		///SmartUGraph is \e not copy constructible. Use UGraphCopy() instead.
		///
		SmartUGraph(const SmartUGraph &) : ExtendedSmartUGraphBase() {};
		
		///\brief Assignment of SmartUGraph to another one is \e not allowed.
		///Use UGraphCopy() instead.
		
		///Assignment of SmartUGraph to another one is \e not allowed.
		///Use UGraphCopy() instead.
		void operator=(const SmartUGraph &) {}
		
	public:
		
		typedef ExtendedSmartUGraphBase Parent;
		typedef Parent::OutEdgeIt IncEdgeIt;
		
		/// Constructor
		
		/// Constructor.
		///
		SmartUGraph() {}
		
		///Add a new node to the graph.
		
		/// \return the new node.
		///
		Node addNode() { return Parent::addNode(); }
		
		///Add a new undirected edge to the graph.
		
		///Add a new undirected edge to the graph with node \c s
		///and \c t.
		///\return the new undirected edge.
		UEdge addEdge(const Node& s, const Node& t) { 
			return Parent::addEdge(s, t); 
		}
		
		///Clear the graph.
		
		///Erase all the nodes and edges from the graph.
		///
		void clear() {
			Parent::clear();
		}
		
	public:
		
		class Snapshot;
		
	protected:
		
		void saveSnapshot(Snapshot &s)
		{
			s.graph = this;
			s.node_num = nodes.size();
			s.edge_num = edges.size();
		}
		
		void restoreSnapshot(const Snapshot &s)
		{
			while(s.edge_num<edges.size()) {
				int n=edges.size()-1;
				UEdge edge=uEdgeFromId(n/2);
				Parent::notifier(UEdge()).erase(edge);
				std::vector<Edge> dir;
				dir.push_back(edgeFromId(n));
				dir.push_back(edgeFromId(n-1));
				Parent::notifier(Edge()).erase(dir);
				nodes[edges[n].target].first_out=edges[n].next_out;
				nodes[edges[n-1].target].first_out=edges[n-1].next_out;
				edges.pop_back();
				edges.pop_back();
			}
			while(s.node_num<nodes.size()) {
				int n=nodes.size()-1;
				Node node = nodeFromId(n);
				Parent::notifier(Node()).erase(node);
				nodes.pop_back();
			}
		}    
		
	public:
		
		///Class to make a snapshot of the graph and to restrore to it later.
		
		///Class to make a snapshot of the graph and to restrore to it later.
		///
		///The newly added nodes and edges can be removed using the
		///restore() function.
		///
		///\note After you restore a state, you cannot restore
		///a later state, in other word you cannot add again the edges deleted
		///by restore() using another one Snapshot instance.
		///
		///\warning If you do not use correctly the snapshot that can cause
		///either broken program, invalid state of the graph, valid but
		///not the restored graph or no change. Because the runtime performance
		///the validity of the snapshot is not stored.
		class Snapshot 
		{
			SmartUGraph *graph;
		protected:
			friend class SmartUGraph;
			unsigned int node_num;
			unsigned int edge_num;
		public:
			///Default constructor.
			
			///Default constructor.
			///To actually make a snapshot you must call save().
			///
			Snapshot() : graph(0) {}
			///Constructor that immediately makes a snapshot
			
			///This constructor immediately makes a snapshot of the graph.
			///\param g The graph we make a snapshot of.
			Snapshot(SmartUGraph &g) {
				g.saveSnapshot(*this);
			}
			
			///Make a snapshot.
			
			///Make a snapshot of the graph.
			///
			///This function can be called more than once. In case of a repeated
			///call, the previous snapshot gets lost.
			///\param g The graph we make the snapshot of.
			void save(SmartUGraph &g) 
			{
				g.saveSnapshot(*this);
			}
			
			///Undo the changes until a snapshot.
			
			///Undo the changes until a snapshot created by save().
			///
			///\note After you restored a state, you cannot restore
			///a later state, in other word you cannot add again the edges deleted
			///by restore().
			void restore()
			{
				graph->restoreSnapshot(*this);
			}
		};
	};
	
	
	class SmartBpUGraphBase {
	public:
		
		class NodeSetError : public LogicError {
		public:
			virtual const char* what() const throw() { 
				return "lemon::SmartBpUGraph::NodeSetError";
			}
		};
		
	protected:
		
		struct NodeT {
			int first;
			NodeT() {}
			NodeT(int _first) : first(_first) {}
		};
		
		struct UEdgeT {
			int aNode, next_out;
			int bNode, next_in;
		};
		
		std::vector<NodeT> aNodes;
		std::vector<NodeT> bNodes;
		
		std::vector<UEdgeT> edges;
		
	public:
		
		class Node {
			friend class SmartBpUGraphBase;
		protected:
			int id;
			
			explicit Node(int _id) : id(_id) {}
		public:
			Node() {}
			Node(Invalid) : id(-1) {}
			bool operator==(const Node i) const {return id==i.id;}
			bool operator!=(const Node i) const {return id!=i.id;}
			bool operator<(const Node i) const {return id<i.id;}
		};
		
		class UEdge {
			friend class SmartBpUGraphBase;
		protected:
			int id;
			
			UEdge(int _id) : id(_id) {}
		public:
			UEdge() {}
			UEdge(Invalid) : id(-1) {}
			bool operator==(const UEdge i) const {return id==i.id;}
			bool operator!=(const UEdge i) const {return id!=i.id;}
			bool operator<(const UEdge i) const {return id<i.id;}
		};
		
		void firstANode(Node& node) const {
			node.id = 2 * aNodes.size() - 2;
			if (node.id < 0) node.id = -1; 
		}
		void nextANode(Node& node) const {
			node.id -= 2;
			if (node.id < 0) node.id = -1; 
		}
		
		void firstBNode(Node& node) const {
			node.id = 2 * bNodes.size() - 1;
		}
		void nextBNode(Node& node) const {
			node.id -= 2;
		}
		
		void first(Node& node) const {
			if (aNodes.size() > 0) {
				node.id = 2 * aNodes.size() - 2;
			} else {
				node.id = 2 * bNodes.size() - 1;
			}
		}
		void next(Node& node) const {
			node.id -= 2;
			if (node.id == -2) {
				node.id = 2 * bNodes.size() - 1;
			}
		}
		
		void first(UEdge& edge) const {
			edge.id = edges.size() - 1;
		}
		void next(UEdge& edge) const {
			--edge.id;
		}
		
		void firstFromANode(UEdge& edge, const Node& node) const {
			LEMON_ASSERT((node.id & 1) == 0, NodeSetError());
			edge.id = aNodes[node.id >> 1].first;
		}
		void nextFromANode(UEdge& edge) const {
			edge.id = edges[edge.id].next_out;
		}
		
		void firstFromBNode(UEdge& edge, const Node& node) const {
			LEMON_ASSERT((node.id & 1) == 1, NodeSetError());
			edge.id = bNodes[node.id >> 1].first;
		}
		void nextFromBNode(UEdge& edge) const {
			edge.id = edges[edge.id].next_in;
		}
		
		static int id(const Node& node) {
			return node.id;
		}
		static Node nodeFromId(int id) {
			return Node(id);
		}
		int maxNodeId() const {
			return aNodes.size() > bNodes.size() ?
			aNodes.size() * 2 - 2 : bNodes.size() * 2 - 1;
		}
		
		static int id(const UEdge& edge) {
			return edge.id;
		}
		static UEdge uEdgeFromId(int id) {
			return UEdge(id);
		}
		int maxUEdgeId() const {
			return edges.size();
		}
		
		static int aNodeId(const Node& node) {
			return node.id >> 1;
		}
		static Node nodeFromANodeId(int id) {
			return Node(id << 1);
		}
		int maxANodeId() const {
			return aNodes.size();
		}
		
		static int bNodeId(const Node& node) {
			return node.id >> 1;
		}
		static Node nodeFromBNodeId(int id) {
			return Node((id << 1) + 1);
		}
		int maxBNodeId() const {
			return bNodes.size();
		}
		
		Node aNode(const UEdge& edge) const {
			return Node(edges[edge.id].aNode);
		}
		Node bNode(const UEdge& edge) const {
			return Node(edges[edge.id].bNode);
		}
		
		static bool aNode(const Node& node) {
			return (node.id & 1) == 0;
		}
		
		static bool bNode(const Node& node) {
			return (node.id & 1) == 1;
		}
		
		Node addANode() {
			NodeT nodeT;
			nodeT.first = -1;
			aNodes.push_back(nodeT);
			return Node(aNodes.size() * 2 - 2);
		}
		
		Node addBNode() {
			NodeT nodeT;
			nodeT.first = -1;
			bNodes.push_back(nodeT);
			return Node(bNodes.size() * 2 - 1);
		}
		
		UEdge addEdge(const Node& source, const Node& target) {
			LEMON_ASSERT(((source.id ^ target.id) & 1) == 1, NodeSetError());
			UEdgeT edgeT;
			if ((source.id & 1) == 0) {
				edgeT.aNode = source.id;
				edgeT.bNode = target.id;
			} else {
				edgeT.aNode = target.id;
				edgeT.bNode = source.id;
			}
			edgeT.next_out = aNodes[edgeT.aNode >> 1].first;
			aNodes[edgeT.aNode >> 1].first = edges.size();
			edgeT.next_in = bNodes[edgeT.bNode >> 1].first;
			bNodes[edgeT.bNode >> 1].first = edges.size();
			edges.push_back(edgeT);
			return UEdge(edges.size() - 1);
		}
		
		void reserveANode(int n) { aNodes.reserve(n); };
		void reserveBNode(int n) { bNodes.reserve(n); };
		
		void reserveEdge(int m) { edges.reserve(m); };
		
		void clear() {
			aNodes.clear();
			bNodes.clear();
			edges.clear();
		}
		
		typedef True NodeNumTag;
		int nodeNum() const { return aNodes.size() + bNodes.size(); }
		int aNodeNum() const { return aNodes.size(); }
		int bNodeNum() const { return bNodes.size(); }
		
		typedef True EdgeNumTag;
		int uEdgeNum() const { return edges.size(); }
		
	};
	
	
	typedef BpUGraphExtender<BidirBpUGraphExtender<SmartBpUGraphBase> >
		ExtendedSmartBpUGraphBase;
	
	/// \ingroup graphs
	///
	/// \brief A smart bipartite undirected graph class.
	///
	/// This is a simple and fast bipartite undirected graph implementation.
	/// It is also quite memory efficient, but at the price
	/// that <b> it does not support node and edge deletions</b>.
	/// Except from this it conforms to 
	/// the \ref concepts::BpUGraph "BpUGraph concept".
	///
	///It also has an
	///important extra feature that
	///its maps are real \ref concepts::ReferenceMap "reference map"s.
	///
	/// \sa concepts::BpUGraph.
	///
	class SmartBpUGraph : public ExtendedSmartBpUGraphBase {
	private:
		
		/// \brief SmartBpUGraph is \e not copy constructible.
		///
		///SmartBpUGraph is \e not copy constructible.
		SmartBpUGraph(const SmartBpUGraph &) : ExtendedSmartBpUGraphBase() {};
		
		/// \brief Assignment of SmartBpUGraph to another one is \e not
		/// allowed.
		///
		/// Assignment of SmartBpUGraph to another one is \e not allowed.
		void operator=(const SmartBpUGraph &) {}
		
	public:
		
		typedef ExtendedSmartBpUGraphBase Parent;
		
		///Constructor
		
		///Constructor.
		///
		SmartBpUGraph() : ExtendedSmartBpUGraphBase() {}
		
		///Add a new ANode to the graph.
		
		/// \return the new node.
		///
		Node addANode() { return Parent::addANode(); }
		
		///Add a new BNode to the graph.
		
		/// \return the new node.
		///
		Node addBNode() { return Parent::addBNode(); }
		
		///Add a new undirected edge to the graph.
		
		///Add a new undirected edge to the graph with node \c s
		///and \c t.
		///\return the new undirected edge.
		UEdge addEdge(const Node& s, const Node& t) { 
			return Parent::addEdge(s, t); 
		}
		
		///Clear the graph.
		
		///Erase all the nodes and edges from the graph.
		///
		void clear() {
			Parent::clear();
		}
		
	public:
		
		class Snapshot;
		
	protected:
		
		void restoreSnapshot(const Snapshot &s)
		{
			while(s.edge_num<edges.size()) {
				UEdge edge = uEdgeFromId(edges.size()-1);
				Parent::notifier(UEdge()).erase(edge);
				std::vector<Edge> dir;
				dir.push_back(Parent::direct(edge, true));
				dir.push_back(Parent::direct(edge, false));
				Parent::notifier(Edge()).erase(dir);
				aNodes[edges.back().aNode >> 1].first=edges.back().next_out;
				bNodes[edges.back().bNode >> 1].first=edges.back().next_in;
				edges.pop_back();
			}
			while(s.anode_num<aNodes.size()) {
				Node node = nodeFromANodeId(aNodes.size() - 1);
				Parent::notifier(ANode()).erase(node);
				Parent::notifier(Node()).erase(node);
				aNodes.pop_back();
			}
			while(s.bnode_num<bNodes.size()) {
				Node node = nodeFromBNodeId(bNodes.size() - 1);
				Parent::notifier(BNode()).erase(node);
				Parent::notifier(Node()).erase(node);
				bNodes.pop_back();
			}
		}    
		
	public:
		
		///Class to make a snapshot of the graph and to restrore to it later.
		
		///Class to make a snapshot of the graph and to restrore to it later.
		///
		///The newly added nodes and edges can be removed using the
		///restore() function.
		///
		///\note After you restore a state, you cannot restore
		///a later state, in other word you cannot add again the edges deleted
		///by restore() using another one Snapshot instance.
		///
		///\warning If you do not use correctly the snapshot that can cause
		///either broken program, invalid state of the graph, valid but
		///not the restored graph or no change. Because the runtime performance
		///the validity of the snapshot is not stored.
		class Snapshot 
		{
			SmartBpUGraph *g;
		protected:
			friend class SmartBpUGraph;
			unsigned int anode_num;
			unsigned int bnode_num;
			unsigned int edge_num;
		public:
			///Default constructor.
			
			///Default constructor.
			///To actually make a snapshot you must call save().
			///
			Snapshot() : g(0) {}
			
			///Constructor that immediately makes a snapshot
			
			///This constructor immediately makes a snapshot of the graph.
			///\param _g The graph we make a snapshot of.
			Snapshot(SmartBpUGraph &_g) : g(&_g) {
				anode_num=g->aNodes.size();
				bnode_num=g->bNodes.size();
				edge_num=g->edges.size();
			}
			
			///Make a snapshot.
			
			///Make a snapshot of the graph.
			///
			///This function can be called more than once. In case of a repeated
			///call, the previous snapshot gets lost.
			///\param _g The graph we make the snapshot of.
			void save(SmartBpUGraph &_g) 
			{
				g=&_g;
				anode_num=g->aNodes.size();
				bnode_num=g->bNodes.size();
				edge_num=g->edges.size();
			}
			
			///Undo the changes until a snapshot.
			
			///Undo the changes until a snapshot created by save().
			///
			///\note After you restored a state, you cannot restore
			///a later state, in other word you cannot add again the edges deleted
			///by restore().
			void restore()
			{
				g->restoreSnapshot(*this);
			}
		};
	};
	
	
	/// @}  
} //namespace lemon


#endif //LEMON_SMART_GRAPH_H

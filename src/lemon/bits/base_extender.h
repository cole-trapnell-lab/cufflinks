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

#ifndef LEMON_BITS_BASE_EXTENDER_H
#define LEMON_BITS_BASE_EXTENDER_H

#include <lemon/bits/invalid.h>
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
	/// \brief BaseGraph to BaseUGraph extender
	template <typename Base>
	class UndirGraphExtender : public Base {
		
	public:
		
		typedef Base Parent;
		typedef typename Parent::Edge UEdge;
		typedef typename Parent::Node Node;
		
		typedef True UndirectedTag;
		
		class Edge : public UEdge {
			friend class UndirGraphExtender;
			
		protected:
			bool forward;
			
			Edge(const UEdge &ue, bool _forward) :
			UEdge(ue), forward(_forward) {}
			
		public:
			Edge() {}
			
			/// Invalid edge constructor
			Edge(Invalid i) : UEdge(i), forward(true) {}
			
			bool operator==(const Edge &that) const {
				return forward==that.forward && UEdge(*this)==UEdge(that);
			}
			bool operator!=(const Edge &that) const {
				return forward!=that.forward || UEdge(*this)!=UEdge(that);
			}
			bool operator<(const Edge &that) const {
				return forward<that.forward ||
				(!(that.forward<forward) && UEdge(*this)<UEdge(that));
			}
		};
		
		
		
		using Parent::source;
		
		/// Source of the given Edge.
		Node source(const Edge &e) const {
			return e.forward ? Parent::source(e) : Parent::target(e);
		}
		
		using Parent::target;
		
		/// Target of the given Edge.
		Node target(const Edge &e) const {
			return e.forward ? Parent::target(e) : Parent::source(e);
		}
		
		/// \brief Directed edge from an undirected edge.
		///
		/// Returns a directed edge corresponding to the specified UEdge.
		/// If the given bool is true the given undirected edge and the
		/// returned edge have the same source node.
		static Edge direct(const UEdge &ue, bool d) {
			return Edge(ue, d);
		}
		
		/// Returns whether the given directed edge is same orientation as the
		/// corresponding undirected edge.
		///
		/// \todo reference to the corresponding point of the undirected graph
		/// concept. "What does the direction of an undirected edge mean?"
		static bool direction(const Edge &e) { return e.forward; }
		
		
		using Parent::first;
		using Parent::next;
		
		void first(Edge &e) const {
			Parent::first(e);
			e.forward=true;
		}
		
		void next(Edge &e) const {
			if( e.forward ) {
				e.forward = false;
			}
			else {
				Parent::next(e);
				e.forward = true;
			}
		}
		
		void firstOut(Edge &e, const Node &n) const {
			Parent::firstIn(e,n);
			if( UEdge(e) != INVALID ) {
				e.forward = false;
			}
			else {
				Parent::firstOut(e,n);
				e.forward = true;
			}
		}
		void nextOut(Edge &e) const {
			if( ! e.forward ) {
				Node n = Parent::target(e);
				Parent::nextIn(e);
				if( UEdge(e) == INVALID ) {
					Parent::firstOut(e, n);
					e.forward = true;
				}
			}
			else {
				Parent::nextOut(e);
			}
		}
		
		void firstIn(Edge &e, const Node &n) const {
			Parent::firstOut(e,n);
			if( UEdge(e) != INVALID ) {
				e.forward = false;
			}
			else {
				Parent::firstIn(e,n);
				e.forward = true;
			}
		}
		void nextIn(Edge &e) const {
			if( ! e.forward ) {
				Node n = Parent::source(e);
				Parent::nextOut(e);
				if( UEdge(e) == INVALID ) {
					Parent::firstIn(e, n);
					e.forward = true;
				}
			}
			else {
				Parent::nextIn(e);
			}
		}
		
		void firstInc(UEdge &e, bool &d, const Node &n) const {
			d = true;
			Parent::firstOut(e, n);
			if (e != INVALID) return;
			d = false;
			Parent::firstIn(e, n);
		}
		
		void nextInc(UEdge &e, bool &d) const {
			if (d) {
				Node s = Parent::source(e);
				Parent::nextOut(e);
				if (e != INVALID) return;
				d = false;
				Parent::firstIn(e, s);
			} else {
				Parent::nextIn(e);
			}
		}
		
		Node nodeFromId(int ix) const {
			return Parent::nodeFromId(ix);
		}
		
		Edge edgeFromId(int ix) const {
			return direct(Parent::edgeFromId(ix >> 1), bool(ix & 1));
		}
		
		UEdge uEdgeFromId(int ix) const {
			return Parent::edgeFromId(ix);
		}
		
		int id(const Node &n) const {
			return Parent::id(n);
		}
		
		int id(const UEdge &e) const {
			return Parent::id(e);
		}
		
		int id(const Edge &e) const {
			return 2 * Parent::id(e) + int(e.forward);
		}
		
		int maxNodeId() const {
			return Parent::maxNodeId();
		}
		
		int maxEdgeId() const {
			return 2 * Parent::maxEdgeId() + 1;
		}
		
		int maxUEdgeId() const {
			return Parent::maxEdgeId();
		}
		
		
		int edgeNum() const {
			return 2 * Parent::edgeNum();
		}
		
		int uEdgeNum() const {
			return Parent::edgeNum();
		}
		
		Edge findEdge(Node s, Node t, Edge p = INVALID) const {
			if (p == INVALID) {
				UEdge edge = Parent::findEdge(s, t);
				if (edge != INVALID) return direct(edge, true);
				edge = Parent::findEdge(t, s);
				if (edge != INVALID) return direct(edge, false);
			} else if (direction(p)) {
				UEdge edge = Parent::findEdge(s, t, p);
				if (edge != INVALID) return direct(edge, true);
				edge = Parent::findEdge(t, s);
				if (edge != INVALID) return direct(edge, false);	
			} else {
				UEdge edge = Parent::findEdge(t, s, p);
				if (edge != INVALID) return direct(edge, false);	      
			}
			return INVALID;
		}
		
		UEdge findUEdge(Node s, Node t, UEdge p = INVALID) const {
			if (s != t) {
				if (p == INVALID) {
					UEdge edge = Parent::findEdge(s, t);
					if (edge != INVALID) return edge;
					edge = Parent::findEdge(t, s);
					if (edge != INVALID) return edge;
				} else if (Parent::s(p) == s) {
					UEdge edge = Parent::findEdge(s, t, p);
					if (edge != INVALID) return edge;
					edge = Parent::findEdge(t, s);
					if (edge != INVALID) return edge;	
				} else {
					UEdge edge = Parent::findEdge(t, s, p);
					if (edge != INVALID) return edge;	      
				}
			} else {
				return Parent::findEdge(s, t, p);
			}
			return INVALID;
		}
	};
	
	template <typename Base>
	class BidirBpUGraphExtender : public Base {
	public:
		typedef Base Parent;
		typedef BidirBpUGraphExtender Graph;
		
		typedef typename Parent::Node Node;
		typedef typename Parent::UEdge UEdge;
		
		
		using Parent::first;
		using Parent::next;
		
		using Parent::id;
		
		class ANode : public Node {
			friend class BidirBpUGraphExtender;
		public:
			ANode() {}
			ANode(const Node& node) : Node(node) {
				LEMON_ASSERT(Parent::aNode(node) || node == INVALID, 
							 typename Parent::NodeSetError());
			}
			ANode& operator=(const Node& node) {
				LEMON_ASSERT(Parent::aNode(node) || node == INVALID, 
							 typename Parent::NodeSetError());
				Node::operator=(node);
				return *this;
			}
			ANode(Invalid) : Node(INVALID) {}
			ANode& operator=(Invalid) {
				Node::operator=(INVALID);
				return *this;
			}
		};
		
		void first(ANode& node) const {
			Parent::firstANode(static_cast<Node&>(node));
		}
		void next(ANode& node) const {
			Parent::nextANode(static_cast<Node&>(node));
		}
		
		int id(const ANode& node) const {
			return Parent::aNodeId(node);
		}
		
		class BNode : public Node {
			friend class BidirBpUGraphExtender;
		public:
			BNode() {}
			BNode(const Node& node) : Node(node) {
				LEMON_ASSERT(Parent::bNode(node) || node == INVALID,
							 typename Parent::NodeSetError());
			}
			BNode& operator=(const Node& node) {
				LEMON_ASSERT(Parent::bNode(node) || node == INVALID, 
							 typename Parent::NodeSetError());
				Node::operator=(node);
				return *this;
			}
			BNode(Invalid) : Node(INVALID) {}
			BNode& operator=(Invalid) {
				Node::operator=(INVALID);
				return *this;
			}
		};
		
		void first(BNode& node) const {
			Parent::firstBNode(static_cast<Node&>(node));
		}
		void next(BNode& node) const {
			Parent::nextBNode(static_cast<Node&>(node));
		}
		
		int id(const BNode& node) const {
			return Parent::aNodeId(node);
		}
		
		Node source(const UEdge& edge) const {
			return aNode(edge);
		}
		Node target(const UEdge& edge) const {
			return bNode(edge);
		}
		
		void firstInc(UEdge& edge, bool& dir, const Node& node) const {
			if (Parent::aNode(node)) {
				Parent::firstFromANode(edge, node);
				dir = true;
			} else {
				Parent::firstFromBNode(edge, node);
				dir = static_cast<UEdge&>(edge) == INVALID;
			}
		}
		void nextInc(UEdge& edge, bool& dir) const {
			if (dir) {
				Parent::nextFromANode(edge);
			} else {
				Parent::nextFromBNode(edge);
				if (edge == INVALID) dir = true;
			}
		}
		
		class Edge : public UEdge {
			friend class BidirBpUGraphExtender;
		protected:
			bool forward;
			
			Edge(const UEdge& edge, bool _forward)
			: UEdge(edge), forward(_forward) {}
			
		public:
			Edge() {}
			Edge (Invalid) : UEdge(INVALID), forward(true) {}
			bool operator==(const Edge& i) const {
				return UEdge::operator==(i) && forward == i.forward;
			}
			bool operator!=(const Edge& i) const {
				return UEdge::operator!=(i) || forward != i.forward;
			}
			bool operator<(const Edge& i) const {
				return UEdge::operator<(i) || 
				(!(i.forward<forward) && UEdge(*this)<UEdge(i));
			}
		};
		
		void first(Edge& edge) const {
			Parent::first(static_cast<UEdge&>(edge));
			edge.forward = true;
		}
		
		void next(Edge& edge) const {
			if (!edge.forward) {
				Parent::next(static_cast<UEdge&>(edge));
			}
			edge.forward = !edge.forward;
		}
		
		void firstOut(Edge& edge, const Node& node) const {
			if (Parent::aNode(node)) {
				Parent::firstFromANode(edge, node);
				edge.forward = true;
			} else {
				Parent::firstFromBNode(edge, node);
				edge.forward = static_cast<UEdge&>(edge) == INVALID;
			}
		}
		void nextOut(Edge& edge) const {
			if (edge.forward) {
				Parent::nextFromANode(edge);
			} else {
				Parent::nextFromBNode(edge);
				edge.forward = static_cast<UEdge&>(edge) == INVALID;
			}
		}
		
		void firstIn(Edge& edge, const Node& node) const {
			if (Parent::bNode(node)) {
				Parent::firstFromBNode(edge, node);
				edge.forward = true;	
			} else {
				Parent::firstFromANode(edge, node);
				edge.forward = static_cast<UEdge&>(edge) == INVALID;
			}
		}
		void nextIn(Edge& edge) const {
			if (edge.forward) {
				Parent::nextFromBNode(edge);
			} else {
				Parent::nextFromANode(edge);
				edge.forward = static_cast<UEdge&>(edge) == INVALID;
			}
		}
		
		Node source(const Edge& edge) const {
			return edge.forward ? Parent::aNode(edge) : Parent::bNode(edge);
		}
		Node target(const Edge& edge) const {
			return edge.forward ? Parent::bNode(edge) : Parent::aNode(edge);
		}
		
		int id(const Edge& edge) const {
			return (Parent::id(static_cast<const UEdge&>(edge)) << 1) + 
			(edge.forward ? 0 : 1);
		}
		Edge edgeFromId(int ix) const {
			return Edge(Parent::fromUEdgeId(ix >> 1), (ix & 1) == 0);
		}
		int maxEdgeId() const {
			return (Parent::maxUEdgeId() << 1) + 1;
		}
		
		bool direction(const Edge& edge) const {
			return edge.forward;
		}
		
		Edge direct(const UEdge& edge, bool dir) const {
			return Edge(edge, dir);
		}
		
		int edgeNum() const {
			return 2 * Parent::uEdgeNum();
		}
		
		int uEdgeNum() const {
			return Parent::uEdgeNum();
		}
		
		
	};
}

#endif

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

#ifndef LEMON_BIPARTITE_MATCHING
#define LEMON_BIPARTITE_MATCHING

#include <functional>

#include <lemon/bin_heap.h>
#include <lemon/fib_heap.h>
#include <lemon/maps.h>

#include <iostream>

///\ingroup matching
///\file
///\brief Maximum matching algorithms in bipartite graphs.
///
///\note The pr_bipartite_matching.h file also contains algorithms to
///solve maximum cardinality bipartite matching problems.

namespace lemon {
	
	/// \ingroup matching
	///
	/// \brief Bipartite Max Cardinality Matching algorithm
	///
	/// Bipartite Max Cardinality Matching algorithm. This class implements
	/// the Hopcroft-Karp algorithm which has \f$ O(e\sqrt{n}) \f$ time
	/// complexity.
	///
	/// \note In several cases the push-relabel based algorithms have
	/// better runtime performance than the augmenting path based ones. 
	///
	/// \see PrBipartiteMatching
	template <typename BpUGraph>
	class MaxBipartiteMatching {
	protected:
		
		typedef BpUGraph Graph;
		
		typedef typename Graph::Node Node;
		typedef typename Graph::ANodeIt ANodeIt;
		typedef typename Graph::BNodeIt BNodeIt;
		typedef typename Graph::UEdge UEdge;
		typedef typename Graph::UEdgeIt UEdgeIt;
		typedef typename Graph::IncEdgeIt IncEdgeIt;
		
		typedef typename BpUGraph::template ANodeMap<UEdge> ANodeMatchingMap;
		typedef typename BpUGraph::template BNodeMap<UEdge> BNodeMatchingMap;
		
		
	public:
		
		/// \brief Constructor.
		///
		/// Constructor of the algorithm. 
		MaxBipartiteMatching(const BpUGraph& graph) 
		: _matching(graph), _rmatching(graph), _reached(graph), _graph(&graph) {}
		
		/// \name Execution control
		/// The simplest way to execute the algorithm is to use
		/// one of the member functions called \c run().
		/// \n
		/// If you need more control on the execution,
		/// first you must call \ref init() or one alternative for it.
		/// Finally \ref start() will perform the matching computation or
		/// with step-by-step execution you can augment the solution.
		
		/// @{
		
		/// \brief Initalize the data structures.
		///
		/// It initalizes the data structures and creates an empty matching.
		void init() {
			for (ANodeIt it(*_graph); it != INVALID; ++it) {
				_matching.set(it, INVALID);
			}
			for (BNodeIt it(*_graph); it != INVALID; ++it) {
				_rmatching.set(it, INVALID);
				_reached.set(it, -1);
			}
			_size = 0;
			_phase = -1;
		}
		
		/// \brief Initalize the data structures.
		///
		/// It initalizes the data structures and creates a greedy
		/// matching.  From this matching sometimes it is faster to get
		/// the matching than from the initial empty matching.
		void greedyInit() {
			_size = 0;
			for (BNodeIt it(*_graph); it != INVALID; ++it) {
				_rmatching.set(it, INVALID);
				_reached.set(it, 0);
			}
			for (ANodeIt it(*_graph); it != INVALID; ++it) {
				_matching[it] = INVALID;
				for (IncEdgeIt jt(*_graph, it); jt != INVALID; ++jt) {
					if (_rmatching[_graph->bNode(jt)] == INVALID) {
						_matching.set(it, jt);
						_rmatching.set(_graph->bNode(jt), jt);
						_reached.set(_graph->bNode(jt), -1);
						++_size;
						break;
					}
				}
			}
			_phase = 0;
		}
		
		/// \brief Initalize the data structures with an initial matching.
		///
		/// It initalizes the data structures with an initial matching.
		template <typename MatchingMap>
		void matchingInit(const MatchingMap& mm) {
			for (ANodeIt it(*_graph); it != INVALID; ++it) {
				_matching.set(it, INVALID);
			}
			for (BNodeIt it(*_graph); it != INVALID; ++it) {
				_rmatching.set(it, INVALID);
				_reached.set(it, 0);
			}
			_size = 0;
			for (UEdgeIt it(*_graph); it != INVALID; ++it) {
				if (mm[it]) {
					++_size;
					_matching.set(_graph->aNode(it), it);
					_rmatching.set(_graph->bNode(it), it);
					_reached.set(_graph->bNode(it), 0);
				}
			}
			_phase = 0;
		}
		
		/// \brief Initalize the data structures with an initial matching.
		///
		/// It initalizes the data structures with an initial matching.
		/// \return %True when the given map contains really a matching.
		template <typename MatchingMap>
		bool checkedMatchingInit(const MatchingMap& mm) {
			for (ANodeIt it(*_graph); it != INVALID; ++it) {
				_matching.set(it, INVALID);
			}
			for (BNodeIt it(*_graph); it != INVALID; ++it) {
				_rmatching.set(it, INVALID);
				_reached.set(it, 0);
			}
			_size = 0;
			for (UEdgeIt it(*_graph); it != INVALID; ++it) {
				if (mm[it]) {
					++_size;
					if (_matching[_graph->aNode(it)] != INVALID) {
						return false;
					}
					_matching.set(_graph->aNode(it), it);
					if (_matching[_graph->bNode(it)] != INVALID) {
						return false;
					}
					_matching.set(_graph->bNode(it), it);
					_reached.set(_graph->bNode(it), -1);
				}
			}
			_phase = 0;
			return true;
		}
		
	private:
		
		bool _find_path(Node anode, int maxlevel,
						typename Graph::template BNodeMap<int>& level) {
			for (IncEdgeIt it(*_graph, anode); it != INVALID; ++it) {
				Node bnode = _graph->bNode(it); 
				if (level[bnode] == maxlevel) {
					level.set(bnode, -1);
					if (maxlevel == 0) {
						_matching.set(anode, it);
						_rmatching.set(bnode, it);
						return true;
					} else {
						Node nnode = _graph->aNode(_rmatching[bnode]);
						if (_find_path(nnode, maxlevel - 1, level)) {
							_matching.set(anode, it);
							_rmatching.set(bnode, it);
							return true;
						}
					}
				}
			}
			return false;
		}
		
	public:
		
		/// \brief An augmenting phase of the Hopcroft-Karp algorithm
		///
		/// It runs an augmenting phase of the Hopcroft-Karp
		/// algorithm. This phase finds maximal edge disjoint augmenting
		/// paths and augments on these paths. The algorithm consists at
		/// most of \f$ O(\sqrt{n}) \f$ phase and one phase is \f$ O(e)
		/// \f$ long.
		bool augment() {
			
			++_phase;
			
			typename Graph::template BNodeMap<int> _level(*_graph, -1);
			//typename Graph::template ANodeMap<int> _found(*_graph, false);
			typename Graph::template ANodeMap<bool> _found(*_graph, false);
			std::vector<Node> queue, aqueue;
			for (BNodeIt it(*_graph); it != INVALID; ++it) {
				if (_rmatching[it] == INVALID) {
					queue.push_back(it);
					_reached.set(it, _phase);
					_level.set(it, 0);
				}
			}
			
			bool success = false;
			
			int level = 0;
			//std::vector<Node> nqueue;
			while (!success && !queue.empty()) {
				//nqueue.clear();
				std::vector<Node> nqueue;
				for (int i = 0; i < int(queue.size()); ++i) {
					Node bnode = queue[i];
					for (IncEdgeIt jt(*_graph, bnode); jt != INVALID; ++jt) {
						Node anode = _graph->aNode(jt);
						if (_matching[anode] == INVALID) {
							
							if (!_found[anode]) {
								if (_find_path(anode, level, _level)) {
									++_size;
								}
								_found.set(anode, true);
							}
							success = true;
						} else {           
							Node nnode = _graph->bNode(_matching[anode]);
							if (_reached[nnode] != _phase) {
								_reached.set(nnode, _phase);
								nqueue.push_back(nnode);
								_level.set(nnode, level + 1);
							}
						}
					}
				}
				++level;
				queue.swap(nqueue);
			}
			
			return success;
		}
	private:
		
		void _find_path_bfs(Node anode,
							typename Graph::template ANodeMap<UEdge>& pred) {
			while (true) {
				UEdge uedge = pred[anode];
				Node bnode = _graph->bNode(uedge);
				
				UEdge nedge = _rmatching[bnode];
				
				_matching.set(anode, uedge);
				_rmatching.set(bnode, uedge);
				
				if (nedge == INVALID) break;
				anode = _graph->aNode(nedge);
			}
		}
		
	public:
		
		/// \brief An augmenting phase with single path augementing
		///
		/// This phase finds only one augmenting paths and augments on
		/// these paths. The algorithm consists at most of \f$ O(n) \f$
		/// phase and one phase is \f$ O(e) \f$ long.
		bool simpleAugment() { 
			++_phase;
			
			typename Graph::template ANodeMap<UEdge> _pred(*_graph);
			
			std::vector<Node> queue, aqueue;
			for (BNodeIt it(*_graph); it != INVALID; ++it) {
				if (_rmatching[it] == INVALID) {
					queue.push_back(it);
					_reached.set(it, _phase);
				}
			}
			
			bool success = false;
			
			int level = 0;
			while (!success && !queue.empty()) {
				std::vector<Node> nqueue;
				for (int i = 0; i < int(queue.size()); ++i) {
					Node bnode = queue[i];
					for (IncEdgeIt jt(*_graph, bnode); jt != INVALID; ++jt) {
						Node anode = _graph->aNode(jt);
						if (_matching[anode] == INVALID) {
							_pred.set(anode, jt);
							_find_path_bfs(anode, _pred);
							++_size;
							return true;
						} else {           
							Node nnode = _graph->bNode(_matching[anode]);
							if (_reached[nnode] != _phase) {
								_pred.set(anode, jt);
								_reached.set(nnode, _phase);
								nqueue.push_back(nnode);
							}
						}
					}
				}
				++level;
				queue.swap(nqueue);
			}
			
			return success;
		}
		
		
		
		/// \brief Starts the algorithm.
		///
		/// Starts the algorithm. It runs augmenting phases until the optimal
		/// solution reached.
		void start() {
			while (augment()) {}
		}
		
		/// \brief Runs the algorithm.
		///
		/// It just initalize the algorithm and then start it.
		void run() {
			greedyInit();
			start();
		}
		
		/// @}
		
		/// \name Query Functions
		/// The result of the %Matching algorithm can be obtained using these
		/// functions.\n
		/// Before the use of these functions,
		/// either run() or start() must be called.
		
		///@{
		
		/// \brief Return true if the given uedge is in the matching.
		/// 
		/// It returns true if the given uedge is in the matching.
		bool matchingEdge(const UEdge& edge) const {
			return _matching[_graph->aNode(edge)] == edge;
		}
		
		/// \brief Returns the matching edge from the node.
		/// 
		/// Returns the matching edge from the node. If there is not such
		/// edge it gives back \c INVALID.
		/// \note If the parameter node is a B-node then the running time is
		/// propotional to the degree of the node.
		UEdge matchingEdge(const Node& node) const {
			if (_graph->aNode(node)) {
				return _matching[node];
			} else {
				return _rmatching[node];
			}
		}
		
		/// \brief Set true all matching uedge in the map.
		/// 
		/// Set true all matching uedge in the map. It does not change the
		/// value mapped to the other uedges.
		/// \return The number of the matching edges.
		template <typename MatchingMap>
		int quickMatching(MatchingMap& mm) const {
			for (ANodeIt it(*_graph); it != INVALID; ++it) {
				if (_matching[it] != INVALID) {
					mm.set(_matching[it], true);
				}
			}
			return _size;
		}
		
		/// \brief Set true all matching uedge in the map and the others to false.
		/// 
		/// Set true all matching uedge in the map and the others to false.
		/// \return The number of the matching edges.
		template <typename MatchingMap>
		int matching(MatchingMap& mm) const {
			for (UEdgeIt it(*_graph); it != INVALID; ++it) {
				mm.set(it, it == _matching[_graph->aNode(it)]);
			}
			return _size;
		}
		
		///Gives back the matching in an ANodeMap.
		
		///Gives back the matching in an ANodeMap. The parameter should
		///be a write ANodeMap of UEdge values.
		///\return The number of the matching edges.
		template<class MatchingMap>
		int aMatching(MatchingMap& mm) const {
			for (ANodeIt it(*_graph); it != INVALID; ++it) {
				mm.set(it, _matching[it]);
			}
			return _size;
		}
		
		///Gives back the matching in a BNodeMap.
		
		///Gives back the matching in a BNodeMap. The parameter should
		///be a write BNodeMap of UEdge values.
		///\return The number of the matching edges.
		template<class MatchingMap>
		int bMatching(MatchingMap& mm) const {
			for (BNodeIt it(*_graph); it != INVALID; ++it) {
				mm.set(it, _rmatching[it]);
			}
			return _size;
		}
		
		/// \brief Returns a minimum covering of the nodes.
		///
		/// The minimum covering set problem is the dual solution of the
		/// maximum bipartite matching. It provides a solution for this
		/// problem what is proof of the optimality of the matching.
		/// \return The size of the cover set.
		template <typename CoverMap>
		int coverSet(CoverMap& covering) const {
			
			int size = 0;
			for (ANodeIt it(*_graph); it != INVALID; ++it) {
				bool cn = _matching[it] != INVALID && 
				_reached[_graph->bNode(_matching[it])] == _phase;
				covering.set(it, cn);
				if (cn) ++size;
			}
			for (BNodeIt it(*_graph); it != INVALID; ++it) {
				bool cn = _reached[it] != _phase;
				covering.set(it, cn);
				if (cn) ++size;
			}
			return size;
		}
		
		/// \brief Gives back a barrier on the A-nodes
		///    
		/// The barrier is s subset of the nodes on the same side of the
		/// graph, which size minus its neighbours is exactly the
		/// unmatched nodes on the A-side.  
		/// \retval barrier A WriteMap on the ANodes with bool value.
		template <typename BarrierMap>
		void aBarrier(BarrierMap& barrier) const {
			
			for (ANodeIt it(*_graph); it != INVALID; ++it) {
				barrier.set(it, _matching[it] == INVALID || 
							_reached[_graph->bNode(_matching[it])] != _phase);
			}
		}
		
		/// \brief Gives back a barrier on the B-nodes
		///    
		/// The barrier is s subset of the nodes on the same side of the
		/// graph, which size minus its neighbours is exactly the
		/// unmatched nodes on the B-side.  
		/// \retval barrier A WriteMap on the BNodes with bool value.
		template <typename BarrierMap>
		void bBarrier(BarrierMap& barrier) const {
			
			for (BNodeIt it(*_graph); it != INVALID; ++it) {
				barrier.set(it, _reached[it] == _phase);
			}
		}
		
		/// \brief Gives back the number of the matching edges.
		///
		/// Gives back the number of the matching edges.
		int matchingSize() const {
			return _size;
		}
		
		/// @}
		
	private:
		
		typename BpUGraph::template ANodeMap<UEdge> _matching;
		typename BpUGraph::template BNodeMap<UEdge> _rmatching;
		
		typename BpUGraph::template BNodeMap<int> _reached;
		
		int _phase;
		const Graph *_graph;
		
		int _size;
		
	};
	
	/// \ingroup matching
	///
	/// \brief Maximum cardinality bipartite matching
	///
	/// This function calculates the maximum cardinality matching
	/// in a bipartite graph. It gives back the matching in an undirected
	/// edge map.
	///
	/// \param graph The bipartite graph.
	/// \return The size of the matching.
	template <typename BpUGraph>
	int maxBipartiteMatching(const BpUGraph& graph) {
		MaxBipartiteMatching<BpUGraph> bpmatching(graph);
		bpmatching.run();
		return bpmatching.matchingSize();
	}
	
	/// \ingroup matching
	///
	/// \brief Maximum cardinality bipartite matching
	///
	/// This function calculates the maximum cardinality matching
	/// in a bipartite graph. It gives back the matching in an undirected
	/// edge map.
	///
	/// \param graph The bipartite graph.
	/// \retval matching The ANodeMap of UEdges which will be set to covered
	/// matching undirected edge.
	/// \return The size of the matching.
	template <typename BpUGraph, typename MatchingMap>
	int maxBipartiteMatching(const BpUGraph& graph, MatchingMap& matching) {
		MaxBipartiteMatching<BpUGraph> bpmatching(graph);
		bpmatching.run();
		bpmatching.aMatching(matching);
		return bpmatching.matchingSize();
	}
	
	/// \ingroup matching
	///
	/// \brief Maximum cardinality bipartite matching
	///
	/// This function calculates the maximum cardinality matching
	/// in a bipartite graph. It gives back the matching in an undirected
	/// edge map.
	///
	/// \param graph The bipartite graph.
	/// \retval matching The ANodeMap of UEdges which will be set to covered
	/// matching undirected edge.
	/// \retval barrier The BNodeMap of bools which will be set to a barrier
	/// of the BNode-set.
	/// \return The size of the matching.
	template <typename BpUGraph, typename MatchingMap, typename BarrierMap>
	int maxBipartiteMatching(const BpUGraph& graph, 
							 MatchingMap& matching, BarrierMap& barrier) {
		MaxBipartiteMatching<BpUGraph> bpmatching(graph);
		bpmatching.run();
		bpmatching.aMatching(matching);
		bpmatching.bBarrier(barrier);
		return bpmatching.matchingSize();
	}
	
	/// \brief Default traits class for weighted bipartite matching algoritms.
	///
	/// Default traits class for weighted bipartite matching algoritms.
	/// \param _BpUGraph The bipartite undirected graph type.
	/// \param _WeightMap Type of weight map.
	template <typename _BpUGraph, typename _WeightMap>
	struct MaxWeightedBipartiteMatchingDefaultTraits {
		/// \brief The type of the weight of the undirected edges.
		typedef typename _WeightMap::Value Value;
		
		/// The undirected bipartite graph type the algorithm runs on. 
		typedef _BpUGraph BpUGraph;
		
		/// The map of the edges weights
		typedef _WeightMap WeightMap;
		
		/// \brief The cross reference type used by heap.
		///
		/// The cross reference type used by heap.
		/// Usually it is \c Graph::ANodeMap<int>.
		typedef typename BpUGraph::template ANodeMap<int> HeapCrossRef;
		
		/// \brief Instantiates a HeapCrossRef.
		///
		/// This function instantiates a \ref HeapCrossRef. 
		/// \param graph is the graph, to which we would like to define the 
		/// HeapCrossRef.
		static HeapCrossRef *createHeapCrossRef(const BpUGraph &graph) {
			return new HeapCrossRef(graph);
		}
		
		/// \brief The heap type used by weighted matching algorithms.
		///
		/// The heap type used by weighted matching algorithms. It should
		/// minimize the priorities and the heap's key type is the graph's
		/// anode graph's node.
		///
		/// \sa BinHeap
		//typedef BinHeap<Value, HeapCrossRef> Heap;
		typedef FibHeap<Value, HeapCrossRef> Heap;

		/// \brief Instantiates a Heap.
		///
		/// This function instantiates a \ref Heap. 
		/// \param crossref The cross reference of the heap.
		static Heap *createHeap(HeapCrossRef& crossref) {
			return new Heap(crossref);
		}
		
	};
	
	
	/// \ingroup matching
	///
	/// \brief Bipartite Max Weighted Matching algorithm
	///
	/// This class implements the bipartite Max Weighted Matching
	/// algorithm.  It uses the successive shortest path algorithm to
	/// calculate the maximum weighted matching in the bipartite
	/// graph. The algorithm can be used also to calculate the maximum
	/// cardinality maximum weighted matching. The time complexity
	/// of the algorithm is \f$ O(ne\log(n)) \f$ with the default binary
	/// heap implementation but this can be improved to 
	/// \f$ O(n^2\log(n)+ne) \f$ if we use fibonacci heaps.
	///
	/// The algorithm also provides a potential function on the nodes
	/// which a dual solution of the matching algorithm and it can be
	/// used to proof the optimality of the given pimal solution.
#ifdef DOXYGEN
	template <typename _BpUGraph, typename _WeightMap, typename _Traits>
#else
	template <typename _BpUGraph, 
	typename _WeightMap = typename _BpUGraph::template UEdgeMap<int>,
	typename _Traits = 
	MaxWeightedBipartiteMatchingDefaultTraits<_BpUGraph, _WeightMap> >
#endif
	class MaxWeightedBipartiteMatching {
public:
	
    typedef _Traits Traits;
    typedef typename Traits::BpUGraph BpUGraph;
    typedef typename Traits::WeightMap WeightMap;
    typedef typename Traits::Value Value;
	
protected:
	
    typedef typename Traits::HeapCrossRef HeapCrossRef;
    typedef typename Traits::Heap Heap; 
	
    
    typedef typename BpUGraph::Node Node;
    typedef typename BpUGraph::ANodeIt ANodeIt;
    typedef typename BpUGraph::BNodeIt BNodeIt;
    typedef typename BpUGraph::UEdge UEdge;
    typedef typename BpUGraph::UEdgeIt UEdgeIt;
    typedef typename BpUGraph::IncEdgeIt IncEdgeIt;
	
    typedef typename BpUGraph::template ANodeMap<UEdge> ANodeMatchingMap;
    typedef typename BpUGraph::template BNodeMap<UEdge> BNodeMatchingMap;
	
    typedef typename BpUGraph::template ANodeMap<Value> ANodePotentialMap;
    typedef typename BpUGraph::template BNodeMap<Value> BNodePotentialMap;
	
	
public:
	
    /// \brief \ref Exception for uninitialized parameters.
    ///
    /// This error represents problems in the initialization
    /// of the parameters of the algorithms.
    class UninitializedParameter : public lemon::UninitializedParameter {
    public:
		virtual const char* what() const throw() {
			return "lemon::MaxWeightedBipartiteMatching::UninitializedParameter";
		}
    };
	
    ///\name Named template parameters
	
    ///@{
	
    template <class H, class CR>
    struct DefHeapTraits : public Traits {
		typedef CR HeapCrossRef;
		typedef H Heap;
		static HeapCrossRef *createHeapCrossRef(const BpUGraph &) {
			throw UninitializedParameter();
		}
		static Heap *createHeap(HeapCrossRef &) {
			throw UninitializedParameter();
		}
    };
	
    /// \brief \ref named-templ-param "Named parameter" for setting heap 
    /// and cross reference type
    ///
    /// \ref named-templ-param "Named parameter" for setting heap and cross 
    /// reference type
    template <class H, class CR = typename BpUGraph::template NodeMap<int> >
    struct DefHeap
	: public MaxWeightedBipartiteMatching<BpUGraph, WeightMap, 
	DefHeapTraits<H, CR> > { 
	typedef MaxWeightedBipartiteMatching<BpUGraph, WeightMap, 
	DefHeapTraits<H, CR> > Create;
};

template <class H, class CR>
struct DefStandardHeapTraits : public Traits {
	typedef CR HeapCrossRef;
	typedef H Heap;
	static HeapCrossRef *createHeapCrossRef(const BpUGraph &graph) {
		return new HeapCrossRef(graph);
	}
	static Heap *createHeap(HeapCrossRef &crossref) {
		return new Heap(crossref);
	}
};

/// \brief \ref named-templ-param "Named parameter" for setting heap and 
/// cross reference type with automatic allocation
///
/// \ref named-templ-param "Named parameter" for setting heap and cross 
/// reference type. It can allocate the heap and the cross reference 
/// object if the cross reference's constructor waits for the graph as 
/// parameter and the heap's constructor waits for the cross reference.
template <class H, class CR = typename BpUGraph::template NodeMap<int> >
struct DefStandardHeap
: public MaxWeightedBipartiteMatching<BpUGraph, WeightMap, 
DefStandardHeapTraits<H, CR> > { 
typedef MaxWeightedBipartiteMatching<BpUGraph, WeightMap, 
DefStandardHeapTraits<H, CR> > 
Create;
};

///@}


/// \brief Constructor.
///
/// Constructor of the algorithm. 
MaxWeightedBipartiteMatching(const BpUGraph& _graph, 
							 const WeightMap& _weight) 
: graph(&_graph), weight(&_weight),
anode_matching(_graph), bnode_matching(_graph),
anode_potential(_graph), bnode_potential(_graph),
_heap_cross_ref(0), local_heap_cross_ref(false),
_heap(0), local_heap(0) {}

/// \brief Destructor.
///
/// Destructor of the algorithm.
~MaxWeightedBipartiteMatching() {
	destroyStructures();
}

/// \brief Sets the heap and the cross reference used by algorithm.
///
/// Sets the heap and the cross reference used by algorithm.
/// If you don't use this function before calling \ref run(),
/// it will allocate one. The destuctor deallocates this
/// automatically allocated map, of course.
/// \return \c (*this)
MaxWeightedBipartiteMatching& heap(Heap& hp, HeapCrossRef &cr) {
	if(local_heap_cross_ref) {
		delete _heap_cross_ref;
		local_heap_cross_ref = false;
	}
	_heap_cross_ref = &cr;
	if(local_heap) {
		delete _heap;
		local_heap = false;
	}
	_heap = &hp;
	return *this;
}

/// \name Execution control
/// The simplest way to execute the algorithm is to use
/// one of the member functions called \c run().
/// \n
/// If you need more control on the execution,
/// first you must call \ref init() or one alternative for it.
/// Finally \ref start() will perform the matching computation or
/// with step-by-step execution you can augment the solution.

/// @{

/// \brief Initalize the data structures.
///
/// It initalizes the data structures and creates an empty matching.
void init() {
	initStructures();
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        anode_matching[it] = INVALID;
        anode_potential[it] = 0;
	}
	for (BNodeIt it(*graph); it != INVALID; ++it) {
        bnode_matching[it] = INVALID;
        bnode_potential[it] = 0;
        for (IncEdgeIt jt(*graph, it); jt != INVALID; ++jt) {
			if ((*weight)[jt] > bnode_potential[it]) {
				bnode_potential[it] = (*weight)[jt];
			}
        }
	}
	matching_value = 0;
	matching_size = 0;
}


/// \brief An augmenting phase of the weighted matching algorithm
///
/// It runs an augmenting phase of the weighted matching 
/// algorithm. This phase finds the best augmenting path and 
/// augments only on this paths. 
///
/// The algorithm consists at most 
/// of \f$ O(n) \f$ phase and one phase is \f$ O(n\log(n)+e) \f$ 
/// long with Fibonacci heap or \f$ O((n+e)\log(n)) \f$ long 
/// with binary heap.
/// \param decrease If the given parameter true the matching value
/// can be decreased in the augmenting phase. If we would like
/// to calculate the maximum cardinality maximum weighted matching
/// then we should let the algorithm to decrease the matching
/// value in order to increase the number of the matching edges.
bool augment(bool decrease = false) {
	
	typename BpUGraph::template BNodeMap<Value> bdist(*graph);
	typename BpUGraph::template BNodeMap<UEdge> bpred(*graph, INVALID);
	
	Node bestNode = INVALID;
	Value bestValue = 0;
	
	_heap->clear();
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        (*_heap_cross_ref)[it] = Heap::PRE_HEAP;
	}
	
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        if (anode_matching[it] == INVALID) {
			_heap->push(it, 0);
        }
	}
	
	Value bdistMax = 0;
	while (!_heap->empty()) {
        Node anode = _heap->top();
        Value avalue = _heap->prio();
        _heap->pop();
        for (IncEdgeIt jt(*graph, anode); jt != INVALID; ++jt) {
			if (jt == anode_matching[anode]) continue;
			Node bnode = graph->bNode(jt);
			Value bvalue = avalue  - (*weight)[jt] +
            anode_potential[anode] + bnode_potential[bnode];
			if (bvalue > bdistMax) {
				bdistMax = bvalue;
			}
			if (bpred[bnode] == INVALID || bvalue < bdist[bnode]) {
				bdist[bnode] = bvalue;
				bpred[bnode] = jt;
			} else continue;
			if (bnode_matching[bnode] != INVALID) {
				Node newanode = graph->aNode(bnode_matching[bnode]);
				switch (_heap->state(newanode)) {
					case Heap::PRE_HEAP:
						_heap->push(newanode, bvalue);
						break;
					case Heap::IN_HEAP:
						if (bvalue < (*_heap)[newanode]) {
							_heap->decrease(newanode, bvalue);
						}
						break;
					case Heap::POST_HEAP:
						break;
				}
			} else {
				if (bestNode == INVALID || 
					bnode_potential[bnode] - bvalue > bestValue) {
					bestValue = bnode_potential[bnode] - bvalue;
					bestNode = bnode;
				}
			}
        }
	}
	
	if (bestNode == INVALID || (!decrease && bestValue < 0)) {
        return false;
	}
	
	matching_value += bestValue;
	++matching_size;
	
	for (BNodeIt it(*graph); it != INVALID; ++it) {
        if (bpred[it] != INVALID) {
			bnode_potential[it] -= bdist[it];
        } else {
			bnode_potential[it] -= bdistMax;
        }
	}
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        if (anode_matching[it] != INVALID) {
			Node bnode = graph->bNode(anode_matching[it]);
			if (bpred[bnode] != INVALID) {
				anode_potential[it] += bdist[bnode];
			} else {
				anode_potential[it] += bdistMax;
			}
        }
	}
	
	while (bestNode != INVALID) {
        UEdge uedge = bpred[bestNode];
        Node anode = graph->aNode(uedge);
        
        bnode_matching[bestNode] = uedge;
        if (anode_matching[anode] != INVALID) {
			bestNode = graph->bNode(anode_matching[anode]);
        } else {
			bestNode = INVALID;
        }
        anode_matching[anode] = uedge;
	}
	
	
	return true;
}

/// \brief Starts the algorithm.
///
/// Starts the algorithm. It runs augmenting phases until the
/// optimal solution reached.
///
/// \param maxCardinality If the given value is true it will
/// calculate the maximum cardinality maximum matching instead of
/// the maximum matching.
void start(bool maxCardinality = false) {
	while (augment(maxCardinality)) {}
}

/// \brief Runs the algorithm.
///
/// It just initalize the algorithm and then start it.
///
/// \param maxCardinality If the given value is true it will
/// calculate the maximum cardinality maximum matching instead of
/// the maximum matching.
void run(bool maxCardinality = false) {
	init();
	start(maxCardinality);
}

/// @}

/// \name Query Functions
/// The result of the %Matching algorithm can be obtained using these
/// functions.\n
/// Before the use of these functions,
/// either run() or start() must be called.

///@{

/// \brief Gives back the potential in the NodeMap
///
/// Gives back the potential in the NodeMap. The matching is optimal
/// with the current number of edges if \f$ \pi(a) + \pi(b) - w(ab) = 0 \f$
/// for each matching edges and \f$ \pi(a) + \pi(b) - w(ab) \ge 0 \f$
/// for each edges. 
template <typename PotentialMap>
void potential(PotentialMap& pt) const {
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        pt.set(it, anode_potential[it]);
	}
	for (BNodeIt it(*graph); it != INVALID; ++it) {
        pt.set(it, bnode_potential[it]);
	}
}

/// \brief Set true all matching uedge in the map.
/// 
/// Set true all matching uedge in the map. It does not change the
/// value mapped to the other uedges.
/// \return The number of the matching edges.
template <typename MatchingMap>
int quickMatching(MatchingMap& mm) const {
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        if (anode_matching[it] != INVALID) {
			mm.set(anode_matching[it], true);
        }
	}
	return matching_size;
}

/// \brief Set true all matching uedge in the map and the others to false.
/// 
/// Set true all matching uedge in the map and the others to false.
/// \return The number of the matching edges.
template <typename MatchingMap>
int matching(MatchingMap& mm) const {
	for (UEdgeIt it(*graph); it != INVALID; ++it) {
        mm.set(it, it == anode_matching[graph->aNode(it)]);
	}
	return matching_size;
}

///Gives back the matching in an ANodeMap.

///Gives back the matching in an ANodeMap. The parameter should
///be a write ANodeMap of UEdge values.
///\return The number of the matching edges.
template<class MatchingMap>
int aMatching(MatchingMap& mm) const {
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        mm.set(it, anode_matching[it]);
	}
	return matching_size;
}

///Gives back the matching in a BNodeMap.

///Gives back the matching in a BNodeMap. The parameter should
///be a write BNodeMap of UEdge values.
///\return The number of the matching edges.
template<class MatchingMap>
int bMatching(MatchingMap& mm) const {
	for (BNodeIt it(*graph); it != INVALID; ++it) {
        mm.set(it, bnode_matching[it]);
	}
	return matching_size;
}


/// \brief Return true if the given uedge is in the matching.
/// 
/// It returns true if the given uedge is in the matching.
bool matchingEdge(const UEdge& edge) const {
	return anode_matching[graph->aNode(edge)] == edge;
}

/// \brief Returns the matching edge from the node.
/// 
/// Returns the matching edge from the node. If there is not such
/// edge it gives back \c INVALID.
UEdge matchingEdge(const Node& node) const {
	if (graph->aNode(node)) {
        return anode_matching[node];
	} else {
        return bnode_matching[node];
	}
}

/// \brief Gives back the sum of weights of the matching edges.
///
/// Gives back the sum of weights of the matching edges.
Value matchingValue() const {
	return matching_value;
}

/// \brief Gives back the number of the matching edges.
///
/// Gives back the number of the matching edges.
int matchingSize() const {
	return matching_size;
}

/// @}

private:

void initStructures() {
	if (!_heap_cross_ref) {
		local_heap_cross_ref = true;
		_heap_cross_ref = Traits::createHeapCrossRef(*graph);
	}
	if (!_heap) {
		local_heap = true;
		_heap = Traits::createHeap(*_heap_cross_ref);
	}
}

void destroyStructures() {
	if (local_heap_cross_ref) delete _heap_cross_ref;
	if (local_heap) delete _heap;
}


private:

const BpUGraph *graph;
const WeightMap* weight;

ANodeMatchingMap anode_matching;
BNodeMatchingMap bnode_matching;

ANodePotentialMap anode_potential;
BNodePotentialMap bnode_potential;

Value matching_value;
int matching_size;

HeapCrossRef *_heap_cross_ref;
bool local_heap_cross_ref;

Heap *_heap;
bool local_heap;

};

/// \ingroup matching
///
/// \brief Maximum weighted bipartite matching
///
/// This function calculates the maximum weighted matching
/// in a bipartite graph. It gives back the matching in an undirected
/// edge map.
///
/// \param graph The bipartite graph.
/// \param weight The undirected edge map which contains the weights.
/// \retval matching The undirected edge map which will be set to 
/// the matching.
/// \return The value of the matching.
template <typename BpUGraph, typename WeightMap, typename MatchingMap>
typename WeightMap::Value 
maxWeightedBipartiteMatching(const BpUGraph& graph, const WeightMap& weight,
							 MatchingMap& matching) {
    MaxWeightedBipartiteMatching<BpUGraph, WeightMap> 
	bpmatching(graph, weight);
    bpmatching.run();
    bpmatching.matching(matching);
    return bpmatching.matchingValue();
}

/// \ingroup matching
///
/// \brief Maximum weighted maximum cardinality bipartite matching
///
/// This function calculates the maximum weighted of the maximum cardinality
/// matchings of a bipartite graph. It gives back the matching in an 
/// undirected edge map.
///
/// \param graph The bipartite graph.
/// \param weight The undirected edge map which contains the weights.
/// \retval matching The undirected edge map which will be set to 
/// the matching.
/// \return The value of the matching.
template <typename BpUGraph, typename WeightMap, typename MatchingMap>
typename WeightMap::Value 
maxWeightedMaxBipartiteMatching(const BpUGraph& graph, 
								const WeightMap& weight,
								MatchingMap& matching) {
    MaxWeightedBipartiteMatching<BpUGraph, WeightMap> 
	bpmatching(graph, weight);
    bpmatching.run(true);
    bpmatching.matching(matching);
    return bpmatching.matchingValue();
}

/// \brief Default traits class for minimum cost bipartite matching
/// algoritms.
///
/// Default traits class for minimum cost bipartite matching
/// algoritms.  
///
/// \param _BpUGraph The bipartite undirected graph
/// type.  
///
/// \param _CostMap Type of cost map.
template <typename _BpUGraph, typename _CostMap>
struct MinCostMaxBipartiteMatchingDefaultTraits {
    /// \brief The type of the cost of the undirected edges.
    typedef typename _CostMap::Value Value;
	
    /// The undirected bipartite graph type the algorithm runs on. 
    typedef _BpUGraph BpUGraph;
	
    /// The map of the edges costs
    typedef _CostMap CostMap;
	
    /// \brief The cross reference type used by heap.
    ///
    /// The cross reference type used by heap.
    /// Usually it is \c Graph::NodeMap<int>.
    typedef typename BpUGraph::template NodeMap<int> HeapCrossRef;
	
    /// \brief Instantiates a HeapCrossRef.
    ///
    /// This function instantiates a \ref HeapCrossRef. 
    /// \param graph is the graph, to which we would like to define the 
    /// HeapCrossRef.
    static HeapCrossRef *createHeapCrossRef(const BpUGraph &graph) {
		return new HeapCrossRef(graph);
    }
    
    /// \brief The heap type used by costed matching algorithms.
    ///
    /// The heap type used by costed matching algorithms. It should
    /// minimize the priorities and the heap's key type is the graph's
    /// anode graph's node.
    ///
    /// \sa BinHeap
    //typedef BinHeap<Value, HeapCrossRef> Heap;
	typedef FibHeap<Value, HeapCrossRef> Heap;
    
    /// \brief Instantiates a Heap.
    ///
    /// This function instantiates a \ref Heap. 
    /// \param crossref The cross reference of the heap.
    static Heap *createHeap(HeapCrossRef& crossref) {
		return new Heap(crossref);
    }
	
};


/// \ingroup matching
///
/// \brief Bipartite Min Cost Matching algorithm
///
/// This class implements the bipartite Min Cost Matching algorithm.
/// It uses the successive shortest path algorithm to calculate the
/// minimum cost maximum matching in the bipartite graph. The time
/// complexity of the algorithm is \f$ O(ne\log(n)) \f$ with the
/// default binary heap implementation but this can be improved to
/// \f$ O(n^2\log(n)+ne) \f$ if we use fibonacci heaps.
///
/// The algorithm also provides a potential function on the nodes
/// which a dual solution of the matching algorithm and it can be
/// used to proof the optimality of the given pimal solution.
#ifdef DOXYGEN
template <typename _BpUGraph, typename _CostMap, typename _Traits>
#else
template <typename _BpUGraph, 
typename _CostMap = typename _BpUGraph::template UEdgeMap<int>,
typename _Traits = 
MinCostMaxBipartiteMatchingDefaultTraits<_BpUGraph, _CostMap> >
#endif
class MinCostMaxBipartiteMatching {
public:

typedef _Traits Traits;
typedef typename Traits::BpUGraph BpUGraph;
typedef typename Traits::CostMap CostMap;
typedef typename Traits::Value Value;

protected:

typedef typename Traits::HeapCrossRef HeapCrossRef;
typedef typename Traits::Heap Heap; 


typedef typename BpUGraph::Node Node;
typedef typename BpUGraph::ANodeIt ANodeIt;
typedef typename BpUGraph::BNodeIt BNodeIt;
typedef typename BpUGraph::UEdge UEdge;
typedef typename BpUGraph::UEdgeIt UEdgeIt;
typedef typename BpUGraph::IncEdgeIt IncEdgeIt;

typedef typename BpUGraph::template ANodeMap<UEdge> ANodeMatchingMap;
typedef typename BpUGraph::template BNodeMap<UEdge> BNodeMatchingMap;

typedef typename BpUGraph::template ANodeMap<Value> ANodePotentialMap;
typedef typename BpUGraph::template BNodeMap<Value> BNodePotentialMap;


public:

/// \brief \ref Exception for uninitialized parameters.
///
/// This error represents problems in the initialization
/// of the parameters of the algorithms.
class UninitializedParameter : public lemon::UninitializedParameter {
public:
	virtual const char* what() const throw() {
		return "lemon::MinCostMaxBipartiteMatching::UninitializedParameter";
	}
};

///\name Named template parameters

///@{

template <class H, class CR>
struct DefHeapTraits : public Traits {
	typedef CR HeapCrossRef;
	typedef H Heap;
	static HeapCrossRef *createHeapCrossRef(const BpUGraph &) {
		throw UninitializedParameter();
	}
	static Heap *createHeap(HeapCrossRef &) {
		throw UninitializedParameter();
	}
};

/// \brief \ref named-templ-param "Named parameter" for setting heap 
/// and cross reference type
///
/// \ref named-templ-param "Named parameter" for setting heap and cross 
/// reference type
template <class H, class CR = typename BpUGraph::template NodeMap<int> >
struct DefHeap
: public MinCostMaxBipartiteMatching<BpUGraph, CostMap, 
DefHeapTraits<H, CR> > { 
typedef MinCostMaxBipartiteMatching<BpUGraph, CostMap, 
DefHeapTraits<H, CR> > Create;
};

template <class H, class CR>
struct DefStandardHeapTraits : public Traits {
	typedef CR HeapCrossRef;
	typedef H Heap;
	static HeapCrossRef *createHeapCrossRef(const BpUGraph &graph) {
		return new HeapCrossRef(graph);
	}
	static Heap *createHeap(HeapCrossRef &crossref) {
		return new Heap(crossref);
	}
};

/// \brief \ref named-templ-param "Named parameter" for setting heap and 
/// cross reference type with automatic allocation
///
/// \ref named-templ-param "Named parameter" for setting heap and cross 
/// reference type. It can allocate the heap and the cross reference 
/// object if the cross reference's constructor waits for the graph as 
/// parameter and the heap's constructor waits for the cross reference.
template <class H, class CR = typename BpUGraph::template NodeMap<int> >
struct DefStandardHeap
: public MinCostMaxBipartiteMatching<BpUGraph, CostMap, 
DefStandardHeapTraits<H, CR> > { 
typedef MinCostMaxBipartiteMatching<BpUGraph, CostMap, 
DefStandardHeapTraits<H, CR> > 
Create;
};

///@}


/// \brief Constructor.
///
/// Constructor of the algorithm. 
MinCostMaxBipartiteMatching(const BpUGraph& _graph, 
							const CostMap& _cost) 
: graph(&_graph), cost(&_cost),
anode_matching(_graph), bnode_matching(_graph),
anode_potential(_graph), bnode_potential(_graph),
_heap_cross_ref(0), local_heap_cross_ref(false),
_heap(0), local_heap(0) {}

/// \brief Destructor.
///
/// Destructor of the algorithm.
~MinCostMaxBipartiteMatching() {
	destroyStructures();
}

/// \brief Sets the heap and the cross reference used by algorithm.
///
/// Sets the heap and the cross reference used by algorithm.
/// If you don't use this function before calling \ref run(),
/// it will allocate one. The destuctor deallocates this
/// automatically allocated map, of course.
/// \return \c (*this)
MinCostMaxBipartiteMatching& heap(Heap& hp, HeapCrossRef &cr) {
	if(local_heap_cross_ref) {
		delete _heap_cross_ref;
		local_heap_cross_ref = false;
	}
	_heap_cross_ref = &cr;
	if(local_heap) {
		delete _heap;
		local_heap = false;
	}
	_heap = &hp;
	return *this;
}

/// \name Execution control
/// The simplest way to execute the algorithm is to use
/// one of the member functions called \c run().
/// \n
/// If you need more control on the execution,
/// first you must call \ref init() or one alternative for it.
/// Finally \ref start() will perform the matching computation or
/// with step-by-step execution you can augment the solution.

/// @{

/// \brief Initalize the data structures.
///
/// It initalizes the data structures and creates an empty matching.
void init() {
	initStructures();
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        anode_matching[it] = INVALID;
        anode_potential[it] = 0;
	}
	for (BNodeIt it(*graph); it != INVALID; ++it) {
        bnode_matching[it] = INVALID;
        bnode_potential[it] = 0;
	}
	matching_cost = 0;
	matching_size = 0;
}


/// \brief An augmenting phase of the costed matching algorithm
///
/// It runs an augmenting phase of the matching algorithm. The
/// phase finds the best augmenting path and augments only on this
/// paths.
///
/// The algorithm consists at most 
/// of \f$ O(n) \f$ phase and one phase is \f$ O(n\log(n)+e) \f$ 
/// long with Fibonacci heap or \f$ O((n+e)\log(n)) \f$ long 
/// with binary heap.
bool augment() {
	
	typename BpUGraph::template BNodeMap<Value> bdist(*graph);
	typename BpUGraph::template BNodeMap<UEdge> bpred(*graph, INVALID);
	
	Node bestNode = INVALID;
	Value bestValue = 0;
	
	_heap->clear();
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        (*_heap_cross_ref)[it] = Heap::PRE_HEAP;
	}
	
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        if (anode_matching[it] == INVALID) {
			_heap->push(it, 0);
        }
	}
	Value bdistMax = 0;
	
	while (!_heap->empty()) {
        Node anode = _heap->top();
        Value avalue = _heap->prio();
        _heap->pop();
        for (IncEdgeIt jt(*graph, anode); jt != INVALID; ++jt) {
			if (jt == anode_matching[anode]) continue;
			Node bnode = graph->bNode(jt);
			Value bvalue = avalue + (*cost)[jt] + 
            anode_potential[anode] - bnode_potential[bnode];
			if (bvalue > bdistMax) {
				bdistMax = bvalue;
			}
			if (bpred[bnode] == INVALID || bvalue < bdist[bnode]) {
				bdist[bnode] = bvalue;
				bpred[bnode] = jt;
			} else continue;
			if (bnode_matching[bnode] != INVALID) {
				Node newanode = graph->aNode(bnode_matching[bnode]);
				switch (_heap->state(newanode)) {
					case Heap::PRE_HEAP:
						_heap->push(newanode, bvalue);
						break;
					case Heap::IN_HEAP:
						if (bvalue < (*_heap)[newanode]) {
							_heap->decrease(newanode, bvalue);
						}
						break;
					case Heap::POST_HEAP:
						break;
				}
			} else {
				if (bestNode == INVALID || 
					bvalue + bnode_potential[bnode] < bestValue) {
					bestValue = bvalue + bnode_potential[bnode];
					bestNode = bnode;
				}
			}
        }
	}
	
	if (bestNode == INVALID) {
        return false;
	}
	
	matching_cost += bestValue;
	++matching_size;
	
	for (BNodeIt it(*graph); it != INVALID; ++it) {
        if (bpred[it] != INVALID) {
			bnode_potential[it] += bdist[it];
        } else {
			bnode_potential[it] += bdistMax;
        }
	}
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        if (anode_matching[it] != INVALID) {
			Node bnode = graph->bNode(anode_matching[it]);
			if (bpred[bnode] != INVALID) {
				anode_potential[it] += bdist[bnode];
			} else {
				anode_potential[it] += bdistMax;
			}
        }
	}
	
	while (bestNode != INVALID) {
        UEdge uedge = bpred[bestNode];
        Node anode = graph->aNode(uedge);
        
        bnode_matching[bestNode] = uedge;
        if (anode_matching[anode] != INVALID) {
			bestNode = graph->bNode(anode_matching[anode]);
        } else {
			bestNode = INVALID;
        }
        anode_matching[anode] = uedge;
	}
	
	
	return true;
}

/// \brief Starts the algorithm.
///
/// Starts the algorithm. It runs augmenting phases until the
/// optimal solution reached.
void start() {
	while (augment()) {}
}

/// \brief Runs the algorithm.
///
/// It just initalize the algorithm and then start it.
void run() {
	init();
	start();
}

/// @}

/// \name Query Functions
/// The result of the %Matching algorithm can be obtained using these
/// functions.\n
/// Before the use of these functions,
/// either run() or start() must be called.

///@{

/// \brief Gives back the potential in the NodeMap
///
/// Gives back the potential in the NodeMap. The matching is optimal
/// with the current number of edges if \f$ \pi(a) + \pi(b) - w(ab) = 0 \f$
/// for each matching edges and \f$ \pi(a) + \pi(b) - w(ab) \ge 0 \f$
/// for each edges. 
template <typename PotentialMap>
void potential(PotentialMap& pt) const {
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        pt.set(it, anode_potential[it]);
	}
	for (BNodeIt it(*graph); it != INVALID; ++it) {
        pt.set(it, bnode_potential[it]);
	}
}

/// \brief Set true all matching uedge in the map.
/// 
/// Set true all matching uedge in the map. It does not change the
/// value mapped to the other uedges.
/// \return The number of the matching edges.
template <typename MatchingMap>
int quickMatching(MatchingMap& mm) const {
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        if (anode_matching[it] != INVALID) {
			mm.set(anode_matching[it], true);
        }
	}
	return matching_size;
}

/// \brief Set true all matching uedge in the map and the others to false.
/// 
/// Set true all matching uedge in the map and the others to false.
/// \return The number of the matching edges.
template <typename MatchingMap>
int matching(MatchingMap& mm) const {
	for (UEdgeIt it(*graph); it != INVALID; ++it) {
        mm.set(it, it == anode_matching[graph->aNode(it)]);
	}
	return matching_size;
}

/// \brief Gives back the matching in an ANodeMap.
///
/// Gives back the matching in an ANodeMap. The parameter should
/// be a write ANodeMap of UEdge values.
/// \return The number of the matching edges.
template<class MatchingMap>
int aMatching(MatchingMap& mm) const {
	for (ANodeIt it(*graph); it != INVALID; ++it) {
        mm.set(it, anode_matching[it]);
	}
	return matching_size;
}

/// \brief Gives back the matching in a BNodeMap.
///
/// Gives back the matching in a BNodeMap. The parameter should
/// be a write BNodeMap of UEdge values.
/// \return The number of the matching edges.
template<class MatchingMap>
int bMatching(MatchingMap& mm) const {
	for (BNodeIt it(*graph); it != INVALID; ++it) {
        mm.set(it, bnode_matching[it]);
	}
	return matching_size;
}

/// \brief Return true if the given uedge is in the matching.
/// 
/// It returns true if the given uedge is in the matching.
bool matchingEdge(const UEdge& edge) const {
	return anode_matching[graph->aNode(edge)] == edge;
}

/// \brief Returns the matching edge from the node.
/// 
/// Returns the matching edge from the node. If there is not such
/// edge it gives back \c INVALID.
UEdge matchingEdge(const Node& node) const {
	if (graph->aNode(node)) {
        return anode_matching[node];
	} else {
        return bnode_matching[node];
	}
}

/// \brief Gives back the sum of costs of the matching edges.
///
/// Gives back the sum of costs of the matching edges.
Value matchingCost() const {
	return matching_cost;
}

/// \brief Gives back the number of the matching edges.
///
/// Gives back the number of the matching edges.
int matchingSize() const {
	return matching_size;
}

/// @}

private:

void initStructures() {
	if (!_heap_cross_ref) {
		local_heap_cross_ref = true;
		_heap_cross_ref = Traits::createHeapCrossRef(*graph);
	}
	if (!_heap) {
		local_heap = true;
		_heap = Traits::createHeap(*_heap_cross_ref);
	}
}

void destroyStructures() {
	if (local_heap_cross_ref) delete _heap_cross_ref;
	if (local_heap) delete _heap;
}


private:

const BpUGraph *graph;
const CostMap* cost;

ANodeMatchingMap anode_matching;
BNodeMatchingMap bnode_matching;

ANodePotentialMap anode_potential;
BNodePotentialMap bnode_potential;

Value matching_cost;
int matching_size;

HeapCrossRef *_heap_cross_ref;
bool local_heap_cross_ref;

Heap *_heap;
bool local_heap;

};

/// \ingroup matching
///
/// \brief Minimum cost maximum cardinality bipartite matching
///
/// This function calculates the maximum cardinality matching with
/// minimum cost of a bipartite graph. It gives back the matching in
/// an undirected edge map.
///
/// \param graph The bipartite graph.
/// \param cost The undirected edge map which contains the costs.
/// \retval matching The undirected edge map which will be set to 
/// the matching.
/// \return The cost of the matching.
template <typename BpUGraph, typename CostMap, typename MatchingMap>
typename CostMap::Value 
minCostMaxBipartiteMatching(const BpUGraph& graph, 
							const CostMap& cost,
							MatchingMap& matching) {
    MinCostMaxBipartiteMatching<BpUGraph, CostMap> 
	bpmatching(graph, cost);
    bpmatching.run();
    bpmatching.matching(matching);
    return bpmatching.matchingCost();
}

}

#endif

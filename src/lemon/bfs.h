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

#ifndef LEMON_BFS_H
#define LEMON_BFS_H

///\ingroup search
///\file
///\brief Bfs algorithm.

#include <lemon/list_graph.h>
#include <lemon/graph_utils.h>
#include <lemon/bits/path_dump.h>
#include <lemon/bits/invalid.h>
#include <lemon/error.h>
#include <lemon/maps.h>

namespace lemon {


  
  ///Default traits class of Bfs class.

  ///Default traits class of Bfs class.
  ///\param GR Graph type.
  template<class GR>
  struct BfsDefaultTraits
  {
    ///The graph type the algorithm runs on. 
    typedef GR Graph;
    ///\brief The type of the map that stores the last
    ///edges of the shortest paths.
    /// 
    ///The type of the map that stores the last
    ///edges of the shortest paths.
    ///It must meet the \ref concepts::WriteMap "WriteMap" concept.
    ///
    typedef typename Graph::template NodeMap<typename GR::Edge> PredMap;
    ///Instantiates a PredMap.
 
    ///This function instantiates a \ref PredMap. 
    ///\param G is the graph, to which we would like to define the PredMap.
    ///\todo The graph alone may be insufficient to initialize
    static PredMap *createPredMap(const GR &G) 
    {
      return new PredMap(G);
    }
    ///The type of the map that indicates which nodes are processed.
 
    ///The type of the map that indicates which nodes are processed.
    ///It must meet the \ref concepts::WriteMap "WriteMap" concept.
    ///\todo named parameter to set this type, function to read and write.
    typedef NullMap<typename Graph::Node,bool> ProcessedMap;
    ///Instantiates a ProcessedMap.
 
    ///This function instantiates a \ref ProcessedMap. 
    ///\param g is the graph, to which
    ///we would like to define the \ref ProcessedMap
#ifdef DOXYGEN
    static ProcessedMap *createProcessedMap(const GR &g)
#else
    static ProcessedMap *createProcessedMap(const GR &)
#endif
    {
      return new ProcessedMap();
    }
    ///The type of the map that indicates which nodes are reached.
 
    ///The type of the map that indicates which nodes are reached.
    ///It must meet the \ref concepts::WriteMap "WriteMap" concept.
    ///\todo named parameter to set this type, function to read and write.
    typedef typename Graph::template NodeMap<bool> ReachedMap;
    ///Instantiates a ReachedMap.
 
    ///This function instantiates a \ref ReachedMap. 
    ///\param G is the graph, to which
    ///we would like to define the \ref ReachedMap.
    static ReachedMap *createReachedMap(const GR &G)
    {
      return new ReachedMap(G);
    }
    ///The type of the map that stores the dists of the nodes.
 
    ///The type of the map that stores the dists of the nodes.
    ///It must meet the \ref concepts::WriteMap "WriteMap" concept.
    ///
    typedef typename Graph::template NodeMap<int> DistMap;
    ///Instantiates a DistMap.
 
    ///This function instantiates a \ref DistMap. 
    ///\param G is the graph, to which we would like to define the \ref DistMap
    static DistMap *createDistMap(const GR &G)
    {
      return new DistMap(G);
    }
  };
  
  ///%BFS algorithm class.
  
  ///\ingroup search
  ///This class provides an efficient implementation of the %BFS algorithm.
  ///
  ///\param GR The graph type the algorithm runs on. The default value is
  ///\ref ListGraph. The value of GR is not used directly by Bfs, it
  ///is only passed to \ref BfsDefaultTraits.
  ///\param TR Traits class to set various data types used by the algorithm.
  ///The default traits class is
  ///\ref BfsDefaultTraits "BfsDefaultTraits<GR>".
  ///See \ref BfsDefaultTraits for the documentation of
  ///a Bfs traits class.
  ///
  ///\author Alpar Juttner

#ifdef DOXYGEN
  template <typename GR,
	    typename TR>
#else
  template <typename GR=ListGraph,
	    typename TR=BfsDefaultTraits<GR> >
#endif
  class Bfs {
  public:
    /**
     * \brief \ref Exception for uninitialized parameters.
     *
     * This error represents problems in the initialization
     * of the parameters of the algorithms.
     */
    class UninitializedParameter : public lemon::UninitializedParameter {
    public:
      virtual const char* what() const throw() {
	return "lemon::Bfs::UninitializedParameter";
      }
    };

    typedef TR Traits;
    ///The type of the underlying graph.
    typedef typename TR::Graph Graph;
    
    ///\brief The type of the map that stores the last
    ///edges of the shortest paths.
    typedef typename TR::PredMap PredMap;
    ///The type of the map indicating which nodes are reached.
    typedef typename TR::ReachedMap ReachedMap;
    ///The type of the map indicating which nodes are processed.
    typedef typename TR::ProcessedMap ProcessedMap;
    ///The type of the map that stores the dists of the nodes.
    typedef typename TR::DistMap DistMap;
  private:

    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::OutEdgeIt OutEdgeIt;

    /// Pointer to the underlying graph.
    const Graph *G;
    ///Pointer to the map of predecessors edges.
    PredMap *_pred;
    ///Indicates if \ref _pred is locally allocated (\c true) or not.
    bool local_pred;
    ///Pointer to the map of distances.
    DistMap *_dist;
    ///Indicates if \ref _dist is locally allocated (\c true) or not.
    bool local_dist;
    ///Pointer to the map of reached status of the nodes.
    ReachedMap *_reached;
    ///Indicates if \ref _reached is locally allocated (\c true) or not.
    bool local_reached;
    ///Pointer to the map of processed status of the nodes.
    ProcessedMap *_processed;
    ///Indicates if \ref _processed is locally allocated (\c true) or not.
    bool local_processed;

    std::vector<typename Graph::Node> _queue;
    int _queue_head,_queue_tail,_queue_next_dist;
    int _curr_dist;

    ///Creates the maps if necessary.
    
    ///\todo Better memory allocation (instead of new).
    void create_maps() 
    {
      if(!_pred) {
	local_pred = true;
	_pred = Traits::createPredMap(*G);
      }
      if(!_dist) {
	local_dist = true;
	_dist = Traits::createDistMap(*G);
      }
      if(!_reached) {
	local_reached = true;
	_reached = Traits::createReachedMap(*G);
      }
      if(!_processed) {
	local_processed = true;
	_processed = Traits::createProcessedMap(*G);
      }
    }

  protected:
    
    Bfs() {}
    
  public:
 
    typedef Bfs Create;

    ///\name Named template parameters

    ///@{

    template <class T>
    struct DefPredMapTraits : public Traits {
      typedef T PredMap;
      static PredMap *createPredMap(const Graph &) 
      {
	throw UninitializedParameter();
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///PredMap type
    ///
    ///\ref named-templ-param "Named parameter" for setting PredMap type
    ///
    template <class T>
    struct DefPredMap : public Bfs< Graph, DefPredMapTraits<T> > { 
      typedef Bfs< Graph, DefPredMapTraits<T> > Create;
    };
    
    template <class T>
    struct DefDistMapTraits : public Traits {
      typedef T DistMap;
      static DistMap *createDistMap(const Graph &) 
      {
	throw UninitializedParameter();
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///DistMap type
    ///
    ///\ref named-templ-param "Named parameter" for setting DistMap type
    ///
    template <class T>
    struct DefDistMap : public Bfs< Graph, DefDistMapTraits<T> > { 
      typedef Bfs< Graph, DefDistMapTraits<T> > Create;
    };
    
    template <class T>
    struct DefReachedMapTraits : public Traits {
      typedef T ReachedMap;
      static ReachedMap *createReachedMap(const Graph &) 
      {
	throw UninitializedParameter();
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///ReachedMap type
    ///
    ///\ref named-templ-param "Named parameter" for setting ReachedMap type
    ///
    template <class T>
    struct DefReachedMap : public Bfs< Graph, DefReachedMapTraits<T> > { 
      typedef Bfs< Graph, DefReachedMapTraits<T> > Create;
    };
    
    template <class T>
    struct DefProcessedMapTraits : public Traits {
      typedef T ProcessedMap;
      static ProcessedMap *createProcessedMap(const Graph &) 
      {
	throw UninitializedParameter();
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///ProcessedMap type
    ///
    ///\ref named-templ-param "Named parameter" for setting ProcessedMap type
    ///
    template <class T>
    struct DefProcessedMap : public Bfs< Graph, DefProcessedMapTraits<T> > {
      typedef Bfs< Graph, DefProcessedMapTraits<T> > Create;
    };
    
    struct DefGraphProcessedMapTraits : public Traits {
      typedef typename Graph::template NodeMap<bool> ProcessedMap;
      static ProcessedMap *createProcessedMap(const Graph &G) 
      {
	return new ProcessedMap(G);
      }
    };
    ///\brief \ref named-templ-param "Named parameter"
    ///for setting the ProcessedMap type to be Graph::NodeMap<bool>.
    ///
    ///\ref named-templ-param "Named parameter"
    ///for setting the ProcessedMap type to be Graph::NodeMap<bool>.
    ///If you don't set it explicitly, it will be automatically allocated.
    template <class T>
    struct DefProcessedMapToBeDefaultMap :
      public Bfs< Graph, DefGraphProcessedMapTraits> { 
      typedef Bfs< Graph, DefGraphProcessedMapTraits> Create;
    };
    
    ///@}

  public:      
    
    ///Constructor.
    
    ///\param _G the graph the algorithm will run on.
    ///
    Bfs(const Graph& _G) :
      G(&_G),
      _pred(NULL), local_pred(false),
      _dist(NULL), local_dist(false),
      _reached(NULL), local_reached(false),
      _processed(NULL), local_processed(false)
    { }
    
    ///Destructor.
    ~Bfs() 
    {
      if(local_pred) delete _pred;
      if(local_dist) delete _dist;
      if(local_reached) delete _reached;
      if(local_processed) delete _processed;
    }

    ///Sets the map storing the predecessor edges.

    ///Sets the map storing the predecessor edges.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destructor deallocates this
    ///automatically allocated map, of course.
    ///\return <tt> (*this) </tt>
    Bfs &predMap(PredMap &m) 
    {
      if(local_pred) {
	delete _pred;
	local_pred=false;
      }
      _pred = &m;
      return *this;
    }

    ///Sets the map indicating the reached nodes.

    ///Sets the map indicating the reached nodes.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destructor deallocates this
    ///automatically allocated map, of course.
    ///\return <tt> (*this) </tt>
    Bfs &reachedMap(ReachedMap &m) 
    {
      if(local_reached) {
	delete _reached;
	local_reached=false;
      }
      _reached = &m;
      return *this;
    }

    ///Sets the map indicating the processed nodes.

    ///Sets the map indicating the processed nodes.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destructor deallocates this
    ///automatically allocated map, of course.
    ///\return <tt> (*this) </tt>
    Bfs &processedMap(ProcessedMap &m) 
    {
      if(local_processed) {
	delete _processed;
	local_processed=false;
      }
      _processed = &m;
      return *this;
    }

    ///Sets the map storing the distances calculated by the algorithm.

    ///Sets the map storing the distances calculated by the algorithm.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destructor deallocates this
    ///automatically allocated map, of course.
    ///\return <tt> (*this) </tt>
    Bfs &distMap(DistMap &m) 
    {
      if(local_dist) {
	delete _dist;
	local_dist=false;
      }
      _dist = &m;
      return *this;
    }

  public:
    ///\name Execution control
    ///The simplest way to execute the algorithm is to use
    ///one of the member functions called \c run(...).
    ///\n
    ///If you need more control on the execution,
    ///first you must call \ref init(), then you can add several source nodes
    ///with \ref addSource().
    ///Finally \ref start() will perform the actual path
    ///computation.

    ///@{

    ///\brief Initializes the internal data structures.
    ///
    ///Initializes the internal data structures.
    ///
    void init()
    {
      create_maps();
      _queue.resize(countNodes(*G));
      _queue_head=_queue_tail=0;
      _curr_dist=1;
      for ( NodeIt u(*G) ; u!=INVALID ; ++u ) {
	_pred->set(u,INVALID);
	_reached->set(u,false);
	_processed->set(u,false);
      }
    }
    
    ///Adds a new source node.

    ///Adds a new source node to the set of nodes to be processed.
    ///
    void addSource(Node s)
    {
      if(!(*_reached)[s])
	{
	  _reached->set(s,true);
	  _pred->set(s,INVALID);
	  _dist->set(s,0);
	  _queue[_queue_head++]=s;
	  _queue_next_dist=_queue_head;
	}
    }
    
    ///Processes the next node.

    ///Processes the next node.
    ///
    ///\return The processed node.
    ///
    ///\warning The queue must not be empty!
    Node processNextNode()
    {
      if(_queue_tail==_queue_next_dist) {
	_curr_dist++;
	_queue_next_dist=_queue_head;
      }
      Node n=_queue[_queue_tail++];
      _processed->set(n,true);
      Node m;
      for(OutEdgeIt e(*G,n);e!=INVALID;++e)
	if(!(*_reached)[m=G->target(e)]) {
	  _queue[_queue_head++]=m;
	  _reached->set(m,true);
	  _pred->set(m,e);
	  _dist->set(m,_curr_dist);
	}
      return n;
    }

    ///Processes the next node.

    ///Processes the next node. And checks that the given target node
    ///is reached. If the target node is reachable from the processed
    ///node then the reached parameter will be set true. The reached
    ///parameter should be initially false.
    ///
    ///\param target The target node.
    ///\retval reach Indicates that the target node is reached.
    ///\return The processed node.
    ///
    ///\warning The queue must not be empty!
    Node processNextNode(Node target, bool& reach)
    {
      if(_queue_tail==_queue_next_dist) {
	_curr_dist++;
	_queue_next_dist=_queue_head;
      }
      Node n=_queue[_queue_tail++];
      _processed->set(n,true);
      Node m;
      for(OutEdgeIt e(*G,n);e!=INVALID;++e)
	if(!(*_reached)[m=G->target(e)]) {
	  _queue[_queue_head++]=m;
	  _reached->set(m,true);
	  _pred->set(m,e);
	  _dist->set(m,_curr_dist);
          reach = reach || (target == m);
	}
      return n;
    }

    ///Processes the next node.

    ///Processes the next node. And checks that at least one of
    ///reached node has true value in the \c nm node map. If one node
    ///with true value is reachable from the processed node then the
    ///rnode parameter will be set to the first of such nodes.
    ///
    ///\param nm The node map of possible targets.
    ///\retval rnode The reached target node.
    ///\return The processed node.
    ///
    ///\warning The queue must not be empty!
    template<class NM>
    Node processNextNode(const NM& nm, Node& rnode)
    {
      if(_queue_tail==_queue_next_dist) {
	_curr_dist++;
	_queue_next_dist=_queue_head;
      }
      Node n=_queue[_queue_tail++];
      _processed->set(n,true);
      Node m;
      for(OutEdgeIt e(*G,n);e!=INVALID;++e)
	if(!(*_reached)[m=G->target(e)]) {
	  _queue[_queue_head++]=m;
	  _reached->set(m,true);
	  _pred->set(m,e);
	  _dist->set(m,_curr_dist);
	  if (nm[m] && rnode == INVALID) rnode = m;
	}
      return n;
    }
      
    ///Next node to be processed.

    ///Next node to be processed.
    ///
    ///\return The next node to be processed or INVALID if the queue is
    /// empty.
    Node nextNode()
    { 
      return _queue_tail<_queue_head?_queue[_queue_tail]:INVALID;
    }
 
    ///\brief Returns \c false if there are nodes
    ///to be processed in the queue
    ///
    ///Returns \c false if there are nodes
    ///to be processed in the queue
    bool emptyQueue() { return _queue_tail==_queue_head; }
    ///Returns the number of the nodes to be processed.
    
    ///Returns the number of the nodes to be processed in the queue.
    int queueSize() { return _queue_head-_queue_tail; }
    
    ///Executes the algorithm.

    ///Executes the algorithm.
    ///
    ///\pre init() must be called and at least one node should be added
    ///with addSource() before using this function.
    ///
    ///This method runs the %BFS algorithm from the root node(s)
    ///in order to
    ///compute the
    ///shortest path to each node. The algorithm computes
    ///- The shortest path tree.
    ///- The distance of each node from the root(s).
    void start()
    {
      while ( !emptyQueue() ) processNextNode();
    }
    
    ///Executes the algorithm until \c dest is reached.

    ///Executes the algorithm until \c dest is reached.
    ///
    ///\pre init() must be called and at least one node should be added
    ///with addSource() before using this function.
    ///
    ///This method runs the %BFS algorithm from the root node(s)
    ///in order to compute the shortest path to \c dest.
    ///The algorithm computes
    ///- The shortest path to \c  dest.
    ///- The distance of \c dest from the root(s).
    void start(Node dest)
    {
      bool reach = false;
      while ( !emptyQueue() && !reach ) processNextNode(dest, reach);
    }
    
    ///Executes the algorithm until a condition is met.

    ///Executes the algorithm until a condition is met.
    ///
    ///\pre init() must be called and at least one node should be added
    ///with addSource() before using this function.
    ///
    ///\param nm must be a bool (or convertible) node map. The
    ///algorithm will stop when it reaches a node \c v with
    /// <tt>nm[v]</tt> true.
    ///
    ///\return The reached node \c v with <tt>nm[v]</tt> true or
    ///\c INVALID if no such node was found.
    template<class NM>
    Node start(const NM &nm)
    {
      Node rnode = INVALID;
      while ( !emptyQueue() && rnode == INVALID ) {
	processNextNode(nm, rnode);
      }
      return rnode;
    }
    
    ///Runs %BFS algorithm from node \c s.
    
    ///This method runs the %BFS algorithm from a root node \c s
    ///in order to
    ///compute the
    ///shortest path to each node. The algorithm computes
    ///- The shortest path tree.
    ///- The distance of each node from the root.
    ///
    ///\note b.run(s) is just a shortcut of the following code.
    ///\code
    ///  b.init();
    ///  b.addSource(s);
    ///  b.start();
    ///\endcode
    void run(Node s) {
      init();
      addSource(s);
      start();
    }
    
    ///Finds the shortest path between \c s and \c t.
    
    ///Finds the shortest path between \c s and \c t.
    ///
    ///\return The length of the shortest s---t path if there exists one,
    ///0 otherwise.
    ///\note Apart from the return value, b.run(s) is
    ///just a shortcut of the following code.
    ///\code
    ///  b.init();
    ///  b.addSource(s);
    ///  b.start(t);
    ///\endcode
    int run(Node s,Node t) {
      init();
      addSource(s);
      start(t);
      return reached(t) ? _curr_dist : 0;
    }
    
    ///@}

    ///\name Query Functions
    ///The result of the %BFS algorithm can be obtained using these
    ///functions.\n
    ///Before the use of these functions,
    ///either run() or start() must be calleb.
    
    ///@{

    typedef PredMapPath<Graph, PredMap> Path;

    ///Gives back the shortest path.
    
    ///Gives back the shortest path.
    ///\pre The \c t should be reachable from the source.
    Path path(Node t) 
    {
      return Path(*G, *_pred, t);
    }

    ///The distance of a node from the root(s).

    ///Returns the distance of a node from the root(s).
    ///\pre \ref run() must be called before using this function.
    ///\warning If node \c v in unreachable from the root(s) the return value
    ///of this function is undefined.
    int dist(Node v) const { return (*_dist)[v]; }

    ///Returns the 'previous edge' of the shortest path tree.

    ///For a node \c v it returns the 'previous edge'
    ///of the shortest path tree,
    ///i.e. it returns the last edge of a shortest path from the root(s) to \c
    ///v. It is \ref INVALID
    ///if \c v is unreachable from the root(s) or \c v is a root. The
    ///shortest path tree used here is equal to the shortest path tree used in
    ///\ref predNode().
    ///\pre Either \ref run() or \ref start() must be called before using
    ///this function.
    Edge predEdge(Node v) const { return (*_pred)[v];}

    ///Returns the 'previous node' of the shortest path tree.

    ///For a node \c v it returns the 'previous node'
    ///of the shortest path tree,
    ///i.e. it returns the last but one node from a shortest path from the
    ///root(a) to \c /v.
    ///It is INVALID if \c v is unreachable from the root(s) or
    ///if \c v itself a root.
    ///The shortest path tree used here is equal to the shortest path
    ///tree used in \ref predEdge().
    ///\pre Either \ref run() or \ref start() must be called before
    ///using this function.
    Node predNode(Node v) const { return (*_pred)[v]==INVALID ? INVALID:
				  G->source((*_pred)[v]); }
    
    ///Returns a reference to the NodeMap of distances.

    ///Returns a reference to the NodeMap of distances.
    ///\pre Either \ref run() or \ref init() must
    ///be called before using this function.
    const DistMap &distMap() const { return *_dist;}
 
    ///Returns a reference to the shortest path tree map.

    ///Returns a reference to the NodeMap of the edges of the
    ///shortest path tree.
    ///\pre Either \ref run() or \ref init()
    ///must be called before using this function.
    const PredMap &predMap() const { return *_pred;}
 
    ///Checks if a node is reachable from the root.

    ///Returns \c true if \c v is reachable from the root.
    ///\warning The source nodes are indicated as unreached.
    ///\pre Either \ref run() or \ref start()
    ///must be called before using this function.
    ///
    bool reached(Node v) { return (*_reached)[v]; }
    
    ///@}
  };

  ///Default traits class of Bfs function.

  ///Default traits class of Bfs function.
  ///\param GR Graph type.
  template<class GR>
  struct BfsWizardDefaultTraits
  {
    ///The graph type the algorithm runs on. 
    typedef GR Graph;
    ///\brief The type of the map that stores the last
    ///edges of the shortest paths.
    /// 
    ///The type of the map that stores the last
    ///edges of the shortest paths.
    ///It must meet the \ref concepts::WriteMap "WriteMap" concept.
    ///
    typedef NullMap<typename Graph::Node,typename GR::Edge> PredMap;
    ///Instantiates a PredMap.
 
    ///This function instantiates a \ref PredMap. 
    ///\param g is the graph, to which we would like to define the PredMap.
    ///\todo The graph alone may be insufficient to initialize
#ifdef DOXYGEN
    static PredMap *createPredMap(const GR &g) 
#else
    static PredMap *createPredMap(const GR &) 
#endif
    {
      return new PredMap();
    }

    ///The type of the map that indicates which nodes are processed.
 
    ///The type of the map that indicates which nodes are processed.
    ///It must meet the \ref concepts::WriteMap "WriteMap" concept.
    ///\todo named parameter to set this type, function to read and write.
    typedef NullMap<typename Graph::Node,bool> ProcessedMap;
    ///Instantiates a ProcessedMap.
 
    ///This function instantiates a \ref ProcessedMap. 
    ///\param g is the graph, to which
    ///we would like to define the \ref ProcessedMap
#ifdef DOXYGEN
    static ProcessedMap *createProcessedMap(const GR &g)
#else
    static ProcessedMap *createProcessedMap(const GR &)
#endif
    {
      return new ProcessedMap();
    }
    ///The type of the map that indicates which nodes are reached.
 
    ///The type of the map that indicates which nodes are reached.
    ///It must meet the \ref concepts::WriteMap "WriteMap" concept.
    ///\todo named parameter to set this type, function to read and write.
    typedef typename Graph::template NodeMap<bool> ReachedMap;
    ///Instantiates a ReachedMap.
 
    ///This function instantiates a \ref ReachedMap. 
    ///\param G is the graph, to which
    ///we would like to define the \ref ReachedMap.
    static ReachedMap *createReachedMap(const GR &G)
    {
      return new ReachedMap(G);
    }
    ///The type of the map that stores the dists of the nodes.
 
    ///The type of the map that stores the dists of the nodes.
    ///It must meet the \ref concepts::WriteMap "WriteMap" concept.
    ///
    typedef NullMap<typename Graph::Node,int> DistMap;
    ///Instantiates a DistMap.
 
    ///This function instantiates a \ref DistMap. 
    ///\param g is the graph, to which we would like to define the \ref DistMap
#ifdef DOXYGEN
    static DistMap *createDistMap(const GR &g)
#else
    static DistMap *createDistMap(const GR &)
#endif
    {
      return new DistMap();
    }
  };
  
  /// Default traits used by \ref BfsWizard

  /// To make it easier to use Bfs algorithm
  ///we have created a wizard class.
  /// This \ref BfsWizard class needs default traits,
  ///as well as the \ref Bfs class.
  /// The \ref BfsWizardBase is a class to be the default traits of the
  /// \ref BfsWizard class.
  template<class GR>
  class BfsWizardBase : public BfsWizardDefaultTraits<GR>
  {

    typedef BfsWizardDefaultTraits<GR> Base;
  protected:
    /// Type of the nodes in the graph.
    typedef typename Base::Graph::Node Node;

    /// Pointer to the underlying graph.
    void *_g;
    ///Pointer to the map of reached nodes.
    void *_reached;
    ///Pointer to the map of processed nodes.
    void *_processed;
    ///Pointer to the map of predecessors edges.
    void *_pred;
    ///Pointer to the map of distances.
    void *_dist;
    ///Pointer to the source node.
    Node _source;
    
    public:
    /// Constructor.
    
    /// This constructor does not require parameters, therefore it initiates
    /// all of the attributes to default values (0, INVALID).
    BfsWizardBase() : _g(0), _reached(0), _processed(0), _pred(0),
			   _dist(0), _source(INVALID) {}

    /// Constructor.
    
    /// This constructor requires some parameters,
    /// listed in the parameters list.
    /// Others are initiated to 0.
    /// \param g is the initial value of  \ref _g
    /// \param s is the initial value of  \ref _source
    BfsWizardBase(const GR &g, Node s=INVALID) :
      _g(reinterpret_cast<void*>(const_cast<GR*>(&g))), 
      _reached(0), _processed(0), _pred(0), _dist(0), _source(s) {}

  };
  
  /// A class to make the usage of Bfs algorithm easier

  /// This class is created to make it easier to use Bfs algorithm.
  /// It uses the functions and features of the plain \ref Bfs,
  /// but it is much simpler to use it.
  ///
  /// Simplicity means that the way to change the types defined
  /// in the traits class is based on functions that returns the new class
  /// and not on templatable built-in classes.
  /// When using the plain \ref Bfs
  /// the new class with the modified type comes from
  /// the original class by using the ::
  /// operator. In the case of \ref BfsWizard only
  /// a function have to be called and it will
  /// return the needed class.
  ///
  /// It does not have own \ref run method. When its \ref run method is called
  /// it initiates a plain \ref Bfs class, and calls the \ref Bfs::run
  /// method of it.
  template<class TR>
  class BfsWizard : public TR
  {
    typedef TR Base;

    ///The type of the underlying graph.
    typedef typename TR::Graph Graph;
    //\e
    typedef typename Graph::Node Node;
    //\e
    typedef typename Graph::NodeIt NodeIt;
    //\e
    typedef typename Graph::Edge Edge;
    //\e
    typedef typename Graph::OutEdgeIt OutEdgeIt;
    
    ///\brief The type of the map that stores
    ///the reached nodes
    typedef typename TR::ReachedMap ReachedMap;
    ///\brief The type of the map that stores
    ///the processed nodes
    typedef typename TR::ProcessedMap ProcessedMap;
    ///\brief The type of the map that stores the last
    ///edges of the shortest paths.
    typedef typename TR::PredMap PredMap;
    ///The type of the map that stores the dists of the nodes.
    typedef typename TR::DistMap DistMap;

  public:
    /// Constructor.
    BfsWizard() : TR() {}

    /// Constructor that requires parameters.

    /// Constructor that requires parameters.
    /// These parameters will be the default values for the traits class.
    BfsWizard(const Graph &g, Node s=INVALID) :
      TR(g,s) {}

    ///Copy constructor
    BfsWizard(const TR &b) : TR(b) {}

    ~BfsWizard() {}

    ///Runs Bfs algorithm from a given node.
    
    ///Runs Bfs algorithm from a given node.
    ///The node can be given by the \ref source function.
    void run()
    {
      if(Base::_source==INVALID) throw UninitializedParameter();
      Bfs<Graph,TR> alg(*reinterpret_cast<const Graph*>(Base::_g));
      if(Base::_reached)
	alg.reachedMap(*reinterpret_cast<ReachedMap*>(Base::_reached));
      if(Base::_processed) 
        alg.processedMap(*reinterpret_cast<ProcessedMap*>(Base::_processed));
      if(Base::_pred) 
        alg.predMap(*reinterpret_cast<PredMap*>(Base::_pred));
      if(Base::_dist) 
        alg.distMap(*reinterpret_cast<DistMap*>(Base::_dist));
      alg.run(Base::_source);
    }

    ///Runs Bfs algorithm from the given node.

    ///Runs Bfs algorithm from the given node.
    ///\param s is the given source.
    void run(Node s)
    {
      Base::_source=s;
      run();
    }

    template<class T>
    struct DefPredMapBase : public Base {
      typedef T PredMap;
      static PredMap *createPredMap(const Graph &) { return 0; };
      DefPredMapBase(const TR &b) : TR(b) {}
    };
    
    ///\brief \ref named-templ-param "Named parameter"
    ///function for setting PredMap
    ///
    /// \ref named-templ-param "Named parameter"
    ///function for setting PredMap
    ///
    template<class T>
    BfsWizard<DefPredMapBase<T> > predMap(const T &t) 
    {
      Base::_pred=reinterpret_cast<void*>(const_cast<T*>(&t));
      return BfsWizard<DefPredMapBase<T> >(*this);
    }
    
 
    template<class T>
    struct DefReachedMapBase : public Base {
      typedef T ReachedMap;
      static ReachedMap *createReachedMap(const Graph &) { return 0; };
      DefReachedMapBase(const TR &b) : TR(b) {}
    };
    
    ///\brief \ref named-templ-param "Named parameter"
    ///function for setting ReachedMap
    ///
    /// \ref named-templ-param "Named parameter"
    ///function for setting ReachedMap
    ///
    template<class T>
    BfsWizard<DefReachedMapBase<T> > reachedMap(const T &t) 
    {
      Base::_pred=reinterpret_cast<void*>(const_cast<T*>(&t));
      return BfsWizard<DefReachedMapBase<T> >(*this);
    }
    

    template<class T>
    struct DefProcessedMapBase : public Base {
      typedef T ProcessedMap;
      static ProcessedMap *createProcessedMap(const Graph &) { return 0; };
      DefProcessedMapBase(const TR &b) : TR(b) {}
    };
    
    ///\brief \ref named-templ-param "Named parameter"
    ///function for setting ProcessedMap
    ///
    /// \ref named-templ-param "Named parameter"
    ///function for setting ProcessedMap
    ///
    template<class T>
    BfsWizard<DefProcessedMapBase<T> > processedMap(const T &t) 
    {
      Base::_pred=reinterpret_cast<void*>(const_cast<T*>(&t));
      return BfsWizard<DefProcessedMapBase<T> >(*this);
    }
    
   
    template<class T>
    struct DefDistMapBase : public Base {
      typedef T DistMap;
      static DistMap *createDistMap(const Graph &) { return 0; };
      DefDistMapBase(const TR &b) : TR(b) {}
    };
    
    ///\brief \ref named-templ-param "Named parameter"
    ///function for setting DistMap type
    ///
    /// \ref named-templ-param "Named parameter"
    ///function for setting DistMap type
    ///
    template<class T>
    BfsWizard<DefDistMapBase<T> > distMap(const T &t) 
    {
      Base::_dist=reinterpret_cast<void*>(const_cast<T*>(&t));
      return BfsWizard<DefDistMapBase<T> >(*this);
    }
    
    /// Sets the source node, from which the Bfs algorithm runs.

    /// Sets the source node, from which the Bfs algorithm runs.
    /// \param s is the source node.
    BfsWizard<TR> &source(Node s) 
    {
      Base::_source=s;
      return *this;
    }
    
  };
  
  ///Function type interface for Bfs algorithm.

  /// \ingroup search
  ///Function type interface for Bfs algorithm.
  ///
  ///This function also has several
  ///\ref named-templ-func-param "named parameters",
  ///they are declared as the members of class \ref BfsWizard.
  ///The following
  ///example shows how to use these parameters.
  ///\code
  ///  bfs(g,source).predMap(preds).run();
  ///\endcode
  ///\warning Don't forget to put the \ref BfsWizard::run() "run()"
  ///to the end of the parameter list.
  ///\sa BfsWizard
  ///\sa Bfs
  template<class GR>
  BfsWizard<BfsWizardBase<GR> >
  bfs(const GR &g,typename GR::Node s=INVALID)
  {
    return BfsWizard<BfsWizardBase<GR> >(g,s);
  }

#ifdef DOXYGEN
  /// \brief Visitor class for bfs.
  ///  
  /// This class defines the interface of the BfsVisit events, and
  /// it could be the base of a real Visitor class.
  template <typename _Graph>
  struct BfsVisitor {
    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::Node Node;
    /// \brief Called when the edge reach a node.
    /// 
    /// It is called when the bfs find an edge which target is not
    /// reached yet.
    void discover(const Edge& edge) {}
    /// \brief Called when the node reached first time.
    /// 
    /// It is Called when the node reached first time.
    void reach(const Node& node) {}
    /// \brief Called when the edge examined but target of the edge 
    /// already discovered.
    /// 
    /// It called when the edge examined but the target of the edge 
    /// already discovered.
    void examine(const Edge& edge) {}
    /// \brief Called for the source node of the bfs.
    /// 
    /// It is called for the source node of the bfs.
    void start(const Node& node) {}
    /// \brief Called when the node processed.
    /// 
    /// It is Called when the node processed.
    void process(const Node& node) {}
  };
#else
  template <typename _Graph>
  struct BfsVisitor {
    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::Node Node;
    void discover(const Edge&) {}
    void reach(const Node&) {}
    void examine(const Edge&) {}
    void start(const Node&) {}
    void process(const Node&) {}

    template <typename _Visitor>
    struct Constraints {
      void constraints() {
	Edge edge;
	Node node;
	visitor.discover(edge);
	visitor.reach(node);
	visitor.examine(edge);
	visitor.start(node);
        visitor.process(node);
      }
      _Visitor& visitor;
    };
  };
#endif

  /// \brief Default traits class of BfsVisit class.
  ///
  /// Default traits class of BfsVisit class.
  /// \param _Graph Graph type.
  template<class _Graph>
  struct BfsVisitDefaultTraits {

    /// \brief The graph type the algorithm runs on. 
    typedef _Graph Graph;

    /// \brief The type of the map that indicates which nodes are reached.
    /// 
    /// The type of the map that indicates which nodes are reached.
    /// It must meet the \ref concepts::WriteMap "WriteMap" concept.
    /// \todo named parameter to set this type, function to read and write.
    typedef typename Graph::template NodeMap<bool> ReachedMap;

    /// \brief Instantiates a ReachedMap.
    ///
    /// This function instantiates a \ref ReachedMap. 
    /// \param graph is the graph, to which
    /// we would like to define the \ref ReachedMap.
    static ReachedMap *createReachedMap(const Graph &graph) {
      return new ReachedMap(graph);
    }

  };

  /// \ingroup search
  ///  
  /// \brief %BFS Visit algorithm class.
  ///  
  /// This class provides an efficient implementation of the %BFS algorithm
  /// with visitor interface.
  ///
  /// The %BfsVisit class provides an alternative interface to the Bfs
  /// class. It works with callback mechanism, the BfsVisit object calls
  /// on every bfs event the \c Visitor class member functions. 
  ///
  /// \param _Graph The graph type the algorithm runs on. The default value is
  /// \ref ListGraph. The value of _Graph is not used directly by Bfs, it
  /// is only passed to \ref BfsDefaultTraits.
  /// \param _Visitor The Visitor object for the algorithm. The 
  /// \ref BfsVisitor "BfsVisitor<_Graph>" is an empty Visitor which
  /// does not observe the Bfs events. If you want to observe the bfs
  /// events you should implement your own Visitor class.
  /// \param _Traits Traits class to set various data types used by the 
  /// algorithm. The default traits class is
  /// \ref BfsVisitDefaultTraits "BfsVisitDefaultTraits<_Graph>".
  /// See \ref BfsVisitDefaultTraits for the documentation of
  /// a Bfs visit traits class.
  ///
  /// \author Jacint Szabo, Alpar Juttner and Balazs Dezso
#ifdef DOXYGEN
  template <typename _Graph, typename _Visitor, typename _Traits>
#else
  template <typename _Graph = ListGraph,
	    typename _Visitor = BfsVisitor<_Graph>,
	    typename _Traits = BfsDefaultTraits<_Graph> >
#endif
  class BfsVisit {
  public:
    
    /// \brief \ref Exception for uninitialized parameters.
    ///
    /// This error represents problems in the initialization
    /// of the parameters of the algorithms.
    class UninitializedParameter : public lemon::UninitializedParameter {
    public:
      virtual const char* what() const throw() 
      {
	return "lemon::BfsVisit::UninitializedParameter";
      }
    };

    typedef _Traits Traits;

    typedef typename Traits::Graph Graph;

    typedef _Visitor Visitor;

    ///The type of the map indicating which nodes are reached.
    typedef typename Traits::ReachedMap ReachedMap;

  private:

    typedef typename Graph::Node Node;
    typedef typename Graph::NodeIt NodeIt;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::OutEdgeIt OutEdgeIt;

    /// Pointer to the underlying graph.
    const Graph *_graph;
    /// Pointer to the visitor object.
    Visitor *_visitor;
    ///Pointer to the map of reached status of the nodes.
    ReachedMap *_reached;
    ///Indicates if \ref _reached is locally allocated (\c true) or not.
    bool local_reached;

    std::vector<typename Graph::Node> _list;
    int _list_front, _list_back;

    /// \brief Creates the maps if necessary.
    ///
    /// Creates the maps if necessary.
    void create_maps() {
      if(!_reached) {
	local_reached = true;
	_reached = Traits::createReachedMap(*_graph);
      }
    }

  protected:

    BfsVisit() {}
    
  public:

    typedef BfsVisit Create;

    /// \name Named template parameters

    ///@{
    template <class T>
    struct DefReachedMapTraits : public Traits {
      typedef T ReachedMap;
      static ReachedMap *createReachedMap(const Graph &graph) {
	throw UninitializedParameter();
      }
    };
    /// \brief \ref named-templ-param "Named parameter" for setting 
    /// ReachedMap type
    ///
    /// \ref named-templ-param "Named parameter" for setting ReachedMap type
    template <class T>
    struct DefReachedMap : public BfsVisit< Graph, Visitor,
					    DefReachedMapTraits<T> > {
      typedef BfsVisit< Graph, Visitor, DefReachedMapTraits<T> > Create;
    };
    ///@}

  public:      
    
    /// \brief Constructor.
    ///
    /// Constructor.
    ///
    /// \param graph the graph the algorithm will run on.
    /// \param visitor The visitor of the algorithm.
    ///
    BfsVisit(const Graph& graph, Visitor& visitor) 
      : _graph(&graph), _visitor(&visitor),
	_reached(0), local_reached(false) {}
    
    /// \brief Destructor.
    ///
    /// Destructor.
    ~BfsVisit() {
      if(local_reached) delete _reached;
    }

    /// \brief Sets the map indicating if a node is reached.
    ///
    /// Sets the map indicating if a node is reached.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated map, of course.
    /// \return <tt> (*this) </tt>
    BfsVisit &reachedMap(ReachedMap &m) {
      if(local_reached) {
	delete _reached;
	local_reached = false;
      }
      _reached = &m;
      return *this;
    }

  public:
    /// \name Execution control
    /// The simplest way to execute the algorithm is to use
    /// one of the member functions called \c run(...).
    /// \n
    /// If you need more control on the execution,
    /// first you must call \ref init(), then you can adda source node
    /// with \ref addSource().
    /// Finally \ref start() will perform the actual path
    /// computation.

    /// @{
    /// \brief Initializes the internal data structures.
    ///
    /// Initializes the internal data structures.
    ///
    void init() {
      create_maps();
      _list.resize(countNodes(*_graph));
      _list_front = _list_back = -1;
      for (NodeIt u(*_graph) ; u != INVALID ; ++u) {
	_reached->set(u, false);
      }
    }
    
    /// \brief Adds a new source node.
    ///
    /// Adds a new source node to the set of nodes to be processed.
    void addSource(Node s) {
      if(!(*_reached)[s]) {
	  _reached->set(s,true);
	  _visitor->start(s);
	  _visitor->reach(s);
          _list[++_list_back] = s;
	}
    }
    
    /// \brief Processes the next node.
    ///
    /// Processes the next node.
    ///
    /// \return The processed node.
    ///
    /// \pre The queue must not be empty!
    Node processNextNode() { 
      Node n = _list[++_list_front];
      _visitor->process(n);
      Edge e;
      for (_graph->firstOut(e, n); e != INVALID; _graph->nextOut(e)) {
        Node m = _graph->target(e);
        if (!(*_reached)[m]) {
          _visitor->discover(e);
          _visitor->reach(m);
          _reached->set(m, true);
          _list[++_list_back] = m;
        } else {
          _visitor->examine(e);
        }
      }
      return n;
    }

    /// \brief Processes the next node.
    ///
    /// Processes the next node. And checks that the given target node
    /// is reached. If the target node is reachable from the processed
    /// node then the reached parameter will be set true. The reached
    /// parameter should be initially false.
    ///
    /// \param target The target node.
    /// \retval reach Indicates that the target node is reached.
    /// \return The processed node.
    ///
    /// \warning The queue must not be empty!
    Node processNextNode(Node target, bool& reach) {
      Node n = _list[++_list_front];
      _visitor->process(n);
      Edge e;
      for (_graph->firstOut(e, n); e != INVALID; _graph->nextOut(e)) {
        Node m = _graph->target(e);
        if (!(*_reached)[m]) {
          _visitor->discover(e);
          _visitor->reach(m);
          _reached->set(m, true);
          _list[++_list_back] = m;
          reach = reach || (target == m);
        } else {
          _visitor->examine(e);
        }
      }
      return n;
    }

    /// \brief Processes the next node.
    ///
    /// Processes the next node. And checks that at least one of
    /// reached node has true value in the \c nm node map. If one node
    /// with true value is reachable from the processed node then the
    /// rnode parameter will be set to the first of such nodes.
    ///
    /// \param nm The node map of possible targets.
    /// \retval rnode The reached target node.
    /// \return The processed node.
    ///
    /// \warning The queue must not be empty!
    template <typename NM>
    Node processNextNode(const NM& nm, Node& rnode) {
      Node n = _list[++_list_front];
      _visitor->process(n);
      Edge e;
      for (_graph->firstOut(e, n); e != INVALID; _graph->nextOut(e)) {
        Node m = _graph->target(e);
        if (!(*_reached)[m]) {
          _visitor->discover(e);
          _visitor->reach(m);
          _reached->set(m, true);
          _list[++_list_back] = m;
          if (nm[m] && rnode == INVALID) rnode = m;
        } else {
          _visitor->examine(e);
        }
      }
      return n;
    }

    /// \brief Next node to be processed.
    ///
    /// Next node to be processed.
    ///
    /// \return The next node to be processed or INVALID if the stack is
    /// empty.
    Node nextNode() { 
      return _list_front != _list_back ? _list[_list_front + 1] : INVALID;
    }

    /// \brief Returns \c false if there are nodes
    /// to be processed in the queue
    ///
    /// Returns \c false if there are nodes
    /// to be processed in the queue
    bool emptyQueue() { return _list_front == _list_back; }

    /// \brief Returns the number of the nodes to be processed.
    ///
    /// Returns the number of the nodes to be processed in the queue.
    int queueSize() { return _list_back - _list_front; }
    
    /// \brief Executes the algorithm.
    ///
    /// Executes the algorithm.
    ///
    /// \pre init() must be called and at least one node should be added
    /// with addSource() before using this function.
    void start() {
      while ( !emptyQueue() ) processNextNode();
    }
    
    /// \brief Executes the algorithm until \c dest is reached.
    ///
    /// Executes the algorithm until \c dest is reached.
    ///
    /// \pre init() must be called and at least one node should be added
    /// with addSource() before using this function.
    void start(Node dest) {
      bool reach = false;
      while ( !emptyQueue() && !reach ) processNextNode(dest, reach);
    }
    
    /// \brief Executes the algorithm until a condition is met.
    ///
    /// Executes the algorithm until a condition is met.
    ///
    /// \pre init() must be called and at least one node should be added
    /// with addSource() before using this function.
    ///
    ///\param nm must be a bool (or convertible) node map. The
    ///algorithm will stop when it reaches a node \c v with
    /// <tt>nm[v]</tt> true.
    ///
    ///\return The reached node \c v with <tt>nm[v]</tt> true or
    ///\c INVALID if no such node was found.
    template <typename NM>
    Node start(const NM &nm) {
      Node rnode = INVALID;
      while ( !emptyQueue() && rnode == INVALID ) {
	processNextNode(nm, rnode);
      }
      return rnode;
    }

    /// \brief Runs %BFSVisit algorithm from node \c s.
    ///
    /// This method runs the %BFS algorithm from a root node \c s.
    /// \note b.run(s) is just a shortcut of the following code.
    ///\code
    ///   b.init();
    ///   b.addSource(s);
    ///   b.start();
    ///\endcode
    void run(Node s) {
      init();
      addSource(s);
      start();
    }

    /// \brief Runs %BFSVisit algorithm to visit all nodes in the graph.
    ///    
    /// This method runs the %BFS algorithm in order to
    /// compute the %BFS path to each node. The algorithm computes
    /// - The %BFS tree.
    /// - The distance of each node from the root in the %BFS tree.
    ///
    ///\note b.run() is just a shortcut of the following code.
    ///\code
    ///  b.init();
    ///  for (NodeIt it(graph); it != INVALID; ++it) {
    ///    if (!b.reached(it)) {
    ///      b.addSource(it);
    ///      b.start();
    ///    }
    ///  }
    ///\endcode
    void run() {
      init();
      for (NodeIt it(*_graph); it != INVALID; ++it) {
        if (!reached(it)) {
          addSource(it);
          start();
        }
      }
    }
    ///@}

    /// \name Query Functions
    /// The result of the %BFS algorithm can be obtained using these
    /// functions.\n
    /// Before the use of these functions,
    /// either run() or start() must be called.
    ///@{

    /// \brief Checks if a node is reachable from the root.
    ///
    /// Returns \c true if \c v is reachable from the root(s).
    /// \warning The source nodes are inditated as unreachable.
    /// \pre Either \ref run() or \ref start()
    /// must be called before using this function.
    ///
    bool reached(Node v) { return (*_reached)[v]; }
    ///@}
  };

} //END OF NAMESPACE LEMON

#endif


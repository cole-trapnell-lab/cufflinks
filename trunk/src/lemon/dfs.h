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

#ifndef LEMON_DFS_H
#define LEMON_DFS_H

///\ingroup search
///\file
///\brief Dfs algorithm.

#include <lemon/list_graph.h>
#include <lemon/graph_utils.h>
#include <lemon/bits/path_dump.h>
#include <lemon/bits/invalid.h>
#include <lemon/error.h>
#include <lemon/maps.h>

#include <lemon/concept_check.h>

namespace lemon {

  
  ///Default traits class of Dfs class.

  ///Default traits class of Dfs class.
  ///\param GR Graph type.
  template<class GR>
  struct DfsDefaultTraits
  {
    ///The graph type the algorithm runs on. 
    typedef GR Graph;
    ///\brief The type of the map that stores the last
    ///edges of the %DFS paths.
    /// 
    ///The type of the map that stores the last
    ///edges of the %DFS paths.
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
  
  ///%DFS algorithm class.
  
  ///\ingroup search
  ///This class provides an efficient implementation of the %DFS algorithm.
  ///
  ///\param GR The graph type the algorithm runs on. The default value is
  ///\ref ListGraph. The value of GR is not used directly by Dfs, it
  ///is only passed to \ref DfsDefaultTraits.
  ///\param TR Traits class to set various data types used by the algorithm.
  ///The default traits class is
  ///\ref DfsDefaultTraits "DfsDefaultTraits<GR>".
  ///See \ref DfsDefaultTraits for the documentation of
  ///a Dfs traits class.
  ///
  ///\author Jacint Szabo and Alpar Juttner
#ifdef DOXYGEN
  template <typename GR,
	    typename TR>
#else
  template <typename GR=ListGraph,
	    typename TR=DfsDefaultTraits<GR> >
#endif
  class Dfs {
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
	return "lemon::Dfs::UninitializedParameter";
      }
    };

    typedef TR Traits;
    ///The type of the underlying graph.
    typedef typename TR::Graph Graph;
    ///\e
    typedef typename Graph::Node Node;
    ///\e
    typedef typename Graph::NodeIt NodeIt;
    ///\e
    typedef typename Graph::Edge Edge;
    ///\e
    typedef typename Graph::OutEdgeIt OutEdgeIt;
    
    ///\brief The type of the map that stores the last
    ///edges of the %DFS paths.
    typedef typename TR::PredMap PredMap;
    ///The type of the map indicating which nodes are reached.
    typedef typename TR::ReachedMap ReachedMap;
    ///The type of the map indicating which nodes are processed.
    typedef typename TR::ProcessedMap ProcessedMap;
    ///The type of the map that stores the dists of the nodes.
    typedef typename TR::DistMap DistMap;
  private:
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

    std::vector<typename Graph::OutEdgeIt> _stack;
    int _stack_head;

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

    Dfs() {}
    
  public:

    typedef Dfs Create;

    ///\name Named template parameters

    ///@{

    template <class T>
    struct DefPredMapTraits : public Traits {
      typedef T PredMap;
      static PredMap *createPredMap(const Graph &G) 
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
    struct DefPredMap : public Dfs<Graph, DefPredMapTraits<T> > {
      typedef Dfs<Graph, DefPredMapTraits<T> > Create;
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
    ///\ref named-templ-param "Named parameter" for setting DistMap
    ///type
    template <class T>
    struct DefDistMap {
      typedef Dfs<Graph, DefDistMapTraits<T> > Create;
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
    struct DefReachedMap : public Dfs< Graph, DefReachedMapTraits<T> > {
      typedef Dfs< Graph, DefReachedMapTraits<T> > Create;
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
    struct DefProcessedMap : public Dfs< Graph, DefProcessedMapTraits<T> > { 
      typedef Dfs< Graph, DefProcessedMapTraits<T> > Create;
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
    ///If you don't set it explicitely, it will be automatically allocated.
    template <class T>
    class DefProcessedMapToBeDefaultMap :
      public Dfs< Graph, DefGraphProcessedMapTraits> { 
      typedef Dfs< Graph, DefGraphProcessedMapTraits> Create;
    };
    
    ///@}

  public:      
    
    ///Constructor.
    
    ///\param _G the graph the algorithm will run on.
    ///
    Dfs(const Graph& _G) :
      G(&_G),
      _pred(NULL), local_pred(false),
      _dist(NULL), local_dist(false),
      _reached(NULL), local_reached(false),
      _processed(NULL), local_processed(false)
    { }
    
    ///Destructor.
    ~Dfs() 
    {
      if(local_pred) delete _pred;
      if(local_dist) delete _dist;
      if(local_reached) delete _reached;
      if(local_processed) delete _processed;
    }

    ///Sets the map storing the predecessor edges.

    ///Sets the map storing the predecessor edges.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destuctor deallocates this
    ///automatically allocated map, of course.
    ///\return <tt> (*this) </tt>
    Dfs &predMap(PredMap &m) 
    {
      if(local_pred) {
	delete _pred;
	local_pred=false;
      }
      _pred = &m;
      return *this;
    }

    ///Sets the map storing the distances calculated by the algorithm.

    ///Sets the map storing the distances calculated by the algorithm.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destuctor deallocates this
    ///automatically allocated map, of course.
    ///\return <tt> (*this) </tt>
    Dfs &distMap(DistMap &m) 
    {
      if(local_dist) {
	delete _dist;
	local_dist=false;
      }
      _dist = &m;
      return *this;
    }

    ///Sets the map indicating if a node is reached.

    ///Sets the map indicating if a node is reached.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destuctor deallocates this
    ///automatically allocated map, of course.
    ///\return <tt> (*this) </tt>
    Dfs &reachedMap(ReachedMap &m) 
    {
      if(local_reached) {
	delete _reached;
	local_reached=false;
      }
      _reached = &m;
      return *this;
    }

    ///Sets the map indicating if a node is processed.

    ///Sets the map indicating if a node is processed.
    ///If you don't use this function before calling \ref run(),
    ///it will allocate one. The destuctor deallocates this
    ///automatically allocated map, of course.
    ///\return <tt> (*this) </tt>
    Dfs &processedMap(ProcessedMap &m) 
    {
      if(local_processed) {
	delete _processed;
	local_processed=false;
      }
      _processed = &m;
      return *this;
    }

  public:
    ///\name Execution control
    ///The simplest way to execute the algorithm is to use
    ///one of the member functions called \c run(...).
    ///\n
    ///If you need more control on the execution,
    ///first you must call \ref init(), then you can add a source node
    ///with \ref addSource().
    ///Finally \ref start() will perform the actual path
    ///computation.

    ///@{

    ///Initializes the internal data structures.

    ///Initializes the internal data structures.
    ///
    void init()
    {
      create_maps();
      _stack.resize(countNodes(*G));
      _stack_head=-1;
      for ( NodeIt u(*G) ; u!=INVALID ; ++u ) {
	_pred->set(u,INVALID);
	// _predNode->set(u,INVALID);
	_reached->set(u,false);
	_processed->set(u,false);
      }
    }
    
    ///Adds a new source node.

    ///Adds a new source node to the set of nodes to be processed.
    ///
    ///\warning dists are wrong (or at least strange)
    ///in case of multiple sources.
    void addSource(Node s)
    {
      if(!(*_reached)[s])
	{
	  _reached->set(s,true);
	  _pred->set(s,INVALID);
	  OutEdgeIt e(*G,s);
	  if(e!=INVALID) {
	    _stack[++_stack_head]=e;
	    _dist->set(s,_stack_head);
	  }
	  else {
	    _processed->set(s,true);
	    _dist->set(s,0);
	  }
	}
    }
    
    ///Processes the next edge.

    ///Processes the next edge.
    ///
    ///\return The processed edge.
    ///
    ///\pre The stack must not be empty!
    Edge processNextEdge()
    { 
      Node m;
      Edge e=_stack[_stack_head];
      if(!(*_reached)[m=G->target(e)]) {
	_pred->set(m,e);
	_reached->set(m,true);
	++_stack_head;
	_stack[_stack_head] = OutEdgeIt(*G, m);
	_dist->set(m,_stack_head);
      }
      else {
	m=G->source(e);
	++_stack[_stack_head];
      }
      while(_stack_head>=0 && _stack[_stack_head]==INVALID) {
	_processed->set(m,true);
	--_stack_head;
	if(_stack_head>=0) {
	  m=G->source(_stack[_stack_head]);
	  ++_stack[_stack_head];
	}
      }
      return e;
    }
    ///Next edge to be processed.

    ///Next edge to be processed.
    ///
    ///\return The next edge to be processed or INVALID if the stack is
    /// empty.
    OutEdgeIt nextEdge()
    { 
      return _stack_head>=0?_stack[_stack_head]:INVALID;
    }

    ///\brief Returns \c false if there are nodes
    ///to be processed in the queue
    ///
    ///Returns \c false if there are nodes
    ///to be processed in the queue
    bool emptyQueue() { return _stack_head<0; }
    ///Returns the number of the nodes to be processed.
    
    ///Returns the number of the nodes to be processed in the queue.
    int queueSize() { return _stack_head+1; }
    
    ///Executes the algorithm.

    ///Executes the algorithm.
    ///
    ///\pre init() must be called and at least one node should be added
    ///with addSource() before using this function.
    ///
    ///This method runs the %DFS algorithm from the root node(s)
    ///in order to
    ///compute the
    ///%DFS path to each node. The algorithm computes
    ///- The %DFS tree.
    ///- The distance of each node from the root(s) in the %DFS tree.
    ///
    void start()
    {
      while ( !emptyQueue() ) processNextEdge();
    }
    
    ///Executes the algorithm until \c dest is reached.

    ///Executes the algorithm until \c dest is reached.
    ///
    ///\pre init() must be called and at least one node should be added
    ///with addSource() before using this function.
    ///
    ///This method runs the %DFS algorithm from the root node(s)
    ///in order to
    ///compute the
    ///%DFS path to \c dest. The algorithm computes
    ///- The %DFS path to \c  dest.
    ///- The distance of \c dest from the root(s) in the %DFS tree.
    ///
    void start(Node dest)
    {
      while ( !emptyQueue() && G->target(_stack[_stack_head])!=dest ) 
	processNextEdge();
    }
    
    ///Executes the algorithm until a condition is met.

    ///Executes the algorithm until a condition is met.
    ///
    ///\pre init() must be called and at least one node should be added
    ///with addSource() before using this function.
    ///
    ///\param em must be a bool (or convertible) edge map. The algorithm
    ///will stop when it reaches an edge \c e with <tt>em[e]</tt> true.
    ///
    ///\return The reached edge \c e with <tt>em[e]</tt> true or
    ///\c INVALID if no such edge was found.
    ///
    ///\warning Contrary to \ref Bfs and \ref Dijkstra, \c em is an edge map,
    ///not a node map.
    template<class EM>
    Edge start(const EM &em)
    {
      while ( !emptyQueue() && !em[_stack[_stack_head]] )
        processNextEdge();
      return emptyQueue() ? INVALID : _stack[_stack_head];
    }

    ///Runs %DFS algorithm to visit all nodes in the graph.
    
    ///This method runs the %DFS algorithm in order to
    ///compute the
    ///%DFS path to each node. The algorithm computes
    ///- The %DFS tree.
    ///- The distance of each node from the root in the %DFS tree.
    ///
    ///\note d.run() is just a shortcut of the following code.
    ///\code
    ///  d.init();
    ///  for (NodeIt it(graph); it != INVALID; ++it) {
    ///    if (!d.reached(it)) {
    ///      d.addSource(it);
    ///      d.start();
    ///    }
    ///  }
    ///\endcode
    void run() {
      init();
      for (NodeIt it(*G); it != INVALID; ++it) {
        if (!reached(it)) {
          addSource(it);
          start();
        }
      }
    }

    ///Runs %DFS algorithm from node \c s.
    
    ///This method runs the %DFS algorithm from a root node \c s
    ///in order to
    ///compute the
    ///%DFS path to each node. The algorithm computes
    ///- The %DFS tree.
    ///- The distance of each node from the root in the %DFS tree.
    ///
    ///\note d.run(s) is just a shortcut of the following code.
    ///\code
    ///  d.init();
    ///  d.addSource(s);
    ///  d.start();
    ///\endcode
    void run(Node s) {
      init();
      addSource(s);
      start();
    }
    
    ///Finds the %DFS path between \c s and \c t.
    
    ///Finds the %DFS path between \c s and \c t.
    ///
    ///\return The length of the %DFS s---t path if there exists one,
    ///0 otherwise.
    ///\note Apart from the return value, d.run(s,t) is
    ///just a shortcut of the following code.
    ///\code
    ///  d.init();
    ///  d.addSource(s);
    ///  d.start(t);
    ///\endcode
    int run(Node s,Node t) {
      init();
      addSource(s);
      start(t);
      return reached(t)?_stack_head+1:0;
    }
    
    ///@}

    ///\name Query Functions
    ///The result of the %DFS algorithm can be obtained using these
    ///functions.\n
    ///Before the use of these functions,
    ///either run() or start() must be called.
    
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
    ///\warning If node \c v is unreachable from the root(s) then the return 
    ///value of this funcion is undefined.
    int dist(Node v) const { return (*_dist)[v]; }

    ///Returns the 'previous edge' of the %DFS tree.

    ///For a node \c v it returns the 'previous edge'
    ///of the %DFS path,
    ///i.e. it returns the last edge of a %DFS path from the root(s) to \c
    ///v. It is \ref INVALID
    ///if \c v is unreachable from the root(s) or \c v is a root. The
    ///%DFS tree used here is equal to the %DFS tree used in
    ///\ref predNode().
    ///\pre Either \ref run() or \ref start() must be called before using
    ///this function.
    Edge predEdge(Node v) const { return (*_pred)[v];}

    ///Returns the 'previous node' of the %DFS tree.

    ///For a node \c v it returns the 'previous node'
    ///of the %DFS tree,
    ///i.e. it returns the last but one node from a %DFS path from the
    ///root(s) to \c v.
    ///It is INVALID if \c v is unreachable from the root(s) or
    ///if \c v itself a root.
    ///The %DFS tree used here is equal to the %DFS
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
 
    ///Returns a reference to the %DFS edge-tree map.

    ///Returns a reference to the NodeMap of the edges of the
    ///%DFS tree.
    ///\pre Either \ref run() or \ref init()
    ///must be called before using this function.
    const PredMap &predMap() const { return *_pred;}
 
    ///Checks if a node is reachable from the root.

    ///Returns \c true if \c v is reachable from the root(s).
    ///\warning The source nodes are inditated as unreachable.
    ///\pre Either \ref run() or \ref start()
    ///must be called before using this function.
    ///
    bool reached(Node v) { return (*_reached)[v]; }
    
    ///@}
  };

  ///Default traits class of Dfs function.

  ///Default traits class of Dfs function.
  ///\param GR Graph type.
  template<class GR>
  struct DfsWizardDefaultTraits
  {
    ///The graph type the algorithm runs on. 
    typedef GR Graph;
    ///\brief The type of the map that stores the last
    ///edges of the %DFS paths.
    /// 
    ///The type of the map that stores the last
    ///edges of the %DFS paths.
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
  
  /// Default traits used by \ref DfsWizard

  /// To make it easier to use Dfs algorithm
  ///we have created a wizard class.
  /// This \ref DfsWizard class needs default traits,
  ///as well as the \ref Dfs class.
  /// The \ref DfsWizardBase is a class to be the default traits of the
  /// \ref DfsWizard class.
  template<class GR>
  class DfsWizardBase : public DfsWizardDefaultTraits<GR>
  {

    typedef DfsWizardDefaultTraits<GR> Base;
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
    DfsWizardBase() : _g(0), _reached(0), _processed(0), _pred(0),
			   _dist(0), _source(INVALID) {}

    /// Constructor.
    
    /// This constructor requires some parameters,
    /// listed in the parameters list.
    /// Others are initiated to 0.
    /// \param g is the initial value of  \ref _g
    /// \param s is the initial value of  \ref _source
    DfsWizardBase(const GR &g, Node s=INVALID) :
      _g(reinterpret_cast<void*>(const_cast<GR*>(&g))), 
      _reached(0), _processed(0), _pred(0), _dist(0), _source(s) {}

  };
  
  /// A class to make the usage of the Dfs algorithm easier

  /// This class is created to make it easier to use the Dfs algorithm.
  /// It uses the functions and features of the plain \ref Dfs,
  /// but it is much simpler to use it.
  ///
  /// Simplicity means that the way to change the types defined
  /// in the traits class is based on functions that returns the new class
  /// and not on templatable built-in classes.
  /// When using the plain \ref Dfs
  /// the new class with the modified type comes from
  /// the original class by using the ::
  /// operator. In the case of \ref DfsWizard only
  /// a function have to be called and it will
  /// return the needed class.
  ///
  /// It does not have own \ref run method. When its \ref run method is called
  /// it initiates a plain \ref Dfs object, and calls the \ref Dfs::run
  /// method of it.
  template<class TR>
  class DfsWizard : public TR
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
    ///edges of the %DFS paths.
    typedef typename TR::PredMap PredMap;
    ///The type of the map that stores the distances of the nodes.
    typedef typename TR::DistMap DistMap;

  public:
    /// Constructor.
    DfsWizard() : TR() {}

    /// Constructor that requires parameters.

    /// Constructor that requires parameters.
    /// These parameters will be the default values for the traits class.
    DfsWizard(const Graph &g, Node s=INVALID) :
      TR(g,s) {}

    ///Copy constructor
    DfsWizard(const TR &b) : TR(b) {}

    ~DfsWizard() {}

    ///Runs Dfs algorithm from a given node.
    
    ///Runs Dfs algorithm from a given node.
    ///The node can be given by the \ref source function.
    void run()
    {
      if(Base::_source==INVALID) throw UninitializedParameter();
      Dfs<Graph,TR> alg(*reinterpret_cast<const Graph*>(Base::_g));
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

    ///Runs Dfs algorithm from the given node.

    ///Runs Dfs algorithm from the given node.
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
    ///function for setting PredMap type
    ///
    /// \ref named-templ-param "Named parameter"
    ///function for setting PredMap type
    ///
    template<class T>
    DfsWizard<DefPredMapBase<T> > predMap(const T &t) 
    {
      Base::_pred=reinterpret_cast<void*>(const_cast<T*>(&t));
      return DfsWizard<DefPredMapBase<T> >(*this);
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
    DfsWizard<DefReachedMapBase<T> > reachedMap(const T &t) 
    {
      Base::_pred=reinterpret_cast<void*>(const_cast<T*>(&t));
      return DfsWizard<DefReachedMapBase<T> >(*this);
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
    DfsWizard<DefProcessedMapBase<T> > processedMap(const T &t) 
    {
      Base::_pred=reinterpret_cast<void*>(const_cast<T*>(&t));
      return DfsWizard<DefProcessedMapBase<T> >(*this);
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
    DfsWizard<DefDistMapBase<T> > distMap(const T &t) 
    {
      Base::_dist=reinterpret_cast<void*>(const_cast<T*>(&t));
      return DfsWizard<DefDistMapBase<T> >(*this);
    }
    
    /// Sets the source node, from which the Dfs algorithm runs.

    /// Sets the source node, from which the Dfs algorithm runs.
    /// \param s is the source node.
    DfsWizard<TR> &source(Node s) 
    {
      Base::_source=s;
      return *this;
    }
    
  };
  
  ///Function type interface for Dfs algorithm.

  ///\ingroup search
  ///Function type interface for Dfs algorithm.
  ///
  ///This function also has several
  ///\ref named-templ-func-param "named parameters",
  ///they are declared as the members of class \ref DfsWizard.
  ///The following
  ///example shows how to use these parameters.
  ///\code
  ///  dfs(g,source).predMap(preds).run();
  ///\endcode
  ///\warning Don't forget to put the \ref DfsWizard::run() "run()"
  ///to the end of the parameter list.
  ///\sa DfsWizard
  ///\sa Dfs
  template<class GR>
  DfsWizard<DfsWizardBase<GR> >
  dfs(const GR &g,typename GR::Node s=INVALID)
  {
    return DfsWizard<DfsWizardBase<GR> >(g,s);
  }

#ifdef DOXYGEN
  /// \brief Visitor class for dfs.
  ///  
  /// It gives a simple interface for a functional interface for dfs 
  /// traversal. The traversal on a linear data structure. 
  template <typename _Graph>
  struct DfsVisitor {
    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::Node Node;
    /// \brief Called when the edge reach a node.
    /// 
    /// It is called when the dfs find an edge which target is not
    /// reached yet.
    void discover(const Edge& edge) {}
    /// \brief Called when the node reached first time.
    /// 
    /// It is Called when the node reached first time.
    void reach(const Node& node) {}
    /// \brief Called when we step back on an edge.
    /// 
    /// It is called when the dfs should step back on the edge.
    void backtrack(const Edge& edge) {}
    /// \brief Called when we step back from the node.
    /// 
    /// It is called when we step back from the node.
    void leave(const Node& node) {}
    /// \brief Called when the edge examined but target of the edge 
    /// already discovered.
    /// 
    /// It called when the edge examined but the target of the edge 
    /// already discovered.
    void examine(const Edge& edge) {}
    /// \brief Called for the source node of the dfs.
    /// 
    /// It is called for the source node of the dfs.
    void start(const Node& node) {}
    /// \brief Called when we leave the source node of the dfs.
    /// 
    /// It is called when we leave the source node of the dfs.
    void stop(const Node& node) {}

  };
#else
  template <typename _Graph>
  struct DfsVisitor {
    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::Node Node;
    void discover(const Edge&) {}
    void reach(const Node&) {}
    void backtrack(const Edge&) {}
    void leave(const Node&) {}
    void examine(const Edge&) {}
    void start(const Node&) {}
    void stop(const Node&) {}

    template <typename _Visitor>
    struct Constraints {
      void constraints() {
	Edge edge;
	Node node;
	visitor.discover(edge);
	visitor.reach(node);
	visitor.backtrack(edge);
	visitor.leave(node);
	visitor.examine(edge);
	visitor.start(node);
	visitor.stop(edge);
      }
      _Visitor& visitor;
    };
  };
#endif

  /// \brief Default traits class of DfsVisit class.
  ///
  /// Default traits class of DfsVisit class.
  /// \param _Graph Graph type.
  template<class _Graph>
  struct DfsVisitDefaultTraits {

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
  
  /// %DFS Visit algorithm class.
  
  /// \ingroup search
  /// This class provides an efficient implementation of the %DFS algorithm
  /// with visitor interface.
  ///
  /// The %DfsVisit class provides an alternative interface to the Dfs
  /// class. It works with callback mechanism, the DfsVisit object calls
  /// on every dfs event the \c Visitor class member functions. 
  ///
  /// \param _Graph The graph type the algorithm runs on. The default value is
  /// \ref ListGraph. The value of _Graph is not used directly by Dfs, it
  /// is only passed to \ref DfsDefaultTraits.
  /// \param _Visitor The Visitor object for the algorithm. The 
  /// \ref DfsVisitor "DfsVisitor<_Graph>" is an empty Visitor which
  /// does not observe the Dfs events. If you want to observe the dfs
  /// events you should implement your own Visitor class.
  /// \param _Traits Traits class to set various data types used by the 
  /// algorithm. The default traits class is
  /// \ref DfsVisitDefaultTraits "DfsVisitDefaultTraits<_Graph>".
  /// See \ref DfsVisitDefaultTraits for the documentation of
  /// a Dfs visit traits class.
  ///
  /// \author Jacint Szabo, Alpar Juttner and Balazs Dezso
#ifdef DOXYGEN
  template <typename _Graph, typename _Visitor, typename _Traits>
#else
  template <typename _Graph = ListGraph,
	    typename _Visitor = DfsVisitor<_Graph>,
	    typename _Traits = DfsDefaultTraits<_Graph> >
#endif
  class DfsVisit {
  public:
    
    /// \brief \ref Exception for uninitialized parameters.
    ///
    /// This error represents problems in the initialization
    /// of the parameters of the algorithms.
    class UninitializedParameter : public lemon::UninitializedParameter {
    public:
      virtual const char* what() const throw() 
      {
	return "lemon::DfsVisit::UninitializedParameter";
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

    std::vector<typename Graph::Edge> _stack;
    int _stack_head;

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

    DfsVisit() {}
    
  public:

    typedef DfsVisit Create;

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
    struct DefReachedMap : public DfsVisit< Graph, Visitor,
					    DefReachedMapTraits<T> > {
      typedef DfsVisit< Graph, Visitor, DefReachedMapTraits<T> > Create;
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
    DfsVisit(const Graph& graph, Visitor& visitor) 
      : _graph(&graph), _visitor(&visitor),
	_reached(0), local_reached(false) {}
    
    /// \brief Destructor.
    ///
    /// Destructor.
    ~DfsVisit() {
      if(local_reached) delete _reached;
    }

    /// \brief Sets the map indicating if a node is reached.
    ///
    /// Sets the map indicating if a node is reached.
    /// If you don't use this function before calling \ref run(),
    /// it will allocate one. The destuctor deallocates this
    /// automatically allocated map, of course.
    /// \return <tt> (*this) </tt>
    DfsVisit &reachedMap(ReachedMap &m) {
      if(local_reached) {
	delete _reached;
	local_reached=false;
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
      _stack.resize(countNodes(*_graph));
      _stack_head = -1;
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
	  Edge e; 
	  _graph->firstOut(e, s);
	  if (e != INVALID) {
	    _stack[++_stack_head] = e;
	  } else {
	    _visitor->leave(s);
	  }
	}
    }
    
    /// \brief Processes the next edge.
    ///
    /// Processes the next edge.
    ///
    /// \return The processed edge.
    ///
    /// \pre The stack must not be empty!
    Edge processNextEdge() { 
      Edge e = _stack[_stack_head];
      Node m = _graph->target(e);
      if(!(*_reached)[m]) {
	_visitor->discover(e);
	_visitor->reach(m);
	_reached->set(m, true);
	_graph->firstOut(_stack[++_stack_head], m);
      } else {
	_visitor->examine(e);
	m = _graph->source(e);
	_graph->nextOut(_stack[_stack_head]);
      }
      while (_stack_head>=0 && _stack[_stack_head] == INVALID) {
	_visitor->leave(m);
	--_stack_head;
	if (_stack_head >= 0) {
	  _visitor->backtrack(_stack[_stack_head]);
	  m = _graph->source(_stack[_stack_head]);
	  _graph->nextOut(_stack[_stack_head]);
	} else {
	  _visitor->stop(m);	  
	}
      }
      return e;
    }

    /// \brief Next edge to be processed.
    ///
    /// Next edge to be processed.
    ///
    /// \return The next edge to be processed or INVALID if the stack is
    /// empty.
    Edge nextEdge() { 
      return _stack_head >= 0 ? _stack[_stack_head] : INVALID;
    }

    /// \brief Returns \c false if there are nodes
    /// to be processed in the queue
    ///
    /// Returns \c false if there are nodes
    /// to be processed in the queue
    bool emptyQueue() { return _stack_head < 0; }

    /// \brief Returns the number of the nodes to be processed.
    ///
    /// Returns the number of the nodes to be processed in the queue.
    int queueSize() { return _stack_head + 1; }
    
    /// \brief Executes the algorithm.
    ///
    /// Executes the algorithm.
    ///
    /// \pre init() must be called and at least one node should be added
    /// with addSource() before using this function.
    void start() {
      while ( !emptyQueue() ) processNextEdge();
    }
    
    /// \brief Executes the algorithm until \c dest is reached.
    ///
    /// Executes the algorithm until \c dest is reached.
    ///
    /// \pre init() must be called and at least one node should be added
    /// with addSource() before using this function.
    void start(Node dest) {
      while ( !emptyQueue() && _graph->target(_stack[_stack_head]) != dest ) 
	processNextEdge();
    }
    
    /// \brief Executes the algorithm until a condition is met.
    ///
    /// Executes the algorithm until a condition is met.
    ///
    /// \pre init() must be called and at least one node should be added
    /// with addSource() before using this function.
    ///
    /// \param em must be a bool (or convertible) edge map. The algorithm
    /// will stop when it reaches an edge \c e with <tt>em[e]</tt> true.
    ///
    ///\return The reached edge \c e with <tt>em[e]</tt> true or
    ///\c INVALID if no such edge was found.
    ///
    /// \warning Contrary to \ref Bfs and \ref Dijkstra, \c em is an edge map,
    /// not a node map.
    template <typename EM>
    Edge start(const EM &em) {
      while ( !emptyQueue() && !em[_stack[_stack_head]] )
        processNextEdge();
      return emptyQueue() ? INVALID : _stack[_stack_head];
    }

    /// \brief Runs %DFSVisit algorithm from node \c s.
    ///
    /// This method runs the %DFS algorithm from a root node \c s.
    /// \note d.run(s) is just a shortcut of the following code.
    ///\code
    ///   d.init();
    ///   d.addSource(s);
    ///   d.start();
    ///\endcode
    void run(Node s) {
      init();
      addSource(s);
      start();
    }

    /// \brief Runs %DFSVisit algorithm to visit all nodes in the graph.
    
    /// This method runs the %DFS algorithm in order to
    /// compute the %DFS path to each node. The algorithm computes
    /// - The %DFS tree.
    /// - The distance of each node from the root in the %DFS tree.
    ///
    ///\note d.run() is just a shortcut of the following code.
    ///\code
    ///  d.init();
    ///  for (NodeIt it(graph); it != INVALID; ++it) {
    ///    if (!d.reached(it)) {
    ///      d.addSource(it);
    ///      d.start();
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
    /// The result of the %DFS algorithm can be obtained using these
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


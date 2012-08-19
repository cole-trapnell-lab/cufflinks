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

#ifndef LEMON_LIST_GRAPH_H
#define LEMON_LIST_GRAPH_H

///\ingroup graphs
///\file
///\brief ListGraph, ListUGraph classes.

#include <lemon/bits/base_extender.h>
#include <lemon/bits/graph_extender.h>

#include <lemon/error.h>

#include <vector>
#include <list>

namespace lemon {

  class ListGraphBase {

  protected:
    struct NodeT {
      int first_in, first_out;
      int prev, next;
    };
 
    struct EdgeT {
      int target, source;
      int prev_in, prev_out;
      int next_in, next_out;
    };

    std::vector<NodeT> nodes;

    int first_node;

    int first_free_node;

    std::vector<EdgeT> edges;

    int first_free_edge;
    
  public:
    
    typedef ListGraphBase Graph;
    
    class Node {
      friend class ListGraphBase;
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

    class Edge {
      friend class ListGraphBase;
    protected:

      int id;
      explicit Edge(int pid) { id = pid;}

    public:
      Edge() {}
      Edge (Invalid) { id = -1; }
      bool operator==(const Edge& edge) const {return id == edge.id;}
      bool operator!=(const Edge& edge) const {return id != edge.id;}
      bool operator<(const Edge& edge) const {return id < edge.id;}
    };



    ListGraphBase()
      : nodes(), first_node(-1),
	first_free_node(-1), edges(), first_free_edge(-1) {}

    
    int maxNodeId() const { return nodes.size()-1; } 
    int maxEdgeId() const { return edges.size()-1; }

    Node source(Edge e) const { return Node(edges[e.id].source); }
    Node target(Edge e) const { return Node(edges[e.id].target); }


    void first(Node& node) const { 
      node.id = first_node;
    }

    void next(Node& node) const {
      node.id = nodes[node.id].next;
    }


    void first(Edge& e) const { 
      int n;
      for(n = first_node; 
	  n!=-1 && nodes[n].first_in == -1; 
	  n = nodes[n].next);
      e.id = (n == -1) ? -1 : nodes[n].first_in;
    }

    void next(Edge& edge) const {
      if (edges[edge.id].next_in != -1) {
	edge.id = edges[edge.id].next_in;
      } else {
	int n;
	for(n = nodes[edges[edge.id].target].next;
	  n!=-1 && nodes[n].first_in == -1; 
	  n = nodes[n].next);
	edge.id = (n == -1) ? -1 : nodes[n].first_in;
      }      
    }

    void firstOut(Edge &e, const Node& v) const {
      e.id = nodes[v.id].first_out;
    }
    void nextOut(Edge &e) const {
      e.id=edges[e.id].next_out;
    }

    void firstIn(Edge &e, const Node& v) const {
      e.id = nodes[v.id].first_in;
    }
    void nextIn(Edge &e) const {
      e.id=edges[e.id].next_in;
    }

    
    static int id(Node v) { return v.id; }
    static int id(Edge e) { return e.id; }

    static Node nodeFromId(int id) { return Node(id);}
    static Edge edgeFromId(int id) { return Edge(id);}

    Node addNode() {     
      int n;
      
      if(first_free_node==-1) {
	n = nodes.size();
	nodes.push_back(NodeT());
      } else {
	n = first_free_node;
	first_free_node = nodes[n].next;
      }
      
      nodes[n].next = first_node;
      if(first_node != -1) nodes[first_node].prev = n;
      first_node = n;
      nodes[n].prev = -1;
      
      nodes[n].first_in = nodes[n].first_out = -1;
      
      return Node(n);
    }
    
    Edge addEdge(Node u, Node v) {
      int n;      

      if (first_free_edge == -1) {
	n = edges.size();
	edges.push_back(EdgeT());
      } else {
	n = first_free_edge;
	first_free_edge = edges[n].next_in;
      }
      
      edges[n].source = u.id; 
      edges[n].target = v.id;

      edges[n].next_out = nodes[u.id].first_out;
      if(nodes[u.id].first_out != -1) {
	edges[nodes[u.id].first_out].prev_out = n;
      }
      
      edges[n].next_in = nodes[v.id].first_in;
      if(nodes[v.id].first_in != -1) {
	edges[nodes[v.id].first_in].prev_in = n;
      }
      
      edges[n].prev_in = edges[n].prev_out = -1;
	
      nodes[u.id].first_out = nodes[v.id].first_in = n;

      return Edge(n);
    }
    
    void erase(const Node& node) {
      int n = node.id;
      
      if(nodes[n].next != -1) {
	nodes[nodes[n].next].prev = nodes[n].prev;
      }
      
      if(nodes[n].prev != -1) {
	nodes[nodes[n].prev].next = nodes[n].next;
      } else {
	first_node = nodes[n].next;
      }
      
      nodes[n].next = first_free_node;
      first_free_node = n;

    }
    
    void erase(const Edge& edge) {
      int n = edge.id;
      
      if(edges[n].next_in!=-1) {
	edges[edges[n].next_in].prev_in = edges[n].prev_in;
      }

      if(edges[n].prev_in!=-1) {
	edges[edges[n].prev_in].next_in = edges[n].next_in;
      } else {
	nodes[edges[n].target].first_in = edges[n].next_in;
      }

      
      if(edges[n].next_out!=-1) {
	edges[edges[n].next_out].prev_out = edges[n].prev_out;
      } 

      if(edges[n].prev_out!=-1) {
	edges[edges[n].prev_out].next_out = edges[n].next_out;
      } else {
	nodes[edges[n].source].first_out = edges[n].next_out;
      }
      
      edges[n].next_in = first_free_edge;
      first_free_edge = n;      

    }

    void clear() {
      edges.clear();
      nodes.clear();
      first_node = first_free_node = first_free_edge = -1;
    }

  protected:
    void changeTarget(Edge e, Node n) 
    {
      if(edges[e.id].next_in != -1)
	edges[edges[e.id].next_in].prev_in = edges[e.id].prev_in;
      if(edges[e.id].prev_in != -1)
	edges[edges[e.id].prev_in].next_in = edges[e.id].next_in;
      else nodes[edges[e.id].target].first_in = edges[e.id].next_in;
      if (nodes[n.id].first_in != -1) {
	edges[nodes[n.id].first_in].prev_in = e.id;
      }
      edges[e.id].target = n.id;
      edges[e.id].prev_in = -1;
      edges[e.id].next_in = nodes[n.id].first_in;
      nodes[n.id].first_in = e.id;
    }
    void changeSource(Edge e, Node n) 
    {
      if(edges[e.id].next_out != -1)
	edges[edges[e.id].next_out].prev_out = edges[e.id].prev_out;
      if(edges[e.id].prev_out != -1)
	edges[edges[e.id].prev_out].next_out = edges[e.id].next_out;
      else nodes[edges[e.id].source].first_out = edges[e.id].next_out;
      if (nodes[n.id].first_out != -1) {
	edges[nodes[n.id].first_out].prev_out = e.id;
      }
      edges[e.id].source = n.id;
      edges[e.id].prev_out = -1;
      edges[e.id].next_out = nodes[n.id].first_out;
      nodes[n.id].first_out = e.id;
    }

  };

  typedef GraphExtender<ListGraphBase> ExtendedListGraphBase;

  /// \addtogroup graphs
  /// @{

  ///A list graph class.

  ///This is a simple and fast graph implementation.
  ///
  ///It conforms to the \ref concepts::Graph "Graph concept" and it
  ///also provides several additional useful extra functionalities.
  ///The most of the member functions and nested classes are
  ///documented only in the concept class.
  ///
  ///An important extra feature of this graph implementation is that
  ///its maps are real \ref concepts::ReferenceMap "reference map"s.
  ///
  ///\sa concepts::Graph.

  class ListGraph : public ExtendedListGraphBase {
  private:
    ///ListGraph is \e not copy constructible. Use GraphCopy() instead.
    
    ///ListGraph is \e not copy constructible. Use GraphCopy() instead.
    ///
    ListGraph(const ListGraph &) :ExtendedListGraphBase() {};
    ///\brief Assignment of ListGraph to another one is \e not allowed.
    ///Use GraphCopy() instead.

    ///Assignment of ListGraph to another one is \e not allowed.
    ///Use GraphCopy() instead.
    void operator=(const ListGraph &) {}
  public:

    typedef ExtendedListGraphBase Parent;

    /// Constructor
    
    /// Constructor.
    ///
    ListGraph() {}

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

    /// Changes the target of \c e to \c n

    /// Changes the target of \c e to \c n
    ///
    ///\note The <tt>EdgeIt</tt>s and <tt>OutEdgeIt</tt>s referencing
    ///the changed edge remain valid. However <tt>InEdgeIt</tt>s are
    ///invalidated.
    ///\warning This functionality cannot be used together with the Snapshot
    ///feature.
    void changeTarget(Edge e, Node n) { 
      Parent::changeTarget(e,n); 
    }
    /// Changes the source of \c e to \c n

    /// Changes the source of \c e to \c n
    ///
    ///\note The <tt>EdgeIt</tt>s and <tt>InEdgeIt</tt>s referencing
    ///the changed edge remain valid. However <tt>OutEdgeIt</tt>s are
    ///invalidated.
    ///\warning This functionality cannot be used together with the Snapshot
    ///feature.
    void changeSource(Edge e, Node n) { 
      Parent::changeSource(e,n);
    }

    /// Invert the direction of an edge.

    ///\note The <tt>EdgeIt</tt>s referencing the changed edge remain
    ///valid. However <tt>OutEdgeIt</tt>s and <tt>InEdgeIt</tt>s are
    ///invalidated.
    ///\warning This functionality cannot be used together with the Snapshot
    ///feature.
    void reverseEdge(Edge e) {
      Node t=target(e);
      changeTarget(e,source(e));
      changeSource(e,t);
    }

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

    ///Contract two nodes.

    ///This function contracts two nodes.
    ///
    ///Node \p b will be removed but instead of deleting
    ///incident edges, they will be joined to \p a.
    ///The last parameter \p r controls whether to remove loops. \c true
    ///means that loops will be removed.
    ///
    ///\note The <tt>EdgeIt</tt>s
    ///referencing a moved edge remain
    ///valid. However <tt>InEdgeIt</tt>s and <tt>OutEdgeIt</tt>s
    ///may be invalidated.
    ///\warning This functionality cannot be used together with the Snapshot
    ///feature.
    void contract(Node a, Node b, bool r = true) 
    {
      for(OutEdgeIt e(*this,b);e!=INVALID;) {
	OutEdgeIt f=e;
	++f;
	if(r && target(e)==a) erase(e);
	else changeSource(e,a);
	e=f;
      }
      for(InEdgeIt e(*this,b);e!=INVALID;) {
	InEdgeIt f=e;
	++f;
	if(r && source(e)==a) erase(e);
	else changeTarget(e,a);
	e=f;
      }
      erase(b);
    }

    ///Split a node.

    ///This function splits a node. First a new node is added to the graph,
    ///then the source of each outgoing edge of \c n is moved to this new node.
    ///If \c connect is \c true (this is the default value), then a new edge
    ///from \c n to the newly created node is also added.
    ///\return The newly created node.
    ///
    ///\note The <tt>EdgeIt</tt>s referencing a moved edge remain
    ///valid. However <tt>InEdgeIt</tt>s and <tt>OutEdgeIt</tt>s may
    ///be invalidated.  
    ///
    ///\warning This functionality cannot be used together with the
    ///Snapshot feature.  \todo It could be implemented in a bit
    ///faster way.
    Node split(Node n, bool connect = true) {
      Node b = addNode();
      for(OutEdgeIt e(*this,n);e!=INVALID;) {
 	OutEdgeIt f=e;
	++f;
	changeSource(e,b);
	e=f;
      }
      if (connect) addEdge(n,b);
      return b;
    }
      
    ///Split an edge.

    ///This function splits an edge. First a new node \c b is added to
    ///the graph, then the original edge is re-targeted to \c
    ///b. Finally an edge from \c b to the original target is added.
    ///\return The newly created node.  
    ///\warning This functionality
    ///cannot be used together with the Snapshot feature.
    Node split(Edge e) {
      Node b = addNode();
      addEdge(b,target(e));
      changeTarget(e,b);
      return b;
    }
      
    /// \brief Class to make a snapshot of the graph and restore
    /// to it later.
    ///
    /// Class to make a snapshot of the graph and to restore it
    /// later.
    ///
    /// The newly added nodes and edges can be removed using the
    /// restore() function.
    ///
    /// \warning Edge and node deletions cannot be restored. This
    /// events invalidate the snapshot. 
    class Snapshot {
    protected:

      typedef Parent::NodeNotifier NodeNotifier;

      class NodeObserverProxy : public NodeNotifier::ObserverBase {
      public:

        NodeObserverProxy(Snapshot& _snapshot)
          : snapshot(_snapshot) {}

        using NodeNotifier::ObserverBase::attach;
        using NodeNotifier::ObserverBase::detach;
        using NodeNotifier::ObserverBase::attached;
        
      protected:
        
        virtual void add(const Node& node) {
          snapshot.addNode(node);
        }
        virtual void add(const std::vector<Node>& nodes) {
          for (int i = nodes.size() - 1; i >= 0; ++i) {
            snapshot.addNode(nodes[i]);
          }
        }
        virtual void erase(const Node& node) {
          snapshot.eraseNode(node);
        }
        virtual void erase(const std::vector<Node>& nodes) {
          for (int i = 0; i < int(nodes.size()); ++i) {
            snapshot.eraseNode(nodes[i]);
          }
        }
        virtual void build() {
          Node node;
          std::vector<Node> nodes;
          for (notifier()->first(node); node != INVALID; 
               notifier()->next(node)) {
            nodes.push_back(node);
          }
          for (int i = nodes.size() - 1; i >= 0; --i) {
            snapshot.addNode(nodes[i]);
          }
        }
        virtual void clear() {
          Node node;
          for (notifier()->first(node); node != INVALID; 
               notifier()->next(node)) {
            snapshot.eraseNode(node);
          }
        }

        Snapshot& snapshot;
      };

      class EdgeObserverProxy : public EdgeNotifier::ObserverBase {
      public:

        EdgeObserverProxy(Snapshot& _snapshot)
          : snapshot(_snapshot) {}

        using EdgeNotifier::ObserverBase::attach;
        using EdgeNotifier::ObserverBase::detach;
        using EdgeNotifier::ObserverBase::attached;
        
      protected:

        virtual void add(const Edge& edge) {
          snapshot.addEdge(edge);
        }
        virtual void add(const std::vector<Edge>& edges) {
          for (int i = edges.size() - 1; i >= 0; ++i) {
            snapshot.addEdge(edges[i]);
          }
        }
        virtual void erase(const Edge& edge) {
          snapshot.eraseEdge(edge);
        }
        virtual void erase(const std::vector<Edge>& edges) {
          for (int i = 0; i < int(edges.size()); ++i) {
            snapshot.eraseEdge(edges[i]);
          }
        }
        virtual void build() {
          Edge edge;
          std::vector<Edge> edges;
          for (notifier()->first(edge); edge != INVALID; 
               notifier()->next(edge)) {
            edges.push_back(edge);
          }
          for (int i = edges.size() - 1; i >= 0; --i) {
            snapshot.addEdge(edges[i]);
          }
        }
        virtual void clear() {
          Edge edge;
          for (notifier()->first(edge); edge != INVALID; 
               notifier()->next(edge)) {
            snapshot.eraseEdge(edge);
          }
        }

        Snapshot& snapshot;
      };
      
      ListGraph *graph;

      NodeObserverProxy node_observer_proxy;
      EdgeObserverProxy edge_observer_proxy;

      std::list<Node> added_nodes;
      std::list<Edge> added_edges;


      void addNode(const Node& node) {
        added_nodes.push_front(node);        
      }
      void eraseNode(const Node& node) {
        std::list<Node>::iterator it = 
          std::find(added_nodes.begin(), added_nodes.end(), node);
        if (it == added_nodes.end()) {
          clear();
          edge_observer_proxy.detach();
          throw NodeNotifier::ImmediateDetach();
        } else {
          added_nodes.erase(it);
        }
      }

      void addEdge(const Edge& edge) {
        added_edges.push_front(edge);        
      }
      void eraseEdge(const Edge& edge) {
        std::list<Edge>::iterator it = 
          std::find(added_edges.begin(), added_edges.end(), edge);
        if (it == added_edges.end()) {
          clear();
          node_observer_proxy.detach(); 
          throw EdgeNotifier::ImmediateDetach();
        } else {
          added_edges.erase(it);
        }        
      }

      void attach(ListGraph &_graph) {
	graph = &_graph;
	node_observer_proxy.attach(graph->notifier(Node()));
        edge_observer_proxy.attach(graph->notifier(Edge()));
      }
            
      void detach() {
	node_observer_proxy.detach();
	edge_observer_proxy.detach();
      }

      bool attached() const {
        return node_observer_proxy.attached();
      }

      void clear() {
        added_nodes.clear();
        added_edges.clear();        
      }

    public:

      /// \brief Default constructor.
      ///
      /// Default constructor.
      /// To actually make a snapshot you must call save().
      Snapshot() 
        : graph(0), node_observer_proxy(*this), 
          edge_observer_proxy(*this) {}
      
      /// \brief Constructor that immediately makes a snapshot.
      ///      
      /// This constructor immediately makes a snapshot of the graph.
      /// \param _graph The graph we make a snapshot of.
      Snapshot(ListGraph &_graph) 
        : node_observer_proxy(*this), 
          edge_observer_proxy(*this) {
	attach(_graph);
      }
      
      /// \brief Make a snapshot.
      ///
      /// Make a snapshot of the graph.
      ///
      /// This function can be called more than once. In case of a repeated
      /// call, the previous snapshot gets lost.
      /// \param _graph The graph we make the snapshot of.
      void save(ListGraph &_graph) {
        if (attached()) {
          detach();
          clear();
        }
        attach(_graph);
      }
      
      /// \brief Undo the changes until the last snapshot.
      // 
      /// Undo the changes until the last snapshot created by save().
      void restore() {
	detach();
	for(std::list<Edge>::iterator it = added_edges.begin(); 
            it != added_edges.end(); ++it) {
	  graph->erase(*it);
	}
	for(std::list<Node>::iterator it = added_nodes.begin(); 
            it != added_nodes.end(); ++it) {
	  graph->erase(*it);
	}
        clear();
      }

      /// \brief Gives back true when the snapshot is valid.
      ///
      /// Gives back true when the snapshot is valid.
      bool valid() const {
        return attached();
      }
    };
    
  };

  ///@}

  class ListUGraphBase {

  protected:

    struct NodeT {
      int first_out;
      int prev, next;
    };
 
    struct EdgeT {
      int target;
      int prev_out, next_out;
    };

    std::vector<NodeT> nodes;

    int first_node;

    int first_free_node;

    std::vector<EdgeT> edges;

    int first_free_edge;
    
  public:
    
    typedef ListUGraphBase Graph;

    class Node;
    class Edge;
    class UEdge;
    
    class Node {
      friend class ListUGraphBase;
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
      friend class ListUGraphBase;
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
      friend class ListUGraphBase;
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



    ListUGraphBase()
      : nodes(), first_node(-1),
	first_free_node(-1), edges(), first_free_edge(-1) {}

    
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
      node.id = first_node;
    }

    void next(Node& node) const {
      node.id = nodes[node.id].next;
    }

    void first(Edge& e) const { 
      int n = first_node;
      while (n != -1 && nodes[n].first_out == -1) {
        n = nodes[n].next;
      }
      e.id = (n == -1) ? -1 : nodes[n].first_out;
    }

    void next(Edge& e) const {
      if (edges[e.id].next_out != -1) {
	e.id = edges[e.id].next_out;
      } else {
	int n = nodes[edges[e.id ^ 1].target].next;
        while(n != -1 && nodes[n].first_out == -1) {
          n = nodes[n].next;
        }
	e.id = (n == -1) ? -1 : nodes[n].first_out;
      }      
    }

    void first(UEdge& e) const { 
      int n = first_node;
      while (n != -1) {
        e.id = nodes[n].first_out;
        while ((e.id & 1) != 1) {
          e.id = edges[e.id].next_out;
        }
        if (e.id != -1) {
          e.id /= 2;
          return;
        } 
        n = nodes[n].next;
      }
      e.id = -1;
    }

    void next(UEdge& e) const {
      int n = edges[e.id * 2].target;
      e.id = edges[(e.id * 2) | 1].next_out;
      while ((e.id & 1) != 1) {
        e.id = edges[e.id].next_out;
      }
      if (e.id != -1) {
        e.id /= 2;
        return;
      } 
      n = nodes[n].next;
      while (n != -1) {
        e.id = nodes[n].first_out;
        while ((e.id & 1) != 1) {
          e.id = edges[e.id].next_out;
        }
        if (e.id != -1) {
          e.id /= 2;
          return;
        } 
        n = nodes[n].next;
      }
      e.id = -1;
    }

    void firstOut(Edge &e, const Node& v) const {
      e.id = nodes[v.id].first_out;
    }
    void nextOut(Edge &e) const {
      e.id = edges[e.id].next_out;
    }

    void firstIn(Edge &e, const Node& v) const {
      e.id = ((nodes[v.id].first_out) ^ 1);
      if (e.id == -2) e.id = -1;
    }
    void nextIn(Edge &e) const {
      e.id = ((edges[e.id ^ 1].next_out) ^ 1);
      if (e.id == -2) e.id = -1;
    }

    void firstInc(UEdge &e, bool& d, const Node& v) const {
      int de = nodes[v.id].first_out;
      if (de != -1 ) {
        e.id = de / 2;
        d = ((de & 1) == 1);
      } else {
        e.id = -1;
        d = true;
      }
    }
    void nextInc(UEdge &e, bool& d) const {
      int de = (edges[(e.id * 2) | (d ? 1 : 0)].next_out);
      if (de != -1 ) {
        e.id = de / 2;
        d = ((de & 1) == 1);
      } else {
        e.id = -1;
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
      int n;
      
      if(first_free_node==-1) {
	n = nodes.size();
	nodes.push_back(NodeT());
      } else {
	n = first_free_node;
	first_free_node = nodes[n].next;
      }
      
      nodes[n].next = first_node;
      if (first_node != -1) nodes[first_node].prev = n;
      first_node = n;
      nodes[n].prev = -1;
      
      nodes[n].first_out = -1;
      
      return Node(n);
    }
    
    UEdge addEdge(Node u, Node v) {
      int n;      

      if (first_free_edge == -1) {
	n = edges.size();
	edges.push_back(EdgeT());
	edges.push_back(EdgeT());
      } else {
	n = first_free_edge;
	first_free_edge = edges[n].next_out;
      }
      
      edges[n].target = u.id;
      edges[n | 1].target = v.id;

      edges[n].next_out = nodes[v.id].first_out;
      if (nodes[v.id].first_out != -1) {
	edges[nodes[v.id].first_out].prev_out = n;
      }      
      edges[n].prev_out = -1;
      nodes[v.id].first_out = n;
      
      edges[n | 1].next_out = nodes[u.id].first_out;
      if (nodes[u.id].first_out != -1) {
	edges[nodes[u.id].first_out].prev_out = (n | 1);
      }
      edges[n | 1].prev_out = -1;      
      nodes[u.id].first_out = (n | 1);

      return UEdge(n / 2);
    }
    
    void erase(const Node& node) {
      int n = node.id;
      
      if(nodes[n].next != -1) {
	nodes[nodes[n].next].prev = nodes[n].prev;
      }
      
      if(nodes[n].prev != -1) {
	nodes[nodes[n].prev].next = nodes[n].next;
      } else {
	first_node = nodes[n].next;
      }
      
      nodes[n].next = first_free_node;
      first_free_node = n;

    }
    
    void erase(const UEdge& edge) {
      int n = edge.id * 2;
      
      if (edges[n].next_out != -1) {
	edges[edges[n].next_out].prev_out = edges[n].prev_out;
      } 

      if (edges[n].prev_out != -1) {
	edges[edges[n].prev_out].next_out = edges[n].next_out;
      } else {
	nodes[edges[n | 1].target].first_out = edges[n].next_out;
      }

      if (edges[n | 1].next_out != -1) {
	edges[edges[n | 1].next_out].prev_out = edges[n | 1].prev_out;
      } 

      if (edges[n | 1].prev_out != -1) {
	edges[edges[n | 1].prev_out].next_out = edges[n | 1].next_out;
      } else {
	nodes[edges[n].target].first_out = edges[n | 1].next_out;
      }
      
      edges[n].next_out = first_free_edge;
      first_free_edge = n;      

    }

    void clear() {
      edges.clear();
      nodes.clear();
      first_node = first_free_node = first_free_edge = -1;
    }

  protected:

    void changeTarget(UEdge e, Node n) {
      if(edges[2 * e.id].next_out != -1) {
	edges[edges[2 * e.id].next_out].prev_out = edges[2 * e.id].prev_out;
      }
      if(edges[2 * e.id].prev_out != -1) {
	edges[edges[2 * e.id].prev_out].next_out = 
          edges[2 * e.id].next_out;
      } else {
        nodes[edges[(2 * e.id) | 1].target].first_out = 
          edges[2 * e.id].next_out;
      }

      if (nodes[n.id].first_out != -1) {
	edges[nodes[n.id].first_out].prev_out = 2 * e.id;
      }
      edges[(2 * e.id) | 1].target = n.id;
      edges[2 * e.id].prev_out = -1;
      edges[2 * e.id].next_out = nodes[n.id].first_out;
      nodes[n.id].first_out = 2 * e.id;
    }

    void changeSource(UEdge e, Node n) {
      if(edges[(2 * e.id) | 1].next_out != -1) {
	edges[edges[(2 * e.id) | 1].next_out].prev_out = 
          edges[(2 * e.id) | 1].prev_out;
      }
      if(edges[(2 * e.id) | 1].prev_out != -1) {
	edges[edges[(2 * e.id) | 1].prev_out].next_out = 
          edges[(2 * e.id) | 1].next_out;
      } else {
        nodes[edges[2 * e.id].target].first_out = 
          edges[(2 * e.id) | 1].next_out;
      }

      if (nodes[n.id].first_out != -1) {
	edges[nodes[n.id].first_out].prev_out = ((2 * e.id) | 1);
      }
      edges[2 * e.id].target = n.id;
      edges[(2 * e.id) | 1].prev_out = -1;
      edges[(2 * e.id) | 1].next_out = nodes[n.id].first_out;
      nodes[n.id].first_out = ((2 * e.id) | 1);
    }

  };

//   typedef UGraphExtender<UndirGraphExtender<ListGraphBase> > 
//   ExtendedListUGraphBase;

  typedef UGraphExtender<ListUGraphBase> ExtendedListUGraphBase;



  /// \addtogroup graphs
  /// @{

  ///An undirected list graph class.

  ///This is a simple and fast undirected graph implementation.
  ///
  ///An important extra feature of this graph implementation is that
  ///its maps are real \ref concepts::ReferenceMap "reference map"s.
  ///
  ///It conforms to the
  ///\ref concepts::UGraph "UGraph concept".
  ///
  ///\sa concepts::UGraph.
  ///
  class ListUGraph : public ExtendedListUGraphBase {
  private:
    ///ListUGraph is \e not copy constructible. Use UGraphCopy() instead.

    ///ListUGraph is \e not copy constructible. Use UGraphCopy() instead.
    ///
    ListUGraph(const ListUGraph &) :ExtendedListUGraphBase()  {};
    ///\brief Assignment of ListUGraph to another one is \e not allowed.
    ///Use UGraphCopy() instead.

    ///Assignment of ListUGraph to another one is \e not allowed.
    ///Use UGraphCopy() instead.
    void operator=(const ListUGraph &) {}
  public:
    /// Constructor
    
    /// Constructor.
    ///
    ListUGraph() {}

    typedef ExtendedListUGraphBase Parent;

    typedef Parent::OutEdgeIt IncEdgeIt;

    /// \brief Add a new node to the graph.
    ///
    /// \return the new node.
    ///
    Node addNode() { return Parent::addNode(); }

    /// \brief Add a new edge to the graph.
    ///
    /// Add a new edge to the graph with source node \c s
    /// and target node \c t.
    /// \return the new undirected edge.
    UEdge addEdge(const Node& s, const Node& t) { 
      return Parent::addEdge(s, t); 
    }
    /// \brief Changes the source of \c e to \c n
    ///
    /// Changes the source of \c e to \c n
    ///
    ///\note The <tt>EdgeIt</tt>s and <tt>InEdgeIt</tt>s
    ///referencing the changed edge remain
    ///valid. However <tt>OutEdgeIt</tt>s are invalidated.
    void changeSource(UEdge e, Node n) { 
      Parent::changeSource(e,n); 
    }    
    /// \brief Changes the target of \c e to \c n
    ///
    /// Changes the target of \c e to \c n
    ///
    /// \note The <tt>EdgeIt</tt>s referencing the changed edge remain
    /// valid. However the other iterators may be invalidated.
    void changeTarget(UEdge e, Node n) { 
      Parent::changeTarget(e,n); 
    }
    /// \brief Changes the source of \c e to \c n
    ///
    /// Changes the source of \c e to \c n. It changes the proper
    /// node of the represented undirected edge.
    ///
    ///\note The <tt>EdgeIt</tt>s and <tt>InEdgeIt</tt>s
    ///referencing the changed edge remain
    ///valid. However <tt>OutEdgeIt</tt>s are invalidated.
    void changeSource(Edge e, Node n) { 
      if (Parent::direction(e)) {
        Parent::changeSource(e,n);
      } else {
        Parent::changeTarget(e,n);
      } 
    }
    /// \brief Changes the target of \c e to \c n
    ///
    /// Changes the target of \c e to \c n. It changes the proper
    /// node of the represented undirected edge.
    ///
    ///\note The <tt>EdgeIt</tt>s and <tt>OutEdgeIt</tt>s
    ///referencing the changed edge remain
    ///valid. However <tt>InEdgeIt</tt>s are invalidated.
    void changeTarget(Edge e, Node n) { 
      if (Parent::direction(e)) {
        Parent::changeTarget(e,n);
      } else {
        Parent::changeSource(e,n);
      } 
    }
    /// \brief Contract two nodes.
    ///
    /// This function contracts two nodes.
    ///
    /// Node \p b will be removed but instead of deleting
    /// its neighboring edges, they will be joined to \p a.
    /// The last parameter \p r controls whether to remove loops. \c true
    /// means that loops will be removed.
    ///
    /// \note The <tt>EdgeIt</tt>s referencing a moved edge remain
    /// valid.
    void contract(Node a, Node b, bool r = true) {
      for(IncEdgeIt e(*this, b); e!=INVALID;) {
	IncEdgeIt f = e; ++f;
	if (r && runningNode(e) == a) {
	  erase(e);
	} else if (source(e) == b) {
	  changeSource(e, a);
	} else {
	  changeTarget(e, a);
	}
	e = f;
      }
      erase(b);
    }


    /// \brief Class to make a snapshot of the graph and restore
    /// to it later.
    ///
    /// Class to make a snapshot of the graph and to restore it
    /// later.
    ///
    /// The newly added nodes and undirected edges can be removed
    /// using the restore() function.
    ///
    /// \warning Edge and node deletions cannot be restored. This
    /// events invalidate the snapshot. 
    class Snapshot {
    protected:

      typedef Parent::NodeNotifier NodeNotifier;

      class NodeObserverProxy : public NodeNotifier::ObserverBase {
      public:

        NodeObserverProxy(Snapshot& _snapshot)
          : snapshot(_snapshot) {}

        using NodeNotifier::ObserverBase::attach;
        using NodeNotifier::ObserverBase::detach;
        using NodeNotifier::ObserverBase::attached;
        
      protected:
        
        virtual void add(const Node& node) {
          snapshot.addNode(node);
        }
        virtual void add(const std::vector<Node>& nodes) {
          for (int i = nodes.size() - 1; i >= 0; ++i) {
            snapshot.addNode(nodes[i]);
          }
        }
        virtual void erase(const Node& node) {
          snapshot.eraseNode(node);
        }
        virtual void erase(const std::vector<Node>& nodes) {
          for (int i = 0; i < int(nodes.size()); ++i) {
            snapshot.eraseNode(nodes[i]);
          }
        }
        virtual void build() {
          Node node;
          std::vector<Node> nodes;
          for (notifier()->first(node); node != INVALID; 
               notifier()->next(node)) {
            nodes.push_back(node);
          }
          for (int i = nodes.size() - 1; i >= 0; --i) {
            snapshot.addNode(nodes[i]);
          }
        }
        virtual void clear() {
          Node node;
          for (notifier()->first(node); node != INVALID; 
               notifier()->next(node)) {
            snapshot.eraseNode(node);
          }
        }

        Snapshot& snapshot;
      };

      class UEdgeObserverProxy : public UEdgeNotifier::ObserverBase {
      public:

        UEdgeObserverProxy(Snapshot& _snapshot)
          : snapshot(_snapshot) {}

        using UEdgeNotifier::ObserverBase::attach;
        using UEdgeNotifier::ObserverBase::detach;
        using UEdgeNotifier::ObserverBase::attached;
        
      protected:

        virtual void add(const UEdge& edge) {
          snapshot.addUEdge(edge);
        }
        virtual void add(const std::vector<UEdge>& edges) {
          for (int i = edges.size() - 1; i >= 0; ++i) {
            snapshot.addUEdge(edges[i]);
          }
        }
        virtual void erase(const UEdge& edge) {
          snapshot.eraseUEdge(edge);
        }
        virtual void erase(const std::vector<UEdge>& edges) {
          for (int i = 0; i < int(edges.size()); ++i) {
            snapshot.eraseUEdge(edges[i]);
          }
        }
        virtual void build() {
          UEdge edge;
          std::vector<UEdge> edges;
          for (notifier()->first(edge); edge != INVALID; 
               notifier()->next(edge)) {
            edges.push_back(edge);
          }
          for (int i = edges.size() - 1; i >= 0; --i) {
            snapshot.addUEdge(edges[i]);
          }
        }
        virtual void clear() {
          UEdge edge;
          for (notifier()->first(edge); edge != INVALID; 
               notifier()->next(edge)) {
            snapshot.eraseUEdge(edge);
          }
        }

        Snapshot& snapshot;
      };
      
      ListUGraph *graph;

      NodeObserverProxy node_observer_proxy;
      UEdgeObserverProxy edge_observer_proxy;

      std::list<Node> added_nodes;
      std::list<UEdge> added_edges;


      void addNode(const Node& node) {
        added_nodes.push_front(node);        
      }
      void eraseNode(const Node& node) {
        std::list<Node>::iterator it = 
          std::find(added_nodes.begin(), added_nodes.end(), node);
        if (it == added_nodes.end()) {
          clear();
          edge_observer_proxy.detach();
          throw NodeNotifier::ImmediateDetach();
        } else {
          added_nodes.erase(it);
        }
      }

      void addUEdge(const UEdge& edge) {
        added_edges.push_front(edge);        
      }
      void eraseUEdge(const UEdge& edge) {
        std::list<UEdge>::iterator it = 
          std::find(added_edges.begin(), added_edges.end(), edge);
        if (it == added_edges.end()) {
          clear();
          node_observer_proxy.detach();
          throw UEdgeNotifier::ImmediateDetach();
        } else {
          added_edges.erase(it);
        }        
      }

      void attach(ListUGraph &_graph) {
	graph = &_graph;
	node_observer_proxy.attach(graph->notifier(Node()));
        edge_observer_proxy.attach(graph->notifier(UEdge()));
      }
            
      void detach() {
	node_observer_proxy.detach();
	edge_observer_proxy.detach();
      }

      bool attached() const {
        return node_observer_proxy.attached();
      }

      void clear() {
        added_nodes.clear();
        added_edges.clear();        
      }

    public:

      /// \brief Default constructor.
      ///
      /// Default constructor.
      /// To actually make a snapshot you must call save().
      Snapshot() 
        : graph(0), node_observer_proxy(*this), 
          edge_observer_proxy(*this) {}
      
      /// \brief Constructor that immediately makes a snapshot.
      ///      
      /// This constructor immediately makes a snapshot of the graph.
      /// \param _graph The graph we make a snapshot of.
      Snapshot(ListUGraph &_graph) 
        : node_observer_proxy(*this), 
          edge_observer_proxy(*this) {
	attach(_graph);
      }
      
      /// \brief Make a snapshot.
      ///
      /// Make a snapshot of the graph.
      ///
      /// This function can be called more than once. In case of a repeated
      /// call, the previous snapshot gets lost.
      /// \param _graph The graph we make the snapshot of.
      void save(ListUGraph &_graph) {
        if (attached()) {
          detach();
          clear();
        }
        attach(_graph);
      }
      
      /// \brief Undo the changes until the last snapshot.
      // 
      /// Undo the changes until the last snapshot created by save().
      void restore() {
	detach();
	for(std::list<UEdge>::iterator it = added_edges.begin(); 
            it != added_edges.end(); ++it) {
	  graph->erase(*it);
	}
	for(std::list<Node>::iterator it = added_nodes.begin(); 
            it != added_nodes.end(); ++it) {
	  graph->erase(*it);
	}
        clear();
      }

      /// \brief Gives back true when the snapshot is valid.
      ///
      /// Gives back true when the snapshot is valid.
      bool valid() const {
        return attached();
      }
    };
  };


  class ListBpUGraphBase {
  public:

    class NodeSetError : public LogicError {
    public:
      virtual const char* what() const throw() { 
	return "lemon::ListBpUGraph::NodeSetError";
      }
    };

  protected:

    struct NodeT {
      int first_edge, prev, next;
    };

    struct UEdgeT {
      int aNode, prev_out, next_out;
      int bNode, prev_in, next_in;
    };

    std::vector<NodeT> aNodes;
    std::vector<NodeT> bNodes;

    std::vector<UEdgeT> edges;

    int first_anode;
    int first_free_anode;

    int first_bnode;
    int first_free_bnode;

    int first_free_edge;

  public:
  
    class Node {
      friend class ListBpUGraphBase;
    protected:
      int id;

      explicit Node(int _id) : id(_id) {}
    public:
      Node() {}
      Node(Invalid) { id = -1; }
      bool operator==(const Node i) const {return id==i.id;}
      bool operator!=(const Node i) const {return id!=i.id;}
      bool operator<(const Node i) const {return id<i.id;}
    };

    class UEdge {
      friend class ListBpUGraphBase;
    protected:
      int id;

      explicit UEdge(int _id) { id = _id;}
    public:
      UEdge() {}
      UEdge (Invalid) { id = -1; }
      bool operator==(const UEdge i) const {return id==i.id;}
      bool operator!=(const UEdge i) const {return id!=i.id;}
      bool operator<(const UEdge i) const {return id<i.id;}
    };

    ListBpUGraphBase()
      : first_anode(-1), first_free_anode(-1),
        first_bnode(-1), first_free_bnode(-1),
        first_free_edge(-1) {}

    void firstANode(Node& node) const {
      node.id = first_anode != -1 ? (first_anode << 1) : -1;
    }
    void nextANode(Node& node) const {
      node.id = aNodes[node.id >> 1].next;
    }

    void firstBNode(Node& node) const {
      node.id = first_bnode != -1 ? (first_bnode << 1) + 1 : -1;
    }
    void nextBNode(Node& node) const {
      node.id = bNodes[node.id >> 1].next;
    }

    void first(Node& node) const {
      if (first_anode != -1) {
        node.id = (first_anode << 1);
      } else if (first_bnode != -1) {
        node.id = (first_bnode << 1) + 1;
      } else {
        node.id = -1;
      }
    }
    void next(Node& node) const {
      if (aNode(node)) {
        node.id = aNodes[node.id >> 1].next;
        if (node.id == -1) {
          if (first_bnode != -1) {
            node.id = (first_bnode << 1) + 1;
          }
        }
      } else {
        node.id = bNodes[node.id >> 1].next;
      }
    }
  
    void first(UEdge& edge) const {
      int aid = first_anode;
      while (aid != -1 && aNodes[aid].first_edge == -1) {
        aid = aNodes[aid].next != -1 ? 
          aNodes[aid].next >> 1 : -1;
      }
      if (aid != -1) {
        edge.id = aNodes[aid].first_edge;
      } else {
        edge.id = -1;
      }
    }
    void next(UEdge& edge) const {
      int aid = edges[edge.id].aNode >> 1;
      edge.id = edges[edge.id].next_out;
      if (edge.id == -1) {
        aid = aNodes[aid].next != -1 ? 
          aNodes[aid].next >> 1 : -1;
        while (aid != -1 && aNodes[aid].first_edge == -1) {
          aid = aNodes[aid].next != -1 ? 
          aNodes[aid].next >> 1 : -1;
        }
        if (aid != -1) {
          edge.id = aNodes[aid].first_edge;
        } else {
          edge.id = -1;
        }
      }
    }

    void firstFromANode(UEdge& edge, const Node& node) const {
      LEMON_ASSERT((node.id & 1) == 0, NodeSetError());
      edge.id = aNodes[node.id >> 1].first_edge;
    }
    void nextFromANode(UEdge& edge) const {
      edge.id = edges[edge.id].next_out;
    }

    void firstFromBNode(UEdge& edge, const Node& node) const {
      LEMON_ASSERT((node.id & 1) == 1, NodeSetError());
      edge.id = bNodes[node.id >> 1].first_edge;
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
      int aid;
      if (first_free_anode == -1) {
        aid = aNodes.size();
        aNodes.push_back(NodeT());
      } else {
        aid = first_free_anode;
        first_free_anode = aNodes[first_free_anode].next;
      }
      if (first_anode != -1) {
        aNodes[aid].next = first_anode << 1;
        aNodes[first_anode].prev = aid << 1;
      } else {
        aNodes[aid].next = -1;
      }
      aNodes[aid].prev = -1;
      first_anode = aid;
      aNodes[aid].first_edge = -1;
      return Node(aid << 1);
    }

    Node addBNode() {
      int bid;
      if (first_free_bnode == -1) {
        bid = bNodes.size();
        bNodes.push_back(NodeT());
      } else {
        bid = first_free_bnode;
        first_free_bnode = bNodes[first_free_bnode].next;
      }
      if (first_bnode != -1) {
        bNodes[bid].next = (first_bnode << 1) + 1;
        bNodes[first_bnode].prev = (bid << 1) + 1;
      } else {
        bNodes[bid].next = -1;
      }
      bNodes[bid].prev = -1;
      first_bnode = bid;
      bNodes[bid].first_edge = -1;
      return Node((bid << 1) + 1);
    }

    UEdge addEdge(const Node& source, const Node& target) {
      LEMON_ASSERT(((source.id ^ target.id) & 1) == 1, NodeSetError());
      int edgeId;
      if (first_free_edge != -1) {
        edgeId = first_free_edge;
        first_free_edge = edges[edgeId].next_out;
      } else {
        edgeId = edges.size();
        edges.push_back(UEdgeT());
      }
      if ((source.id & 1) == 0) {
	edges[edgeId].aNode = source.id;
	edges[edgeId].bNode = target.id;
      } else {
	edges[edgeId].aNode = target.id;
	edges[edgeId].bNode = source.id;
      }
      edges[edgeId].next_out = aNodes[edges[edgeId].aNode >> 1].first_edge;
      edges[edgeId].prev_out = -1;
      if (aNodes[edges[edgeId].aNode >> 1].first_edge != -1) {
        edges[aNodes[edges[edgeId].aNode >> 1].first_edge].prev_out = edgeId;
      }
      aNodes[edges[edgeId].aNode >> 1].first_edge = edgeId;
      edges[edgeId].next_in = bNodes[edges[edgeId].bNode >> 1].first_edge;
      edges[edgeId].prev_in = -1;
      if (bNodes[edges[edgeId].bNode >> 1].first_edge != -1) {
        edges[bNodes[edges[edgeId].bNode >> 1].first_edge].prev_in = edgeId;
      }
      bNodes[edges[edgeId].bNode >> 1].first_edge = edgeId;
      return UEdge(edgeId);
    }

    void erase(const Node& node) {
      if (aNode(node)) {
        int aid = node.id >> 1;
        if (aNodes[aid].prev != -1) {
          aNodes[aNodes[aid].prev >> 1].next = aNodes[aid].next;
        } else {
          first_anode = 
            aNodes[aid].next != -1 ? aNodes[aid].next >> 1 : -1;
        }
        if (aNodes[aid].next != -1) {
          aNodes[aNodes[aid].next >> 1].prev = aNodes[aid].prev;
        }
        aNodes[aid].next = first_free_anode;
        first_free_anode = aid;
      } else {
        int bid = node.id >> 1;
        if (bNodes[bid].prev != -1) {
          bNodes[bNodes[bid].prev >> 1].next = bNodes[bid].next;
        } else {
          first_bnode = 
            bNodes[bid].next != -1 ? bNodes[bid].next >> 1 : -1;
        }
        if (bNodes[bid].next != -1) {
          bNodes[bNodes[bid].next >> 1].prev = bNodes[bid].prev;
        }
        bNodes[bid].next = first_free_bnode;
        first_free_bnode = bid;
      }
    }

    void erase(const UEdge& edge) {

      if (edges[edge.id].prev_out != -1) {
        edges[edges[edge.id].prev_out].next_out = edges[edge.id].next_out;
      } else {
        aNodes[edges[edge.id].aNode >> 1].first_edge = edges[edge.id].next_out;
      }
      if (edges[edge.id].next_out != -1) {
        edges[edges[edge.id].next_out].prev_out = edges[edge.id].prev_out;
      }

      if (edges[edge.id].prev_in != -1) {
        edges[edges[edge.id].prev_in].next_in = edges[edge.id].next_in;
      } else {
        bNodes[edges[edge.id].bNode >> 1].first_edge = edges[edge.id].next_in;
      }
      if (edges[edge.id].next_in != -1) {
        edges[edges[edge.id].next_in].prev_in = edges[edge.id].prev_in;
      }

      edges[edge.id].next_out = first_free_edge;
      first_free_edge = edge.id;
    }
 
    void clear() {
      aNodes.clear();
      bNodes.clear();
      edges.clear();
      first_anode = -1;
      first_free_anode = -1;
      first_bnode = -1;
      first_free_bnode = -1;
      first_free_edge = -1;
    }

    void changeANode(const UEdge& edge, const Node& node) {
      LEMON_ASSERT((node.id & 1) == 0, NodeSetError());
      if (edges[edge.id].prev_out != -1) {
        edges[edges[edge.id].prev_out].next_out = edges[edge.id].next_out;
      } else {
        aNodes[edges[edge.id].aNode >> 1].first_edge = edges[edge.id].next_out;
      }
      if (edges[edge.id].next_out != -1) {
        edges[edges[edge.id].next_out].prev_out = edges[edge.id].prev_out;  
      }
      if (aNodes[node.id >> 1].first_edge != -1) {
        edges[aNodes[node.id >> 1].first_edge].prev_out = edge.id;
      }
      edges[edge.id].prev_out = -1;
      edges[edge.id].next_out = aNodes[node.id >> 1].first_edge;
      aNodes[node.id >> 1].first_edge = edge.id;
      edges[edge.id].aNode = node.id;
    } 

    void changeBNode(const UEdge& edge, const Node& node) {
      LEMON_ASSERT((node.id & 1) == 1, NodeSetError());
      if (edges[edge.id].prev_in != -1) {
        edges[edges[edge.id].prev_in].next_in = edges[edge.id].next_in;
      } else {
        bNodes[edges[edge.id].bNode >> 1].first_edge = edges[edge.id].next_in;
      }
      if (edges[edge.id].next_in != -1) {
        edges[edges[edge.id].next_in].prev_in = edges[edge.id].prev_in;  
      }
      if (bNodes[node.id >> 1].first_edge != -1) {
        edges[bNodes[node.id >> 1].first_edge].prev_in = edge.id;
      }
      edges[edge.id].prev_in = -1;
      edges[edge.id].next_in = bNodes[node.id >> 1].first_edge;
      bNodes[node.id >> 1].first_edge = edge.id;
      edges[edge.id].bNode = node.id;
    } 

  };


  typedef BpUGraphExtender<BidirBpUGraphExtender<ListBpUGraphBase> > 
  ExtendedListBpUGraphBase;

  /// \ingroup graphs
  ///
  /// \brief A smart bipartite undirected graph class.
  ///
  /// This is a bipartite undirected graph implementation.
  /// It is conforms to the \ref concepts::BpUGraph "BpUGraph concept".
  ///
  ///An important extra feature of this graph implementation is that
  ///its maps are real \ref concepts::ReferenceMap "reference map"s.
  ///
  /// \sa concepts::BpUGraph.
  ///
  class ListBpUGraph : public ExtendedListBpUGraphBase {
    /// \brief ListBpUGraph is \e not copy constructible.
    ///
    ///ListBpUGraph is \e not copy constructible.
    ListBpUGraph(const ListBpUGraph &) :ExtendedListBpUGraphBase()  {};
    /// \brief Assignment of ListBpUGraph to another one is \e not
    /// allowed.
    ///
    /// Assignment of ListBpUGraph to another one is \e not allowed.
    void operator=(const ListBpUGraph &) {}
  public:
    /// \brief Constructor
    ///    
    /// Constructor.
    ///
    ListBpUGraph() {}

    typedef ExtendedListBpUGraphBase Parent;
    /// \brief Add a new ANode to the graph.
    ///
    /// \return the new node.
    ///
    Node addANode() { return Parent::addANode(); }

    /// \brief Add a new BNode to the graph.
    ///
    /// \return the new node.
    ///
    Node addBNode() { return Parent::addBNode(); }

    /// \brief Add a new edge to the graph.
    ///
    /// Add a new edge to the graph with an ANode and a BNode.
    /// \return the new undirected edge.
    UEdge addEdge(const Node& s, const Node& t) { 
      return Parent::addEdge(s, t); 
    }

    /// \brief Changes the ANode of \c e to \c n
    ///
    /// Changes the ANode of \c e to \c n
    ///
    ///\note The <tt>EdgeIt</tt>s and <tt>InEdgeIt</tt>s referencing
    ///the changed edge remain valid. However <tt>OutEdgeIt</tt>s are
    ///invalidated.
    void changeANode(UEdge e, Node n) { 
      Parent::changeANode(e,n); 
    }

    /// \brief Changes the BNode of \c e to \c n
    ///
    /// Changes the BNode of \c e to \c n
    ///
    /// \note The <tt>EdgeIt</tt>s and <tt>OutEdgeIt</tt>s
    /// referencing the changed edge remain
    /// valid. However <tt>InEdgeIt</tt>s are invalidated.
    void changeBNode(UEdge e, Node n) { 
      Parent::changeBNode(e,n); 
    }

    /// \brief Changes the source(ANode) of \c e to \c n
    ///
    /// Changes the source(ANode) of \c e to \c n
    ///
    ///\note The <tt>EdgeIt</tt>s and <tt>InEdgeIt</tt>s referencing
    ///the changed edge remain valid. However <tt>OutEdgeIt</tt>s are
    ///invalidated.
    void changeSource(UEdge e, Node n) { 
      Parent::changeANode(e,n); 
    }

    /// \brief Changes the target(BNode) of \c e to \c n
    ///
    /// Changes the target(BNode) of \c e to \c n
    ///
    /// \note The <tt>EdgeIt</tt>s and <tt>OutEdgeIt</tt>s
    /// referencing the changed edge remain
    /// valid. However <tt>InEdgeIt</tt>s are invalidated.
    void changeTarget(UEdge e, Node n) { 
      Parent::changeBNode(e,n); 
    }

    /// \brief Changes the source of \c e to \c n
    ///
    /// Changes the source of \c e to \c n. It changes the proper
    /// node of the represented undirected edge.
    ///
    ///\note The <tt>EdgeIt</tt>s and <tt>InEdgeIt</tt>s
    ///referencing the changed edge remain
    ///valid. However <tt>OutEdgeIt</tt>s are invalidated.
    void changeSource(Edge e, Node n) { 
      if (Parent::direction(e)) {
        Parent::changeANode(e,n);
      } else {
        Parent::changeBNode(e,n);
      } 
    }
    /// \brief Changes the target of \c e to \c n
    ///
    /// Changes the target of \c e to \c n. It changes the proper
    /// node of the represented undirected edge.
    ///
    ///\note The <tt>EdgeIt</tt>s and <tt>OutEdgeIt</tt>s
    ///referencing the changed edge remain
    ///valid. However <tt>InEdgeIt</tt>s are invalidated.
    void changeTarget(Edge e, Node n) { 
      if (Parent::direction(e)) {
        Parent::changeBNode(e,n);
      } else {
        Parent::changeANode(e,n);
      } 
    }
    /// \brief Contract two nodes.
    ///
    /// This function contracts two nodes.
    ///
    /// Node \p b will be removed but instead of deleting its
    /// neighboring edges, they will be joined to \p a.  The two nodes
    /// should be from the same nodeset, of course.
    ///
    /// \note The <tt>EdgeIt</tt>s referencing a moved edge remain
    /// valid.
    void contract(const Node& a, const Node& b) {
      LEMON_ASSERT(Parent::aNode(a) == Parent::aNode(b), NodeSetError());
      if (Parent::aNode(a)) {
        for (IncEdgeIt e(*this, b); e!=INVALID;) {
          IncEdgeIt f = e; ++f;
          changeSource(e, a);
          e = f;
        }
      } else {
        for (IncEdgeIt e(*this, b); e!=INVALID;) {
          IncEdgeIt f = e; ++f;
          changeTarget(e, a);
          e = f;
        }
      }
      erase(b);
    }

    /// \brief Class to make a snapshot of the graph and restore
    /// to it later.
    ///
    /// Class to make a snapshot of the graph and to restore it
    /// later.
    ///
    /// The newly added nodes and undirected edges can be removed
    /// using the restore() function.
    ///
    /// \warning Edge and node deletions cannot be restored. This
    /// events invalidate the snapshot. 
    class Snapshot {
    protected:

      typedef Parent::NodeNotifier NodeNotifier;

      class NodeObserverProxy : public NodeNotifier::ObserverBase {
      public:

        NodeObserverProxy(Snapshot& _snapshot)
          : snapshot(_snapshot) {}

        using NodeNotifier::ObserverBase::attach;
        using NodeNotifier::ObserverBase::detach;
        using NodeNotifier::ObserverBase::attached;
        
      protected:
        
        virtual void add(const Node& node) {
          snapshot.addNode(node);
        }
        virtual void add(const std::vector<Node>& nodes) {
          for (int i = nodes.size() - 1; i >= 0; ++i) {
            snapshot.addNode(nodes[i]);
          }
        }
        virtual void erase(const Node& node) {
          snapshot.eraseNode(node);
        }
        virtual void erase(const std::vector<Node>& nodes) {
          for (int i = 0; i < int(nodes.size()); ++i) {
            snapshot.eraseNode(nodes[i]);
          }
        }
        virtual void build() {
          Node node;
          std::vector<Node> nodes;
          for (notifier()->first(node); node != INVALID; 
               notifier()->next(node)) {
            nodes.push_back(node);
          }
          for (int i = nodes.size() - 1; i >= 0; --i) {
            snapshot.addNode(nodes[i]);
          }
        }
        virtual void clear() {
          Node node;
          for (notifier()->first(node); node != INVALID; 
               notifier()->next(node)) {
            snapshot.eraseNode(node);
          }
        }

        Snapshot& snapshot;
      };

      class UEdgeObserverProxy : public UEdgeNotifier::ObserverBase {
      public:

        UEdgeObserverProxy(Snapshot& _snapshot)
          : snapshot(_snapshot) {}

        using UEdgeNotifier::ObserverBase::attach;
        using UEdgeNotifier::ObserverBase::detach;
        using UEdgeNotifier::ObserverBase::attached;
        
      protected:

        virtual void add(const UEdge& edge) {
          snapshot.addUEdge(edge);
        }
        virtual void add(const std::vector<UEdge>& edges) {
          for (int i = edges.size() - 1; i >= 0; ++i) {
            snapshot.addUEdge(edges[i]);
          }
        }
        virtual void erase(const UEdge& edge) {
          snapshot.eraseUEdge(edge);
        }
        virtual void erase(const std::vector<UEdge>& edges) {
          for (int i = 0; i < int(edges.size()); ++i) {
            snapshot.eraseUEdge(edges[i]);
          }
        }
        virtual void build() {
          UEdge edge;
          std::vector<UEdge> edges;
          for (notifier()->first(edge); edge != INVALID; 
               notifier()->next(edge)) {
            edges.push_back(edge);
          }
          for (int i = edges.size() - 1; i >= 0; --i) {
            snapshot.addUEdge(edges[i]);
          }
        }
        virtual void clear() {
          UEdge edge;
          for (notifier()->first(edge); edge != INVALID; 
               notifier()->next(edge)) {
            snapshot.eraseUEdge(edge);
          }
        }

        Snapshot& snapshot;
      };
      
      ListBpUGraph *graph;

      NodeObserverProxy node_observer_proxy;
      UEdgeObserverProxy edge_observer_proxy;

      std::list<Node> added_nodes;
      std::list<UEdge> added_edges;


      void addNode(const Node& node) {
        added_nodes.push_front(node);        
      }
      void eraseNode(const Node& node) {
        std::list<Node>::iterator it = 
          std::find(added_nodes.begin(), added_nodes.end(), node);
        if (it == added_nodes.end()) {
          clear();
          edge_observer_proxy.detach();
          throw NodeNotifier::ImmediateDetach();
        } else {
          added_nodes.erase(it);
        }
      }

      void addUEdge(const UEdge& edge) {
        added_edges.push_front(edge);        
      }
      void eraseUEdge(const UEdge& edge) {
        std::list<UEdge>::iterator it = 
          std::find(added_edges.begin(), added_edges.end(), edge);
        if (it == added_edges.end()) {
          clear();
          node_observer_proxy.detach();
          throw UEdgeNotifier::ImmediateDetach();
        } else {
          added_edges.erase(it);
        }        
      }

      void attach(ListBpUGraph &_graph) {
	graph = &_graph;
	node_observer_proxy.attach(graph->notifier(Node()));
        edge_observer_proxy.attach(graph->notifier(UEdge()));
      }
            
      void detach() {
	node_observer_proxy.detach();
	edge_observer_proxy.detach();
      }

      bool attached() const {
        return node_observer_proxy.attached();
      }

      void clear() {
        added_nodes.clear();
        added_edges.clear();        
      }

    public:

      /// \brief Default constructor.
      ///
      /// Default constructor.
      /// To actually make a snapshot you must call save().
      Snapshot() 
        : graph(0), node_observer_proxy(*this), 
          edge_observer_proxy(*this) {}
      
      /// \brief Constructor that immediately makes a snapshot.
      ///      
      /// This constructor immediately makes a snapshot of the graph.
      /// \param _graph The graph we make a snapshot of.
      Snapshot(ListBpUGraph &_graph) 
        : node_observer_proxy(*this), 
          edge_observer_proxy(*this) {
	attach(_graph);
      }
      
      /// \brief Make a snapshot.
      ///
      /// Make a snapshot of the graph.
      ///
      /// This function can be called more than once. In case of a repeated
      /// call, the previous snapshot gets lost.
      /// \param _graph The graph we make the snapshot of.
      void save(ListBpUGraph &_graph) {
        if (attached()) {
          detach();
          clear();
        }
        attach(_graph);
      }
      
      /// \brief Undo the changes until the last snapshot.
      // 
      /// Undo the changes until the last snapshot created by save().
      void restore() {
	detach();
	for(std::list<UEdge>::iterator it = added_edges.begin(); 
            it != added_edges.end(); ++it) {
	  graph->erase(*it);
	}
	for(std::list<Node>::iterator it = added_nodes.begin(); 
            it != added_nodes.end(); ++it) {
	  graph->erase(*it);
	}
        clear();
      }

      /// \brief Gives back true when the snapshot is valid.
      ///
      /// Gives back true when the snapshot is valid.
      bool valid() const {
        return attached();
      }
    };
  };

  
  /// @}  
} //namespace lemon
  

#endif

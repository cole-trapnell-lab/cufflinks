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

#ifndef LEMON_FIB_HEAP_H
#define LEMON_FIB_HEAP_H

///\file
///\ingroup auxdat
///\brief Fibonacci Heap implementation.

#include <vector>
#include <functional>
#include <lemon/math.h>

namespace lemon {
  
  /// \ingroup auxdat
  ///
  ///\brief Fibonacci Heap.
  ///
  ///This class implements the \e Fibonacci \e heap data structure. A \e heap
  ///is a data structure for storing items with specified values called \e
  ///priorities in such a way that finding the item with minimum priority is
  ///efficient. \c Compare specifies the ordering of the priorities. In a heap
  ///one can change the priority of an item, add or erase an item, etc.
  ///
  ///The methods \ref increase and \ref erase are not efficient in a Fibonacci
  ///heap. In case of many calls to these operations, it is better to use a
  ///\ref BinHeap "binary heap".
  ///
  ///\param _Prio Type of the priority of the items.
  ///\param _ItemIntMap A read and writable Item int map, used internally
  ///to handle the cross references.
  ///\param _Compare A class for the ordering of the priorities. The
  ///default is \c std::less<_Prio>.
  ///
  ///\sa BinHeap
  ///\sa Dijkstra
  ///\author Jacint Szabo 
 
#ifdef DOXYGEN
  template <typename _Prio, 
	    typename _ItemIntMap, 
	    typename _Compare>
#else
  template <typename _Prio, 
	    typename _ItemIntMap, 
	    typename _Compare = std::less<_Prio> >
#endif
  class FibHeap {
  public:
    typedef _ItemIntMap ItemIntMap;
    typedef _Prio Prio;
    typedef typename ItemIntMap::Key Item;
    typedef std::pair<Item,Prio> Pair;
    typedef _Compare Compare;
    
  private:
    class store;
    
    std::vector<store> container;
    int minimum;
    ItemIntMap &iimap;
    Compare comp;
    int num_items;
    
  public:
    ///Status of the nodes
    enum State {
      ///The node is in the heap
      IN_HEAP = 0,
      ///The node has never been in the heap
      PRE_HEAP = -1,
      ///The node was in the heap but it got out of it
      POST_HEAP = -2
    };
    
    /// \brief The constructor
    ///
    /// \c _iimap should be given to the constructor, since it is
    ///   used internally to handle the cross references.
    explicit FibHeap(ItemIntMap &_iimap) 
      : minimum(0), iimap(_iimap), num_items() {} 
 
    /// \brief The constructor
    ///
    /// \c _iimap should be given to the constructor, since it is used
    /// internally to handle the cross references. \c _comp is an
    /// object for ordering of the priorities. 
    FibHeap(ItemIntMap &_iimap, const Compare &_comp) 
      : minimum(0), iimap(_iimap), comp(_comp), num_items() {}
    
    /// \brief The number of items stored in the heap.
    ///
    /// Returns the number of items stored in the heap.
    int size() const { return num_items; }

    /// \brief Checks if the heap stores no items.
    ///
    ///   Returns \c true if and only if the heap stores no items.
    bool empty() const { return num_items==0; }

    /// \brief Make empty this heap.
    /// 
    /// Make empty this heap. It does not change the cross reference
    /// map.  If you want to reuse a heap what is not surely empty you
    /// should first clear the heap and after that you should set the
    /// cross reference map for each item to \c PRE_HEAP.
    void clear() {
      container.clear(); minimum = 0; num_items = 0;
    }

    /// \brief \c item gets to the heap with priority \c value independently 
    /// if \c item was already there.
    ///
    /// This method calls \ref push(\c item, \c value) if \c item is not
    /// stored in the heap and it calls \ref decrease(\c item, \c value) or
    /// \ref increase(\c item, \c value) otherwise.
    void set (const Item& item, const Prio& value) {
      int i=iimap[item];
      if ( i >= 0 && container[i].in ) {
	if ( comp(value, container[i].prio) ) decrease(item, value); 
	if ( comp(container[i].prio, value) ) increase(item, value); 
      } else push(item, value);
    }
    
    /// \brief Adds \c item to the heap with priority \c value. 
    ///    
    /// Adds \c item to the heap with priority \c value. 
    /// \pre \c item must not be stored in the heap. 
    void push (const Item& item, const Prio& value) {
      int i=iimap[item];      
      if ( i < 0 ) {
	int s=container.size();
	iimap.set( item, s );	
	store st;
	st.name=item;
	container.push_back(st);
	i=s;
      } else {
	container[i].parent=container[i].child=-1;
	container[i].degree=0;
	container[i].in=true;
	container[i].marked=false;
      }

      if ( num_items ) {
	container[container[minimum].right_neighbor].left_neighbor=i;
	container[i].right_neighbor=container[minimum].right_neighbor;
	container[minimum].right_neighbor=i;
	container[i].left_neighbor=minimum;
	if ( comp( value, container[minimum].prio) ) minimum=i; 
      } else {
	container[i].right_neighbor=container[i].left_neighbor=i;
	minimum=i;	
      }
      container[i].prio=value;
      ++num_items;
    }
    
    /// \brief Returns the item with minimum priority relative to \c Compare.
    ///
    /// This method returns the item with minimum priority relative to \c
    /// Compare.  
    /// \pre The heap must be nonempty.  
    Item top() const { return container[minimum].name; }

    /// \brief Returns the minimum priority relative to \c Compare.
    ///
    /// It returns the minimum priority relative to \c Compare.
    /// \pre The heap must be nonempty.
    const Prio& prio() const { return container[minimum].prio; }
        
    /// \brief Returns the priority of \c item.
    ///
    /// It returns the priority of \c item.
    /// \pre \c item must be in the heap.
    const Prio& operator[](const Item& item) const { 
      return container[iimap[item]].prio; 
    }

    /// \brief Deletes the item with minimum priority relative to \c Compare.
    ///
    /// This method deletes the item with minimum priority relative to \c
    /// Compare from the heap.  
    /// \pre The heap must be non-empty.  
    void pop() {
      /*The first case is that there are only one root.*/
      if ( container[minimum].left_neighbor==minimum ) {
	container[minimum].in=false;
	if ( container[minimum].degree!=0 ) { 
	  makeroot(container[minimum].child);
	  minimum=container[minimum].child;
	  balance();
	}
      } else {
	int right=container[minimum].right_neighbor;
	unlace(minimum);
	container[minimum].in=false;
	if ( container[minimum].degree > 0 ) {
	  int left=container[minimum].left_neighbor;
	  int child=container[minimum].child;
	  int last_child=container[child].left_neighbor;
	  
	  makeroot(child);
	  
	  container[left].right_neighbor=child;
	  container[child].left_neighbor=left;
	  container[right].left_neighbor=last_child;
	  container[last_child].right_neighbor=right;
	}
	minimum=right;
	balance();
      } // the case where there are more roots
      --num_items;   
    }

    /// \brief Deletes \c item from the heap.
    ///
    /// This method deletes \c item from the heap, if \c item was already
    /// stored in the heap. It is quite inefficient in Fibonacci heaps.
    void erase (const Item& item) {
      int i=iimap[item];
      
      if ( i >= 0 && container[i].in ) { 	
	if ( container[i].parent!=-1 ) {
	  int p=container[i].parent;
	  cut(i,p);	    
	  cascade(p);
	}
	minimum=i;     //As if its prio would be -infinity
	pop();
      }
    }

    /// \brief Decreases the priority of \c item to \c value.
    ///
    /// This method decreases the priority of \c item to \c value.
    /// \pre \c item must be stored in the heap with priority at least \c
    ///   value relative to \c Compare.
    void decrease (Item item, const Prio& value) {
      int i=iimap[item];
      container[i].prio=value;
      int p=container[i].parent;
      
      if ( p!=-1 && comp(value, container[p].prio) ) {
	cut(i,p);	    
	cascade(p);
      }      
      if ( comp(value, container[minimum].prio) ) minimum=i; 
    }

    /// \brief Increases the priority of \c item to \c value.
    ///
    /// This method sets the priority of \c item to \c value. Though
    /// there is no precondition on the priority of \c item, this
    /// method should be used only if it is indeed necessary to increase
    /// (relative to \c Compare) the priority of \c item, because this
    /// method is inefficient.
    void increase (Item item, const Prio& value) {
      erase(item);
      push(item, value);
    }


    /// \brief Returns if \c item is in, has already been in, or has never 
    /// been in the heap.
    ///
    /// This method returns PRE_HEAP if \c item has never been in the
    /// heap, IN_HEAP if it is in the heap at the moment, and POST_HEAP
    /// otherwise. In the latter case it is possible that \c item will
    /// get back to the heap again.
    State state(const Item &item) const {
      int i=iimap[item];
      if( i>=0 ) {
	if ( container[i].in ) i=0;
	else i=-2; 
      }
      return State(i);
    }    

    /// \brief Sets the state of the \c item in the heap.
    ///
    /// Sets the state of the \c item in the heap. It can be used to
    /// manually clear the heap when it is important to achive the
    /// better time complexity.
    /// \param i The item.
    /// \param st The state. It should not be \c IN_HEAP. 
    void state(const Item& i, State st) {
      switch (st) {
      case POST_HEAP:
      case PRE_HEAP:
        if (state(i) == IN_HEAP) {
          erase(i);
        }
        iimap[i] = st;
        break;
      case IN_HEAP:
        break;
      }
    }
    
  private:
    
    void balance() {

      int maxdeg=int( std::floor( 2.08*log(double(container.size()))))+1;
  
      std::vector<int> A(maxdeg,-1); 
    
      /*
       *Recall that now minimum does not point to the minimum prio element.
       *We set minimum to this during balance().
       */
      int anchor=container[minimum].left_neighbor; 
      int next=minimum; 
      bool end=false; 
    	
      do {
	int active=next;
	if ( anchor==active ) end=true;
	int d=container[active].degree;
	next=container[active].right_neighbor;

	while (A[d]!=-1) {	  
	  if( comp(container[active].prio, container[A[d]].prio) ) {
	    fuse(active,A[d]); 
	  } else { 
	    fuse(A[d],active);
	    active=A[d];
	  } 
	  A[d]=-1;
	  ++d;
	}	
	A[d]=active;
      } while ( !end );


      while ( container[minimum].parent >=0 ) 
	minimum=container[minimum].parent;
      int s=minimum;
      int m=minimum;
      do {  
	if ( comp(container[s].prio, container[minimum].prio) ) minimum=s;
	s=container[s].right_neighbor;
      } while ( s != m );
    }

    void makeroot(int c) {
      int s=c;
      do {  
	container[s].parent=-1;
	s=container[s].right_neighbor;
      } while ( s != c );
    }

    void cut(int a, int b) {
      /*
       *Replacing a from the children of b.
       */
      --container[b].degree;
    
      if ( container[b].degree !=0 ) {
	int child=container[b].child;
	if ( child==a ) 
	  container[b].child=container[child].right_neighbor;
	unlace(a);
      }
    
    
      /*Lacing a to the roots.*/
      int right=container[minimum].right_neighbor;
      container[minimum].right_neighbor=a;
      container[a].left_neighbor=minimum;
      container[a].right_neighbor=right;
      container[right].left_neighbor=a;
    
      container[a].parent=-1;
      container[a].marked=false;
    }

    void cascade(int a) {
      if ( container[a].parent!=-1 ) {
	int p=container[a].parent;
	
	if ( container[a].marked==false ) container[a].marked=true;
	else {
	  cut(a,p);
	  cascade(p);
	}
      }
    }

    void fuse(int a, int b) {
      unlace(b);
      
      /*Lacing b under a.*/
      container[b].parent=a;

      if (container[a].degree==0) {
	container[b].left_neighbor=b;
	container[b].right_neighbor=b;
	container[a].child=b;	
      } else {
	int child=container[a].child;
	int last_child=container[child].left_neighbor;
	container[child].left_neighbor=b;
	container[b].right_neighbor=child;
	container[last_child].right_neighbor=b;
	container[b].left_neighbor=last_child;
      }

      ++container[a].degree;
      
      container[b].marked=false;
    }

    /*
     *It is invoked only if a has siblings.
     */
    void unlace(int a) {
      int leftn=container[a].left_neighbor;
      int rightn=container[a].right_neighbor;
      container[leftn].right_neighbor=rightn;
      container[rightn].left_neighbor=leftn;
    }


    class store {
      friend class FibHeap;
      
      Item name;
      int parent;
      int left_neighbor;
      int right_neighbor;
      int child;
      int degree;  
      bool marked;
      bool in;
      Prio prio;
      
      store() : parent(-1), child(-1), degree(), marked(false), in(true) {} 
    };
  };    

} //namespace lemon

#endif //LEMON_FIB_HEAP_H


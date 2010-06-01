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

///\ingroup concept
///\file
///\brief Classes for representing heaps.
///

#ifndef LEMON_CONCEPT_HEAP_H
#define LEMON_CONCEPT_HEAP_H

#include <lemon/bits/invalid.h>

namespace lemon {
  namespace concepts {
    /// \addtogroup concept
    /// @{


    /// \brief A concept structure describes the main interface of heaps.
    ///
    /// A concept structure describes the main interface of heaps.
    ///
    template <typename Prio, typename ItemIntMap>
    class Heap {
    public:

      ///\brief Type of the items stored in the heap.
      typedef typename ItemIntMap::Key  Item;
  

      /// \brief Type to represent the items states.
      ///
      /// Each Item element have a state associated to it. It may be "in heap",
      /// "pre heap" or "post heap". The later two are indifferent from the
      /// heap's point of view, but may be useful to the user.
      ///
      /// The ItemIntMap _should_ be initialized in such way, that it maps
      /// PRE_HEAP (-1) to any element to be put in the heap...
      enum State {
	IN_HEAP = 0,
	PRE_HEAP = -1,
	POST_HEAP = -2
      };
      
      /// \brief The constructor.
      ///
      /// The constructor.
      /// \param _iim should be given to the constructor, since it is used
      /// internally to handle the cross references. The value of the map
      /// should be PRE_HEAP (-1) for each element.
      explicit Heap(ItemIntMap &_iim) {}

      /// \brief The number of items stored in the heap.
      ///
      /// Returns the number of items stored in the heap.
      int size() const { return 0; }

      /// \brief Checks if the heap stores no items.
      ///
      /// Returns \c true if and only if the heap stores no items.
      bool empty() const { return false; }

      /// \brief Makes empty this heap.
      ///
      /// Makes this heap empty.
      void clear();

      /// \brief Insert an item into the heap with the given heap.
      ///    
      /// Adds \c i to the heap with priority \c p. 
      /// \param i The item to insert.
      /// \param p The priority of the item.
      void push(const Item &i, const Prio &p) {}

      /// \brief Returns the item with minimum priority.
      ///
      /// This method returns the item with minimum priority.  
      /// \pre The heap must be nonempty.  
      Item top() const {}

      /// \brief Returns the minimum priority.
      ///
      /// It returns the minimum priority.
      /// \pre The heap must be nonempty.
      Prio prio() const {}

      /// \brief Deletes the item with minimum priority.
      ///
      /// This method deletes the item with minimum priority.
      /// \pre The heap must be non-empty.  
      void pop() {}

      /// \brief Deletes \c i from the heap.
      ///
      /// This method deletes item \c i from the heap, if \c i was
      /// already stored in the heap.
      /// \param i The item to erase. 
      void erase(const Item &i) {}

      /// \brief Returns the priority of \c i.
      ///
      /// This function returns the priority of item \c i.  
      /// \pre \c i must be in the heap.
      /// \param i The item.
      Prio operator[](const Item &i) const {}

      /// \brief \c i gets to the heap with priority \c p independently 
      /// if \c i was already there.
      ///
      /// This method calls \ref push(\c i, \c p) if \c i is not stored
      /// in the heap and sets the priority of \c i to \c p otherwise.
      /// It may throw an \e UnderFlowPriorityException. 
      /// \param i The item.
      /// \param p The priority.
      void set(const Item &i, const Prio &p) {}
      
      /// \brief Decreases the priority of \c i to \c p.
      ///
      /// This method decreases the priority of item \c i to \c p.
      /// \pre \c i must be stored in the heap with priority at least \c p.
      /// \param i The item.
      /// \param p The priority.
      void decrease(const Item &i, const Prio &p) {}

      /// \brief Increases the priority of \c i to \c p.
      ///
      /// This method sets the priority of item \c i to \c p. 
      /// \pre \c i must be stored in the heap with priority at most \c
      /// p relative to \c Compare.
      /// \param i The item.
      /// \param p The priority.
      void increase(const Item &i, const Prio &p) {}

      /// \brief Returns if \c item is in, has already been in, or has 
      /// never been in the heap.
      ///
      /// This method returns PRE_HEAP if \c item has never been in the
      /// heap, IN_HEAP if it is in the heap at the moment, and POST_HEAP
      /// otherwise. In the latter case it is possible that \c item will
      /// get back to the heap again.
      /// \param i The item.
      State state(const Item &i) const {}

      /// \brief Sets the state of the \c item in the heap.
      ///
      /// Sets the state of the \c item in the heap. It can be used to
      /// manually clear the heap when it is important to achive the
      /// better time complexity.
      /// \param i The item.
      /// \param st The state. It should not be \c IN_HEAP. 
      void state(const Item& i, State st) {}


      template <typename _Heap>
      struct Constraints {
      public:
    
	void constraints() {
	  Item item;
	  Prio prio;

	  item=Item();
	  prio=Prio();

	  ignore_unused_variable_warning(item);
	  ignore_unused_variable_warning(prio);

	  typedef typename _Heap::State State;
	  State state;

	  ignore_unused_variable_warning(state);
      
	  _Heap heap1 = _Heap(map);

	  ignore_unused_variable_warning(heap1);
      
	  heap.push(item, prio);

	  prio = heap.prio();
	  item = heap.top();

	  heap.pop();

	  heap.set(item, prio);
	  heap.decrease(item, prio);
	  heap.increase(item, prio);
	  prio = heap[item];

	  heap.erase(item);

	  state = heap.state(item);

	  state = _Heap::PRE_HEAP;
	  state = _Heap::IN_HEAP;
	  state = _Heap::POST_HEAP;

	  heap.clear();
	}
    
	_Heap& heap;
	ItemIntMap& map;

	Constraints() : heap(0), map(0) {}
      };
    };

    /// @}
  } // namespace lemon
}
#endif // LEMON_CONCEPT_PATH_H

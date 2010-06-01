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

#ifndef LEMON_BUCKET_HEAP_H
#define LEMON_BUCKET_HEAP_H

///\ingroup auxdat
///\file
///\brief Bucket Heap implementation.

#include <vector>
#include <utility>
#include <functional>

namespace lemon {

  /// \ingroup auxdat
  ///
  /// \brief A Bucket Heap implementation.
  ///
  /// This class implements the \e bucket \e heap data structure. A \e heap
  /// is a data structure for storing items with specified values called \e
  /// priorities in such a way that finding the item with minimum priority is
  /// efficient. The bucket heap is very simple implementation, it can store
  /// only integer priorities and it stores for each priority in the 
  /// \f$ [0..C) \f$ range a list of items. So it should be used only when 
  /// the priorities are small. It is not intended to use as dijkstra heap.
  ///
  /// \param _ItemIntMap A read and writable Item int map, used internally
  /// to handle the cross references.
  /// \param minimize If the given parameter is true then the heap gives back
  /// the lowest priority. 
  template <typename _ItemIntMap, bool minimize = true >
  class BucketHeap {

  public:
    /// \e
    typedef typename _ItemIntMap::Key Item;
    /// \e
    typedef int Prio;
    /// \e
    typedef std::pair<Item, Prio> Pair;
    /// \e
    typedef _ItemIntMap ItemIntMap;

    /// \brief Type to represent the items states.
    ///
    /// Each Item element have a state associated to it. It may be "in heap",
    /// "pre heap" or "post heap". The latter two are indifferent from the
    /// heap's point of view, but may be useful to the user.
    ///
    /// The ItemIntMap \e should be initialized in such way that it maps
    /// PRE_HEAP (-1) to any element to be put in the heap...
    enum State {
      IN_HEAP = 0,
      PRE_HEAP = -1,
      POST_HEAP = -2
    };

  public:
    /// \brief The constructor.
    ///
    /// The constructor.
    /// \param _index should be given to the constructor, since it is used
    /// internally to handle the cross references. The value of the map
    /// should be PRE_HEAP (-1) for each element.
    explicit BucketHeap(ItemIntMap &_index) : index(_index), minimal(0) {}
    
    /// The number of items stored in the heap.
    ///
    /// \brief Returns the number of items stored in the heap.
    int size() const { return data.size(); }
    
    /// \brief Checks if the heap stores no items.
    ///
    /// Returns \c true if and only if the heap stores no items.
    bool empty() const { return data.empty(); }

    /// \brief Make empty this heap.
    /// 
    /// Make empty this heap. It does not change the cross reference
    /// map.  If you want to reuse a heap what is not surely empty you
    /// should first clear the heap and after that you should set the
    /// cross reference map for each item to \c PRE_HEAP.
    void clear() { 
      data.clear(); first.clear(); minimal = 0;
    }

  private:

    void relocate_last(int idx) {
      if (idx + 1 < int(data.size())) {
	data[idx] = data.back();
	if (data[idx].prev != -1) {
	  data[data[idx].prev].next = idx;
	} else {
	  first[data[idx].value] = idx;
	}
	if (data[idx].next != -1) {
	  data[data[idx].next].prev = idx;
	}
	index[data[idx].item] = idx;
      }
      data.pop_back();
    }

    void unlace(int idx) {
      if (data[idx].prev != -1) {
	data[data[idx].prev].next = data[idx].next;
      } else {
	first[data[idx].value] = data[idx].next;
      }
      if (data[idx].next != -1) {
	data[data[idx].next].prev = data[idx].prev;
      }
    }

    void lace(int idx) {
      if (int(first.size()) <= data[idx].value) {
	first.resize(data[idx].value + 1, -1);
      }
      data[idx].next = first[data[idx].value];
      if (data[idx].next != -1) {
	data[data[idx].next].prev = idx;
      }
      first[data[idx].value] = idx;
      data[idx].prev = -1;
    }

  public:
    /// \brief Insert a pair of item and priority into the heap.
    ///
    /// Adds \c p.first to the heap with priority \c p.second.
    /// \param p The pair to insert.
    void push(const Pair& p) {
      push(p.first, p.second);
    }

    /// \brief Insert an item into the heap with the given priority.
    ///    
    /// Adds \c i to the heap with priority \c p. 
    /// \param i The item to insert.
    /// \param p The priority of the item.
    void push(const Item &i, const Prio &p) { 
      int idx = data.size();
      index[i] = idx;
      data.push_back(BucketItem(i, p));
      lace(idx);
      if (p < minimal) {
	minimal = p;
      }
    }

    /// \brief Returns the item with minimum priority.
    ///
    /// This method returns the item with minimum priority.
    /// \pre The heap must be nonempty.  
    Item top() const {
      while (first[minimal] == -1) {
	++minimal;
      }
      return data[first[minimal]].item;
    }

    /// \brief Returns the minimum priority.
    ///
    /// It returns the minimum priority.
    /// \pre The heap must be nonempty.
    Prio prio() const {
      while (first[minimal] == -1) {
	++minimal;
      }
      return minimal;
    }

    /// \brief Deletes the item with minimum priority.
    ///
    /// This method deletes the item with minimum priority from the heap.  
    /// \pre The heap must be non-empty.  
    void pop() {
      while (first[minimal] == -1) {
	++minimal;
      }
      int idx = first[minimal];
      index[data[idx].item] = -2;
      unlace(idx);
      relocate_last(idx);
    }

    /// \brief Deletes \c i from the heap.
    ///
    /// This method deletes item \c i from the heap, if \c i was
    /// already stored in the heap.
    /// \param i The item to erase. 
    void erase(const Item &i) {
      int idx = index[i];
      index[data[idx].item] = -2;
      unlace(idx);
      relocate_last(idx);
    }

    
    /// \brief Returns the priority of \c i.
    ///
    /// This function returns the priority of item \c i.  
    /// \pre \c i must be in the heap.
    /// \param i The item.
    Prio operator[](const Item &i) const {
      int idx = index[i];
      return data[idx].value;
    }

    /// \brief \c i gets to the heap with priority \c p independently 
    /// if \c i was already there.
    ///
    /// This method calls \ref push(\c i, \c p) if \c i is not stored
    /// in the heap and sets the priority of \c i to \c p otherwise.
    /// \param i The item.
    /// \param p The priority.
    void set(const Item &i, const Prio &p) {
      int idx = index[i];
      if (idx < 0) {
	push(i,p);
      } else if (p > data[idx].value) {
	increase(i, p);
      } else {
	decrease(i, p);
      }
    }

    /// \brief Decreases the priority of \c i to \c p.
    ///
    /// This method decreases the priority of item \c i to \c p.
    /// \pre \c i must be stored in the heap with priority at least \c
    /// p relative to \c Compare.
    /// \param i The item.
    /// \param p The priority.
    void decrease(const Item &i, const Prio &p) {
      int idx = index[i];
      unlace(idx);
      data[idx].value = p;
      if (p < minimal) {
	minimal = p;
      }
      lace(idx);
    }
    
    /// \brief Increases the priority of \c i to \c p.
    ///
    /// This method sets the priority of item \c i to \c p. 
    /// \pre \c i must be stored in the heap with priority at most \c
    /// p relative to \c Compare.
    /// \param i The item.
    /// \param p The priority.
    void increase(const Item &i, const Prio &p) {
      int idx = index[i];
      unlace(idx);
      data[idx].value = p;
      lace(idx);
    }

    /// \brief Returns if \c item is in, has already been in, or has 
    /// never been in the heap.
    ///
    /// This method returns PRE_HEAP if \c item has never been in the
    /// heap, IN_HEAP if it is in the heap at the moment, and POST_HEAP
    /// otherwise. In the latter case it is possible that \c item will
    /// get back to the heap again.
    /// \param i The item.
    State state(const Item &i) const {
      int idx = index[i];
      if (idx >= 0) idx = 0;
      return State(idx);
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
        index[i] = st;
        break;
      case IN_HEAP:
        break;
      }
    }

  private:

    struct BucketItem {
      BucketItem(const Item& _item, int _value) 
	: item(_item), value(_value) {}

      Item item;
      int value;

      int prev, next;
    };

    ItemIntMap& index;
    std::vector<int> first;
    std::vector<BucketItem> data;
    mutable int minimal;

  }; // class BucketHeap


  template <typename _ItemIntMap>
  class BucketHeap<_ItemIntMap, false> {

  public:
    typedef typename _ItemIntMap::Key Item;
    typedef int Prio;
    typedef std::pair<Item, Prio> Pair;
    typedef _ItemIntMap ItemIntMap;

    enum State {
      IN_HEAP = 0,
      PRE_HEAP = -1,
      POST_HEAP = -2
    };

  public:

    explicit BucketHeap(ItemIntMap &_index) : index(_index), maximal(-1) {}

    int size() const { return data.size(); }
    bool empty() const { return data.empty(); }

    void clear() { 
      data.clear(); first.clear(); maximal = -1; 
    }

  private:

    void relocate_last(int idx) {
      if (idx + 1 != int(data.size())) {
	data[idx] = data.back();
	if (data[idx].prev != -1) {
	  data[data[idx].prev].next = idx;
	} else {
	  first[data[idx].value] = idx;
	}
	if (data[idx].next != -1) {
	  data[data[idx].next].prev = idx;
	}
	index[data[idx].item] = idx;
      }
      data.pop_back();
    }

    void unlace(int idx) {
      if (data[idx].prev != -1) {
	data[data[idx].prev].next = data[idx].next;
      } else {
	first[data[idx].value] = data[idx].next;
      }
      if (data[idx].next != -1) {
	data[data[idx].next].prev = data[idx].prev;
      }
    }

    void lace(int idx) {
      if (int(first.size()) <= data[idx].value) {
	first.resize(data[idx].value + 1, -1);
      }
      data[idx].next = first[data[idx].value];
      if (data[idx].next != -1) {
	data[data[idx].next].prev = idx;
      }
      first[data[idx].value] = idx;
      data[idx].prev = -1;
    }

  public:

    void push(const Pair& p) {
      push(p.first, p.second);
    }

    void push(const Item &i, const Prio &p) { 
      int idx = data.size();
      index[i] = idx;
      data.push_back(BucketItem(i, p));
      lace(idx);
      if (data[idx].value > maximal) {
	maximal = data[idx].value;
      }
    }

    Item top() const {
      while (first[maximal] == -1) {
	--maximal;
      }
      return data[first[maximal]].item;
    }

    Prio prio() const {
      while (first[maximal] == -1) {
	--maximal;
      }
      return maximal;
    }

    void pop() {
      while (first[maximal] == -1) {
	--maximal;
      }
      int idx = first[maximal];
      index[data[idx].item] = -2;
      unlace(idx);
      relocate_last(idx);
    }

    void erase(const Item &i) {
      int idx = index[i];
      index[data[idx].item] = -2;
      unlace(idx);
      relocate_last(idx);
    }

    Prio operator[](const Item &i) const {
      int idx = index[i];
      return data[idx].value;
    }

    void set(const Item &i, const Prio &p) {
      int idx = index[i];
      if (idx < 0) {
	push(i,p);
      } else if (p > data[idx].value) {
	decrease(i, p);
      } else {
	increase(i, p);
      }
    }

    void decrease(const Item &i, const Prio &p) {
      int idx = index[i];
      unlace(idx);
      data[idx].value = p;
      if (p > maximal) {
	maximal = p;
      }
      lace(idx);
    }
    
    void increase(const Item &i, const Prio &p) {
      int idx = index[i];
      unlace(idx);
      data[idx].value = p;
      lace(idx);
    }

    State state(const Item &i) const {
      int idx = index[i];
      if (idx >= 0) idx = 0;
      return State(idx);
    }

    void state(const Item& i, State st) {
      switch (st) {
      case POST_HEAP:
      case PRE_HEAP:
        if (state(i) == IN_HEAP) {
          erase(i);
        }
        index[i] = st;
        break;
      case IN_HEAP:
        break;
      }
    }

  private:

    struct BucketItem {
      BucketItem(const Item& _item, int _value) 
	: item(_item), value(_value) {}

      Item item;
      int value;

      int prev, next;
    };

    ItemIntMap& index;
    std::vector<int> first;
    std::vector<BucketItem> data;
    mutable int maximal;

  }; // class BucketHeap

  /// \ingroup auxdat
  ///
  /// \brief A Simplified Bucket Heap implementation.
  ///
  /// This class implements a simplified \e bucket \e heap data
  /// structure.  It does not provide some functionality but it faster
  /// and simplier data structure than the BucketHeap. The main
  /// difference is that the BucketHeap stores for every key a double
  /// linked list while this class stores just simple lists. In the
  /// other way it does not supports erasing each elements just the
  /// minimal and it does not supports key increasing, decreasing.
  ///
  /// \param _ItemIntMap A read and writable Item int map, used internally
  /// to handle the cross references.
  /// \param minimize If the given parameter is true then the heap gives back
  /// the lowest priority.
  ///
  /// \sa BucketHeap 
  template <typename _ItemIntMap, bool minimize = true >
  class SimpleBucketHeap {

  public:
    typedef typename _ItemIntMap::Key Item;
    typedef int Prio;
    typedef std::pair<Item, Prio> Pair;
    typedef _ItemIntMap ItemIntMap;

    /// \brief Type to represent the items states.
    ///
    /// Each Item element have a state associated to it. It may be "in heap",
    /// "pre heap" or "post heap". The latter two are indifferent from the
    /// heap's point of view, but may be useful to the user.
    ///
    /// The ItemIntMap \e should be initialized in such way that it maps
    /// PRE_HEAP (-1) to any element to be put in the heap...
    enum State {
      IN_HEAP = 0,
      PRE_HEAP = -1,
      POST_HEAP = -2
    };

  public:

    /// \brief The constructor.
    ///
    /// The constructor.
    /// \param _index should be given to the constructor, since it is used
    /// internally to handle the cross references. The value of the map
    /// should be PRE_HEAP (-1) for each element.
    explicit SimpleBucketHeap(ItemIntMap &_index) 
      : index(_index), free(-1), num(0), minimal(0) {}
    
    /// \brief Returns the number of items stored in the heap.
    ///
    /// The number of items stored in the heap.
    int size() const { return num; }
    
    /// \brief Checks if the heap stores no items.
    ///
    /// Returns \c true if and only if the heap stores no items.
    bool empty() const { return num == 0; }

    /// \brief Make empty this heap.
    /// 
    /// Make empty this heap. It does not change the cross reference
    /// map.  If you want to reuse a heap what is not surely empty you
    /// should first clear the heap and after that you should set the
    /// cross reference map for each item to \c PRE_HEAP.
    void clear() { 
      data.clear(); first.clear(); free = -1; num = 0; minimal = 0;
    }

    /// \brief Insert a pair of item and priority into the heap.
    ///
    /// Adds \c p.first to the heap with priority \c p.second.
    /// \param p The pair to insert.
    void push(const Pair& p) {
      push(p.first, p.second);
    }

    /// \brief Insert an item into the heap with the given priority.
    ///    
    /// Adds \c i to the heap with priority \c p. 
    /// \param i The item to insert.
    /// \param p The priority of the item.
    void push(const Item &i, const Prio &p) {
      int idx;
      if (free == -1) {
        idx = data.size();
        data.push_back(BucketItem(i));
      } else {
        idx = free;
        free = data[idx].next;
        data[idx].item = i;
      }
      index[i] = idx;
      if (p >= int(first.size())) first.resize(p + 1, -1);
      data[idx].next = first[p];
      first[p] = idx;
      if (p < minimal) {
	minimal = p;
      }
      ++num;
    }

    /// \brief Returns the item with minimum priority.
    ///
    /// This method returns the item with minimum priority.
    /// \pre The heap must be nonempty.  
    Item top() const {
      while (first[minimal] == -1) {
	++minimal;
      }
      return data[first[minimal]].item;
    }

    /// \brief Returns the minimum priority.
    ///
    /// It returns the minimum priority.
    /// \pre The heap must be nonempty.
    Prio prio() const {
      while (first[minimal] == -1) {
	++minimal;
      }
      return minimal;
    }

    /// \brief Deletes the item with minimum priority.
    ///
    /// This method deletes the item with minimum priority from the heap.  
    /// \pre The heap must be non-empty.  
    void pop() {
      while (first[minimal] == -1) {
	++minimal;
      }
      int idx = first[minimal];
      index[data[idx].item] = -2;
      first[minimal] = data[idx].next;
      data[idx].next = free;
      free = idx;
      --num;
    }
    
    /// \brief Returns the priority of \c i.
    ///
    /// This function returns the priority of item \c i.
    /// \warning This operator is not a constant time function
    /// because it scans the whole data structure to find the proper
    /// value.  
    /// \pre \c i must be in the heap.
    /// \param i The item.
    Prio operator[](const Item &i) const {
      for (int k = 0; k < first.size(); ++k) {
        int idx = first[k];
        while (idx != -1) {
          if (data[idx].item == i) {
            return k;
          }
          idx = data[idx].next;
        }
      }
      return -1;
    }

    /// \brief Returns if \c item is in, has already been in, or has 
    /// never been in the heap.
    ///
    /// This method returns PRE_HEAP if \c item has never been in the
    /// heap, IN_HEAP if it is in the heap at the moment, and POST_HEAP
    /// otherwise. In the latter case it is possible that \c item will
    /// get back to the heap again.
    /// \param i The item.
    State state(const Item &i) const {
      int idx = index[i];
      if (idx >= 0) idx = 0;
      return State(idx);
    }

  private:

    struct BucketItem {
      BucketItem(const Item& _item) 
	: item(_item) {}

      Item item;
      int next;
    };

    ItemIntMap& index;
    std::vector<int> first;
    std::vector<BucketItem> data;
    int free, num;
    mutable int minimal;

  }; // class SimpleBucketHeap

  template <typename _ItemIntMap>
  class SimpleBucketHeap<_ItemIntMap, false> {

  public:
    typedef typename _ItemIntMap::Key Item;
    typedef int Prio;
    typedef std::pair<Item, Prio> Pair;
    typedef _ItemIntMap ItemIntMap;

    enum State {
      IN_HEAP = 0,
      PRE_HEAP = -1,
      POST_HEAP = -2
    };

  public:

    explicit SimpleBucketHeap(ItemIntMap &_index) 
      : index(_index), free(-1), num(0), maximal(0) {}
    
    int size() const { return num; }
    
    bool empty() const { return num == 0; }

    void clear() { 
      data.clear(); first.clear(); free = -1; num = 0; maximal = 0;
    }

    void push(const Pair& p) {
      push(p.first, p.second);
    }

    void push(const Item &i, const Prio &p) {
      int idx;
      if (free == -1) {
        idx = data.size();
        data.push_back(BucketItem(i));
      } else {
        idx = free;
        free = data[idx].next;
        data[idx].item = i;
      }
      index[i] = idx;
      if (p >= int(first.size())) first.resize(p + 1, -1);
      data[idx].next = first[p];
      first[p] = idx;
      if (p > maximal) {
	maximal = p;
      }
      ++num;
    }

    Item top() const {
      while (first[maximal] == -1) {
	--maximal;
      }
      return data[first[maximal]].item;
    }

    Prio prio() const {
      while (first[maximal] == -1) {
	--maximal;
      }
      return maximal;
    }

    void pop() {
      while (first[maximal] == -1) {
	--maximal;
      }
      int idx = first[maximal];
      index[data[idx].item] = -2;
      first[maximal] = data[idx].next;
      data[idx].next = free;
      free = idx;
      --num;
    }
    
    Prio operator[](const Item &i) const {
      for (int k = 0; k < first.size(); ++k) {
        int idx = first[k];
        while (idx != -1) {
          if (data[idx].item == i) {
            return k;
          }
          idx = data[idx].next;
        }
      }
      return -1;
    }

    State state(const Item &i) const {
      int idx = index[i];
      if (idx >= 0) idx = 0;
      return State(idx);
    }

  private:

    struct BucketItem {
      BucketItem(const Item& _item) : item(_item) {}

      Item item;

      int next;
    };

    ItemIntMap& index;
    std::vector<int> first;
    std::vector<BucketItem> data;
    int free, num;
    mutable int maximal;

  };

}
  
#endif

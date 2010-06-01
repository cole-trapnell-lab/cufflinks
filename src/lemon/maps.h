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

#ifndef LEMON_MAPS_H
#define LEMON_MAPS_H

#include <iterator>
#include <functional>
#include <vector>

#include <lemon/bits/utility.h>
#include <lemon/bits/traits.h>

///\file
///\ingroup maps
///\brief Miscellaneous property maps
///
#include <map>

namespace lemon {

  /// \addtogroup maps
  /// @{

  /// Base class of maps.

  /// Base class of maps.
  /// It provides the necessary <tt>typedef</tt>s required by the map concept.
  template<typename K, typename T>
  class MapBase {
  public:
    /// The key type of the map.
    typedef K Key;
    /// The value type of the map. (The type of objects associated with the keys).
    typedef T Value;
  };

  /// Null map. (a.k.a. DoNothingMap)

  /// This map can be used if you have to provide a map only for
  /// its type definitions, or if you have to provide a writable map, 
  /// but data written to it is not required (i.e. it will be sent to 
  /// <tt>/dev/null</tt>).
  template<typename K, typename T>
  class NullMap : public MapBase<K, T> {
  public:
    typedef MapBase<K, T> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;
    
    /// Gives back a default constructed element.
    T operator[](const K&) const { return T(); }
    /// Absorbs the value.
    void set(const K&, const T&) {}
  };

  ///Returns a \c NullMap class

  ///This function just returns a \c NullMap class.
  ///\relates NullMap
  template <typename K, typename V> 
  NullMap<K, V> nullMap() {
    return NullMap<K, V>();
  }


  /// Constant map.

  /// This is a \ref concepts::ReadMap "readable" map which assigns a 
  /// specified value to each key.
  /// In other aspects it is equivalent to \c NullMap.
  template<typename K, typename T>
  class ConstMap : public MapBase<K, T> {
  private:
    T v;
  public:

    typedef MapBase<K, T> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    /// Default constructor

    /// Default constructor.
    /// The value of the map will be uninitialized. 
    /// (More exactly it will be default constructed.)
    ConstMap() {}
    
    /// Constructor with specified initial value

    /// Constructor with specified initial value.
    /// \param _v is the initial value of the map.
    ConstMap(const T &_v) : v(_v) {}
    
    ///\e
    T operator[](const K&) const { return v; }

    ///\e
    void setAll(const T &t) {
      v = t;
    }    

    template<typename T1>
    struct rebind {
      typedef ConstMap<K, T1> other;
    };

    template<typename T1>
    ConstMap(const ConstMap<K, T1> &, const T &_v) : v(_v) {}
  };

  ///Returns a \c ConstMap class

  ///This function just returns a \c ConstMap class.
  ///\relates ConstMap
  template<typename K, typename V> 
  inline ConstMap<K, V> constMap(const V &v) {
    return ConstMap<K, V>(v);
  }


  template<typename T, T v>
  struct Const { };

  /// Constant map with inlined constant value.

  /// This is a \ref concepts::ReadMap "readable" map which assigns a 
  /// specified value to each key.
  /// In other aspects it is equivalent to \c NullMap.
  template<typename K, typename V, V v>
  class ConstMap<K, Const<V, v> > : public MapBase<K, V> {
  public:
    typedef MapBase<K, V> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ConstMap() { }
    ///\e
    V operator[](const K&) const { return v; }
    ///\e
    void set(const K&, const V&) { }
  };

  ///Returns a \c ConstMap class with inlined value

  ///This function just returns a \c ConstMap class with inlined value.
  ///\relates ConstMap
  template<typename K, typename V, V v> 
  inline ConstMap<K, Const<V, v> > constMap() {
    return ConstMap<K, Const<V, v> >();
  }

  ///Map based on \c std::map

  ///This is essentially a wrapper for \c std::map with addition that
  ///you can specify a default value different from \c Value() .
  ///It meets the \ref concepts::ReferenceMap "ReferenceMap" concept.
  template <typename K, typename T, typename Compare = std::less<K> >
  class StdMap : public MapBase<K, T> {
    template <typename K1, typename T1, typename C1>
    friend class StdMap;
  public:

    typedef MapBase<K, T> Parent;
    ///Key type
    typedef typename Parent::Key Key;
    ///Value type
    typedef typename Parent::Value Value;
    ///Reference Type
    typedef T& Reference;
    ///Const reference type
    typedef const T& ConstReference;

    typedef True ReferenceMapTag;

  private:
    
    typedef std::map<K, T, Compare> Map;
    Value _value;
    Map _map;

  public:

    /// Constructor with specified default value
    StdMap(const T& value = T()) : _value(value) {}
    /// \brief Constructs the map from an appropriate \c std::map, and 
    /// explicitly specifies a default value.
    template <typename T1, typename Comp1>
    StdMap(const std::map<Key, T1, Comp1> &map, const T& value = T()) 
      : _map(map.begin(), map.end()), _value(value) {}
    
    /// \brief Constructs a map from an other \ref StdMap.
    template<typename T1, typename Comp1>
    StdMap(const StdMap<Key, T1, Comp1> &c) 
      : _map(c._map.begin(), c._map.end()), _value(c._value) {}

  private:

    StdMap& operator=(const StdMap&);

  public:

    ///\e
    Reference operator[](const Key &k) {
      typename Map::iterator it = _map.lower_bound(k);
      if (it != _map.end() && !_map.key_comp()(k, it->first))
	return it->second;
      else
	return _map.insert(it, std::make_pair(k, _value))->second;
    }

    /// \e 
    ConstReference operator[](const Key &k) const {
      typename Map::const_iterator it = _map.find(k);
      if (it != _map.end())
	return it->second;
      else
	return _value;
    }

    /// \e 
    void set(const Key &k, const T &t) {
      typename Map::iterator it = _map.lower_bound(k);
      if (it != _map.end() && !_map.key_comp()(k, it->first))
	it->second = t;
      else
	_map.insert(it, std::make_pair(k, t));
    }

    /// \e
    void setAll(const T &t) {
      _value = t;
      _map.clear();
    }    

    template <typename T1, typename C1 = std::less<T1> >
    struct rebind {
      typedef StdMap<Key, T1, C1> other;
    };
  };

  ///Returns a \c StdMap class

  ///This function just returns a \c StdMap class with specified 
  ///default value.
  ///\relates StdMap
  template<typename K, typename V, typename Compare> 
  inline StdMap<K, V, Compare> stdMap(const V& value = V()) {
    return StdMap<K, V, Compare>(value);
  }
  
  template<typename K, typename V> 
  inline StdMap<K, V, std::less<K> > stdMap(const V& value = V()) {
    return StdMap<K, V, std::less<K> >(value);
  }
  
  ///Returns a \c StdMap class created from an appropriate \c std::map

  ///This function just returns a \c StdMap class created from an 
  ///appropriate \c std::map.
  ///\relates StdMap
  template<typename K, typename V, typename Compare> 
  inline StdMap<K, V, Compare> stdMap( const std::map<K, V, Compare> &map, 
                                       const V& value = V() ) {
    return StdMap<K, V, Compare>(map, value);
  }

  /// \brief Map for storing values for keys from the range <tt>[0..size-1]</tt>
  ///
  /// This map has the <tt>[0..size-1]</tt> keyset and the values
  /// are stored in a \c std::vector<T>  container. It can be used with
  /// some data structures, for example \c UnionFind, \c BinHeap, when 
  /// the used items are small integer numbers.
  template <typename T>
  class IntegerMap : public MapBase<int, T> {

    template <typename T1>
    friend class IntegerMap;

  public:

    typedef MapBase<int, T> Parent;
    ///\e
    typedef typename Parent::Key Key;
    ///\e
    typedef typename Parent::Value Value;
    ///\e
    typedef T& Reference;
    ///\e
    typedef const T& ConstReference;

    typedef True ReferenceMapTag;

  private:
    
    typedef std::vector<T> Vector;
    Vector _vector;

  public:

    /// Constructor with specified default value
    IntegerMap(int size = 0, const T& value = T()) : _vector(size, value) {}

    /// \brief Constructs the map from an appropriate \c std::vector.
    template <typename T1>
    IntegerMap(const std::vector<T1>& vector) 
      : _vector(vector.begin(), vector.end()) {}
    
    /// \brief Constructs a map from an other \ref IntegerMap.
    template <typename T1>
    IntegerMap(const IntegerMap<T1> &c) 
      : _vector(c._vector.begin(), c._vector.end()) {}

    /// \brief Resize the container
    void resize(int size, const T& value = T()) {
      _vector.resize(size, value);
    }

  private:

    IntegerMap& operator=(const IntegerMap&);

  public:

    ///\e
    Reference operator[](Key k) {
      return _vector[k];
    }

    /// \e 
    ConstReference operator[](Key k) const {
      return _vector[k];
    }

    /// \e 
    void set(const Key &k, const T& t) {
      _vector[k] = t;
    }

  };

  ///Returns an \c IntegerMap class

  ///This function just returns an \c IntegerMap class.
  ///\relates IntegerMap
  template<typename T>
  inline IntegerMap<T> integerMap(int size = 0, const T& value = T()) {
    return IntegerMap<T>(size, value);
  }

  /// @}

  /// \addtogroup map_adaptors
  /// @{

  /// \brief Identity map.
  ///
  /// This map gives back the given key as value without any
  /// modification. 
  template <typename T>
  class IdentityMap : public MapBase<T, T> {
  public:
    typedef MapBase<T, T> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    /// \e
    const T& operator[](const T& t) const {
      return t;
    }
  };

  ///Returns an \c IdentityMap class

  ///This function just returns an \c IdentityMap class.
  ///\relates IdentityMap
  template<typename T>
  inline IdentityMap<T> identityMap() {
    return IdentityMap<T>();
  }
  

  ///\brief Convert the \c Value of a map to another type using
  ///the default conversion.
  ///
  ///This \ref concepts::ReadMap "read only map"
  ///converts the \c Value of a map to type \c T.
  ///Its \c Key is inherited from \c M.
  template <typename M, typename T> 
  class ConvertMap : public MapBase<typename M::Key, T> {
    const M& m;
  public:
    typedef MapBase<typename M::Key, T> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor

    ///Constructor
    ///\param _m is the underlying map
    ConvertMap(const M &_m) : m(_m) {};

    ///\e
    Value operator[](const Key& k) const {return m[k];}
  };
  
  ///Returns a \c ConvertMap class

  ///This function just returns a \c ConvertMap class.
  ///\relates ConvertMap
  template<typename T, typename M>
  inline ConvertMap<M, T> convertMap(const M &m) {
    return ConvertMap<M, T>(m);
  }

  ///Simple wrapping of a map

  ///This \ref concepts::ReadMap "read only map" returns the simple
  ///wrapping of the given map. Sometimes the reference maps cannot be
  ///combined with simple read maps. This map adaptor wraps the given
  ///map to simple read map.
  ///
  ///\sa SimpleWriteMap
  ///
  /// \todo Revise the misleading name 
  template<typename M> 
  class SimpleMap : public MapBase<typename M::Key, typename M::Value> {
    const M& m;

  public:
    typedef MapBase<typename M::Key, typename M::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    SimpleMap(const M &_m) : m(_m) {};
    ///\e
    Value operator[](Key k) const {return m[k];}
  };

  ///Returns a \c SimpleMap class

  ///This function just returns a \c SimpleMap class.
  ///\relates SimpleMap
  template<typename M>
  inline SimpleMap<M> simpleMap(const M &m) {
    return SimpleMap<M>(m);
  }

  ///Simple writable wrapping of a map

  ///This \ref concepts::ReadWriteMap "read-write map" returns the simple
  ///wrapping of the given map. Sometimes the reference maps cannot be
  ///combined with simple read-write maps. This map adaptor wraps the
  ///given map to simple read-write map.
  ///
  ///\sa SimpleMap
  ///
  /// \todo Revise the misleading name
  template<typename M> 
  class SimpleWriteMap : public MapBase<typename M::Key, typename M::Value> {
    M& m;

  public:
    typedef MapBase<typename M::Key, typename M::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    SimpleWriteMap(M &_m) : m(_m) {};
    ///\e
    Value operator[](Key k) const {return m[k];}
    ///\e
    void set(Key k, const Value& c) { m.set(k, c); }
  };

  ///Returns a \c SimpleWriteMap class

  ///This function just returns a \c SimpleWriteMap class.
  ///\relates SimpleWriteMap
  template<typename M>
  inline SimpleWriteMap<M> simpleWriteMap(M &m) {
    return SimpleWriteMap<M>(m);
  }

  ///Sum of two maps

  ///This \ref concepts::ReadMap "read only map" returns the sum of the two
  ///given maps.
  ///Its \c Key and \c Value are inherited from \c M1.
  ///The \c Key and \c Value of \c M2 must be convertible to those of \c M1.
  template<typename M1, typename M2> 
  class AddMap : public MapBase<typename M1::Key, typename M1::Value> {
    const M1& m1;
    const M2& m2;

  public:
    typedef MapBase<typename M1::Key, typename M1::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    AddMap(const M1 &_m1,const M2 &_m2) : m1(_m1), m2(_m2) {};
    ///\e
    Value operator[](Key k) const {return m1[k]+m2[k];}
  };
  
  ///Returns an \c AddMap class

  ///This function just returns an \c AddMap class.
  ///\todo How to call these type of functions?
  ///
  ///\relates AddMap
  template<typename M1, typename M2> 
  inline AddMap<M1, M2> addMap(const M1 &m1,const M2 &m2) {
    return AddMap<M1, M2>(m1,m2);
  }

  ///Shift a map with a constant.

  ///This \ref concepts::ReadMap "read only map" returns the sum of the
  ///given map and a constant value.
  ///Its \c Key and \c Value is inherited from \c M.
  ///
  ///Actually,
  ///\code
  ///  ShiftMap<X> sh(x,v);
  ///\endcode
  ///is equivalent to
  ///\code
  ///  ConstMap<X::Key, X::Value> c_tmp(v);
  ///  AddMap<X, ConstMap<X::Key, X::Value> > sh(x,v);
  ///\endcode
  ///
  ///\sa ShiftWriteMap
  template<typename M, typename C = typename M::Value> 
  class ShiftMap : public MapBase<typename M::Key, typename M::Value> {
    const M& m;
    C v;
  public:
    typedef MapBase<typename M::Key, typename M::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor

    ///Constructor
    ///\param _m is the undelying map
    ///\param _v is the shift value
    ShiftMap(const M &_m, const C &_v ) : m(_m), v(_v) {};
    ///\e
    Value operator[](Key k) const {return m[k] + v;}
  };

  ///Shift a map with a constant (ReadWrite version).

  ///This \ref concepts::ReadWriteMap "read-write map" returns the sum of the
  ///given map and a constant value. It makes also possible to write the map.
  ///Its \c Key and \c Value are inherited from \c M.
  ///
  ///\sa ShiftMap
  template<typename M, typename C = typename M::Value> 
  class ShiftWriteMap : public MapBase<typename M::Key, typename M::Value> {
    M& m;
    C v;
  public:
    typedef MapBase<typename M::Key, typename M::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor

    ///Constructor
    ///\param _m is the undelying map
    ///\param _v is the shift value
    ShiftWriteMap(M &_m, const C &_v ) : m(_m), v(_v) {};
    /// \e
    Value operator[](Key k) const {return m[k] + v;}
    /// \e
    void set(Key k, const Value& c) { m.set(k, c - v); }
  };
  
  ///Returns a \c ShiftMap class

  ///This function just returns an \c ShiftMap class.
  ///\relates ShiftMap
  template<typename M, typename C> 
  inline ShiftMap<M, C> shiftMap(const M &m,const C &v) {
    return ShiftMap<M, C>(m,v);
  }

  ///Returns a \c ShiftWriteMap class

  ///This function just returns a \c ShiftWriteMap class.
  ///\relates ShiftWriteMap
  template<typename M, typename C> 
  inline ShiftWriteMap<M, C> shiftMap(M &m,const C &v) {
    return ShiftWriteMap<M, C>(m,v);
  }

  ///Difference of two maps

  ///This \ref concepts::ReadMap "read only map" returns the difference
  ///of the values of the two given maps.
  ///Its \c Key and \c Value are inherited from \c M1.
  ///The \c Key and \c Value of \c M2 must be convertible to those of \c M1.

  template<typename M1, typename M2> 
  class SubMap : public MapBase<typename M1::Key, typename M1::Value> {
    const M1& m1;
    const M2& m2;
  public:
    typedef MapBase<typename M1::Key, typename M1::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    SubMap(const M1 &_m1,const M2 &_m2) : m1(_m1), m2(_m2) {};
    /// \e
    Value operator[](Key k) const {return m1[k]-m2[k];}
  };
  
  ///Returns a \c SubMap class

  ///This function just returns a \c SubMap class.
  ///
  ///\relates SubMap
  template<typename M1, typename M2> 
  inline SubMap<M1, M2> subMap(const M1 &m1, const M2 &m2) {
    return SubMap<M1, M2>(m1, m2);
  }

  ///Product of two maps

  ///This \ref concepts::ReadMap "read only map" returns the product of the
  ///values of the two given maps.
  ///Its \c Key and \c Value are inherited from \c M1.
  ///The \c Key and \c Value of \c M2 must be convertible to those of \c M1.
  template<typename M1, typename M2> 
  class MulMap : public MapBase<typename M1::Key, typename M1::Value> {
    const M1& m1;
    const M2& m2;
  public:
    typedef MapBase<typename M1::Key, typename M1::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    MulMap(const M1 &_m1,const M2 &_m2) : m1(_m1), m2(_m2) {};
    /// \e
    Value operator[](Key k) const {return m1[k]*m2[k];}
  };
  
  ///Returns a \c MulMap class

  ///This function just returns a \c MulMap class.
  ///\relates MulMap
  template<typename M1, typename M2> 
  inline MulMap<M1, M2> mulMap(const M1 &m1,const M2 &m2) {
    return MulMap<M1, M2>(m1,m2);
  }
 
  ///Scales a map with a constant.

  ///This \ref concepts::ReadMap "read only map" returns the value of the
  ///given map multiplied from the left side with a constant value.
  ///Its \c Key and \c Value are inherited from \c M.
  ///
  ///Actually,
  ///\code
  ///  ScaleMap<X> sc(x,v);
  ///\endcode
  ///is equivalent to
  ///\code
  ///  ConstMap<X::Key, X::Value> c_tmp(v);
  ///  MulMap<X, ConstMap<X::Key, X::Value> > sc(x,v);
  ///\endcode
  ///
  ///\sa ScaleWriteMap
  template<typename M, typename C = typename M::Value> 
  class ScaleMap : public MapBase<typename M::Key, typename M::Value> {
    const M& m;
    C v;
  public:
    typedef MapBase<typename M::Key, typename M::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor

    ///Constructor
    ///\param _m is the undelying map
    ///\param _v is the scaling value
    ScaleMap(const M &_m, const C &_v ) : m(_m), v(_v) {};
    /// \e
    Value operator[](Key k) const {return v * m[k];}
  };

  ///Scales a map with a constant (ReadWrite version).

  ///This \ref concepts::ReadWriteMap "read-write map" returns the value of the
  ///given map multiplied from the left side with a constant value. It can
  ///also be used as write map if the \c / operator is defined between
  ///\c Value and \c C and the given multiplier is not zero.
  ///Its \c Key and \c Value are inherited from \c M.
  ///
  ///\sa ScaleMap
  template<typename M, typename C = typename M::Value> 
  class ScaleWriteMap : public MapBase<typename M::Key, typename M::Value> {
    M& m;
    C v;
  public:
    typedef MapBase<typename M::Key, typename M::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor

    ///Constructor
    ///\param _m is the undelying map
    ///\param _v is the scaling value
    ScaleWriteMap(M &_m, const C &_v ) : m(_m), v(_v) {};
    /// \e
    Value operator[](Key k) const {return v * m[k];}
    /// \e
    void set(Key k, const Value& c) { m.set(k, c / v);}
  };
  
  ///Returns a \c ScaleMap class

  ///This function just returns a \c ScaleMap class.
  ///\relates ScaleMap
  template<typename M, typename C> 
  inline ScaleMap<M, C> scaleMap(const M &m,const C &v) {
    return ScaleMap<M, C>(m,v);
  }

  ///Returns a \c ScaleWriteMap class

  ///This function just returns a \c ScaleWriteMap class.
  ///\relates ScaleWriteMap
  template<typename M, typename C> 
  inline ScaleWriteMap<M, C> scaleMap(M &m,const C &v) {
    return ScaleWriteMap<M, C>(m,v);
  }

  ///Quotient of two maps

  ///This \ref concepts::ReadMap "read only map" returns the quotient of the
  ///values of the two given maps.
  ///Its \c Key and \c Value are inherited from \c M1.
  ///The \c Key and \c Value of \c M2 must be convertible to those of \c M1.
  template<typename M1, typename M2> 
  class DivMap : public MapBase<typename M1::Key, typename M1::Value> {
    const M1& m1;
    const M2& m2;
  public:
    typedef MapBase<typename M1::Key, typename M1::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    DivMap(const M1 &_m1,const M2 &_m2) : m1(_m1), m2(_m2) {};
    /// \e
    Value operator[](Key k) const {return m1[k]/m2[k];}
  };
  
  ///Returns a \c DivMap class

  ///This function just returns a \c DivMap class.
  ///\relates DivMap
  template<typename M1, typename M2> 
  inline DivMap<M1, M2> divMap(const M1 &m1,const M2 &m2) {
    return DivMap<M1, M2>(m1,m2);
  }
  
  ///Composition of two maps

  ///This \ref concepts::ReadMap "read only map" returns the composition of
  ///two given maps.
  ///That is to say, if \c m1 is of type \c M1 and \c m2 is of \c M2,
  ///then for
  ///\code
  ///  ComposeMap<M1, M2> cm(m1,m2);
  ///\endcode
  /// <tt>cm[x]</tt> will be equal to <tt>m1[m2[x]]</tt>.
  ///
  ///Its \c Key is inherited from \c M2 and its \c Value is from \c M1.
  ///\c M2::Value must be convertible to \c M1::Key.
  ///
  ///\sa CombineMap
  ///
  ///\todo Check the requirements.
  template <typename M1, typename M2> 
  class ComposeMap : public MapBase<typename M2::Key, typename M1::Value> {
    const M1& m1;
    const M2& m2;
  public:
    typedef MapBase<typename M2::Key, typename M1::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    ComposeMap(const M1 &_m1,const M2 &_m2) : m1(_m1), m2(_m2) {};
    
    typename MapTraits<M1>::ConstReturnValue
    /// \e
    operator[](Key k) const {return m1[m2[k]];}
  };
  ///Returns a \c ComposeMap class

  ///This function just returns a \c ComposeMap class.
  ///
  ///\relates ComposeMap
  template <typename M1, typename M2> 
  inline ComposeMap<M1, M2> composeMap(const M1 &m1,const M2 &m2) {
    return ComposeMap<M1, M2>(m1,m2);
  }
  
  ///Combine of two maps using an STL (binary) functor.

  ///Combine of two maps using an STL (binary) functor.
  ///
  ///This \ref concepts::ReadMap "read only map" takes two maps and a
  ///binary functor and returns the composition of the two
  ///given maps unsing the functor. 
  ///That is to say, if \c m1 and \c m2 is of type \c M1 and \c M2
  ///and \c f is of \c F, then for
  ///\code
  ///  CombineMap<M1, M2,F,V> cm(m1,m2,f);
  ///\endcode
  /// <tt>cm[x]</tt> will be equal to <tt>f(m1[x],m2[x])</tt>
  ///
  ///Its \c Key is inherited from \c M1 and its \c Value is \c V.
  ///\c M2::Value and \c M1::Value must be convertible to the corresponding
  ///input parameter of \c F and the return type of \c F must be convertible
  ///to \c V.
  ///
  ///\sa ComposeMap
  ///
  ///\todo Check the requirements.
  template<typename M1, typename M2, typename F,
	   typename V = typename F::result_type> 
  class CombineMap : public MapBase<typename M1::Key, V> {
    const M1& m1;
    const M2& m2;
    F f;
  public:
    typedef MapBase<typename M1::Key, V> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    CombineMap(const M1 &_m1,const M2 &_m2,const F &_f = F())
      : m1(_m1), m2(_m2), f(_f) {};
    /// \e
    Value operator[](Key k) const {return f(m1[k],m2[k]);}
  };
  
  ///Returns a \c CombineMap class

  ///This function just returns a \c CombineMap class.
  ///
  ///For example if \c m1 and \c m2 are both \c double valued maps, then 
  ///\code
  ///combineMap(m1,m2,std::plus<double>())
  ///\endcode
  ///is equivalent to
  ///\code
  ///addMap(m1,m2)
  ///\endcode
  ///
  ///This function is specialized for adaptable binary function
  ///classes and C++ functions.
  ///
  ///\relates CombineMap
  template<typename M1, typename M2, typename F, typename V> 
  inline CombineMap<M1, M2, F, V> 
  combineMap(const M1& m1,const M2& m2, const F& f) {
    return CombineMap<M1, M2, F, V>(m1,m2,f);
  }

  template<typename M1, typename M2, typename F> 
  inline CombineMap<M1, M2, F, typename F::result_type> 
  combineMap(const M1& m1, const M2& m2, const F& f) {
    return combineMap<M1, M2, F, typename F::result_type>(m1,m2,f);
  }

  template<typename M1, typename M2, typename K1, typename K2, typename V> 
  inline CombineMap<M1, M2, V (*)(K1, K2), V> 
  combineMap(const M1 &m1, const M2 &m2, V (*f)(K1, K2)) {
    return combineMap<M1, M2, V (*)(K1, K2), V>(m1,m2,f);
  }

  ///Negative value of a map

  ///This \ref concepts::ReadMap "read only map" returns the negative
  ///value of the value returned by the given map.
  ///Its \c Key and \c Value are inherited from \c M.
  ///The unary \c - operator must be defined for \c Value, of course.
  ///
  ///\sa NegWriteMap
  template<typename M> 
  class NegMap : public MapBase<typename M::Key, typename M::Value> {
    const M& m;
  public:
    typedef MapBase<typename M::Key, typename M::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    NegMap(const M &_m) : m(_m) {};
    /// \e
    Value operator[](Key k) const {return -m[k];}
  };
  
  ///Negative value of a map (ReadWrite version)

  ///This \ref concepts::ReadWriteMap "read-write map" returns the negative
  ///value of the value returned by the given map.
  ///Its \c Key and \c Value are inherited from \c M.
  ///The unary \c - operator must be defined for \c Value, of course.
  ///
  /// \sa NegMap
  template<typename M> 
  class NegWriteMap : public MapBase<typename M::Key, typename M::Value> {
    M& m;
  public:
    typedef MapBase<typename M::Key, typename M::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    NegWriteMap(M &_m) : m(_m) {};
    /// \e
    Value operator[](Key k) const {return -m[k];}
    /// \e
    void set(Key k, const Value& v) { m.set(k, -v); }
  };

  ///Returns a \c NegMap class

  ///This function just returns a \c NegMap class.
  ///\relates NegMap
  template <typename M> 
  inline NegMap<M> negMap(const M &m) {
    return NegMap<M>(m);
  }

  ///Returns a \c NegWriteMap class

  ///This function just returns a \c NegWriteMap class.
  ///\relates NegWriteMap
  template <typename M> 
  inline NegWriteMap<M> negMap(M &m) {
    return NegWriteMap<M>(m);
  }

  ///Absolute value of a map

  ///This \ref concepts::ReadMap "read only map" returns the absolute value
  ///of the value returned by the given map.
  ///Its \c Key and \c Value are inherited from \c M. 
  ///\c Value must be comparable to \c 0 and the unary \c -
  ///operator must be defined for it, of course.
  ///
  ///\bug We need a unified way to handle the situation below:
  ///\code
  ///  struct _UnConvertible {};
  ///  template<class A> inline A t_abs(A a) {return _UnConvertible();}
  ///  template<> inline int t_abs<>(int n) {return abs(n);}
  ///  template<> inline long int t_abs<>(long int n) {return labs(n);}
  ///  template<> inline long long int t_abs<>(long long int n) {return ::llabs(n);}
  ///  template<> inline float t_abs<>(float n) {return fabsf(n);}
  ///  template<> inline double t_abs<>(double n) {return fabs(n);}
  ///  template<> inline long double t_abs<>(long double n) {return fabsl(n);}
  ///\endcode
  

  template<typename M> 
  class AbsMap : public MapBase<typename M::Key, typename M::Value> {
    const M& m;
  public:
    typedef MapBase<typename M::Key, typename M::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    AbsMap(const M &_m) : m(_m) {};
    /// \e
    Value operator[](Key k) const {
      Value tmp = m[k]; 
      return tmp >= 0 ? tmp : -tmp;
    }

  };
  
  ///Returns an \c AbsMap class

  ///This function just returns an \c AbsMap class.
  ///\relates AbsMap
  template<typename M> 
  inline AbsMap<M> absMap(const M &m) {
    return AbsMap<M>(m);
  }

  ///Converts an STL style functor to a map

  ///This \ref concepts::ReadMap "read only map" returns the value
  ///of a given functor.
  ///
  ///Template parameters \c K and \c V will become its
  ///\c Key and \c Value. 
  ///In most cases they have to be given explicitly because a 
  ///functor typically does not provide \c argument_type and 
  ///\c result_type typedefs.
  ///
  ///Parameter \c F is the type of the used functor.
  ///
  ///\sa MapFunctor
  template<typename F, 
	   typename K = typename F::argument_type, 
	   typename V = typename F::result_type> 
  class FunctorMap : public MapBase<K, V> {
    F f;
  public:
    typedef MapBase<K, V> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    FunctorMap(const F &_f = F()) : f(_f) {}
    /// \e
    Value operator[](Key k) const { return f(k);}
  };
  
  ///Returns a \c FunctorMap class

  ///This function just returns a \c FunctorMap class.
  ///
  ///This function is specialized for adaptable binary function
  ///classes and C++ functions.
  ///
  ///\relates FunctorMap
  template<typename K, typename V, typename F> inline 
  FunctorMap<F, K, V> functorMap(const F &f) {
    return FunctorMap<F, K, V>(f);
  }

  template <typename F> inline 
  FunctorMap<F, typename F::argument_type, typename F::result_type> 
  functorMap(const F &f) {
    return FunctorMap<F, typename F::argument_type, 
      typename F::result_type>(f);
  }

  template <typename K, typename V> inline 
  FunctorMap<V (*)(K), K, V> functorMap(V (*f)(K)) {
    return FunctorMap<V (*)(K), K, V>(f);
  }


  ///Converts a map to an STL style (unary) functor

  ///This class Converts a map to an STL style (unary) functor.
  ///That is it provides an <tt>operator()</tt> to read its values.
  ///
  ///For the sake of convenience it also works as
  ///a ususal \ref concepts::ReadMap "readable map",
  ///i.e. <tt>operator[]</tt> and the \c Key and \c Value typedefs also exist.
  ///
  ///\sa FunctorMap
  template <typename M> 
  class MapFunctor : public MapBase<typename M::Key, typename M::Value> {
    const M& m;
  public:
    typedef MapBase<typename M::Key, typename M::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    typedef typename M::Key argument_type;
    typedef typename M::Value result_type;

    ///Constructor
    MapFunctor(const M &_m) : m(_m) {};
    ///\e
    Value operator()(Key k) const {return m[k];}
    ///\e
    Value operator[](Key k) const {return m[k];}
  };
  
  ///Returns a \c MapFunctor class

  ///This function just returns a \c MapFunctor class.
  ///\relates MapFunctor
  template<typename M> 
  inline MapFunctor<M> mapFunctor(const M &m) {
    return MapFunctor<M>(m);
  }

  ///Just readable version of \ref ForkWriteMap

  ///This map has two \ref concepts::ReadMap "readable map"
  ///parameters and each read request will be passed just to the
  ///first map. This class is the just readable map type of \c ForkWriteMap.
  ///
  ///The \c Key and \c Value are inherited from \c M1.
  ///The \c Key and \c Value of \c M2 must be convertible from those of \c M1.
  ///
  ///\sa ForkWriteMap

  template<typename  M1, typename M2> 
  class ForkMap : public MapBase<typename M1::Key, typename M1::Value> {
    const M1& m1;
    const M2& m2;
  public:
    typedef MapBase<typename M1::Key, typename M1::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    ForkMap(const M1 &_m1, const M2 &_m2) : m1(_m1), m2(_m2) {};
    /// \e
    Value operator[](Key k) const {return m1[k];}
  };


  ///Applies all map setting operations to two maps

  ///This map has two \ref concepts::WriteMap "writable map"
  ///parameters and each write request will be passed to both of them.
  ///If \c M1 is also \ref concepts::ReadMap "readable",
  ///then the read operations will return the
  ///corresponding values of \c M1.
  ///
  ///The \c Key and \c Value are inherited from \c M1.
  ///The \c Key and \c Value of \c M2 must be convertible from those of \c M1.
  ///
  ///\sa ForkMap
  template<typename  M1, typename M2> 
  class ForkWriteMap : public MapBase<typename M1::Key, typename M1::Value> {
    M1& m1;
    M2& m2;
  public:
    typedef MapBase<typename M1::Key, typename M1::Value> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    ///Constructor
    ForkWriteMap(M1 &_m1, M2 &_m2) : m1(_m1), m2(_m2) {};
    ///\e
    Value operator[](Key k) const {return m1[k];}
    ///\e
    void set(Key k, const Value &v) {m1.set(k,v); m2.set(k,v);}
  };
  
  ///Returns a \c ForkMap class

  ///This function just returns a \c ForkMap class.
  ///\relates ForkMap
  template <typename M1, typename M2> 
  inline ForkMap<M1, M2> forkMap(const M1 &m1, const M2 &m2) {
    return ForkMap<M1, M2>(m1,m2);
  }

  ///Returns a \c ForkWriteMap class

  ///This function just returns a \c ForkWriteMap class.
  ///\relates ForkWriteMap
  template <typename M1, typename M2> 
  inline ForkWriteMap<M1, M2> forkMap(M1 &m1, M2 &m2) {
    return ForkWriteMap<M1, M2>(m1,m2);
  }


  
  /* ************* BOOL MAPS ******************* */
  
  ///Logical 'not' of a map
  
  ///This bool \ref concepts::ReadMap "read only map" returns the 
  ///logical negation of the value returned by the given map.
  ///Its \c Key is inherited from \c M, its \c Value is \c bool.
  ///
  ///\sa NotWriteMap
  template <typename M> 
  class NotMap : public MapBase<typename M::Key, bool> {
    const M& m;
  public:
    typedef MapBase<typename M::Key, bool> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    /// Constructor
    NotMap(const M &_m) : m(_m) {};
    ///\e
    Value operator[](Key k) const {return !m[k];}
  };

  ///Logical 'not' of a map (ReadWrie version)
  
  ///This bool \ref concepts::ReadWriteMap "read-write map" returns the 
  ///logical negation of the value returned by the given map. When it is set,
  ///the opposite value is set to the original map.
  ///Its \c Key is inherited from \c M, its \c Value is \c bool.
  ///
  ///\sa NotMap
  template <typename M> 
  class NotWriteMap : public MapBase<typename M::Key, bool> {
    M& m;
  public:
    typedef MapBase<typename M::Key, bool> Parent;
    typedef typename Parent::Key Key;
    typedef typename Parent::Value Value;

    /// Constructor
    NotWriteMap(M &_m) : m(_m) {};
    ///\e
    Value operator[](Key k) const {return !m[k];}
    ///\e
    void set(Key k, bool v) { m.set(k, !v); }
  };
  
  ///Returns a \c NotMap class
  
  ///This function just returns a \c NotMap class.
  ///\relates NotMap
  template <typename M> 
  inline NotMap<M> notMap(const M &m) {
    return NotMap<M>(m);
  }
  
  ///Returns a \c NotWriteMap class
  
  ///This function just returns a \c NotWriteMap class.
  ///\relates NotWriteMap
  template <typename M> 
  inline NotWriteMap<M> notMap(M &m) {
    return NotWriteMap<M>(m);
  }

  namespace _maps_bits {

    template <typename Value>
    struct Identity {
      typedef Value argument_type;
      typedef Value result_type;
      Value operator()(const Value& val) const {
	return val;
      }
    };

    template <typename _Iterator, typename Enable = void>
    struct IteratorTraits {
      typedef typename std::iterator_traits<_Iterator>::value_type Value;
    };

    template <typename _Iterator>
    struct IteratorTraits<_Iterator,
      typename exists<typename _Iterator::container_type>::type> 
    {
      typedef typename _Iterator::container_type::value_type Value;
    };

  }
  

  /// \brief Writable bool map for logging each \c true assigned element
  ///
  /// A \ref concepts::ReadWriteMap "read-write" bool map for logging 
  /// each \c true assigned element, i.e it copies all the keys set 
  /// to \c true to the given iterator.
  ///
  /// \note The container of the iterator should contain space 
  /// for each element.
  ///
  /// The following example shows how you can write the edges found by 
  /// the \ref Prim algorithm directly to the standard output.
  ///\code
  /// typedef IdMap<UGraph, UEdge> UEdgeIdMap;
  /// UEdgeIdMap uedgeId(ugraph);
  ///
  /// typedef MapFunctor<UEdgeIdMap> UEdgeIdFunctor;
  /// UEdgeIdFunctor uedgeIdFunctor(uedgeId);
  ///
  /// StoreBoolMap<ostream_iterator<int>, UEdgeIdFunctor> 
  ///   writerMap(ostream_iterator<int>(cout, " "), uedgeIdFunctor);
  ///
  /// prim(ugraph, cost, writerMap);
  ///\endcode
  ///
  ///\sa BackInserterBoolMap 
  ///\sa FrontInserterBoolMap 
  ///\sa InserterBoolMap 
  template <typename _Iterator, 
            typename _Functor =
            _maps_bits::Identity<typename _maps_bits::
                                 IteratorTraits<_Iterator>::Value> >
  class StoreBoolMap {
  public:
    typedef _Iterator Iterator;

    typedef typename _Functor::argument_type Key;
    typedef bool Value;

    typedef _Functor Functor;

    /// Constructor
    StoreBoolMap(Iterator it, const Functor& functor = Functor()) 
      : _begin(it), _end(it), _functor(functor) {}

    /// Gives back the given iterator set for the first key
    Iterator begin() const {
      return _begin;
    }
 
    /// Gives back the the 'after the last' iterator
    Iterator end() const {
      return _end;
    }

    /// The \c set function of the map
    void set(const Key& key, Value value) const {
      if (value) {
	*_end++ = _functor(key);
      }
    }
    
  private:
    Iterator _begin;
    mutable Iterator _end;
    Functor _functor;
  };

  /// \brief Writable bool map for logging each \c true assigned element in 
  /// a back insertable container.
  ///
  /// Writable bool map for logging each \c true assigned element by pushing
  /// them into a back insertable container.
  /// It can be used to retrieve the items into a standard
  /// container. The next example shows how you can store the
  /// edges found by the Prim algorithm in a vector.
  ///
  ///\code
  /// vector<UEdge> span_tree_uedges;
  /// BackInserterBoolMap<vector<UEdge> > inserter_map(span_tree_uedges);
  /// prim(ugraph, cost, inserter_map);
  ///\endcode
  ///
  ///\sa StoreBoolMap
  ///\sa FrontInserterBoolMap
  ///\sa InserterBoolMap
  template <typename Container,
            typename Functor =
            _maps_bits::Identity<typename Container::value_type> >
  class BackInserterBoolMap {
  public:
    typedef typename Functor::argument_type Key;
    typedef bool Value;

    /// Constructor
    BackInserterBoolMap(Container& _container, 
                        const Functor& _functor = Functor()) 
      : container(_container), functor(_functor) {}

    /// The \c set function of the map
    void set(const Key& key, Value value) {
      if (value) {
	container.push_back(functor(key));
      }
    }
    
  private:
    Container& container;
    Functor functor;
  };

  /// \brief Writable bool map for logging each \c true assigned element in 
  /// a front insertable container.
  ///
  /// Writable bool map for logging each \c true assigned element by pushing
  /// them into a front insertable container.
  /// It can be used to retrieve the items into a standard
  /// container. For example see \ref BackInserterBoolMap.
  ///
  ///\sa BackInserterBoolMap
  ///\sa InserterBoolMap
  template <typename Container,
            typename Functor =
            _maps_bits::Identity<typename Container::value_type> >
  class FrontInserterBoolMap {
  public:
    typedef typename Functor::argument_type Key;
    typedef bool Value;

    /// Constructor
    FrontInserterBoolMap(Container& _container,
                         const Functor& _functor = Functor()) 
      : container(_container), functor(_functor) {}

    /// The \c set function of the map
    void set(const Key& key, Value value) {
      if (value) {
	container.push_front(functor(key));
      }
    }
    
  private:
    Container& container;    
    Functor functor;
  };

  /// \brief Writable bool map for storing each \c true assigned element in 
  /// an insertable container.
  ///
  /// Writable bool map for storing each \c true assigned element in an 
  /// insertable container. It will insert all the keys set to \c true into
  /// the container.
  ///
  /// For example, if you want to store the cut arcs of the strongly
  /// connected components in a set you can use the next code:
  ///
  ///\code
  /// set<Edge> cut_edges;
  /// InserterBoolMap<set<Edge> > inserter_map(cut_edges);
  /// stronglyConnectedCutEdges(graph, cost, inserter_map);
  ///\endcode
  ///
  ///\sa BackInserterBoolMap
  ///\sa FrontInserterBoolMap
  template <typename Container,
            typename Functor =
            _maps_bits::Identity<typename Container::value_type> >
  class InserterBoolMap {
  public:
    typedef typename Container::value_type Key;
    typedef bool Value;

    /// Constructor with specified iterator
    
    /// Constructor with specified iterator.
    /// \param _container The container for storing the elements.
    /// \param _it The elements will be inserted before this iterator.
    /// \param _functor The functor that is used when an element is stored.
    InserterBoolMap(Container& _container, typename Container::iterator _it,
                    const Functor& _functor = Functor()) 
      : container(_container), it(_it), functor(_functor) {}

    /// Constructor

    /// Constructor without specified iterator.
    /// The elements will be inserted before <tt>_container.end()</tt>.
    /// \param _container The container for storing the elements.
    /// \param _functor The functor that is used when an element is stored.
    InserterBoolMap(Container& _container, const Functor& _functor = Functor())
      : container(_container), it(_container.end()), functor(_functor) {}

    /// The \c set function of the map
    void set(const Key& key, Value value) {
      if (value) {
	it = container.insert(it, functor(key));
        ++it;
      }
    }
    
  private:
    Container& container;
    typename Container::iterator it;
    Functor functor;
  };

  /// \brief Writable bool map for filling each \c true assigned element with a 
  /// given value.
  ///
  /// Writable bool map for filling each \c true assigned element with a 
  /// given value. The value can set the container.
  ///
  /// The following code finds the connected components of a graph
  /// and stores it in the \c comp map:
  ///\code
  /// typedef UGraph::NodeMap<int> ComponentMap;
  /// ComponentMap comp(ugraph);
  /// typedef FillBoolMap<UGraph::NodeMap<int> > ComponentFillerMap;
  /// ComponentFillerMap filler(comp, 0);
  ///
  /// Dfs<UGraph>::DefProcessedMap<ComponentFillerMap>::Create dfs(ugraph);
  /// dfs.processedMap(filler);
  /// dfs.init();
  /// for (NodeIt it(ugraph); it != INVALID; ++it) {
  ///   if (!dfs.reached(it)) {
  ///     dfs.addSource(it);
  ///     dfs.start();
  ///     ++filler.fillValue();
  ///   }
  /// }
  ///\endcode
  template <typename Map>
  class FillBoolMap {
  public:
    typedef typename Map::Key Key;
    typedef bool Value;

    /// Constructor
    FillBoolMap(Map& _map, const typename Map::Value& _fill) 
      : map(_map), fill(_fill) {}

    /// Constructor
    FillBoolMap(Map& _map) 
      : map(_map), fill() {}

    /// Gives back the current fill value
    const typename Map::Value& fillValue() const {
      return fill;
    } 

    /// Gives back the current fill value
    typename Map::Value& fillValue() {
      return fill;
    } 

    /// Sets the current fill value
    void fillValue(const typename Map::Value& _fill) {
      fill = _fill;
    } 

    /// The \c set function of the map
    void set(const Key& key, Value value) {
      if (value) {
	map.set(key, fill);
      }
    }
    
  private:
    Map& map;
    typename Map::Value fill;
  };


  /// \brief Writable bool map for storing the sequence number of 
  /// \c true assignments.  
  ///
  /// Writable bool map that stores for each \c true assigned elements  
  /// the sequence number of this setting.
  /// It makes it easy to calculate the leaving
  /// order of the nodes in the \c Dfs algorithm.
  ///
  ///\code
  /// typedef Graph::NodeMap<int> OrderMap;
  /// OrderMap order(graph);
  /// typedef SettingOrderBoolMap<OrderMap> OrderSetterMap;
  /// OrderSetterMap setter(order);
  /// Dfs<Graph>::DefProcessedMap<OrderSetterMap>::Create dfs(graph);
  /// dfs.processedMap(setter);
  /// dfs.init();
  /// for (NodeIt it(graph); it != INVALID; ++it) {
  ///   if (!dfs.reached(it)) {
  ///     dfs.addSource(it);
  ///     dfs.start();
  ///   }
  /// }
  ///\endcode
  ///
  /// The storing of the discovering order is more difficult because the
  /// ReachedMap should be readable in the dfs algorithm but the setting
  /// order map is not readable. Thus we must use the fork map:
  ///
  ///\code
  /// typedef Graph::NodeMap<int> OrderMap;
  /// OrderMap order(graph);
  /// typedef SettingOrderBoolMap<OrderMap> OrderSetterMap;
  /// OrderSetterMap setter(order);
  /// typedef Graph::NodeMap<bool> StoreMap;
  /// StoreMap store(graph);
  ///
  /// typedef ForkWriteMap<StoreMap, OrderSetterMap> ReachedMap;
  /// ReachedMap reached(store, setter);
  ///
  /// Dfs<Graph>::DefReachedMap<ReachedMap>::Create dfs(graph);
  /// dfs.reachedMap(reached);
  /// dfs.init();
  /// for (NodeIt it(graph); it != INVALID; ++it) {
  ///   if (!dfs.reached(it)) {
  ///     dfs.addSource(it);
  ///     dfs.start();
  ///   }
  /// }
  ///\endcode
  template <typename Map>
  class SettingOrderBoolMap {
  public:
    typedef typename Map::Key Key;
    typedef bool Value;

    /// Constructor
    SettingOrderBoolMap(Map& _map) 
      : map(_map), counter(0) {}

    /// Number of set operations.
    int num() const {
      return counter;
    }

    /// The \c set function of the map
    void set(const Key& key, Value value) {
      if (value) {
	map.set(key, counter++);
      }
    }
    
  private:
    Map& map;
    int counter;
  };

  /// @}
}

#endif // LEMON_MAPS_H

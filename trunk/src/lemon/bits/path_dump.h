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

#ifndef LEMON_BITS_PRED_MAP_PATH_H
#define LEMON_BITS_PRED_MAP_PATH_H

namespace lemon {

  template <typename _Graph, typename _PredMap>
  class PredMapPath {
  public:
    typedef True RevPathTag;

    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
    typedef _PredMap PredMap;

    PredMapPath(const Graph& _graph, const PredMap& _predMap,
                typename Graph::Node _target)
      : graph(_graph), predMap(_predMap), target(_target) {}

    int length() const {
      int len = 0;
      typename Graph::Node node = target;
      typename Graph::Edge edge;
      while ((edge = predMap[node]) != INVALID) {
        node = graph.source(edge);
        ++len;
      }
      return len;
    }

    bool empty() const {
      return predMap[target] != INVALID;
    }

    class RevEdgeIt {
    public:
      RevEdgeIt() {}
      RevEdgeIt(Invalid) : path(0), current(INVALID) {}
      RevEdgeIt(const PredMapPath& _path) 
        : path(&_path), current(_path.target) {
        if (path->predMap[current] == INVALID) current = INVALID;
      }

      operator const typename Graph::Edge() const {
        return path->predMap[current];
      }

      RevEdgeIt& operator++() {
        current = path->graph.source(path->predMap[current]);
        if (path->predMap[current] == INVALID) current = INVALID;
        return *this;
      }

      bool operator==(const RevEdgeIt& e) const { 
        return current == e.current; 
      }

      bool operator!=(const RevEdgeIt& e) const {
        return current != e.current; 
      }

      bool operator<(const RevEdgeIt& e) const { 
        return current < e.current; 
      }
      
    private:
      const PredMapPath* path;
      typename Graph::Node current;
    };

  private:
    const Graph& graph;
    const PredMap& predMap;
    typename Graph::Node target;
  };


  template <typename _Graph, typename _PredMatrixMap>
  class PredMatrixMapPath {
  public:
    typedef True RevPathTag;

    typedef _Graph Graph;
    typedef typename Graph::Edge Edge;
    typedef _PredMatrixMap PredMatrixMap;

    PredMatrixMapPath(const Graph& _graph, 
                      const PredMatrixMap& _predMatrixMap,
                      typename Graph::Node _source, 
                      typename Graph::Node _target)
      : graph(_graph), predMatrixMap(_predMatrixMap), 
        source(_source), target(_target) {}

    int length() const {
      int len = 0;
      typename Graph::Node node = target;
      typename Graph::Edge edge;
      while ((edge = predMatrixMap(source, node)) != INVALID) {
        node = graph.source(edge);
        ++len;
      }
      return len;
    }

    bool empty() const {
      return source != target;
    }

    class RevEdgeIt {
    public:
      RevEdgeIt() {}
      RevEdgeIt(Invalid) : path(0), current(INVALID) {}
      RevEdgeIt(const PredMatrixMapPath& _path) 
        : path(&_path), current(_path.target) {
        if (path->predMatrixMap(path->source, current) == INVALID) 
          current = INVALID;
      }

      operator const typename Graph::Edge() const {
        return path->predMatrixMap(path->source, current);
      }

      RevEdgeIt& operator++() {
        current = 
          path->graph.source(path->predMatrixMap(path->source, current));
        if (path->predMatrixMap(path->source, current) == INVALID) 
          current = INVALID;
        return *this;
      }

      bool operator==(const RevEdgeIt& e) const { 
        return current == e.current; 
      }

      bool operator!=(const RevEdgeIt& e) const {
        return current != e.current; 
      }

      bool operator<(const RevEdgeIt& e) const { 
        return current < e.current; 
      }
      
    private:
      const PredMatrixMapPath* path;
      typename Graph::Node current;
    };

  private:
    const Graph& graph;
    const PredMatrixMap& predMatrixMap;
    typename Graph::Node source;
    typename Graph::Node target;
  };

}

#endif

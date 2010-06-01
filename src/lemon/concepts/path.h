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
///\brief Classes for representing paths in graphs.
///
///\todo Iterators have obsolete style

#ifndef LEMON_CONCEPT_PATH_H
#define LEMON_CONCEPT_PATH_H

#include <lemon/bits/invalid.h>
#include <lemon/bits/utility.h>
#include <lemon/concept_check.h>

namespace lemon {
  namespace concepts {

    /// \addtogroup concept
    /// @{

    /// \brief A skeleton structure for representing directed paths in
    /// a graph.
    ///
    /// A skeleton structure for representing directed paths in a
    /// graph.  
    /// \param _Graph The graph type in which the path is.
    ///
    /// In a sense, the path can be treated as a list of edges. The
    /// lemon path type stores just this list. As a consequence it
    /// cannot enumerate the nodes in the path and the zero length
    /// paths cannot store the source.
    ///
    template <typename _Graph>
    class Path {
    public:

      /// Type of the underlying graph.
      typedef _Graph Graph;
      /// Edge type of the underlying graph.
      typedef typename Graph::Edge Edge;

      class EdgeIt;

      /// \brief Default constructor
      Path() {}

      /// \brief Template constructor
      template <typename CPath>
      Path(const CPath& cpath) {}

      /// \brief Template assigment
      template <typename CPath>
      Path& operator=(const CPath& cpath) {}

      /// Length of the path ie. the number of edges in the path.
      int length() const { return 0;}

      /// Returns whether the path is empty.
      bool empty() const { return true;}

      /// Resets the path to an empty path.
      void clear() {}

      /// \brief Lemon style iterator for path edges
      ///
      /// This class is used to iterate on the edges of the paths.
      class EdgeIt {
      public:
	/// Default constructor
	EdgeIt() {}
	/// Invalid constructor
	EdgeIt(Invalid) {}
	/// Constructor for first edge
	EdgeIt(const Path &) {}

        /// Conversion to Edge
	operator Edge() const { return INVALID; }

	/// Next edge
	EdgeIt& operator++() {return *this;}

	/// Comparison operator
	bool operator==(const EdgeIt&) const {return true;}
	/// Comparison operator
	bool operator!=(const EdgeIt&) const {return true;}
 	/// Comparison operator
 	bool operator<(const EdgeIt&) const {return false;}

      };

      template <typename _Path>
      struct Constraints {
        void constraints() {
          Path<Graph> pc;
          _Path p, pp(pc);
          int l = p.length();
          int e = p.empty();
          p.clear();

          p = pc;

          typename _Path::EdgeIt id, ii(INVALID), i(p);

          ++i;
          typename Graph::Edge ed = i;

          e = (i == ii);
          e = (i != ii);
          e = (i < ii);

          ignore_unused_variable_warning(l);
          ignore_unused_variable_warning(pp);
          ignore_unused_variable_warning(e);
          ignore_unused_variable_warning(id);
          ignore_unused_variable_warning(ii);
          ignore_unused_variable_warning(ed);
        }
      };

    };

    namespace _path_bits {
      
      template <typename _Graph, typename _Path, typename RevPathTag = void>
      struct PathDumperConstraints {
        void constraints() {
          int l = p.length();
          int e = p.empty();

          typename _Path::EdgeIt id, i(p);

          ++i;
          typename _Graph::Edge ed = i;

          e = (i == INVALID);
          e = (i != INVALID);

          ignore_unused_variable_warning(l);
          ignore_unused_variable_warning(e);
          ignore_unused_variable_warning(id);
          ignore_unused_variable_warning(ed);
        }
        _Path& p;
      };

      template <typename _Graph, typename _Path>
      struct PathDumperConstraints<
        _Graph, _Path, 
        typename enable_if<typename _Path::RevPathTag, void>::type
      > {
        void constraints() {
          int l = p.length();
          int e = p.empty();

          typename _Path::RevEdgeIt id, i(p);

          ++i;
          typename _Graph::Edge ed = i;

          e = (i == INVALID);
          e = (i != INVALID);

          ignore_unused_variable_warning(l);
          ignore_unused_variable_warning(e);
          ignore_unused_variable_warning(id);
          ignore_unused_variable_warning(ed);
        }
        _Path& p;
      };
    
    }


    /// \brief A skeleton structure for path dumpers.
    ///
    /// A skeleton structure for path dumpers. The path dumpers are
    /// the generalization of the paths. The path dumpers can
    /// enumerate the edges of the path wheter in forward or in
    /// backward order.  In most time these classes are not used
    /// directly rather it used to assign a dumped class to a real
    /// path type.
    ///
    /// The main purpose of this concept is that the shortest path
    /// algorithms can enumerate easily the edges in reverse order.
    /// If we would like to give back a real path from these
    /// algorithms then we should create a temporarly path object. In
    /// Lemon such algorithms gives back a path dumper what can
    /// assigned to a real path and the dumpers can be implemented as
    /// an adaptor class to the predecessor map.

    /// \param _Graph  The graph type in which the path is.
    ///
    /// The paths can be constructed from any path type by a
    /// template constructor or a template assignment operator.
    /// 
    template <typename _Graph>
    class PathDumper {
    public:

      /// Type of the underlying graph.
      typedef _Graph Graph;
      /// Edge type of the underlying graph.
      typedef typename Graph::Edge Edge;

      /// Length of the path ie. the number of edges in the path.
      int length() const { return 0;}

      /// Returns whether the path is empty.
      bool empty() const { return true;}

      /// \brief Forward or reverse dumping
      ///
      /// If the RevPathTag is defined and true then reverse dumping
      /// is provided in the path dumper. In this case instead of the
      /// EdgeIt the RevEdgeIt iterator should be implemented in the
      /// dumper.
      typedef False RevPathTag;

      /// \brief Lemon style iterator for path edges
      ///
      /// This class is used to iterate on the edges of the paths.
      class EdgeIt {
      public:
	/// Default constructor
	EdgeIt() {}
	/// Invalid constructor
	EdgeIt(Invalid) {}
	/// Constructor for first edge
	EdgeIt(const PathDumper&) {}

        /// Conversion to Edge
	operator Edge() const { return INVALID; }

	/// Next edge
	EdgeIt& operator++() {return *this;}

	/// Comparison operator
	bool operator==(const EdgeIt&) const {return true;}
	/// Comparison operator
	bool operator!=(const EdgeIt&) const {return true;}
 	/// Comparison operator
 	bool operator<(const EdgeIt&) const {return false;}

      };

      /// \brief Lemon style iterator for path edges
      ///
      /// This class is used to iterate on the edges of the paths in
      /// reverse direction.
      class RevEdgeIt {
      public:
	/// Default constructor
	RevEdgeIt() {}
	/// Invalid constructor
	RevEdgeIt(Invalid) {}
	/// Constructor for first edge
	RevEdgeIt(const PathDumper &) {}

        /// Conversion to Edge
	operator Edge() const { return INVALID; }

	/// Next edge
	RevEdgeIt& operator++() {return *this;}

	/// Comparison operator
	bool operator==(const RevEdgeIt&) const {return true;}
	/// Comparison operator
	bool operator!=(const RevEdgeIt&) const {return true;}
 	/// Comparison operator
 	bool operator<(const RevEdgeIt&) const {return false;}

      };

      template <typename _Path>
      struct Constraints {
        void constraints() {
          function_requires<_path_bits::
            PathDumperConstraints<Graph, _Path> >();
        }
      };

    };


    ///@}
  }

} // namespace lemon

#endif // LEMON_CONCEPT_PATH_H

//-------------------------------------------------------------------------------------
//  Copyright 2014 Michael Peeri
//
//  This file is part of hmmdsl.
//  hmmdsl is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  hmmdsl is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with hmmdsl.  If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------------------
#pragma once
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/arithmetic.hpp>
#include <boost/mpl/assert.hpp>
#include "v2.hpp"
#include "sequential_processing.hpp"


namespace v2
{

namespace tag
{
struct algorithm {};
} // namespace tag

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename ValueType, typename Args, typename Model>
class Algo
{

public:
  typedef Args arg_type;
  typedef double probability_t; // TODO : infer from model
  typedef ValueType result_type;
  typedef Model model_type;

protected:
    /*
     * Provide const access to the model
     * (for write-access, derive from affects_model)
     */
    template <typename Comp>
    const Model& model(Comp& comp) const
    {
        return comp.get_model();
    }

protected:
    /*
     * Provide const access to the model
     * (for write-access, derive from affects_model)
     */
    template <typename Comp>
    const Model& model_c(Comp& comp) const
    {
        return comp.get_model();
    }
};




////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Provide access to input data (at the algorithm's level)
 */
template <typename Level>
class applies_to_data
{
protected:
  template <typename Comp>
  typename Comp::template data<Level>::type
  data(Comp& comp)
  {
    return comp._get_data(Level());
  }

protected:
  template <typename Comp>
  const typename Comp::template data<Level>::type&
  data_c(Comp& comp) const
  {
    return comp._get_data_c(Level());
  }

public:
  typedef Level data_level;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Provide write access to the model
 */
template <typename Model>
class affects_model
{
protected:
  template <typename Comp>
  Model& model(Comp& comp)
  {
    return comp.get_model();
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace detail
{
struct has_dependency_list {};  // class 'has_dependency_list' is a base used to easily determine if a class derives from depends_on (since depends on is a template)
} // namespace detail

template <typename DependencyList, typename Level, typename SelfTag=void /* TODO: this shouldn't have a default */ /* TODO: add a registry metafunction to obtain this automatically? */>
class depends_on : public detail::has_dependency_list
{
public:
  typedef DependencyList dependencies;
  BOOST_MPL_ASSERT((
		    boost::is_same<
		        typename boost::mpl::sequence_tag<dependencies>::type
		      , boost::mpl::sequence_tag< boost::mpl::set0<>
		    >::type > ));

protected:
  // TODO: currently, an algorithm can simply bypass this by directly calling comp.apply(...)
  template <class Comp, class Tag>
  typename boost::enable_if<
  boost::mpl::and_<
      boost::mpl::or_<
	typename boost::mpl::has_key<dependencies, Tag>::type // Only allow calling declared dependencies
	, boost::is_same< Tag, SelfTag> // Or reflexive calls
      >
    , boost::is_same< // Only allow calling apply() within the same level of the data hierarchy
        typename Comp::template algo<Tag>::type::data_level
        , Level
      >
    >
  , typename Comp::template algo<Tag>::type::result_type
  >::type
  apply(Comp& comp, Tag, const typename Comp::template algo<Tag>::type::arg_type& args )
  {
    return comp._apply(Tag(), args); // make friend?
  }


protected:
  // TODO: currently, an algorithm can simply bypass this by directly calling comp.apply(...)
  template <class Comp, class Tag>
  typename boost::enable_if<
  boost::mpl::and_<
      boost::mpl::or_<
	typename boost::mpl::has_key<dependencies, Tag>::type // Only allow calling declared dependencies
	, boost::is_same< Tag, SelfTag> // Or reflexive calls
      >
      , boost::mpl::equal_to< // Only allow calling level down the data hierarchy
	  boost::mpl::plus<
              typename Comp::template algo<Tag>::type::data_level::value
	    , boost::mpl::int_<1>
          >
        , typename Level::value
        >
    >
  , typename Comp::template algo<Tag>::type::result_type
  >::type
  apply(Comp& comp, Tag, const typename Comp::template algo<Tag>::type::arg_type& args )
  {
    // Slow version - do checks
    return comp._apply(Tag(), args); // make friend?
  }

protected:
  template <class Comp, class LevelToMove>
  typename boost::enable_if<
  //      boost::mpl::and_<
	//        typename boost::mpl::has_key<dependencies, Tag>::type // Only allow calling declared dependencies
  //   , boost::mpl::equal_to< // Only allow calling apply_on_segment() one level down the data hierarchy
  //       boost::mpl::plus<
  //            typename Comp::template algo<Tag>::type::data_level::value // The level of the requested algo
  //          , boost::mpl::int_<1>
  //        >
  //      , typename Level::value  // The level of this algo
  //      >
      boost::mpl::not_<boost::is_same<
          LevelToMove
        , typename LevelToMove::info::max_level
        > >
  //    >
  ,  void
  >::type
  move_segment(Comp& comp, LevelToMove, typename Comp::segment_id id )
  {
    comp._move_segment(LevelToMove(), id); // make friend?
  }

};

/*
 * Utility meta-function to get a set of direct dependencies for an algorithm (even if it isn't derived from depends_on<> at all...)
 * Note: the 'has_dependency_list' base is used since depends_on is a template (with variadic argument)
 */
// Base template (never used)
template <class Algo, typename Enable=void>
struct direct_dependencies_of
{
};
// Partial specialization for algorithms derived from depends_on -- simply extract the dependecies set
template <typename Algo>
struct direct_dependencies_of<Algo, typename boost::enable_if< typename boost::is_base_of<detail::has_dependency_list, Algo> >::type >
{
  typedef typename Algo::dependencies type;
};
// Partial specialization for algorithms *not* derived from depends_on -- provide an empty set
template <typename Algo>
struct direct_dependencies_of<Algo, typename boost::disable_if< typename boost::is_base_of<detail::has_dependency_list, Algo> >::type >
{
  typedef boost::mpl::set0<> type;
};

/*
template <typename Algo>
struct direct_dependencies_of
{
  typedef boost::mpl::set<> type;
};
*/


 /*
  * Metafunction to select the suitable implementation for each algorithm.
  * Specialize over different attributes of the Model (possible using
  * enable_if) to provide the best implementation available.
  */
template <typename Tag, typename Model, typename Enable=void>
struct implementations
{
//    BOOST_MPL_ASSERT((boost::is_base_of<v2::tag::algorithm, Tag>));  // TODO - This check is not effective here. Can this be tested in a different way?

    struct missing_impl {};
    typedef missing_impl impl_type; // Default impl; TODO: find a way to remove this!
};


} // namespace v2

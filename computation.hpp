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
#include <boost/utility/enable_if.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/fusion/container/map.hpp>
#include <boost/fusion/view/filter_view.hpp>
#include "v2.hpp"
#include "read_fasta.hpp"

namespace v2
{



namespace detail
{

/*
 * Join two mpl::set<>s to obtain a new mpl::set<>
 */
struct join_sets
{
  template< typename Set1, typename Set2>
  struct apply
  {
    BOOST_MPL_ASSERT(( boost::mpl::is_sequence<Set1> ));
    BOOST_MPL_ASSERT(( boost::mpl::is_sequence<Set2> ));
  
  
    typedef typename boost::mpl::joint_view<Set1, Set2> joined_view;
    BOOST_MPL_ASSERT(( boost::mpl::is_sequence<joined_view > ));
  
  
    typedef typename boost::mpl::fold<
      joined_view
      , boost::mpl::set0<>
      , boost::mpl::insert<boost::mpl::_1, boost::mpl::_2>
      >::type type;
    BOOST_MPL_ASSERT(( boost::mpl::is_sequence<type> ));
  };
};
// TEST only
typedef join_sets::apply< boost::mpl::set<int,double>, boost::mpl::set<double, std::vector<int> > >::type x44;
BOOST_MPL_ASSERT(( boost::mpl::is_sequence<x44> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<x44>::value, ==, 3 );
// TEST only

//////////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the transitive closure of a set.
 * The closure is specified by a unary function  x -> set<x>
 * The container used throughout is mpl::set
 */
// Base template (never used)
template< typename Set, typename UnaryFunc, typename Enable=void>
struct closure
{
};

// Partial specialization 1: empty Set
template< typename Set, typename UnaryFunc>
struct closure<Set,UnaryFunc, typename boost::enable_if< typename boost::mpl::empty<Set> >::type >
{
  BOOST_MPL_ASSERT_RELATION( boost::mpl::size<Set>::value, ==, 0 );
  typedef boost::mpl::set0<> type;
};

// Partial specialization 2: Set contains a single value
template<
    typename Set
  , typename UnaryFunc>
struct closure<
    Set
  , UnaryFunc
  , typename boost::enable_if<
      boost::mpl::equal_to<
	boost::mpl::int_<1>
      , typename boost::mpl::size<Set>
      >
    >::type
  >
{
  BOOST_MPL_ASSERT_RELATION( boost::mpl::size<Set>::value, ==, 1 );
  typedef typename boost::mpl::front<Set>::type element;
  typedef typename boost::mpl::apply1<UnaryFunc, element>::type direct_results;
  typedef typename closure<direct_results, UnaryFunc>::type recursive_results;
  
  typedef typename join_sets::apply< Set, recursive_results>::type type;
};

// Partial specialization 3: Set contains more than a single value
template<
    typename Set
  , typename UnaryFunc>
struct closure<
  Set
  , UnaryFunc
  , typename boost::enable_if<
      boost::mpl::less<
	boost::mpl::int_<1>
      , typename boost::mpl::size<Set>
      >
    >::type
  >
{
  BOOST_MPL_ASSERT_RELATION( boost::mpl::size<Set>::value, >=, 1 );

  typedef typename boost::mpl::begin<Set>::type itfirst;
  typedef typename boost::mpl::deref<itfirst>::type first_element;

  typedef typename boost::mpl::erase_key< Set, first_element >::type rest_of_elements;
  enum { size_of_all  = boost::mpl::size<Set>::value};
  enum { size_of_rest = boost::mpl::size<rest_of_elements>::value };

  BOOST_MPL_ASSERT_RELATION( size_of_all-1, ==, size_of_rest ); // The 'rest' set contains one item less than Set

  typedef typename closure< boost::mpl::set1<first_element>, UnaryFunc >::type closure_of_first; // closure of the first item
  typedef typename closure< rest_of_elements, UnaryFunc >::type closure_of_rest; // recursive closure of the rest of the items
  
  typedef typename join_sets::apply< closure_of_first, closure_of_rest >::type type;
};




  /*
   * sub - get the type referenced by an array
   * TODO - change to use decltype (C++11), boost::result_of or something similar
   * http://www.boost.org/doc/libs/1_53_0/libs/utility/utility.htm#result_of
   */
  template <typename T>
  struct sub
  {
  };
  template <typename T>
  struct sub<T[]>
  {
    typedef T type;
  };
  template <typename T>
  struct sub<std::vector<T> >
  {
    typedef T type;
  };
  template <>
  struct sub<fasta::seq_cont_t>
  {
    typedef fasta::seq_cont_t::value_type type;
  };
  template <typename T, size_t N>
  struct sub<boost::multi_array<T,N> >
  {
    typedef typename boost::multi_array<T,N>::value_type type;
  };

} // namespace detail



//////////////////////////////////////////////////////////////////////////////////////////

struct EIncorrectAddress {};
/*
 * Rules for multi-level computation:
 * Each algorithm belongs to a level.
 * A level represents a set of non-overlapping parts of the data (e.g. partitions of the sequences).
 * An algorithm on level N may depend on algorithms on levels 0..N (but not on level N+1 or higher). ('depending' means requesting values).
 * The implementation of each algorithm may memoize values related to the current part (but not other parts).
 * So before moving to a different part on level N-1, an algorithm on level N must clear all lower levels.
 * To avoid repeated calculations, beforing moving to a different part on level N-1, the algorithm must make sure all work is done on the current part
 *
 * 0 [0]----][1]--------][2]------][3]-------][0]----][1]---------][2]-----]   scope: sequence       elements: chars      algos: e.g. delta_psi, viterbi_decoding
 * 1 [0]-------------------------------------][1]--------------------------]   scope: sequence-list  elements: seqs       algos: e.g. multiseq_viterbi
 * 2 [0]-------------------------------------------------------------------]   scope: seq-seq-list   elements: seq-list   algos: ?
 *
 */
template<typename Model, typename Algorithms, typename Data>
class Computation
{

  typedef Computation<Model,Algorithms,Data> self_type;

protected:
  /*
   * Use the implementations<> mechanism to determine the correct implementation class for a tag
   */
  template< typename Tag>
  struct determine_algo_impl_by_tag
  {
    typedef typename implementations<Tag,Model>::impl_type type;
  };

protected:
  // Helper metafunction for creating tag -> implementation pairs.
  // The implementation is selected using the implementations<> mechanism.
  struct create_impl_pairs
  {
    template< class Tag>
    struct apply
    {
      typedef typename boost::fusion::pair< Tag,
                                            typename determine_algo_impl_by_tag<Tag>::type > type;
    };
  };


public:
  /*
   * Helper metafunction to obtain an mpl::set containing the direct dependencies of an algorithm (specified by its tag)
   */
  template< class Tag>
  struct dependencies_of
  {
    typedef typename determine_algo_impl_by_tag<Tag>::type impl;

    typedef typename direct_dependencies_of<impl>::type direct_dependencies;
    BOOST_MPL_ASSERT((
		      boost::is_same<
		        typename boost::mpl::sequence_tag<direct_dependencies>::type
			, boost::mpl::sequence_tag< boost::mpl::set0<>
		      >::type > ));

    typedef direct_dependencies type;

    //typedef typename impl::dependencies direct_dependencies;
    typedef typename boost::mpl::size<direct_dependencies>::type depcount;
  };

  /*
   * MPL unary meta-function object facade for 'dependencies_of<>'
   */
private:
  struct direct_dependencies_set
  {
    template< class Tag>
    struct apply
    {
      typedef typename dependencies_of<Tag>::type type;
    };
  };

protected:

  // Calculate the dependency closure of all algorithms specifies explicitly
  typedef typename detail::closure<Algorithms, direct_dependencies_set>::type dependency_closure_of_algorithms;

  // Convert the dependencies set to a vector
  typedef typename boost::mpl::insert_range<
    boost::mpl::vector0<>
    , boost::mpl::end< boost::mpl::vector0<> >::type
    , dependency_closure_of_algorithms
    >::type impl_tags;

  // fusion::map storing tag -> implementation_class
  typedef typename boost::fusion::result_of::as_map<
    typename boost::mpl::transform<
        impl_tags
      , create_impl_pairs
    >::type
  >::type impls_t;

protected:
  // Store the actual instances of all contained algorithms
  impls_t _impls;


public:
  typedef size_t segment_id;

protected:
  const Data& _data;
  typedef std::vector<segment_id> data_pos_t;
  data_pos_t _data_pos;
  Model& _model;


  /*
   * Determine the type of an expression using operator[] 'Depth' times on 'Associative'
   * e.g.:
   * data_helper<std::vector<std::vector<int> >, 0> == std::vector<std::vector<int> >
   * data_helper<std::vector<std::vector<int> >, 1> == std::vector<int>
   * data_helper<std::vector<std::vector<int> >, 2> == int
   */
protected:
  template<typename Associative, size_t Depth>
  struct data_helper
  {
    /*
    typedef typename data_helper<
      typename boost::result_of<
	associative_access_functor(Associative)
      >::type
      , Depth-1
      >::result result;
    */
    typedef typename data_helper<
      typename detail::sub<Associative>::type
      , Depth-1
    >::result result;
  };
  template<typename Associative>
  struct data_helper<Associative,0>
  {
    typedef Associative result;
  };

public:
  template <class Level>
  struct data
  {
    typedef typename data_helper<
        Data
      , Level::info::max_level::value - Level::value::value
    >::result type;
  };



public:
  /*
   * Helper metafunction to allow the user access to types exposed by each contained algorithm.
   * e.g. Comp::algo<some_func>::type::return_type
   */
  template <class Tag>
  struct algo
  {
    typedef typename boost::fusion::result_of::value_at_key<impls_t,Tag>::type type;
  };

public:
  template <class Tag>
  struct level_of
  {
    typedef typename boost::fusion::result_of::value_at_key<impls_t,Tag>::type type;
  };
    
  /*
   * Provide access to one of the underlying algorithms (in the context of the
   * complete computation).
   * Apply the given arguments to the algorithm specified by the tag, and return the result.
   *
   * TODO: Change Tag from an actual argument (used to assist template resolution) to an
   *       explicit template argument (i.e. apply<Tag>(args) instead of apply(Tag(), args) ).
   * 
   */
public:  
  template <class Tag>
  typename boost::fusion::result_of::value_at_key<impls_t,Tag>::type::result_type
  _apply(Tag, const typename boost::fusion::result_of::value_at_key<impls_t,Tag>::type::arg_type& args )
  {
    //std::cout<< "apply<" << get_debug_id(Tag())<< ">("<< args<< ")"<< std::endl;

    //boost::fusion::at_key<Tag>(_impls).clear(*this); // resize (to new extents) and clear; move this!

    return boost::fusion::at_key<Tag>(_impls).apply_impl(*this, args);
  }

public:
  template <class Tag>
  const char* get_debug_id(Tag) const
  {
    return boost::fusion::at_key<Tag>(_impls).get_debug_id();
  }

private:
  /*
    Helper meta-function to compare the levels of Algo<>'s stored as fusion::mpl elements
   */
  template <class Level>
  struct get_level
  {
    template <typename T>
    struct apply
    {
      typedef typename boost::fusion::result_of::second<T>::type second;

      typedef boost::is_same<
	  Level
	, typename second::data_level
	> type;
    };
  };

  /*
   * Helper functor for clearing the algorithm impl. conatined in a fusion::map element
   */
private:
  struct clear_algo_storage
  {
    clear_algo_storage( self_type& comp ) : _comp(comp) {}

    template <class Item>
    void operator()(Item& item) const // TODO - refactor this; inelegant, and requires a unique template instance for each element
    {
#ifdef EXCESSIVE_DEBUG_PRINTING
      std::cout<< "Clearing algo '"<< item.second.get_debug_id()<< "'"<< std::endl;
#endif
      item.second.clear(_comp);
    }
  private:
    self_type& _comp;
  };

protected:
  template <class Level>
  typename boost::enable_if<
    boost::mpl::not_<boost::is_same< Level, sequential_processing::empty_scope> >
    , void
  >::type
  reset_storage_for_level()
  {
#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "Reseting storage for level "<< Level::value::value << std::endl;
#endif
    
    typedef boost::fusion::filter_view<
        impls_t
      , get_level<Level>
    > impls_at_level;
    BOOST_MPL_ASSERT(( boost::fusion::traits::is_view< impls_at_level > ));
    BOOST_MPL_ASSERT(( boost::fusion::traits::is_sequence<impls_at_level > ));

    impls_at_level ial(_impls);
    boost::fusion::for_each( impls_at_level(_impls)
                           , clear_algo_storage(*this) );

  }
protected:
  template <class Level>
  typename boost::enable_if<
  boost::is_same< Level, sequential_processing::empty_scope>
  , void
  >::type
  reset_storage_for_level()
  {
  }

public:  
  template <class Level>
  typename boost::enable_if<
  //boost::mpl::and_< 
      boost::is_base_of< // Level must be a sequential_processing_level
        sequential_processing::sequential_processing_level
      , Level
      >
      //    , boost::mpl::not_<boost::is_same< // does not apply to the outermost layer (because it has no index)
      //   Level
      //, typename Level::info::max_level
      //> >
  //>
  , void
  >::type
  _move_segment(Level, segment_id id )
  {
    // Release _data_pos stack when exiting C++ stack; TODO - use standard RAII mechanism
    //struct auto_release
    //{
    //  auto_release( data_pos_t& dp, segment_id id ) : _data_pos(dp) { _data_pos[ algo_data_level::value::value]= id; }
    //  ~auto_release() { /*_data_pos.pop_back();*/ }
    //private:
    //  data_pos_t& _data_pos;
    //} _auto_release(_data_pos, id);

#ifdef EXCESSIVE_DEBUG_PRINTING
    {
      std::cout<< "apply_on_segment(level: "<< Level::value::value<< ", segment: "<< id<< ", _data_pos: [";
      data_pos_t::const_iterator it, it_end;
      it = _data_pos.begin();
      it_end = _data_pos.end();
      for( ; it != it_end; ++it )
        {
          std::cout<< *it<< ", ";
        }
      std::cout<< "])"<< std::endl;
    }
#endif

    const size_t index_to_change = Level::info::max_level::value - Level::value::value;
    assert( index_to_change <= Level::info::max_level::value );

    // If the requested segment id is different than the existing one
    if( id != _data_pos[index_to_change] )
    { 
      // Clear storage for all algorithms at this level
      // TODO - what about lower levels? how to make sure each algorithm is cleared
      // when needed (but not more than once?)

      //reset_storage_for_level<
      //typename sequential_processing::get_previous_level<Level>::type
      //>();
#ifdef EXCESSIVE_DEBUG_PRINTING
      std::cout<< "Algorithm on level "<< Level::value::value<< " moved to segment "<< id<< ". clearing storage on this level."<< std::endl;
#endif

      // Set the index on the current level

#ifdef EXCESSIVE_DEBUG_PRINTING
      std::cout<< "(0)"<< std::flush;
#endif
      //_data_pos[algo_data_level::info::max_level::value - algo_data_level::value::value] = id;
      _data_pos[index_to_change] = id;
      
#ifdef EXCESSIVE_DEBUG_PRINTING
      std::cout<< "(1)"<< std::flush;
#endif
      
      std::fill(
                _data_pos.begin() + index_to_change + 1
              , _data_pos.end()
              , std::numeric_limits<segment_id>::max() );

#ifdef EXCESSIVE_DEBUG_PRINTING
      std::cout<< "(2)"<< std::flush;
#endif

#ifdef EXCESSIVE_DEBUG_PRINTING
      {
        std::cout<< "_data_pos (after): [";
        data_pos_t::const_iterator it, it_end;
        it = _data_pos.begin();
        it_end = _data_pos.end();
        for( ; it != it_end; ++it )
        {
          std::cout<< *it<< ", ";
        }
        std::cout<< "]"<< std::endl;
      }

      std::cout<< "(3)"<< std::flush;
#endif

      reset_storage_for_level<Level>();

#ifdef EXCESSIVE_DEBUG_PRINTING
      std::cout<< "(4)"<< std::flush;
#endif

    }

  }

protected:
  /*
   * helper functor to recursively use consecutive indexes
   * Pos - the index (in idx) of the next index (in u) to use
   * Depth - How many levels of indices remain (after the current one). For the final index, Depth==0
   */
  template <typename U, size_t Pos, size_t Depth>
  struct get_data_element_by_index
  {
    typedef typename data_helper<U,Depth>::result V;

    const V&
    operator()(const U& u, const std::vector<size_t>& idx) const
    {
      return get_data_element_by_index<
	typename data_helper<U,1>::result
	,Pos+1
	,Depth-1
      >()(u[idx[Pos]], idx);
    }
  };
  template<typename U,size_t Pos>
  struct get_data_element_by_index<U,Pos,0>
  {
    const U&
    operator()(const U& u, const std::vector<size_t>&) const
    {
      return u;
    }
  };
   
  
public:
  template <typename Level>
  typename data_helper<
    Data
  , boost::mpl::minus<
        boost::mpl::int_<Level::info::max_level::value>
      , boost::mpl::int_<Level::value::value>
    >::type::value
  >::result
  _get_data(const Level&)
  {
    std::vector<Level> v;

    // TEST ONLY
    /*
    typedef typename data_helper<int,0>::result type1;
    BOOST_MPL_ASSERT(( boost::is_same< type1, int > ));
    typedef typename data_helper<std::vector<int>, 1>::result type2;
    BOOST_MPL_ASSERT(( boost::is_same< type2, int > ));
    typedef typename data_helper<std::vector<std::vector<int> >, 1>::result type3;
    BOOST_MPL_ASSERT(( boost::is_same< type3, std::vector<int> > ));
    typedef typename data_helper<std::vector<std::vector<int> >, 2>::result type4;
    BOOST_MPL_ASSERT(( boost::is_same< type4, int > ));
    */
    // TEST ONLY


#ifdef EXCESSIVE_DEBUG_PRINTING
    const int i = int(Level::value::value);
    std::cout<< "get_data [level: "<< i<< "] pos: [";
    {
      data_pos_t::const_iterator it, it_end;
      it = _data_pos.begin();
      it_end = _data_pos.end();
      for( ; it != it_end; ++it )
        {
          std::cout<< *it<< ", ";
        }
    }
    std::cout<< "]"<< std::endl;
#endif

    //if( _data_pos.size() != Level::value )
    /*
    {
      std::cout<< "get_data<"<< Level::value<< "> called, current address: [";
      {
	data_pos_t::const_iterator it, it_end;
	it = _data_pos.begin();
	it_end = _data_pos.end();
	for( ; it != it_end; ++it )
	{
	  std::cout<< *it<< ", ";
	}
      }
      std::cout<< "]"<< std::endl;
      //throw EIncorrectAddress();
    }
    */

    return get_data_element_by_index<
         Data
      , 1 // Skip index 0 (which always equals 0)
      , Level::info::max_level::value //- Level::value::value // TEST THIS
      >()(_data, _data_pos);
  }


public:
  template <typename Level>
  const typename data_helper<
    Data
  , boost::mpl::minus<
        boost::mpl::int_<Level::info::max_level::value>
      , boost::mpl::int_<Level::value::value>
    >::type::value
  >::result&
  _get_data_c(const Level&) const
  {
    std::vector<Level> v;

    // TEST ONLY
    /*
    typedef typename data_helper<int,0>::result type1;
    BOOST_MPL_ASSERT(( boost::is_same< type1, int > ));
    typedef typename data_helper<std::vector<int>, 1>::result type2;
    BOOST_MPL_ASSERT(( boost::is_same< type2, int > ));
    typedef typename data_helper<std::vector<std::vector<int> >, 1>::result type3;
    BOOST_MPL_ASSERT(( boost::is_same< type3, std::vector<int> > ));
    typedef typename data_helper<std::vector<std::vector<int> >, 2>::result type4;
    BOOST_MPL_ASSERT(( boost::is_same< type4, int > ));
    */
    // TEST ONLY


    
#ifdef EXCESSIVE_DEBUG_PRINTING
    const int i = int(Level::value::value);
    std::cout<< "get_data <const> [level: "<< i<< "] pos: [";
    {
      data_pos_t::const_iterator it, it_end;
      it = _data_pos.begin();
      it_end = _data_pos.end();
      for( ; it != it_end; ++it )
        {
          std::cout<< *it<< ", ";
        }
    }
    std::cout<< "]"<< std::endl;
#endif    

    //if( _data_pos.size() != Level::value )
    /*
    {
      std::cout<< "get_data<"<< Level::value<< "> called, current address: [";
      {
	data_pos_t::const_iterator it, it_end;
	it = _data_pos.begin();
	it_end = _data_pos.end();
	for( ; it != it_end; ++it )
	{
	  std::cout<< *it<< ", ";
	}
      }
      std::cout<< "]"<< std::endl;
      //throw EIncorrectAddress();
    }
    */

    return get_data_element_by_index<
         Data
      , 1
      , Level::info::max_level::value - Level::value::value
    >()(_data, _data_pos);
  }


public:
  typedef typename boost::fusion::result_of::value_of_data<
   typename boost::fusion::result_of::begin<impls_t>::type
  >::type::data_level::info data_levels_info;

public:
  Computation( const Data& data, Model& model )
    : _data(data)
    , _data_pos(data_levels_info::max_level::value + 1, std::numeric_limits<segment_id>::max())
    , _model(model) 
  {
  }

public:
  const Model& get_model() const
  {
    return _model;
  }

public:
  Model& get_model()
  {
    return _model;
  }


public:

  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  //
  /*
  typedef boost::mpl::set<v2::tag::run_all> set766;
  typedef typename detail::closure<set766, direct_dependencies_set>::type closure_of_set766;
  BOOST_MPL_ASSERT((
		    boost::is_same<
		      typename boost::mpl::sequence_tag<closure_of_set766>::type
		      , boost::mpl::sequence_tag< boost::mpl::set0<>
						  >::type > ));
  BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set766>::value, ==, 4 );

  typedef boost::mpl::set<v2::tag::run_all, v2::tag::delta_psi> set767;
  typedef typename detail::closure<set767, direct_dependencies_set>::type closure_of_set767;
  BOOST_MPL_ASSERT((
		    boost::is_same<
		      typename boost::mpl::sequence_tag<closure_of_set767>::type
		      , boost::mpl::sequence_tag< boost::mpl::set0<>
						  >::type > ));
  BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set767>::value, ==, 4 );


  typedef boost::mpl::set<v2::tag::delta_psi> set768;
  typedef typename detail::closure<set768, direct_dependencies_set>::type closure_of_set768;
  BOOST_MPL_ASSERT((
		    boost::is_same<
		      typename boost::mpl::sequence_tag<closure_of_set768>::type
		      , boost::mpl::sequence_tag< boost::mpl::set0<>
						  >::type > ));
  BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set768>::value, ==, 1 );

  typedef boost::mpl::set<v2::tag::viterbi_decoding> set769;
  typedef typename detail::closure<set769, direct_dependencies_set>::type closure_of_set769;
  BOOST_MPL_ASSERT((
		    boost::is_same<
		      typename boost::mpl::sequence_tag<closure_of_set769>::type
		      , boost::mpl::sequence_tag< boost::mpl::set0<>
						  >::type > ));
  BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set769>::value, ==, 2 );

  typedef boost::mpl::set<v2::tag::multiseq_viterbi> set770;
  typedef typename detail::closure<set770, direct_dependencies_set>::type closure_of_set770;
  BOOST_MPL_ASSERT((
		    boost::is_same<
		      typename boost::mpl::sequence_tag<closure_of_set770>::type
		      , boost::mpl::sequence_tag< boost::mpl::set0<>
						  >::type > ));
  BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set770>::value, ==, 3 );
  */
  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  //

  friend struct detail::has_dependency_list;

};


} // namespace v2


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
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#undef BOOST_DISABLE_ASSERTS
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <boost/tuple/tuple.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/front.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/concept_check.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/utility/result_of.hpp>
#include "common.hpp"
#include "v2.hpp"
#include "algo.hpp"
#include "sequential_processing.hpp"
#include "memoized.hpp"
#include "computation.hpp"
#include "read_fasta.hpp"


typedef size_t prev_ptr_t;

template<class Algo>
class ViterbiAlgorithm
{
  // Notation:
  // _v[l][i]  (Durbin: v_l(i) )  probability of the most probable path ending in state l with observation i

protected:
  typedef typename Algo::probability_type P;
  typedef boost::tuple<P, prev_ptr_t> extval_t;
  typedef std::vector<extval_t> alternatives_t;
  enum {_empty=1        };

public:
  P vli( size_t l, size_t i )
  {
    using boost::get;

    // Is the value already cached?
    const P val = _v[l][i];
    if( val != _empty )
      return val;

    extval_t newext;
    //if( i == 0 )
    //{
    //  newext = boost::make_tuple( log(1. / num_states(*_model) ),
    //                    std::numeric_limits<prev_ptr_t>::max() );
    //}
    //else
    //{
      newext = _produce_vli( l, i );
      //}

    // Value not cached - need to calculate it
    _v[l][i]   = get<0>(newext);
    _ptr[l][i] = get<1>(newext);
    return get<0>(newext);
  }
protected:

  /*
    v(l,i) = The probability of reaching state l at time i,
             while taking the most likely path.
   */
  extval_t _produce_vli( size_t l, size_t i )
  {
    using boost::get; using boost::make_tuple;

    if( i==0 )
      {
        return make_tuple( (l==0)?0:-std::numeric_limits<typename Algo::probability_type>::max() , std::numeric_limits<prev_ptr_t>::max() );
      }
    // Note: The constraint that the terminal state must be occupied at time L is not enforced here

    alternatives_t alternatives;
    for( size_t k = 0; k< num_states(*_model); ++k )
    {
      const P vk_prev = vli(k, i-1);
      const P p = vk_prev + _model->a(k, l);
      alternatives.push_back( make_tuple( p, k ) );
    }

    extval_t maxk = getmax( alternatives );
    if( i < length(*_seq) )
      return make_tuple( get<0>(maxk) + _model->e(l, (*_seq)[i] ),
                         get<1>(maxk) );
    else
      return maxk;
  }

  static extval_t getmax( const alternatives_t& alts )
  {
    using boost::get;
    assert( !alts.empty() );
    typename alternatives_t::const_iterator it, it_best;

    if( alts.size() == 1 )
      return alts[0];

    it_best = it = alts.begin(); 
    P Pbest = get<0>(*it);

    for( ++it; it != alts.end(); ++it )
    {
      if( get<0>(*it) > Pbest )
      {
        Pbest = get<0>(*it);
        it_best = it;
      }
    }
    return *it_best;
  }



public:
  P calc()
  {
    using boost::get;
    alternatives_t alternatives;
    //for( size_t k = 0; k< num_states(*_model); ++k )
    //{
    // const P p = vli(k, length(*_seq));
    //  alternatives.push_back( boost::make_tuple( p, k ) );
    // }
    //extval_t maxk = getmax( alternatives );
    const P p = vli(1, length(*_seq));
    
    return p;
  }

public:
  void get_backtrace(std::vector<std::string>& bt)
  {
    using boost::get;
    alternatives_t alternatives;
    
        // Decide what should be the last state
    const size_t last = length(*_seq);
    //for( size_t k = 0; k< num_states(*_model); ++k )
    //{
    //  const P p = vli(k, last);
    //  alternatives.push_back( boost::make_tuple( p, k ) );
    //}
    //extval_t maxk = getmax( alternatives );
    

    // Backtrack from the last state
    int pos = last;
    size_t state = 1;//maxk.get<1>();
    while(pos>=0)
    {
      bt.push_back( _model->GetStateName(state) );
      state = _ptr[state][(size_t)pos];
      --pos;
    }

    std::reverse( bt.begin(), bt.end() );
  }


public:
  ViterbiAlgorithm( boost::shared_ptr<typename Algo::sequence_type> seq, boost::shared_ptr<typename Algo::model_type> model )
    :   _v(boost::extents[num_states(*model)][length(*seq)+1])
    , _ptr(boost::extents[num_states(*model)][length(*seq)+1])
    , _model(model)
    , _seq(seq)
  {
    // TODO Make this work
    //    for( typename v_type::iterator it = _v.begin(); it != _v.end(); ++it )
    //    *it = 0.0;
    const size_t m_end = _v.shape()[0];
    const size_t n_end = _v.shape()[1];
    
    for( size_t m=0; m< m_end; ++m )
        for( size_t n=0; n< n_end; ++n )
        {
          _v[m][n] = _empty;
          _ptr[m][n] = std::numeric_limits<prev_ptr_t>::max();
        }
    


    //std::fill(_v.begin(), _v.end(), P(_empty) );
    //std::fill(_ptr.begin(), _ptr.end(), prev_ptr_t(_empty) );
  }
protected:
  
  typedef boost::multi_array<P,2> v_type;
  v_type _v;
  boost::multi_array<prev_ptr_t,2> _ptr;

  const boost::shared_ptr<typename Algo::model_type> _model;
  const boost::shared_ptr<typename Algo::sequence_type> _seq;

public:
  void debug_print() const
  {
    std::cout<< "ViterbiAlgorithm<Algo>"<< std::endl;
    std::cout<< " v: ";
    ::debug_print(_v);

    std::cout<< " ptr: ";
    ::debug_print(_ptr);
  }

};

namespace v2
{



namespace tag {
  struct viterbi_decoding {};
}


//////////////////////////////////////////////////////////////////////////////////////////

namespace tag {
  struct delta_psi {};
}

namespace detail
{

/*
 * Helper functor for max_elements: compare the first elements of a pair of fusion::vector<>s
 */
struct compare_first_elements
{
  template<typename T>
  bool operator()(const T& a, const T& b) const { return boost::fusion::at_c<0>(a) < boost::fusion::at_c<0>(b); }
};


template<typename Model>
class delta_psi_hsmm :
    public Algo<
      boost::fusion::vector<probability_t, state_t, time_t>
    , boost::fusion::vector<state_t, time_t>
    , Model
  >
  , public memoize<
      delta_psi_hsmm<Model>
    , boost::fusion::vector<state_t, time_t>
    , boost::fusion::vector<probability_t, state_t, time_t>
    , sequential_processing::single_sequence_scope
  >
  , public applies_to_data<
      sequential_processing::single_sequence_scope
  >
  , public depends_on<
      boost::mpl::set<>  // delta_psi only depends on itself; currently, this means we need to derive from depends_on<>, but the dependecies set may be empty
    , sequential_processing::single_sequence_scope
    , tag::delta_psi
  >
{
  BOOST_MPL_ASSERT((boost::is_same<typename model::model_traits<Model>::hidden_process, model::semi_markov_chain>));

public:
  typedef Algo<
      boost::fusion::vector<probability_t, state_t, time_t>
    , boost::fusion::vector<state_t, time_t>
    , Model
  > base;

  typedef typename base::result_type result_type;
  typedef typename base::arg_type arg_type;
  typedef tag::delta_psi tag;

  // Define types for storing alternatives for maximization

  // First, for the inner (j) loop
  typedef std::vector<result_type> alternatives_j_t;

  // Second, for the outer (tau) loop
  //typedef typename boost::fusion::result_of::as_vector<
  //  typename boost::fusion::result_of::push_back<
  //    result_type, time_t
  //  >::type
  //>::type p_tau_j_t; // Define an element type appending an additional time_t at the end of result_type
  // TODO REMOVE THIS
  typedef result_type p_tau_j_t;

  typedef std::vector<p_tau_j_t> alternatives_tau_j_t;

public:
  /*
    Delta -- The probability of the best path in which state l ends at position i
    Psi -- The state last state occupied before position i on the best path
    Tau -- The duration occupied by state l when ending at position i on the best path

    Reference for Viterbi algorithm for HSMMs:
    Datta, Hu, Ray, "On Efficient Viterbi Decoding for Hidden semi-Markov Models"

    @inproceedings{datta2008efficient,
      title={On efficient Viterbi decoding for hidden semi-Markov models},
      author={Datta, Ritendra and Hu, Jianying and Ray, Bonnie},
      booktitle={Pattern Recognition, 2008. ICPR 2008. 19th International Conference on},
      pages={1--4},
      year={2008},
      organization={IEEE}
    }
   */
  template<typename Comp>
  result_type produce(Comp& comp, const arg_type& args)
  {
    typedef typename Comp::template algo<v2::tag::delta_psi>::type::result_type inferred_result_type;
    BOOST_MPL_ASSERT(( boost::is_same< inferred_result_type, result_type > ));

    typedef typename boost::fusion::result_of::size<result_type>::type size_of_delta_psi_tau_t;
    BOOST_MPL_ASSERT_RELATION( size_of_delta_psi_tau_t::value, ==, 3 );


    // DEBUG ONLY
    // std::cout<< "delta_psi_hsmm*(";
    // boost::fusion::for_each<>( args, print_with_comma() );
    // std::cout<< ")"<< std::endl;
    // DEBUG ONLY

    const state_t initial_state = begin_state(this->model(comp));
    const sequence_t& seq = fasta::get_seq(data_c(comp));
    //std::cout<< "*"<< std::flush;

    const state_t l = boost::fusion::at_c<0>(args);
    const time_t  i = boost::fusion::at_c<1>(args);

    // Termination: 
    if( i==0 )
    {
      // At position 0, only the initial state is allowed to end
      return result_type(  (l == initial_state) ? 0.0 : -std::numeric_limits<probability_t>::max()
                           , 0
                           , 0 );
    }
    // Note: The constraint that the terminal state must be occupied at time L is not enforced here

    const size_t tau_max = this->model(comp).get_max_duration();

    /// Pre-calculate all products of emission probabilities for the range [i-tau, i)
    const size_t max_s = std::min( i , tau_max );
    std::vector<probability_t> pi_bl( max_s + 1, -std::numeric_limits<probability_t>::max() );

    // Calculate pi_bl = \lambda \tau -> e(l, seq(i)) * e(l, seq(i-1)) * e(l, seq(i-2)) * ... * e(l, seq(i-tau+1))
    pi_bl[0] = -std::numeric_limits<probability_t>::max(); // invalid value - used just so pi_bl indices correspond with tau values
    if( i<seq.length() )
      pi_bl[1] = this->model(comp).e( l, seq[ i-1 ] );
    else
      pi_bl[1] = 0;
    for( size_t tau=2; tau<= max_s; ++tau )
    {
      //probability_t bl = 0;
      //if(( i-tau ) < seq.length() )
	probability_t bl = this->model(comp).e( l, seq[i-tau] );
	//else
	//	std::cout<< "*()*"<< std::flush;

      pi_bl[tau] = pi_bl[tau-1] + bl;
    }
    //std::cout<< "(b)"<< std::flush;
    //std::cout<< std::endl;

#ifndef NO_TESTS
    { // TEST 

      // Reference (i.e. trivial) implementation of pi_bl calculation - for testing correctness
      for( size_t tau=1; tau <= max_s; ++tau )
      {
        probability_t p = 0.0;
        size_t num_factors = 0; // verify the number of sequence elements summed equals tau

        for( size_t u=i-tau+1; u<=i; ++u, ++num_factors )
        {
	  if(u < seq.length() )
	    p += this->model(comp).e(l, seq[u-1] );
        }
        assert( num_factors == tau );

        if( fabs( p - pi_bl[tau] ) > 1e-10 )
        {
          std::cout<< "(X) tau=\t"<< tau<<  "\tpi_bl[s]=\t"<< pi_bl[tau]<< "\tp_ref(tau)=\t"<< p<< std::endl;
        }
      }
    }
#endif

    alternatives_tau_j_t alternatives_tau_j; alternatives_tau_j.reserve(tau_max+1);
    //--- max tau loop ---------
    for( size_t tau=1; tau <= max_s ; ++tau )
    {
      assert( i-tau >= 0 );
      
      alternatives_j_t alternatives_j; alternatives_j.reserve( num_states(this->model(comp)) );
      //--- max j loop ---------
      for( size_t j = 0; j< num_states(this->model(comp)); ++j )
      {
        const result_type vj_prev = this->apply(comp
                                              , v2::tag::delta_psi()
                                              , typename Comp::template algo<v2::tag::delta_psi>::type::arg_type(j, i-tau) );

        const probability_t p = boost::fusion::at_c<0>(vj_prev) + this->model(comp).a(j, l);

        // Store this alternative
        alternatives_j.push_back( result_type( p, j, std::numeric_limits<time_t>::max() ) );
      }

      // Get the best alternative (for j's) for this tau
      typename alternatives_j_t::const_iterator maxj = std::max_element( alternatives_j.begin(), alternatives_j.end(), compare_first_elements() );
    
      const probability_t p_tau = this->model(comp).p(l, tau);

      const probability_t this_pi_bl = pi_bl[tau]; // Use one of the pre-calculated values of pi_bl

      // Store the best alternative for this tau
      alternatives_tau_j.push_back(
                                   p_tau_j_t( 
                                             boost::fusion::at_c<0>(*maxj) + p_tau + this_pi_bl
                                           , boost::fusion::at_c<1>(*maxj)
                                           , tau ) );
    }
    // Get the best alternative (for any tau and j)
    typename alternatives_tau_j_t::const_iterator max_tau_j = std::max_element( alternatives_tau_j.begin(), alternatives_tau_j.end(), compare_first_elements() );

    return *max_tau_j;
      ////if( i < length(seq) )
      //  return result_type( boost::fusion::at_c<0>(*max_tau_j)
      //                  ,  boost::fusion::at_c<1>(*max_tau_j) );
      //  //else
      // // return result_type( -std::numeric_limits<probability_t>::max(), 0 );
  }

public:
  template<typename Comp>
  boost::fusion::vector<size_t,size_t> get_extents(Comp& comp)
  {
#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "get_extents "<< fasta::get_len(data_c(comp))+1<< "  "<< num_states( this->model(comp) )<< std::endl;
#endif
    return boost::fusion::vector<state_t,time_t>(
                                                 num_states( this->model(comp) )
                                               , fasta::get_len(data_c(comp))+1
                                                );
  }

public:
  static result_type get_unassigned_value() { return result_type(0.0 /* Assumption: probability_t(0.0) can safely be compared exactly (without requiring a nonzero tolerance). */ , std::numeric_limits<state_t>::max(), std::numeric_limits<time_t>::max() ); }

public:
  const char* get_debug_id() const { return "delta_psi[hsmm]"; }

public:
  void reset(size_t level=0) {};

};

//////////////////////////////////////////////////////////////////////////////////////////

template<typename Model>
class delta_psi_hmm :
    public Algo<
      boost::fusion::vector<probability_t,state_t>
    , boost::fusion::vector<state_t, time_t>
    , Model
  >
  , public memoize<
      delta_psi_hmm<Model>
    , boost::fusion::vector<state_t, time_t>
    , boost::fusion::vector<probability_t,state_t>
    , sequential_processing::single_sequence_scope
  >
  , public applies_to_data<
      sequential_processing::single_sequence_scope
  >
  , public depends_on<
      boost::mpl::set<>  // delta_psi only depends on itself; currently, this means we need to derive from depends_on<>, but the dependecies set may be empty
    , sequential_processing::single_sequence_scope
    , tag::delta_psi
  >
{
  BOOST_MPL_ASSERT((boost::is_same<typename model::model_traits<Model>::hidden_process, model::markov_chain>));

public:
  typedef Algo<
      boost::fusion::vector<probability_t, state_t>
    , boost::fusion::vector<state_t, time_t>
    , Model
  > base;

  typedef typename base::result_type result_type;
  typedef typename base::arg_type arg_type;
  typedef tag::delta_psi tag;

  typedef std::vector<result_type> alternatives_t;

public:
  template<typename Comp>
  result_type produce(Comp& comp, const arg_type& args)
  {
#ifdef EXCESSIVE_DEBUG_PRINTING
    //std::cout<< "delta_psi_hmm(";
    boost::fusion::for_each<>( args, print_with_comma() );
    //std::cout<< ")"<< std::endl;
#endif

    const sequence_t& seq = fasta::get_seq(data_c(comp));
    //std::cout<< "*";

    const state_t l = boost::fusion::at_c<0>(args);
    const time_t  i = boost::fusion::at_c<1>(args);

    if( i==0 )
    {
      return result_type(  (l==0) ? 0 : -std::numeric_limits<probability_t>::max()
                           , 0 );
    }
    // Note: The constraint that the terminal state must be occupied at time L is not enforced here

    alternatives_t alternatives; alternatives.reserve(num_states(this->model(comp)));

    for( size_t k = 0; k< num_states(this->model(comp)); ++k )
    {
      const result_type vk_prev = this->apply(comp, v2::tag::delta_psi(), typename Comp::template algo<v2::tag::delta_psi>::type::arg_type(k, i-1) );
      const probability_t p = boost::fusion::at_c<0>(vk_prev) + this->model(comp).a(k, l);
      alternatives.push_back( result_type( p, k ) );
    }

    typename alternatives_t::const_iterator maxk = std::max_element( alternatives.begin(), alternatives.end(), compare_first_elements() );

    if( i < length(seq) )
      return result_type( boost::fusion::at_c<0>(*maxk) + this->model(comp).e( l, seq[i] ),
                          boost::fusion::at_c<1>(*maxk) );
    else
      return *maxk;
  }

public:
  template<typename Comp>
  boost::fusion::vector<size_t,size_t> get_extents(Comp& comp)
  {
    //std::cout<< "get_extents "<< length(data(comp))+1<< "  "<< num_states( this->model(comp) )<< std::endl;
    return boost::fusion::vector<size_t,size_t>(
                                                  num_states( this->model(comp) )
                                                , fasta::get_len(data(comp))+1
                                                );
  }

public:
  static result_type get_unassigned_value() { return result_type(0.0, std::numeric_limits<state_t>::max() ); }

public:
  const char* get_debug_id() const { return "delta_psi[hmm]"; }

public:
  void reset(size_t level=0) {};
};


} // namespace detail 


//////////////////////////////////////////////////////////////////////////////////////////


template <class Model>
struct implementations<
    tag::delta_psi
  , Model
  , typename boost::enable_if<
      typename boost::is_same<
          typename model::model_traits<Model>::hidden_process
        , model::semi_markov_chain
        >::type
    >::type
  >
{
  //  typedef detail::delta_psi<Model> impl_type;
  typedef typename detail::delta_psi_hsmm<Model> impl_type;
};


template <class Model>
struct implementations<
    tag::delta_psi
  , Model
  , typename boost::enable_if<
      typename boost::is_same<
          typename model::model_traits<Model>::hidden_process
        , model::markov_chain
        >::type
    >::type
  >
{
  typedef typename detail::delta_psi_hmm<Model> impl_type;
};
//////////////////////////////////////////////////////////////////////////////////////////


/*
namespace tag {
  struct viterbi_probability {};
}
*/


namespace detail
{

  // TODO - Fix the mess of repeated arguments (required in Barton-Nackman because of incomplete instantiation)
  /*
template<typename Model>
class viterbi_probability :
  public Algo<
      probability_t
    , boost::fusion::vector<>
    , Model
  >
  , public depends_on<
      boost::mpl::set< tag::delta_psi>
    , sequential_processing::single_sequence_scope
  >
  , public memoize<
      viterbi_probability<Model>
    , boost::fusion::vector<>
    , probability_t
    , sequential_processing::single_sequence_scope
  >
{
public:
  typedef Algo<
      probability_t
    , boost::fusion::vector<>
    , Model
  > base;
  typedef tag::viterbi_probability tag;

  typedef typename base::arg_type arg_type;
  typedef typename base::result_type result_type;
  typedef typename sequential_processing::single_sequence_scope data_level;
  
public:
  template <typename Comp>
  result_type produce(Comp& comp, const arg_type& args)
  {
#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "viterbi_probability ONE (";
    boost::fusion::for_each<>( args, print_with_comma() );
    std::cout<< ")"<< std::endl;
#endif

#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "...   calling delta_psi:"<< std::endl;
#endif
    
apply(comp, v2::tag::delta_psi(), typename Comp::template algo<v2::tag::delta_psi>::type::arg_type(2,3) ); // NOTE: Code not in use!
    assert(false);
    throw; // NOTE: Code not in use!

#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "...   done!"<< std::endl;
#endif

    return result_type();
  }

public:
  template<typename Comp>
  boost::fusion::vector<> get_extents(Comp& comp)
  {
    return boost::fusion::vector<>();
  }

public:
  static result_type get_unassigned_value() { return result_type(-1); }

public:
  const char* get_debug_id() const { return "viterbi_probability ONE"; }


};
  */

} //namespace detail

/*
template <class Model>
struct implementations<tag::viterbi_probability, Model >
{
  typedef typename detail::viterbi_probability<Model> impl_type;
};
*/


//////////////////////////////////////////////////////////////////////////////////////////
struct NullDecoding {};

struct decoding_t
{
public:
  decoding_t() : _uninitialized(false), _logp( -std::numeric_limits<double>::max() ) {}
  decoding_t( const NullDecoding& ) : _uninitialized(true), _logp( -std::numeric_limits<double>::max() ) {}

public:
  typedef std::vector<size_t> positions_t;

public: // TODO: consider making this a real class...
  positions_t _p1,_p2;
  sequence_t _s1, _s2;
  std::string _id1, _id2;
  bool _uninitialized;
  double _logp;

public:
  bool operator==(const decoding_t& other) const
  {
    if( _uninitialized ) return other._uninitialized;
    // TODO Optimize this
    // TODO Verify this
    if( _p1 != other._p1 )
      return false;
    if( _p2 != other._p2 )
      return false;
    if( _s1 != other._s1 )
      return false;
    if( _s2 != other._s2 )
      return false;
    if( fabs( _logp - other._logp) > 1e-6 )
      return false;
    return true;
  }
public:
  bool operator!=(const decoding_t& other) const
  {
    return !(this->operator==(other));
  }
public:
  decoding_t& operator=(const decoding_t& other)
  {
    _p1 = other._p1;
    _p2 = other._p2;
    _s1 = other._s1;
    _s2 = other._s2;
    _id1 = other._id1;
    _id2 = other._id2;
    _logp = other._logp;
    _uninitialized = other._uninitialized;
    return *this;
  }
  friend std::ostream& operator<< (std::ostream& stream, const decoding_t& decoding);
};


std::ostream& operator<< (std::ostream& stream, const decoding_t& decoding);

//////////////////////////////////////////////////////////////////////////////////////////

namespace tag {
  struct multiseq_viterbi {};
}

typedef size_t seq_id_t;

namespace detail
{

template<typename Model>
class multiseq_viterbi :
  public Algo<
      decoding_t
    , boost::fusion::vector<seq_id_t>
    , Model
  >
  , public depends_on<
      boost::mpl::set< tag::viterbi_decoding>
    , sequential_processing::multi_sequence_scope
  >
  , public memoize<
      multiseq_viterbi<Model>
    , boost::fusion::vector<seq_id_t>
    , decoding_t
    , sequential_processing::multi_sequence_scope
  >
  , public applies_to_data<
      sequential_processing::multi_sequence_scope
  >
{
public:
  typedef Algo<
      decoding_t
    , boost::fusion::vector<seq_id_t>
    , Model
  > base;
  typedef typename base::result_type result_type;
  typedef typename base::arg_type arg_type;
  typedef typename sequential_processing::multi_sequence_scope data_level;
  
public:
  template <typename Comp>
  result_type produce(Comp& comp, const arg_type& arg)
  {
#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "multiseq_viterbi ("<< arg<< ")"<< std::endl;
#endif
    
    const size_t i = boost::fusion::at_c<0>(arg);

    /*
    typedef const typename Comp::template data<data_level>::type data_t;
    data_t d = data(comp);
    std::cout<< "   data: "<< d.size()<< " seqs"<< std::endl;
    assert( i < d.size() );
    */

    //typename data_t::const_iterator it,it_end;
    //it = d.begin();
    //it_end = d.end();
    //size_t i = 0;
    //for( ; it != it_end; ++it, ++i )
    //{

    move_segment( comp,
                  typename Comp::template algo<v2::tag::viterbi_decoding>::type::data_level(),
                  i );

#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "     multiseq_viterbi: calling viterbi_decoding()..."<< std::endl;
#endif
    result_type ret = 
      apply( comp,
             v2::tag::viterbi_decoding(),
             typename Comp::template algo<v2::tag::viterbi_decoding>::type::arg_type() );
    //}
#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "     multiseq_viterbi: done (got result: "<< ret<< ")"<< std::endl;
#endif

    return ret;
  }

public:
  template<typename Comp>
  boost::fusion::vector<size_t> get_extents(Comp& comp)
  {
    BOOST_MPL_ASSERT((boost::is_same<typename Comp::template data<data_level>::type, fasta::seq_cont_t >));

    // DEBUG ONLY
#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "multiseq_viterbi::get_extents()"<< std::endl;
    const fasta::seq_cont_t& d = data_c(comp);
    const size_t len = fasta::num_sequences(d);
    std::cout<< len<< std::endl;
#endif
    // DEBUG ONLY

    return boost::fusion::vector<size_t>(fasta::num_sequences(data_c(comp)));
  }

public:
  static result_type get_unassigned_value() { return result_type(NullDecoding()); }

public:
  const char* get_debug_id() const { return "multiseq_viterbi"; }

};

} // namespace detail

template <class Model>
struct implementations<tag::multiseq_viterbi, Model >
{
  typedef typename detail::multiseq_viterbi<Model> impl_type;
};


//////////////////////////////////////////////////////////////////////////////////////////

namespace tag {
  struct run_all {};
}

namespace detail
{

template<typename Model>
class run_all :
  public Algo<
      boost::fusion::vector<bool>
    , boost::fusion::vector<>
    , Model
  >
  , public depends_on<
      boost::mpl::set< tag::multiseq_viterbi>
    , sequential_processing::everything
  >
  , public memoize<
      run_all<Model>
    , boost::fusion::vector<>
    , boost::fusion::vector<bool>
    , sequential_processing::everything
  >
  , public applies_to_data<
      sequential_processing::everything
  >
{
public:
  typedef Algo<
      boost::fusion::vector<bool>
    , boost::fusion::vector<>
    , Model
  > base;
  typedef typename base::result_type result_type;
  typedef typename base::arg_type arg_type;
  typedef typename sequential_processing::everything data_level;  

public:
  template <typename Comp>
  result_type produce(Comp& comp, const arg_type&)
  {
#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "run_all ONE (<ignoring args>)"<< std::endl;
#endif

    /*
    const typename Comp::template data<data_level>::type d = data_c(comp);
    std::cout<< "   data: "<< d.size() << " items"<< std::endl;
    assert( d.size() == 1 );

    for(size_t i = 0; i < d[0].size(); ++i )
    {
      std::cout<< "      run_all: calling multiseq("<< i<< ")..."<< std::endl;
      apply_on_segment(comp,
                       v2::tag::multiseq_viterbi(),
                       typename Comp::template algo<v2::tag::multiseq_viterbi>::type::arg_type(i),
                       0 );

    }
    std::cout<< "      run_all: done!"<< std::endl;
    */
    apply_on_segment(comp,
                     v2::tag::multiseq_viterbi(),
                     typename Comp::template algo<v2::tag::multiseq_viterbi>::type::arg_type(0),
                     0 );

    return result_type(true);
  }

public:
  template<typename Comp>
  boost::fusion::vector<size_t> get_extents(Comp& comp)
  {
#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "run_all::get_extents()"<< std::endl;
#endif
    return boost::fusion::vector<size_t>(0);
  }

public:
  static result_type get_unassigned_value() { return result_type(false); }


public:
  const char* get_debug_id() const { return "run_all ONE"; }

};

} // namespace detail

template <class Model>
struct implementations<tag::run_all, Model >
{
  typedef typename detail::run_all<Model> impl_type;
};


//////////////////////////////////////////////////////////////////////////////////////////

namespace detail
{

template <typename Record>
struct adapt_delta_psi_record
{
  typedef boost::fusion::vector<probability_t, state_t, time_t> result_type;
};

template <>
struct adapt_delta_psi_record<boost::fusion::vector<probability_t, state_t, time_t> >
{
  typedef boost::fusion::vector<probability_t, state_t, time_t> result_type;
  boost::fusion::vector<probability_t, state_t, time_t> operator()(const boost::fusion::vector<probability_t, state_t, time_t>& v) { return v; }
};

template <>
struct adapt_delta_psi_record<boost::fusion::vector<probability_t, state_t> >
{
  typedef boost::fusion::vector<probability_t, state_t, time_t> result_type;
  boost::fusion::vector<probability_t, state_t, time_t> operator()(const boost::fusion::vector<probability_t, state_t>& v) { return boost::fusion::vector<probability_t, state_t, time_t>( boost::fusion::at_c<0>(v), boost::fusion::at_c<1>(v), 1); }
};



template<typename Model>
class viterbi_decoding2 :
  public Algo<
      decoding_t
    , boost::fusion::vector<>
    , Model
  >
  , public depends_on<
      boost::mpl::set< v2::tag::delta_psi>
    , sequential_processing::single_sequence_scope
  >
  , public memoize<
      viterbi_decoding2<Model>
    , boost::fusion::vector<>
    , decoding_t
    , sequential_processing::single_sequence_scope
  >
  , public applies_to_data<
      sequential_processing::single_sequence_scope
  >
{
public:
  typedef Algo<
    decoding_t,
    boost::fusion::vector<>,
    Model
  > base;
  typedef typename base::result_type result_type;
  typedef typename base::arg_type arg_type;
  typedef typename sequential_processing::single_sequence_scope data_level;  

public:
  template <typename Comp>
  result_type produce(Comp& comp, const arg_type&)
  {
#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "viterbi_decoding(<ignoring args>)"<< std::endl;
#endif

    typedef typename Comp::template algo<v2::tag::delta_psi>::type::result_type native_delta_psi_tau_t;
    typedef adapt_delta_psi_record<native_delta_psi_tau_t> adaptor_t;
    typedef typename adaptor_t::result_type delta_psi_tau_t;
    typedef typename boost::fusion::result_of::size<delta_psi_tau_t>::type size_of_delta_psi_tau_t;
    BOOST_MPL_ASSERT_RELATION( size_of_delta_psi_tau_t::value, ==, 3 );

    decoding_t decoding;

    const typename Comp::template data<data_level>::type seq = data_c(comp);
#ifdef EXCESSIVE_DEBUG_PRINTING
    std::cout<< "   data: "<< fasta::get_seq(seq)<< std::endl;
#endif

    decoding._s2 = std::string( fasta::get_len(seq), ' ' );

    // Initialization
    //time_t t_next = length(seq);
    state_t psi_next = end_state( this->model(comp) );

    adaptor_t adapt;

    for( time_t t = fasta::get_len(seq)+1; t>0; )
    {
      const delta_psi_tau_t delta_psi_tau = adapt(apply(comp, v2::tag::delta_psi(), typename Comp::template algo<v2::tag::delta_psi>::type::arg_type(psi_next,t) ));
      const state_t psi = boost::fusion::at_c<1>(delta_psi_tau);
      const time_t  tau = boost::fusion::at_c<2>(delta_psi_tau);
      //const probability_t delta = boost::fusion::at_c<0>(delta_psi_tau);
      //std::cout<< "("<< delta<< ", "<< psi<< ", "<< tau<< ")"<< std::endl;
      assert( tau <= t );

      const std::string state_name = this->model(comp).GetStateName(psi_next);
      decoding._s2[t-tau] = char(state_name[0]);

      psi_next = psi;
      t -= tau;

#ifdef EXCESSIVE_DEBUG_PRINTING
      std::cout<< "decoding: t:\t"<< t+tau<< "\t->\t"<< t<< "\ttau:\t"<< tau<< std::endl;
#endif
    }
    
    decoding._s1 = fasta::get_seq(seq);

    decoding._id1 = fasta::get_acc(seq);

    decoding._logp = boost::fusion::at_c<0>( adapt(apply(comp, v2::tag::delta_psi(), typename Comp::template algo<v2::tag::delta_psi>::type::arg_type(
																		      end_state( this->model(comp) )
																		      ,fasta::get_len(seq)+1) )) );

    return decoding;
  }

public:
  template<typename Comp>
  boost::fusion::vector<> get_extents(Comp& comp)
  {
    return boost::fusion::vector<>();
  }

public:
  static result_type get_unassigned_value() { return result_type(NullDecoding()); }

public:
  const char* get_debug_id() const { return "viterbi_decoding"; }

};

} //namespace detail


//////////////////////////////////////////////////////////////////////////////////////////

template <class Model>
struct implementations<tag::viterbi_decoding, Model >
{
  typedef typename detail::viterbi_decoding2<Model> impl_type;
};


} // namespace v2


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class Algo>
void debug_print( const ViterbiAlgorithm<Algo>& obj )
{
  obj.debug_print();
}


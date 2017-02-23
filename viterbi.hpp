#pragma once
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#undef BOOST_DISABLE_ASSERTS
#include <vector>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include "common.hpp"



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
  enum {_empty=1	};

public:
  P vli( size_t l, size_t i )
  {
    using boost::get;

    // Is the value already cached?
    const P val = _v[l][i];
    if( val != _empty )
      return val;

    extval_t newext;
    if( i == 0 )
    {
      newext = boost::make_tuple( log(1. / num_states(*_model) ),
				  std::numeric_limits<prev_ptr_t>::max() );
    }
    else
    {
      newext = _produce_vli( l, i );
    }

    // Value not cached - need to calculate it
    //_v[l][i]   = newext.get<0>();
    _v[l][i]   = get<0>(newext);
    _ptr[l][i] = get<1>(newext);
    return get<0>(newext);
  }
protected:

  extval_t _produce_vli( size_t l, size_t i )
  {
    using boost::get; using boost::make_tuple;

    alternatives_t alternatives;
    for( size_t k = 0; k< num_states(*_model); ++k )
    {
      const P vk_prev = vli(k, i-1);
      const P p = vk_prev + _model->a(k, l);
      alternatives.push_back( make_tuple( p, k ) );
    }

    extval_t maxk = getmax( alternatives );
    return make_tuple( get<0>(maxk) + _model->e(l, (*_seq)[i] ),
		       get<1>(maxk) );
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
    for( size_t k = 0; k< num_states(*_model); ++k )
    {
      const P p = vli(k, length(*_seq)-1);
      alternatives.push_back( boost::make_tuple( p, k ) );
    }
    extval_t maxk = getmax( alternatives );
    
    return get<0>(maxk);
  }

public:
  void get_backtrace(std::vector<std::string>& bt)
  {
    using boost::get;
    alternatives_t alternatives;
    
	// Decide what should be the last state
    const size_t last = length(*_seq)-1;
    for( size_t k = 0; k< num_states(*_model); ++k )
    {
      const P p = vli(k, last);
      alternatives.push_back( boost::make_tuple( p, k ) );
    }
    extval_t maxk = getmax( alternatives );

	// Backtrack from the last state
    int pos = last;
    size_t state = get<1>(maxk);
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


template<class Algo>
void debug_print( const ViterbiAlgorithm<Algo>& obj )
{
  obj.debug_print();
}



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
#include <iostream>
#include <limits>
#include <utility>
#include <cassert>
#include <math.h>
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#undef BOOST_DISABLE_ASSERTS
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/at_c.hpp> 
#include <boost/utility/enable_if.hpp>
#include "logspace.hpp"


//typedef LogspaceDouble<double> double_logscale;
typedef double double_logscale;


template<typename Ret, typename Arg, typename Impl>
struct UnaryFunc
{
	typedef Ret return_type;
	typedef Arg arg_type;
	return_type operator()(const arg_type& arg) {	return return_type();	}
};

template<
  typename Sequence
  , typename P
  , template <typename A> class Model
  , template <typename A> class Viterbi
  //  , template <typename A, typename Imp> class Viterbi
  //  , template <typename A> class ViterbiImpl
  //  , template <typename A, typename Imp> class ViterbiBegin
  // , template <typename A> class ViterbiBeginImpl
  , template <typename A, typename Imp> class Forward
  , template <typename A> class ForwardImpl
  , template <typename A, typename Imp> class ForwardBegin
  , template <typename A> class ForwardBeginImpl
  , template <typename A, typename Imp> class Backward
  , template <typename A> class BackwardImpl
  , template <typename A, typename Imp> class BackwardBegin
  , template <typename A> class BackwardBeginImpl
  , template <typename A> class BaumWelch
  , bool UseEmissionsClasses=true
  , size_t D=30
>
struct HiddenMarkovAlgorithm
{
  typedef HiddenMarkovAlgorithm<Sequence,P,Model,Viterbi,/*ViterbiImpl,ViterbiBegin,ViterbiBeginImpl,*/Forward,ForwardImpl,ForwardBegin,ForwardBeginImpl,Backward,BackwardImpl,BackwardBegin,BackwardBeginImpl,BaumWelch,UseEmissionsClasses,D> self_type;
  typedef Sequence sequence_type;
  typedef P probability_type;
  typedef Model<self_type> model_type;
  typedef Viterbi<self_type> viterbi_algo_type;
  //  typedef Viterbi<self_type, ViterbiImpl<self_type> > viterbi_algo_type;
  // typedef ViterbiBegin<self_type, ViterbiBeginImpl<self_type> > viterbi_begin_algo_type;
  typedef Forward<self_type, ForwardImpl<self_type> > forward_algo_type;
  typedef ForwardBegin<self_type, ForwardBeginImpl<self_type> > forward_begin_algo_type;
  typedef Backward<self_type, BackwardImpl<self_type> > backward_algo_type;
  typedef Backward<self_type, BackwardBeginImpl<self_type> > backward_begin_algo_type;
  typedef BaumWelch<self_type> baumwelch_algo_type;
  typedef typename Sequence::value_type symbol_type;
  typedef boost::mpl::bool_<UseEmissionsClasses> use_emissions_classes;
  enum { max_duration = D };
};


template<typename T>
size_t length(const std::vector<T>& vec) { return vec.size(); }

size_t length(const std::string& t);



template<typename C>
std::string stringify( const C& a, size_t from, size_t to );


template<class T>
void debug_print( const T& obj );

template<class T>
void debug_print( const boost::multi_array<T, 2>& obj )
{
  const size_t m_end = obj.shape()[0];
  const size_t n_end = obj.shape()[1];
  std::cout<< boost::format("multi_array([%d][%d])") % m_end % n_end<< std::endl;

  for( size_t m=0; m< m_end; ++m )
  {
    for( size_t n=0; n< n_end; ++n )
    {
      std::cout<< obj[m][n]<< " ";
    }
    std::cout<< std::endl;
  }
}

template<class T>
void debug_print( const boost::multi_array<T, 1>& obj )
{
  const size_t m_end = obj.shape()[0];
  std::cout<< boost::format("multi_array([%d])") % m_end << std::endl;

  for( size_t m=0; m< m_end; ++m )
  {
      std::cout<< obj[m]<< " "<< std::endl;
  }
}


template<class T>
void debug_print( const boost::multi_array<T, 3>& obj )
{
  const size_t m_end = obj.shape()[0];
  const size_t n_end = obj.shape()[1];
  const size_t o_end = obj.shape()[2];
  std::cout<< boost::format("multi_array([%d][%d][%d])") % m_end % n_end % o_end<< std::endl;

  for( size_t m=0; m< m_end; ++m )
  {
	  std::cout<< " -- "<< m<< " -- "<< std::endl;
	  for( size_t n=0; n< n_end; ++n )
	  {
		  for( size_t o=0; o< o_end; ++o )
		  {
			  std::cout<< obj[m][n][o]<< " ";
		  }
		  std::cout<< std::endl;
	  }
  }
}


/*
template<class T>
void debug_print( const std::pair<T,T>& iterators )
{
  std::cout<< "iterator_range<T>"<< std::endl;
  for( T it = iterators.first; it != iterators.second; ++it )
  {
    debug_print(*it);
  }
  std::cout<< std::endl;
}

void debug_print( const boost::tuple<double,size_t>& alternative )
{
  std::cout<< "<"<< alternative.get<0>()<< ", "<< alternative.get<1>()<< ">"<< std::endl;
}

*/

/*
 * TODO - Fix all manual logspace conversion by using a new class (instead of double)
 */

template<class T>
void debug_print_exp( const boost::multi_array<T, 2>& obj )
{
  const size_t m_end = obj.shape()[0];
  const size_t n_end = obj.shape()[1];
  std::cout<< boost::format("multi_array([%d][%d])") % m_end % n_end<< std::endl;

  for( size_t m=0; m< m_end; ++m )
  {
    for( size_t n=0; n< n_end; ++n )
    {
		std::cout<< exp(obj[m][n])<< " ";
    }
    std::cout<< std::endl;
  }
}

template<class T>
void debug_print_exp( const boost::multi_array<T, 3>& obj )
{
  const size_t m_end = obj.shape()[0];
  const size_t n_end = obj.shape()[1];
  const size_t o_end = obj.shape()[2];
  std::cout<< boost::format("multi_array([%d][%d][%d])") % m_end % n_end % o_end<< std::endl;

  for( size_t m=0; m< m_end; ++m )
  {
	  std::cout<< " -- "<< m<< " -- "<< std::endl;
	  for( size_t n=0; n< n_end; ++n )
	  {
		  for( size_t o=0; o< o_end; ++o )
		  {
			  std::cout<< exp(obj[m][n][o])<< " ";
		  }
		  std::cout<< std::endl;
	  }
  }
}

template<typename T>
bool isvalidprobability(const T val) { return ( (val >= 0.0) && (val <= 1.00001) ); } // TODO FIX THIS!!!!

template<typename C1, typename C2>
void print_seqs( const C1& a, const C2& b )
{
  const size_t width = 60;
  size_t pos = 0;
  //std::cout<< "A("<< a.size()<< ")\nB("<< b.size()<< ")"<< std::endl;

  while( pos<a.size() )
  {
    std::cout<< stringify( a, pos, pos+width )<< std::endl;
    std::cout<< stringify( b, pos, pos+width )<< std::endl;
    std::cout<< std::endl;
    pos += width;
  }
}



struct null_deleter
{
    void operator()(void const *) const {}
};

template<class EndAlgo, class BeginAlgo>
void bind_begin_end_algorithms( EndAlgo& end, BeginAlgo& begin )
{
	end.set_begin_func( begin );
	begin.set_end_func( end );
}

template<class EndAlgo, class BeginAlgo>
void bind_begin_end_algorithms( boost::shared_ptr<EndAlgo> end, boost::shared_ptr<BeginAlgo> begin )
{
	end->set_begin_func( begin );
	begin->set_end_func( end );
}


bool not_suspected_negative_subscript(size_t k);

// Convert boost::fusion sequence to standard container
template<typename Seq>
boost::array<size_t,2> get_standard_container(const Seq& seq)
{
	boost::array<size_t,2> container;
	container[0] = boost::fusion::at_c<0>(seq);
	container[1] = boost::fusion::at_c<1>(seq);
	return container;
}


template< typename P>
void assert_range_probability_logspace(P p) { assert( p <= 0.0 ); }

template< typename P>
void assert_range_probability_realspace(P p) { assert( p <= 1.0 ); assert ( p >= 0.0 ); }

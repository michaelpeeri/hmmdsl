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
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include "common.hpp"

namespace tests {
   
/*
 * Test 28 (Rabiner p. 263)
 * The probabilities of being in any state (as calculated by the forward and backward algorithms) must sum to 1
 * Tests for: forward, backward
 */
template<typename Algo>
void test_28(typename Algo::forward_algo_type& f, typename Algo::backward_begin_algo_type& b, const typename Algo::model_type& model, const typename Algo::sequence_type& seq)
{
  using namespace boost::accumulators;
  const typename Algo::probability_type Pseq = f.calc();
  //std::cout<< "test (28)"<< std::endl;

  if( Pseq <= -std::numeric_limits<typename Algo::probability_type>::max() )
  {
	  // Pseq is used as the denominator in 
	  std::cout<< "test_28: P(seq|model) = 0.0"<< std::endl;
	  assert( Pseq > -std::numeric_limits<typename Algo::probability_type>::max() );
  }

  

  const typename Algo::probability_type epsilon = 1e-10;

  // Verify that Pfor == Pback
  {
	  const typename Algo::probability_type Pfor = f.calc();
	  const typename Algo::probability_type Pback = b.calc();
	  if( fabs(Pfor-Pback) >= epsilon)
	  {
		  std::cout<< "ERROR: Pforw(O|model) != Pback(O|model)"<< std::endl;
		  std::cout<< "ERROR: Pforw(O|model) = "<< Pfor << std::endl;
		  std::cout<< "ERROR: Pback(O|model) = "<< Pback<< std::endl;
		  assert(fabs(Pfor-Pback) < epsilon);
	  }
  }
  
	  
  accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;

  for( size_t i=1; i< length(seq)+2; ++i )
  {
    typename Algo::probability_type sum = -std::numeric_limits<typename Algo::probability_type>::max();

    for( size_t k=0; k< num_states(model); ++k )
    {
		
		//boost::array<size_t,2> args = {k,i};
		const typename Algo::probability_type fki = f( typename Algo::forward_algo_type::arg_type(k,i) );
		if( fki > 1e-10 )
		{		  std::cout<< "test_28: fki>1.0 for i="<< i<< ", k="<< k<< ", f(k,i)="<< fki<< std::endl;	  }
		
		const typename Algo::probability_type bki = b( typename Algo::backward_begin_algo_type::arg_type(k,i) );
		if( bki > 1e-10 )
		{		  std::cout<< "test_28: bki>1.0 for i="<< i<< ", k="<< k<< ", b(k,i)="<< bki<< std::endl;	  }
		
		const typename Algo::probability_type gamma_ti = fki + bki - Pseq;
		if( gamma_ti > 1e-10 )
		{		  std::cout<< "test_28: gamma_ti>1.0 for i="<< i<< ", k="<< k<< ", fki="<< fki<< ", bki="<< bki<< ", Pseq="<< Pseq<< ", gamma_ti = fki+bki-Pseq = "<< gamma_ti<< std::endl;	  }
		
		sum = logspace_sum( gamma_ti, sum );
    }
	if( sum > 1e-10 )
	{
		std::cout<< "sigma_i_gamma_ti > 1.0 for i="<< i<< ", Pseq="<< Pseq<< " (sigma_i_gamma_ti="<< sum<< ")"<< std::endl;
		for( size_t k=0; k< num_states(model); ++k )
		{
			const typename Algo::probability_type fki = f(typename Algo::forward_algo_type::arg_type(k, i) );
			const typename Algo::probability_type bki = b(typename Algo::backward_algo_type::arg_type(k, i) );
			const typename Algo::probability_type gamma_ti = fki + bki - Pseq;
			std::cout<< "k="<< k<< ", fki="<< fki<< ", bki="<< bki<< ", gamma_ti="<< gamma_ti<< std::endl;
		}
	}

	std::cout<< "i="<< i<< ", sum="<< sum<< std::endl;
	
    acc(sum);
  }
#ifndef NO_TESTS
  const double min_err = min(acc);
  const double max_err = max(acc);
  
  // All results are in logspace;
  assert( fabs(min_err) < epsilon );
  assert( fabs(max_err) < epsilon );
#endif
  //std::cout<< min(acc)<< " <= Error <= "<< max(acc)<< "\t(average = "<< mean(acc)<< ")"<< std::endl;
}

/*
 * Pforw(seq)==Pback(seq)
 * Reference: Durbin p.59
 */
template<typename Algo>
void test_D59_forw_back(typename Algo::forward_algo_type& f, typename Algo::backward_begin_algo_type& b)
{
	const double pf = f.calc();
	const double pb = b.calc();
	if( fabs(pf-pb) > 1e-10 )
	{
		std::cout<< "Error: P(seq) yields different results on forward and backard paths!"
				 << "\tPforward(seq) = "<< pf
				 << "\tPbackward(seq) = "<< pb
				 << std::endl;
		//assert( fabs(pf-pb) <= 1e-10 );
	}
	
}

template<typename Algo>
typename boost::enable_if<typename Algo::model_type::is_explicit_duration_model,void>::type
test_duration_valid_pdf(const typename Algo::model_type& model)
{
	using namespace boost::accumulators;
	static const double epsilon = 1e-3; // TODO - consider raising this?

	//std::cout<< "test_duration_valid_pdf: ";
	bool bSuccess = true;
	accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;

	for( size_t k=0; k< num_states(model); ++k )
	{
		//if(!model.IsReservedState(k)) //
		//	std::cout<< "k="<< k<< ", mu="<< model.GetMu(k)<< ", eta="<< model.GetEta(k)<< std::endl;//
		
		typename Algo::probability_type sum = 0.0;
		for( size_t d=1; d<=Algo::max_duration; ++d )
		{
			//if(!model.IsReservedState(k))//
			//	std::cout<< d<< " "<< exp(model.p(k, d))<< std::endl;//
			sum += exp(model.p(k, d));
		}

		if( fabs(sum - 1.0) > epsilon )
		{
			std::cout<< "Error: model.p() for state="<< k<< " is not a valid PDF! (with D="<< Algo::max_duration<< ")"<< std::endl;
			std::cout<< "Error: Sigma_i P(i) = "<< sum<< std::endl;
			for( size_t d=1; d<=Algo::max_duration; ++d )
				std::cout<< "P(d="<< d<< ") = "<< exp(model.p(k, d))<< std::endl;
			bSuccess = false;
		}
		acc(sum - 1.0);
	}
#ifndef NO_TESTS
	const double min_err = min(acc);
	const double max_err = max(acc);
	assert( fabs(min_err) < epsilon );
	assert( fabs(max_err) < epsilon );
#endif
	//if( bSuccess ) std::cout<< "OK (mean err="<< mean(acc)<< ")"<< std::endl;
}


template<typename Algo>
void test_transitions_valid_pdf(const typename Algo::model_type& model)
{
	using namespace boost::accumulators;
	static const double epsilon = 1e-6;

	//std::cout<< "test_transitions_valid_pdf: ";
	bool bSuccess = true;
	accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;

	for( size_t k=0; k< num_states(model); ++k )
	{
		if( k == model.GetTerminalState()) continue; // Don't test the terminal (since it has no transitions!)
		typename Algo::probability_type sum = 0.0;
		for( size_t l=0; l< num_states(model); ++l )
			sum += exp(model.a(k, l));

		if( fabs(sum - 1.0) > epsilon )
		{
			std::cout<< "Error: model.a() for state="<< k<< " is not a valid PDF!"<< std::endl;
			std::cout<< "Error: Sigma_i a(i) = "<< sum<< std::endl;
			for( size_t l=0; l< num_states(model); ++l )
				std::cout<< "a("<< k<< "=>"<< l<< ") = "<< exp(model.a(k, l))<< std::endl;
			bSuccess = false;
		}
		acc(sum - 1.0);
	}
#ifndef NO_TESTS
	const double min_err = min(acc);
	const double max_err = max(acc);
	assert( fabs(min_err) < epsilon );
	assert( fabs(max_err) < epsilon );
#endif
	//if( bSuccess ) std::cout<< "OK (mean err="<< mean(acc)<< ")"<< std::endl;
}


template<typename Algo>
void test_emissions_valid_pdf(const typename Algo::model_type& model)
{
	using namespace boost::accumulators;
	static const double epsilon = 1e-6;

	//std::cout<< "test_emissions_valid_pdf: ";
	bool bSuccess = true;
	accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;

	for( size_t k=0; k< num_states(model); ++k )
	{
		if( model.is_silent(k)) continue; // Silent state don't have emission probs

		typename Algo::probability_type sum = 0.0;
		for( size_t s=0; s< num_symbols(model); ++s )
			sum += exp(model.e(k, s));

		if( fabs(sum - 1.0) > epsilon )
		{
			std::cout<< "Error: model.e() for state="<< k<< " is not a valid PDF!"<< std::endl;
			std::cout<< "Error: Sigma_i e(i) = "<< sum<< std::endl;
			for( size_t s=0; s< num_symbols(model); ++s )
				std::cout<< "e(k="<< k<< ", s="<< s<< ") = "<< exp(model.e(k, s))<< std::endl;
			bSuccess = false;
		}
		acc(sum - 1.0);
	}
#ifndef NO_TESTS
	const double min_err = min(acc);
	const double max_err = max(acc);
	assert( fabs(min_err) < epsilon );
	assert( fabs(max_err) < epsilon );
#endif
	//if( bSuccess ) std::cout<< "OK (mean err="<< mean(acc)<< ")"<< std::endl;
}


template<typename Algo>
typename boost::enable_if<typename Algo::model_type::is_explicit_duration_model,void>::type
test_no_self_transitions(const typename Algo::model_type& model)
{
	using namespace boost::accumulators;
	//static const double epsilon = 1e-6;

	//std::cout<< "test_no_self_transitions: ";
	bool bSuccess = true;

	for( size_t k=0; k< num_states(model); ++k )
	{
		if( exp(model.a(k, k)) >1e-50 )
		{
			std::cout<< "Error: model.a() for state="<< k<< " has a self-transition!"<< std::endl;
			std::cout<< "a("<< k<< "=>"<< k<< ") = "<< exp(model.a(k, k))<< std::endl;
			bSuccess = false;
		}
	}
	assert(bSuccess);
}
	


template<typename Algo>
typename boost::enable_if<typename Algo::model_type::is_explicit_duration_model,void>::type
test_model_valid(const typename Algo::model_type& model)
{
	// Standard tests
	test_transitions_valid_pdf<Algo>(model);
	test_emissions_valid_pdf<Algo>(model);
	// Explicit duration-related tests
	test_duration_valid_pdf<Algo>(model);
	test_no_self_transitions<Algo>(model);
}


template<typename Algo>
typename boost::disable_if<typename Algo::model_type::is_explicit_duration_model,void>::type
test_model_valid(const typename Algo::model_type& model)
{
	test_transitions_valid_pdf<Algo>(model);
	test_emissions_valid_pdf<Algo>(model);
}
	
	
} //  namespace tests

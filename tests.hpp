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
	

/*
 * Test 28 (Rabiner p. 263)
 * The probabilities of being in any state (as calculated by the forward and backward algorithms) must sum to 1
 * Test for: BaumWelchAlgorithm::gamma()
 */

template<typename Algo>
typename boost::enable_if<typename Algo::baumwelch_algo_type::gamma_in_logspace,void>::type
test_28(typename Algo::baumwelch_algo_type& bw, const typename Algo::model_type& model, const typename Algo::sequence_type& seq)
{
    using namespace boost::accumulators;
    //std::cout<< "test (28)"<< std::endl;
    
    const typename Algo::probability_type epsilon = 1e-10;
    
    accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;
    bw.calc_all();
    debug_print(bw);
    
    const size_t len = length(seq)+2;
    for( size_t i=0; i< len; ++i )
    {
        typename Algo::probability_type sigmak_gamma = -std::numeric_limits<typename Algo::probability_type>::max();
        
        for( size_t k=0; k< num_states(model); ++k )
        {
            //if( model.IsReservedState(k) ) continue;
            
            const typename Algo::probability_type gamma_ti = bw.gamma(k, i);
            if( gamma_ti > 1e-10 )
            {			std::cout<< "test_28: gamma out of range for state="<< k<< ", position="<< i<< ", sum="<< exp(gamma_ti)<< std::endl;		}
            
            //std::cout<< "gamma(k="<< k<< ", i="<< i<< ") = "<< gamma_ti<< std::endl;
            
            sigmak_gamma = logspace_sum( gamma_ti, sigmak_gamma );
        }
	if( (sigmak_gamma > 1e-10 ) || (sigmak_gamma < -0.01) )
	{
            std::cout<< "test_28: sum out of range for position "<< i<< ", sum="<< exp(sigmak_gamma)<< std::endl;
            for( size_t k=0; k< num_states(model); ++k )
            {
                if( model.IsReservedState(k) ) continue;
                const typename Algo::probability_type gamma_ti = bw.gamma(k, i);
                std::cout<< gamma_ti<< " ";
            }
            std::cout<< std::endl;
	}
        
	//std::cout<< "i="<< i<< ", sigma_gamma="<< sigmak_gamma<< std::endl;
	// Update statistics
        acc(sigmak_gamma);
    }
    // All results are in logspace;
    const typename Algo::probability_type min_err = min(acc);
    const typename Algo::probability_type max_err = max(acc);
    
    debug_print(bw);
    
    assert( fabs(min_err) < epsilon );
    assert( fabs(max_err) < epsilon );
    //std::cout<< min(acc)<< " <= Error <= "<< max(acc)<< "\t(average = "<< mean(acc)<< ")"<< std::endl;
}

template<typename Algo>
typename boost::disable_if<typename Algo::baumwelch_algo_type::gamma_in_logspace,void>::type
test_28(typename Algo::baumwelch_algo_type& bw, const typename Algo::model_type& model, const typename Algo::sequence_type& seq)
{
    using namespace boost::accumulators;
    //std::cout<< "test (28)"<< std::endl;
    
    const typename Algo::probability_type epsilon = 1e-10;
    
    accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;
    bw.calc_all();
    debug_print(bw);
    
    const size_t len = length(seq)+2;
    for( size_t i=0; i< len; ++i )
    {
        typename Algo::probability_type sigmak_gamma = -std::numeric_limits<typename Algo::probability_type>::max();
        
        for( size_t k=0; k< num_states(model); ++k )
        {
            //if( model.IsReservedState(k) ) continue;
            
            const typename Algo::probability_type gamma_ti = bw.gamma(k, i);
            if( gamma_ti > 1e-10 )
            {			std::cout<< "test_28: gamma out of range for state="<< k<< ", position="<< i<< ", sum="<< exp(gamma_ti)<< std::endl;		}
            
            //std::cout<< "gamma(k="<< k<< ", i="<< i<< ") = "<< gamma_ti<< std::endl;
            
            sigmak_gamma += gamma_ti;
        }
	if( (sigmak_gamma > 1e-10 ) || (sigmak_gamma < -0.01) )
	{
            std::cout<< "test_28: sum out of range for position "<< i<< ", sum="<< exp(sigmak_gamma)<< std::endl;
            for( size_t k=0; k< num_states(model); ++k )
            {
                if( model.IsReservedState(k) ) continue;
                const typename Algo::probability_type gamma_ti = bw.gamma(k, i);
                std::cout<< gamma_ti<< " ";
            }
            std::cout<< std::endl;
	}
        
	//std::cout<< "i="<< i<< ", sigma_gamma="<< sigmak_gamma<< std::endl;
	// Update statistics
        acc(sigmak_gamma);
    }
    // All results are in logspace;
    const typename Algo::probability_type min_err = min(acc);
    const typename Algo::probability_type max_err = max(acc);
    
    debug_print(bw);
    
    assert( fabs(min_err) < epsilon );
    assert( fabs(max_err) < epsilon );
    //std::cout<< min(acc)<< " <= Error <= "<< max(acc)<< "\t(average = "<< mean(acc)<< ")"<< std::endl;
}
	

template<typename Algo> //, typename BW>
void test_38(typename Algo::baumwelch_algo_type& bw, const typename Algo::model_type& model, const typename Algo::sequence_type& seq)
//void test_38(BW& bw, const typename Algo::model_type& model, const typename Algo::sequence_type& seq)
{
    using namespace boost::accumulators;
    //std::cout<< "test (38)"<< std::endl;
    return;// TODO - Restore this test using the new Sigma_il_xi()
    
    const typename Algo::probability_type epsilon = 1e-10;
    
    accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;
    
    for( size_t i=0; i< length(seq)+1; ++i )
    {
        for( size_t k=0; k< num_states(model); ++k )
	{
            const typename Algo::probability_type gammaki = bw.gamma( k, i );
            /*  // TODO - Sigmal_xi with the original signature is no longer available
            const typename Algo::probability_type sigmalxi = bw.Sigmal_xi( k, i );
            const typename Algo::probability_type delta = exp(gammaki) - sigmalxi;
            if( fabs(delta) > 0.001 )
            {			std::cout<< "test_38: delta>0.001: state="<< k<< "\ti="<< i<< ":\tgammaki="<< exp(gammaki)<< "\tsigmal_xi:"<< sigmalxi<< std::endl;		}
            
            // Update statistics
            acc( delta );
            */
	}
	
    }
#ifndef NO_TESTS
    // All results are in logspace;
  const typename Algo::probability_type min_delta = min(acc);
  assert( fabs(min_delta) < epsilon );
  const typename Algo::probability_type max_delta = max(acc);
  assert( fabs(max_delta) < epsilon );
#endif
  //std::cout<< min(acc)<< " <= Error <= "<< max(acc)<< "\t(average = "<< mean(acc)<< ")"<< std::endl;
}


template<typename Algo, typename BW>
void test_43b(BW& bw, const typename Algo::model_type& model)
{
    using namespace boost::accumulators;
    //std::cout<< "test (43b)"<< std::endl;
    
    const typename Algo::probability_type epsilon = 1e-12;
    bool valid = true;
    
    accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;
    
    for( size_t k=0; k< num_states(model); ++k )
    {
        typename Algo::probability_type sigma = 0.0;
        if( k == model.GetTerminalState()) continue;
	
        for( size_t l=0; l< num_states(model); ++l )
        {
            const typename Algo::probability_type p = exp(bw.nexta(k, l));
            assert( (p >= 0.0) && (p <= 1.000000001) );
            sigma += p;
        }
        if( ( sigma > 1.000000001 ) || ( sigma < 0.9 ) )
        {
            std::cout<< "test_43b: Error: invalid sigma_l_nexta(k,l) for state "<< k<< " (="<< sigma<< ")"<< std::endl;	  
            for( size_t l=0; l< num_states(model); ++l )
            {			  std::cout<< exp(bw.nexta(k, l))<< " ";		  }
            std::cout<< std::endl;
            valid = false;
        }
        
        // Update statistics
        acc( 1.0 - sigma );
    }
    // All results are in logspace;
    const typename Algo::probability_type min_err = min(acc);
    const typename Algo::probability_type max_err = max(acc);
    assert( fabs(min_err) < epsilon );
    assert( fabs(max_err) < epsilon );
    assert( valid );
    //std::cout<< min(acc)<< " <= Error <= "<< max(acc)<< "\t(average = "<< mean(acc)<< ")"<< std::endl;
}

template<typename Algo>
void test_40c(typename Algo::baumwelch_algo_type& bw, const typename Algo::model_type& model )
{
    using namespace boost::accumulators;
    //std::cout<< "test (40c)"<< std::endl;
    
    const typename Algo::probability_type epsilon = 1e-10;
    
    accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;
    
    for( size_t k=0; k< num_states(model); ++k )
    {
        const typename Algo::probability_type sigmai_gamma = bw.Sigmai_gamma(k);
        typename Algo::probability_type sigmas_gamma = 0.0;
	
        for( size_t s=0; s< num_symbols(model); ++s )
        {
            //const typename Algo::symbol_type sym = model.get_symbol(s);
            const typename Algo::probability_type sigmai_gamma_obs_s = bw.Sigmai_gamma_observe_s(k, s);
            assert( ( sigmai_gamma_obs_s == 0.0 ) || (!model.is_silent(k) ) );
            sigmas_gamma += sigmai_gamma_obs_s;
        }
        if( model.is_silent(k) ) continue;
	
        const typename Algo::probability_type delta = sigmai_gamma - sigmas_gamma;
        if( fabs(delta) > epsilon )
        {
            std::cout<< "test_40c: delta for symbol "<< k<< " exceeds limit (sigmai_gamma="<< sigmai_gamma<< ", sigmas_gamma="<< sigmas_gamma<< ")"<< std::endl;
        }
	
        // Update statistics
        acc( delta );
    }
    
    // All results are in logspace;
#ifndef NO_TESTS
    const typename Algo::probability_type min_err = min(acc);
    assert( fabs(min_err) < epsilon );
    const typename Algo::probability_type max_err = max(acc);
    assert( fabs(max_err) < epsilon );
#endif
    //std::cout<< min(acc)<< " <= Error <= "<< max(acc)<< "\t(average = "<< mean(acc)<< ")"<< std::endl;
}

} // namespace tests


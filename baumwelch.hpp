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
#include <vector>
#include <boost/core/is_same.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/python.hpp> // for BPL List
#include <boost/fusion/container/vector.hpp>
#include "asa121.hpp"
#include "common.hpp"
#include "v2.hpp"
#include "read_fasta.hpp"
#include "tests.hpp"
#include "python_helpers.hpp"
#include "forward.hpp"
#include "backward.hpp"
#include "gamma.hpp"
#include "xi.hpp"


namespace baumwelch_algo
{
	
template<typename Algo>
class BaumWelchAlgorithm
{
protected:
    typedef typename Algo::probability_type P;
    enum { _empty=1};
    
public:

    BaumWelchAlgorithm(boost::shared_ptr<const typename Algo::model_type> model, 
                       boost::shared_ptr<const typename Algo::sequence_type> seq )
        : _model(model)
	, _seq(seq)
	, _nextA(boost::extents[num_states(*model)][num_states(*model)] )
	, _nextB(boost::extents[num_states(*model)][num_symbols(*model)] )
    {
        _reset();
    }
    
    /*
     * Memoization support wrapper for gamma (HMM)
     */
public:
    P gamma( size_t k, size_t i)
    {
        assert(false);  // TODO: re-implement using v2 (based on HSMMBaumWelch)
        throw;
    }

    /*
     * Calculate Sigma_i_gamma(k) (HMM)
     * Note: Individial gamma values are memoized, but the sigma is not (TODO - check if this is worthwhile optimizing)
     */
public:
    P Sigmai_gamma( size_t k )
    {
        assert(false);  // TODO: re-implement using v2 (based on HSMMBaumWelch)
        throw;
    }

    /*
     * Calculate Sigma_i_gamma_obs_s(k) (HMM)
     * Note: Individial gamma values are memoized, but the sigma is not (TODO - check if this is worthwhile optimizing)
     */
public:
    P Sigmai_gamma_observe_s( size_t k, size_t s )
    {
        assert(false);  // TODO: re-implement using v2 (based on HSMMBaumWelch)
        throw;
    }

public:
    void calc_all()
    {
        // TODO Impl. this
        assert(false);  // TODO: re-implement using v2 (based on HSMMBaumWelch)
        throw;
    }
	
    /*
     * The probability of being in state k at position i (given the sequence and model) (HMM)
     * Rabiner eq. (27)
     */
protected:
    P _produce_gamma( size_t k, size_t i )
    {
        assert(false);  // TODO: re-implement using v2 (based on HSMMBaumWelch)
        throw;
    }
	
			

    /*
     * Memoization support wrapper for xi (HMM)
     */
public:
    P xi( size_t k, size_t l, size_t i )
    {
        assert(false);  // TODO: re-implement using v2 (based on HSMMBaumWelch)
        throw;
    }

    /*
     * Calculate Sigma_i_xi(k,l) (HMM)
     * Note: Individial xi values are memoized, but the sigma is not (TODO - check if this is worthwhile optimizing)
     */
public:
    P Sigmai_xi( size_t k, size_t l)
    {
        assert(false);  // TODO: re-implement using v2 (based on HSMMBaumWelch)
        throw;
    }
	
    /*
     * Calculate Sigma_l_xi(k,l) (HMM)
     * Note: Individial xi values are memoized, but the sigma is not (TODO - check if this is worthwhile optimizing)
     */
public:
    P Sigmal_xi( size_t k)
    {
        // TODO - impl. this properly (see new HSMM impl.)
        assert(false);  // TODO: re-implement using v2 (based on HSMMBaumWelch)
        throw;
    }
	
			
			
    /*
     * Calculate xi(k,l,i) - The frequency of k->l transitions in position i (HMM)
     * Rabiner eq. (37)
     */
protected:
    P _produce_xi( size_t k, size_t l, size_t i )
    {
        assert(false);  // TODO: re-implement using v2 (based on HSMMBaumWelch)
        throw;
    }
    
    /*
     * Memoization support wrapper for a (HMM)
     */
public:
    P nexta( size_t k, size_t l)
     {
         P val =  _nextA[k][l];
         if( val != _empty ) return val;
         
         _state_reset = false;
         _nextA[k][l] = val = _produce_nexta(k, l);
         return val;		
     }
    
    /*
     * Memoization support wrapper for e (HMM)
     */
public:
    P nextb( size_t k, size_t s)
    {
        assert( !_model->is_silent(k) );
        P val = _nextB[k][s];
        if( val != _empty ) return val;
	
        _state_reset = false;
        
        if( !Algo::use_emissions_classes::value )
        {
            _nextB[k][s] = val = _produce_nextb(k,
                                                _model->get_symbol(s) /* _produce_nextb expects a symbol, not an index */ );
        }
        else
        {
            const size_t emclass = _model->GetEmissionsClass(k);
            
            val = _produce_nextb_class(emclass,
                                       _model->get_symbol(s) /* _produce_nextb expects a symbol, not an index */ );
            
            
            // Store the value for all members of this class
            typename Algo::model_type::emissions_class_members_iterator mem_it, mem_end;
            boost::tie( mem_it, mem_end ) = _model->GetEmissionsClassMembersRange(emclass);
            
            for( ; mem_it != mem_end; ++mem_it )
                _nextB[*mem_it][s] = val;
        }
	
        return val;		
    }

    // Unsupported operations for HMMs
    // TODO - Refactor (using v2 infrastructure) so this becomes unnecessary...
public:
	struct UnsupportedOperation {};

public:	
    P nextscale( size_t k )
    {
        // TODO - Refactor to make this unnecessary
        assert(false);
        throw UnsupportedOperation();
    }

public:
    P nextshape( size_t k )
    {
        // TODO - Refactor to make this unnecessary
        assert(false);
        throw UnsupportedOperation();
    }
	

    /*
     * Calculate P(seq|model) (HMM)
     * Rabiner eq. (21)
     */
public:
    P Pseq()
    {
        assert(false);
        throw;
    }

protected:
    void _reset()
    {
        assert(false);
        throw;
    }
	
public:
    void reset()
    {
        //std::cout<< "BW::reset()"<< std::endl;
        if( _state_reset ) return;
        _reset();
    }
	
			
    /*
     * Baum-Welch 'a' (transitions) estimation, HMM, no emissions classes
     * Rabiner eq. (40b) - Overall equation
     */
protected:
    P _produce_nexta( size_t k, size_t l)
    {
        assert(false);  // Deprecated
        throw;
    }
    
    /*
     * Baum-Welch 'e' (emissions) estimation, HMM, no emissions classes
     * Rabiner eq. (40c) - Overall equation
     */
protected:
    P _produce_nextb( size_t k, size_t s )
    {
        assert(false);
        throw;
    }
    
protected:
    /*
     * Baum-Welch 'e' (emissions) estimation, HMM, with emissions classes
     * Rabiner eq. (40c) - Overall equation
     */
    P _produce_nextb_class( size_t k, size_t s )
    {
        assert(false);
        throw;
    }
    
protected:
    const boost::shared_ptr<const typename Algo::model_type> _model;
    const boost::shared_ptr<const typename Algo::sequence_type> _seq;
    
    typedef boost::multi_array<P,2> nextA_type;
    nextA_type _nextA;
    
    typedef boost::multi_array<P,2> nextB_type;
    nextB_type _nextB;
    
    bool _state_reset;

public:
    void debug_print() const
    {
        std::cout<< "BaumWelchAlgorithm<Algo>"<< std::endl;
    }
    
public:
    typedef boost::mpl::bool_<true> gamma_in_logspace;
};

template<class Algo>
void debug_print( const BaumWelchAlgorithm<Algo>& obj )
{
    obj.debug_print();
}


namespace detail
{

template<typename Algo>
struct gamma_environment
{
    typedef boost::tuple<
        boost::shared_ptr<const typename Algo::model_type>
        , boost::shared_ptr<const typename Algo::sequence_type>
        , boost::shared_ptr<typename Algo::forward_algo_type>
        , boost::shared_ptr<typename Algo::forward_begin_algo_type>
        , boost::shared_ptr<typename Algo::backward_algo_type>
        , boost::shared_ptr<typename Algo::backward_begin_algo_type>
        > type;
    typedef boost::tuple<
        boost::shared_ptr<const typename Algo::sequence_type>
        > type2;
    enum {model, seq, forw_end, forw_begin, back_end, back_begin };
};

/*
 * Calculate gamma(l,i) values (expected number of transitions) values for HSMMs
 * Defined for HMMs in Rabiner eq. (26).
 * (In Rabiner, gamma is not defined explicitly for HSMMs)
 * This is useful for doing Baum-Welch calculation of 'e' (emissions) parameters.
 * Calculation is memoized (TODO - verify in practice)
 */
// template<class Algo>
// class GammaImpl : public UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, GammaImpl<Algo> >
// {
// protected:
//     typedef GammaImpl<Algo> self_type;
//     typedef typename Algo::probability_type P;
//     typedef UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, GammaImpl<Algo> > base_type;
//     typedef gamma_environment<Algo> env_;
    
// public:
//     typedef typename base_type::return_type val_type;
//     typedef typename base_type::arg_type arg_type;	
//     typedef typename gamma_environment<Algo>::type env_type;
    
// protected:
//     typedef MemoizedFunc<self_type, typename gamma_environment<Algo>::type, arg_type> host_type;
//     host_type& _host;
//     const boost::shared_ptr<const typename Algo::sequence_type> _seq;
//     const boost::shared_ptr<const typename Algo::model_type> _model;
//     const boost::shared_ptr<typename Algo::forward_algo_type> _forward;
//     const boost::shared_ptr<typename Algo::forward_begin_algo_type> _forward_begin;
//     const boost::shared_ptr<typename Algo::backward_algo_type> _backward;
//     const boost::shared_ptr<typename Algo::backward_begin_algo_type> _backward_begin;
    
// protected:
//     GammaImpl( host_type& host, const env_type& env)
//         : _host(host) 
//         , _seq(boost::get<env_::seq>(env))
//         , _model(boost::get<env_::model>(env))
//         , _forward(boost::get<env_::forw_end>(env))
//         , _forward_begin(boost::get<env_::forw_begin>(env))
//         , _backward(boost::get<env_::back_end>(env))
//         , _backward_begin(boost::get<env_::back_begin>(env))
//     {}
    
//     /*
//      * Calculate gamma(l, i), HSMM
//      * i.e. the chance (expectation) of occupying state l at time i.
//      * This is calculated as in Rabiner eq. (80).
//      * TODO - Verify memoization is working
//      */
// protected:
//     P _produce_fli( const arg_type& args )
//     {
//         const size_t l = boost::fusion::at_c<0>(args);
//         const size_t i = boost::fusion::at_c<1>(args);
        
//         // Currently, gamma is not supported for initial and terminal states.
//         if( i == 0 )
//             return (l==_model->GetInitialState() ) ? 1.0 : 0.0;
//         if( i == length(*_seq)+1 )
//             return (l==_model->GetTerminalState() ) ? 1.0 : 0.0;

//         // If we reached this point, occupying an initial or terminal state is impossible
//         if( _model->IsReservedState(l) ) return 0.0;
        
//         //const size_t lastpos = length(*_seq)+1;
//         //if( i == lastpos )
//         //	return (l==_model->GetTerminalState()) ? 0 : -std::numeric_limits<P>::max();

//         //std::cout<< "gam("<< l<< ", "<< i<< ") ";
//         //const arg_type prev(l, i-1);
        
//         P P_l_begins_before_i = 0.0;
//         for( size_t tau=0; tau< i; ++tau)
//         {
//             const arg_type l_tau(l, tau);

//             //const P p1 = (*_forward_begin)(l_tau);
//             //if( p1 > -std::numeric_limits<P>::max() ) std::cout<< l_tau<< "[1] ";
//             //const P p2 = (*_backward_begin)(l_tau);
//             //if( p2 > -std::numeric_limits<P>::max() ) std::cout<< l_tau<< "[2] ";
            
//             const P P_l_begins_before_tau = exp( (*_forward_begin)(l_tau) + (*_backward_begin)(l_tau) );
//             P_l_begins_before_i += P_l_begins_before_tau;
//         }

//         P P_l_ends_before_i = 0.0;
//         for( size_t tau=0; tau< i; ++tau)
//         {
//             const arg_type l_tau(l, tau);
//             const P P_l_ends_before_tau = exp( (*_forward      )(l_tau) + (*_backward      )(l_tau) );
//             P_l_ends_before_i += P_l_ends_before_tau;
//         }

//         const P P_l_occurs_at_i = P_l_begins_before_i - P_l_ends_before_i;

//         //std::cout<< "gamma(l="<< l<< ", i="<< i<< ") = "<<P_l_occurs_at_i << " (P_l_before_i="<< P_l_begins_before_i<< ", P_l_ends_before_i="<< P_l_ends_before_i<< ")"<< std::endl;
        
//         if( P_l_occurs_at_i< -1e-10 )
//         {
//             std::cout<< "Error: Got gamma(l,i)<0 for l="<< l<< ", i="<< i<< ", P_l_before_i="<< P_l_begins_before_i<< ", P_l_ends_before_i="<< P_l_ends_before_i<< std::endl;
//             return 0.0;
//         }
//         return P_l_occurs_at_i >= 0.0 ? P_l_occurs_at_i : 0.0;
//     }	
    
// public:
//     inline val_type operator()( const arg_type& args ) 		{  return _host->operator()(args); }
    
// private:
//     GammaImpl(const GammaImpl&);
//     const GammaImpl& operator=(const GammaImpl&);
//     GammaImpl& operator=(GammaImpl&);
// };


typedef boost::fusion::vector<size_t,size_t> gamma_algo_arg_t;

template<typename Algo, typename GammaImpl>
class Gamma : public MemoizedFunc<GammaImpl, typename gamma_environment<Algo>::type, gamma_algo_arg_t >
{
public:
    typedef boost::fusion::vector<size_t,size_t> arg_type;
    typedef typename gamma_environment<Algo>::type env_type;
    
protected:
    typedef MemoizedFunc<GammaImpl,typename gamma_environment<Algo>::type, arg_type> base_type;
    const boost::shared_ptr<const typename Algo::sequence_type> _seq;
    const boost::shared_ptr<const typename Algo::model_type> _model; // NOTE: not really necessary
    const boost::shared_ptr<typename Algo::forward_algo_type> _forward;
    const boost::shared_ptr<typename Algo::forward_begin_algo_type> _forward_begin;
    const boost::shared_ptr<typename Algo::backward_algo_type> _backward;
    const boost::shared_ptr<typename Algo::backward_begin_algo_type> _backward_begin;
    typedef gamma_environment<Algo> env_;
    
public:
    typedef typename base_type::val_type P;
    
public:
    using base_type::debug_print;
    using base_type::operator();
    using base_type::reset;
    using base_type::_reset;
    using base_type::fli;
    
public:
    Gamma(const env_type& env)
    : base_type( env, get_standard_container(boost::fusion::vector<size_t,size_t>(
                                                 num_states(*(boost::get<env_::model>(env)))
                                                 , length(*(boost::get<env_::seq>(env)))+2     ) ))
    , _seq(boost::get<env_::seq>(env))
    , _model(boost::get<env_::model>(env))
    , _forward(boost::get<env_::forw_end>(env))
    , _forward_begin(boost::get<env_::forw_begin>(env))
    , _backward(boost::get<env_::back_end>(env))
    , _backward_begin(boost::get<env_::back_begin>(env))
    {
        _reset();
    }
    
public:
    void calc_all()
    {
        const size_t lastpos = length(*_seq)+1;
        for( size_t i = 0; i<= lastpos; ++i )
        {
            for( typename Algo::model_type::StateId k = 0; k< num_states(*_model); ++k )
            {
                fli( arg_type(k, i) );
            }
        }
    }
    
private:
    Gamma(const Gamma&);
    const Gamma& operator=(const Gamma&);
    Gamma& operator=(Gamma&);
    
public:
    typedef boost::mpl::bool_<false> result_in_logspace;
};


template<class Algo, class Impl>
void debug_print( const Gamma<Algo, Impl>& obj )
{
	obj.debug_print();
}	

	
} // namespace detail


template<typename Algo>
class HSMMBaumWelchAlgorithm
{
protected:
    typedef typename Algo::probability_type P;
    using state_t = typename Algo::model_type::StateId;
    typedef LogspaceDouble<> logspace_t;
        

public:

    HSMMBaumWelchAlgorithm(boost::shared_ptr<const typename Algo::model_type> model, 
                           //boost::shared_ptr<typename Algo::forward_algo_type> forward, 
                           //boost::shared_ptr<typename Algo::forward_begin_algo_type> forward_begin, 
                           //boost::shared_ptr<typename Algo::backward_algo_type> backward, 
                           //boost::shared_ptr<typename Algo::backward_begin_algo_type> backward_begin, 
                           boost::shared_ptr<const typename Algo::sequence_type> seq )
        : _comp_xi(*seq, *model)
        , _model(model)
//        , _forward(forward)
//        , _forward_begin(forward_begin)
//        , _backward(backward)
//        , _backward_begin(backward_begin)
        , _seq(seq)
          //  , _gamma( typename gamma_type::env_type( model, seq, forward, forward_begin, backward, backward_begin ) )
        , _nextA(boost::extents[num_states(*model)][num_states(*model)] )
        , _nextB(boost::extents[num_states(*model)][num_symbols(*model)] )
        , _nextScale(boost::extents[num_states(*model)] )
        , _nextShape(boost::extents[num_states(*model)] )
	{
            _reset();
	}
    enum {_empty=1}; // Remove this when all functions use MemoizedFunc
    
			
protected:
    typedef 
    v2::Computation<
      const typename Algo::model_type
      , boost::mpl::set<
            v2::tag::forward
          , v2::tag::forward_begin
          , v2::tag::backward
          , v2::tag::backward_begin
          , v2::tag::gamma
          , v2::tag::xi // TODO: is xi used directly?
          , v2::tag::sigma_i_xi
          , v2::tag::sigma_li_xi
	>
    , fasta::seq_t
      > comp_xi_t;


    comp_xi_t _comp_xi;


public:
    /*
     * Call upon outside implementation of gamma (HSMM)
     */
    inline P gamma( size_t k, size_t i)
    {
        //return _gamma(typename gamma_type::arg_type(k,i));
	return boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::gamma(), typename comp_xi_t::template algo<v2::tag::gamma>::type::arg_type(k, i) ) );
    }
    
public:
    /*
     * Calculate sigma_i_gamma(k) (HSMM)
     * Note: Individial gamma values are memoized, but the sigma is not (TODO - check if this is worthwhile optimizing)
     */
    P Sigmai_gamma( size_t k )
    {
        P sigmai = 0.0;
        //if( k==0 )
        //{
        //    const typename Algo::probability_type gammaki = gamma( k, 0);
        //    std::cout<< "Sigmai_gamma(k="<< k<< ") = "<< gammaki<< std::endl;
        //    return gammaki;
        //}

        for( size_t i = 1; i< length(*_seq)+1; ++i )
        {
            const typename Algo::probability_type gammaki = gamma( k, i);
            //if( i == length(*_seq) )
            //    std::cout<< "Sigmai_gamma(k="<< k<< ", i="<< i<< ") = "<<gammaki<< std::endl;
            sigmai += gammaki;
        }
        //std::cout<< "Sigmai_gamma(k="<< k<< ") = "<< sigmai<< std::endl;
        return sigmai;
    }
    
public:
    /*
     * Calculate sigma_i_gamma_obs_s(k,s) (HSMM)
     * Note: Individial gamma values are memoized, but the sigma is not (TODO - check if this is worthwhile optimizing)
     */
    P Sigmai_gamma_observe_s( size_t k, size_t s )
    {
        P sigmai = 0.0;
        const typename Algo::symbol_type sym = _model->get_symbol(s);
        for( size_t i = 1; i< length(*_seq)+1; ++i ) // Only iterate over positions where symbols are emitted
        {
            if( (*_seq)[i-1] == sym )
            {
                const typename Algo::probability_type gammaki = gamma( k, i);
                sigmai += gammaki;
            }
        }
        return sigmai;
    }

public:
    /*
     * Calculate xi(k,l,i) (HSMM)
     * This is the probability of a k->l transition occuring at time i.
     * Defined (for HMMs) in Rabiner eq. (36)
     * Implemented as in Rabiner eq. (79)
     */
    typename comp_xi_t::template algo<v2::tag::xi>::type::result_type
    xi( size_t k, size_t l, size_t i )
    {
	return _comp_xi._apply(v2::tag::xi(), typename comp_xi_t::template algo<v2::tag::xi>::type::arg_type(k, l, i) );
    }
    /*
     * Calculate Sigma_i_xi(k,l) (HSMM)
     * This is the sum of probabilities of a k->l transition occuring at any time.
     * Described (for HMMs) in Rabiner eq. (39b), and for HSMMs in eq. (79)
     */
public:
    P Sigmai_xi( size_t k, size_t l)
    {
        //return 0.0;
        //std::cout<< "Sigmai_xi("<< k<< "->"<< l<< "): "<< std::endl;
	return boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::sigma_i_xi(), typename comp_xi_t::template algo<v2::tag::sigma_i_xi>::type::arg_type(k, l) ) );

        /*
        if( k==l )
            return 0.0;

        const size_t N = length(*_seq);
        assert(N>0);

        P sigma_i_xi = 0.0;

        //const P a_k_l = _model->a(k,l);

        for( size_t i = 1; i< N+1; ++i )
        {
            //const typename Algo::forward_algo_type::arg_type k_i(k, i);
            //const typename Algo::backward_begin_algo_type::arg_type l_i_1(l, i);
            //const logspace_t f((*_forward)(j_t), true);

            const P xi_k_l_i = boost::fusion::at_c<0>(xi(k, l, i));
            std::cout<< xi_k_l_i<< std::endl;

            //const P p1 = (*_forward)(k_i);
            //const P p2 = a_k_l;
            //const P p3 = (*_backward_begin)(l_i_1);

            //std::cout<< "xi("<< k<< "->"<< l<< ", i="<< i<< ") = "<< p1<< " + "<< p2<< " + "<< p3<< " = "<< p1+p2+p3<< std::endl;

            // Rabiner eq. (79)
            //const P xi_k_l_i =
            //    (*_forward)(k_i)
            //    + a_k_l
            //    + (*_backward_begin)(l_i_1);

            sigma_i_xi += exp(xi_k_l_i);
        }
        std::cout<< "sigma_i_xi("<< k<< "->"<< l<< ") = "<< sigma_i_xi<< std::endl;

        return sigma_i_xi;
        */
    }
    
public:
    /*
     * Calculate Sigma_l_xi(k,i) (HSMM)
     * This is the probablity of any transition occuring when leaving state k at time i,
     * in other words sigma_t_gamma_k
     * Described (for HMMs) in Rabiner eq. (39a)
     */
    P Sigmal_xi( size_t k)
    {
	return boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::sigma_li_xi(), typename comp_xi_t::template algo<v2::tag::sigma_li_xi>::type::arg_type(k) ) );
    }
    
    /*
     * Memoization support wrapper for a (HSMM)
     */
public:
    P nexta( size_t k, size_t l)
    {
        P val =  _nextA[k][l];
        if( val != _empty ) return val;
	
        _state_reset = false;
        _nextA[k][l] = val = _produce_nexta(k, l);
        return val;		
    }

    /*
     * Memoization support wrapper for e (HSMM)
     */
public:
    P nextb( size_t k, size_t s)
    {
        assert( !_model->is_silent(k) );
        P val = _nextB[k][s];
        if( val != _empty ) return val;
        
        _state_reset = false;
        
        if( !Algo::use_emissions_classes::value )
        {
            _nextB[k][s] = val = _produce_nextb(k,
                                                _model->get_symbol(s) /* _produce_nextb expects a symbol, not an index */ );
        }
        else
        {
            const size_t emclass = _model->GetEmissionsClass(k);
            
            val = _produce_nextb_class(emclass,
                                       _model->get_symbol(s) /* _produce_nextb expects a symbol, not an index */ );
            
            
            // Store the value for all members of this class
            typename Algo::model_type::emissions_class_members_iterator mem_it, mem_end;
            boost::tie( mem_it, mem_end ) = _model->GetEmissionsClassMembersRange(emclass);
            
            for( ; mem_it != mem_end; ++mem_it )
                _nextB[*mem_it][s] = val;
        }
	
        return val;		
    }
    
    /*
     * Memoization support wrapper for scale (HSMM)
     */
public:
    P nextscale( size_t k )
    {
        P val =  _nextScale[k];
        if( val != _empty ) return val;
	
        _state_reset = false;
        _nextScale[k] = val = _produce_nextscale(k);
        return val;		
    }

    /*
     * Memoization support wrapper for shape (HSMM)
     */
public:
    P nextshape( size_t k )
    {
        P val =  _nextShape[k];
        if( val != _empty ) return val;
	
        _state_reset = false;
        _nextShape[k] = val = _produce_nextshape(k);
        return val;		
    }

    /*
     * Calculate P(seq|model) (HSMM)
     * Rabiner eq. (70)
     */
public:
    P Pseq()
    {
        const P val = boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::forward(), typename comp_xi_t::template algo<v2::tag::forward>::type::arg_type(
                                                                 _model->GetTerminalState()
                                                                 , length(*_seq)+1
                                                                 ) ) );
        //std::cout<< "Pseq: "<< val<< std::endl;
        assert_range_probability_logspace(val);
        return val;

        //return _forward->calc();
    }
    
protected:
    void _reset()
    {
        {
            const size_t m_end = _nextA.shape()[0];
            const size_t n_end = _nextA.shape()[1];
            
            for( size_t m=0; m< m_end; ++m )
                for( size_t n=0; n< n_end; ++n )
                {	  _nextA[m][n] = _empty;	}
        }
        
        {
            const size_t m_end = _nextB.shape()[0];
            const size_t n_end = _nextB.shape()[1];
            
            for( size_t m=0; m< m_end; ++m )
                for( size_t n=0; n< n_end; ++n )
                {	  _nextB[m][n] = _empty;	}
        }
        
        {
            const size_t m_end = _nextScale.shape()[0];
            
            for( size_t m=0; m< m_end; ++m )
            {	  _nextScale[m] = _empty;	}
        }
        
        {
            const size_t m_end = _nextShape.shape()[0];
            
            for( size_t m=0; m< m_end; ++m )
            {	  _nextShape[m] = _empty;	}
        }
        
        //_forward->reset();
        //_forward_begin->reset();
        //_backward->reset();
        //_backward_begin->reset();
        //_gamma.reset();
        _state_reset = true;

        // TODO - Check this - isn't this redundant?
        _comp_xi._move_segment(typename comp_xi_t::template algo<v2::tag::forward>::type::data_level(), 0 );
        _comp_xi._move_segment(typename comp_xi_t::template algo<v2::tag::forward_begin>::type::data_level(), 0 );
        _comp_xi._move_segment(typename comp_xi_t::template algo<v2::tag::backward>::type::data_level(), 0 );
        _comp_xi._move_segment(typename comp_xi_t::template algo<v2::tag::backward_begin>::type::data_level(), 0 );
        _comp_xi._move_segment(typename comp_xi_t::template algo<v2::tag::gamma>::type::data_level(), 0 );
        _comp_xi._move_segment(typename comp_xi_t::template algo<v2::tag::xi>::type::data_level(), 0 );
    }
    
public:
    void reset()
    {
        //std::cout<< "HSMMBW::reset()"<< std::endl;
        //	if( _state_reset ) return; // TODO FIX THIS!!!
        //std::cout<< "HSMMBW::reset() *"<< std::endl;
        _reset();
    }
    
public:
    void calc_all()
    {
        // TODO - Re-implement or deprecate?
        assert(false);
        throw;
    }
    
    /*
     * Baum-Welch 'a' (transitions) estimation, HSMM, no emissions classes, single sequence
     * Rabiner eq. (79) - Overall equation
     * Note: This part of the impl. is bypassed by MutipleSequenceBW when multiple sequences are used, so is not used anymore in practice.
     */
protected:
    P _produce_nexta( size_t k, size_t l)
    {
        // probabilities from/to reserved states are fixed
        //if( _model->IsReservedState(k) ||
        //	_model->IsReservedState(l)   )
        //	return _model->a(k, l);
        const P sigmai_xi = Sigmai_xi( k, l );
        const P sigmai_gamma = Sigmal_xi( k );
        assert( sigmai_gamma - sigmai_xi > -1e-10 );
        
        /*
          if( k==0 )
          {
          std::cout<< "Debug: nexta(0, "<< l<< "): sigmai_xi="<< sigmai_xi<< ", sigmai_gamma="<< sigmai_gamma<< std::endl;
          }
        */
	
        if( sigmai_gamma < 1e-200 )
        {
            std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate 'a' for state="<< k<< ", t="<<  l<<std::endl;
            return _model->a(k, l);
        }
	
        const P nexta = log(sigmai_xi / sigmai_gamma);
        assert(nexta < 0.001 );
        return (nexta < 0.0) ? nexta : 0.0;
    }
    
protected:
    /*
     * Baum-Welch 'e' (emissions) estimation, HSMM, no emissions classes, single sequence
     * Rabiner eq. (80) - Overall equation
     * Note: This part of the impl. is bypassed by MutipleSequenceBW when multiple sequences are used, so is not used anymore in practice.
     */
    P _produce_nextb( size_t k, size_t s )
    {
        const P sigmai_gamma_obs_s = Sigmai_gamma_observe_s( k, s );
        const P sigmai_gamma = Sigmai_gamma( k );
        assert( sigmai_gamma - sigmai_gamma_obs_s > -1e-10  );
        assert( !Algo::use_emissions_classes::value );

        if( sigmai_gamma < 1e-200 )
        {
            std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate 'e' for state="<< k<< ", sym="<< s<< std::endl;
            return _model->e(k, s);
        }
        
        const P nextb = log(sigmai_gamma_obs_s / sigmai_gamma);
        assert(nextb < 0.001);
        
        if( nextb > -1e-3)
        {
            std::cout<< "Warning: About to return nextb= "<< exp(nextb)<< " for state="<< k<< ", symbol="<< s<< " (Sigma i gamma obs s="<< sigmai_gamma_obs_s<< ", sigma n gamma="<< sigmai_gamma<< ")"<< std::endl;
        }
        
        return (nextb < 0.0) ? nextb : 0.0;
    }
    
protected:
    /*
     * Baum-Welch 'e' (emissions) estimation, HSMM, with emissions classes, single sequence
     * Rabiner eq. (80) - Overall equation
     * Note: This part of the impl. is bypassed by MutipleSequenceBW when multiple sequences are used, so is not used anymore in practice.
     */
    P _produce_nextb_class( size_t k, size_t s )
    {
        assert( Algo::use_emissions_classes::value );
	
        P sigman_gamma_obs_s = 0.0;
        P sigman_gamma = 0.0;
        
        typename Algo::model_type::emissions_class_members_iterator mem_it, mem_end;
        boost::tie( mem_it, mem_end ) = _model->GetEmissionsClassMembersRange(k);
	
        for( ; mem_it != mem_end; ++mem_it )
        {
            const typename Algo::model_type::StateId mem = *mem_it;
            const P sigmai_gamma_obs_s = Sigmai_gamma_observe_s( mem, s );
            const P sigmai_gamma = Sigmai_gamma( mem );
            
            
            assert( sigmai_gamma - sigmai_gamma_obs_s > -1e-10 );
            
            sigman_gamma_obs_s += sigmai_gamma_obs_s;
            sigman_gamma += sigmai_gamma;			
        }
        
        //sigman_gamma_obs_s /= double(sigman_count);
        //sigman_gamma /= double(sigman_count);
        
        const P nextb = log(sigman_gamma_obs_s / sigman_gamma);
        assert( nextb < 0.001 );
        if( nextb > -1e-6)
        {
            std::cout<< "Warning: About to return nextb= "<< exp(nextb)<< " for state="<< k<< ", symbol="<< s<< " (Sigma n gamma obs s="<< sigman_gamma_obs_s<< ", sigma n gamma="<< sigman_gamma<< ")"<< std::endl;
        }
        return (nextb <= 0.0) ? nextb : 0.0;
    }
    
    /*
     * Baum-Welch scale estimation (HSMM), single sequence
     * TODO - Add ref
     */
protected:
    P _produce_nextscale( size_t j )
    {
        // Numerator
        typedef LogspaceDouble<> logspace_t;
        //std::cout<< "nextscale j="<< j << std::endl;
	
        P sigmat_alpha_beta = 0.0;
        for( size_t t = 1; t< length(*_seq)+1; ++t ) // Only iterate over positions where symbols are emitted
        {
            //const typename Algo::forward_algo_type::arg_type j_t(j, t);
            //const typename Algo::forward_algo_type::arg_type j_t_1(j, t-1);

            const logspace_t f( boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::forward(),  typename comp_xi_t::template algo<v2::tag::forward >::type::arg_type(j, t) ) ), true);
            const logspace_t b( boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::backward(), typename comp_xi_t::template algo<v2::tag::backward>::type::arg_type(j, t) ) ), true);
            // v1 impl::
            //const logspace_t f((*_forward)(j_t), true);
            //const logspace_t b((*_backward)(j_t), true);
            
            //std::cout<< boost::format("t=%d, f=%e, b=%e, f*b=%e") % t % f % b % (f*b)<< std::endl;
            
            sigmat_alpha_beta += f * b;
        }
        const P exp1 = sigmat_alpha_beta * _model->GetMu(j);
        //std::cout<< boost::format("sigmat_alpha_beta=%e, Mu=%e, numerator=%e") % sigmat_alpha_beta % _model->GetMu(j) % exp1 << std::endl;
        const P numerator = exp1;
        
        // Denominator
        P sigma_t = 0.0;
        for( size_t t = 1; t< length(*_seq)+1; ++t ) // Only iterate over positions where symbols are emitted
        {
            const typename Algo::forward_algo_type::arg_type j_t(j, t);
            const typename Algo::forward_algo_type::arg_type j_t_1(j, t-1);

            const logspace_t backward( boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::backward(),  typename comp_xi_t::template algo<v2::tag::backward >::type::arg_type(j, t) ) ), true);
            //const logspace_t backward( (*_backward)(j_t), true);
            
            P sigma_tau = 0.0;
            for( size_t tau = 1; tau <= t && tau <= Algo::max_duration; ++tau )
            {
                const logspace_t d_j_tau( _model->p(j, tau), true);
                
                P sigma_i = 0.0;
                for( size_t i=0; i< num_states(*_model); ++i )
                {
                    if( i == j ) continue;
                    
                    const logspace_t forward( boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::forward(),  typename comp_xi_t::template algo<v2::tag::forward >::type::arg_type(i, t-tau) ) ), true);
                    //const logspace_t forward( (*_forward)(typename Algo::forward_algo_type::arg_type( i, t-tau ) ), true );
                    const logspace_t aij( _model->a(i, j), true);
                    
                    logspace_t pi_b = 1.0;
                    for( size_t theta = 1; theta <= tau; ++theta )
                    {
                        pi_b *= logspace_t(_model->e(j,
                                                     (*_seq)[t - tau + theta - 1] ), true);
                    }
                    
                    const logspace_t exp2 = forward * aij * d_j_tau * pi_b * backward ;
                    //std::cout<< boost::format("[%d=>%d tau=%d t=%d] forward=%e * aij=%e * d=%e * pi_b=%e * backward=%e  = %e") % j % i % tau % t % forward % aij % d_j_tau % pi_b % backward % exp2<< std::endl;
                    sigma_i += exp2;
                }
                //std::cout<< boost::format("sigma_tau += sigma_i=%8.8e * tau=%d") % sigma_i % tau<< std::endl;
                sigma_tau += sigma_i * tau;
            }
            sigma_t += sigma_tau;
        }
        const P denominator = sigma_t;
	
        //std::cout<< boost::format("%d %e %e %e") % j % numerator % denominator % (numerator/denominator)<< std::endl;
        return numerator/denominator;
    }
    
    /*
     * Baum-Welch shape estimation (HSMM), single sequence
     * TODO - Add ref
     */
protected:
    P _produce_nextshape( size_t j )
    {
        // Calculate the value of the b constant (Levinson eq. (21))
        // Numerator
        P sigma_t = 0.0;
        for( size_t t = 1; t < length(*_seq)+1; ++t )
        {
            //const typename Algo::forward_algo_type::arg_type j_t(j, t);
            const logspace_t backward( boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::backward(),  typename comp_xi_t::template algo<v2::tag::backward >::type::arg_type(j, t) ) ), true);
            //const logspace_t backward( (*_backward)(j_t), true);
            
            P sigma_tau = 0.0;
            for( size_t tau = 1; tau <= t && tau <= Algo::max_duration; ++tau )
            {
                const P log_eta_tau = log(_model->GetEta(j) * tau);
                const logspace_t d_j_tau( _model->p(j, tau), true );
		
                P sigma_i = 0.0;
                for( size_t i=0; i< num_states(*_model); ++i )
                {
                    if( i == j ) continue;
                    const logspace_t forward( boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::forward(),  typename comp_xi_t::template algo<v2::tag::forward >::type::arg_type(i, t-tau) ) ), true);
                    //const logspace_t forward( (*_forward)(typename Algo::forward_algo_type::arg_type(i, t-tau) ), true );
                    const logspace_t aij( _model->a(i, j), true );
                    
                    logspace_t pi_b = 1.0;
                    for( size_t theta = 1; theta <= tau; ++theta )
                    {
                        pi_b *= logspace_t(_model->e(j,
                                                     (*_seq)[t - tau + theta - 1] ), true);
                    }
                    const logspace_t exp1 = forward * aij * d_j_tau * pi_b * backward ;
                    sigma_i += exp1;
                }
		
                const P exp2 = log_eta_tau * sigma_i;
                sigma_tau += exp2;
            }
            sigma_t += sigma_tau;
        }
        const P numerator = sigma_t;
	
        // Denominator
        P sigmat_alpha_beta = 0.0;
        for( size_t t = 1; t< length(*_seq)+1; ++t ) // Only iterate over positions where symbols are emitted
        {
            //const typename Algo::forward_algo_type::arg_type j_t(j, t);

            const logspace_t f( boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::forward(),   typename comp_xi_t::template algo<v2::tag::forward >::type::arg_type(j, t) ) ), true);
            const logspace_t b( boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::backward(),  typename comp_xi_t::template algo<v2::tag::backward>::type::arg_type(j, t) ) ), true);
            //const logspace_t f((*_forward)(j_t), true);
            //const logspace_t b((*_backward)(j_t), true);
            
            //std::cout<< boost::format("t=%d, f=%e, b=%e, f*b=%e") % t % f % b % (f*b)<< std::endl;
            
            sigmat_alpha_beta += f * b;
        }
        //std::cout<< boost::format("sigmat_alpha_beta=%e, Mu=%e, numerator=%e") % sigmat_alpha_beta % _model->GetMu(j) % exp1 << std::endl;
        const P denominator = sigmat_alpha_beta;
        
        const P b = numerator/denominator;
        
        // Estimate mu using Newton-Raphson method:
        P mu = _model->GetMu(j);
        P delta = 1.0;
        size_t iter = 1;
        int err = 0;
	
        while( (delta > 1e-8) && (iter < 100) )
        {
            const P digamma_mu  = boost::math::digamma(mu);
            const P trigamma_mu = trigam(mu, &err);			assert(err == 0);
            
            const P nextmu = mu - ((digamma_mu - b) / trigamma_mu);
            
            delta = fabs( mu - nextmu );
            //std::cout<< boost::format("iter=%d, mu=%e, nextmu=%e, delta=%e") % iter % mu % nextmu % delta << std::endl;
            
            mu = nextmu;
            if( mu < 1e-10 )
            {
                std::cout<< "Warning: mu<min at _produce_nextshape(j="<< j<< "); iter="<< iter<< ", mu="<< mu<< ", nextmu="<< nextmu<< ", digamma(mu)="<< digamma_mu<< ", trigamma(mu)="<< trigamma_mu<< std::endl;
                mu = 1e-10; // TODO - check if this can be removed
                // NOTE - This probably doesn't matter much, since it appears to only happen for few sequences and so has very little
                //        effect on the average across all sequences.
                //        To investigate why this happens in the first place, I saved a test-case under ./test_case_negative_mu
            }
            
            ++iter;
        }
	
        return mu;
    }
    
protected:
    const boost::shared_ptr<const typename Algo::model_type> _model;
    const boost::shared_ptr<const typename Algo::sequence_type> _seq;
        
    typedef boost::multi_array<P,2> nextA_type;
    nextA_type _nextA;
    
    typedef boost::multi_array<P,2> nextB_type;
    nextB_type _nextB;
    
    typedef boost::multi_array<P,1> nextScale_type;
    nextScale_type _nextScale;
    
    typedef boost::multi_array<P,1> nextShape_type;
    nextShape_type _nextShape;
    
    
    bool _state_reset;

    
public:
    void debug_print() const
    {
        std::cout<< "HSMMBaumWelchAlgorithm<Algo>"<< std::endl;

        std::cout<< " seq: (len="<< length(*_seq)<< ")"<< std::endl;
        std::cout<< *_seq<< std::endl;

        //std::cout<< " gamma: ";
        //detail::debug_print(_gamma);
        std::cout<< " nextA: ";
        ::debug_print(_nextA);
        std::cout<< " nextB: ";
        ::debug_print(_nextB);
        std::cout<< " nextScale: ";
        ::debug_print(_nextScale);
        std::cout<< " nextShape: ";
        ::debug_print(_nextShape);

    }

    bool test()
    {
        tests::test_38<Algo>(*this, *_model, *_seq);

        return true;
    }

    P peek_forward( state_t j, time_t t )
    { 
        return boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::forward(),   typename comp_xi_t::template algo<v2::tag::forward >::type::arg_type(j, t) ) );
    }
    P peek_forward_begin( state_t j, time_t t )
    { 
        return boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::forward_begin(),   typename comp_xi_t::template algo<v2::tag::forward_begin >::type::arg_type(j, t) ) );
    }
    P peek_backward( state_t j, time_t t )
    { 
        return boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::backward(),   typename comp_xi_t::template algo<v2::tag::backward >::type::arg_type(j, t) ) );
    }
    P peek_backward_begin( state_t j, time_t t )
    { 
        return boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::backward_begin(),   typename comp_xi_t::template algo<v2::tag::backward_begin >::type::arg_type(j, t) ) );
    }
    P peek_gamma( state_t j, time_t t )
    { 
        return boost::fusion::at_c<0>(_comp_xi._apply(v2::tag::gamma(),   typename comp_xi_t::template algo<v2::tag::gamma >::type::arg_type(j, t) ) );
    }
    P peek_sigma_t_xi( state_t i, state_t j )
    { 
        return Sigmai_xi( i, j );
    }

    
public:
    typedef typename boost::mpl::false_ gamma_in_logspace;
    
};

template<class Algo>
void debug_print( const HSMMBaumWelchAlgorithm<Algo>& obj )
{
    obj.debug_print();
}

} // namespace baumwelch_algo


namespace baumwelch_algo
{
	
template<typename Algo>
class MultipleSequenceBW
{
protected:
    typedef typename Algo::probability_type P;
    enum {_empty=-1}; // Remove this when all memoized funcs use MemoizedFunc
    using state_t = typename Algo::model_type::StateId;
	
public:
    
    template<typename It>
    MultipleSequenceBW(boost::shared_ptr<const typename Algo::model_type> model,
                       It first, It last)
        : _model(model)
	, _nextA(boost::extents[num_states(*model)][num_states(*model)] )
	, _nextB(boost::extents[num_states(*model)][num_symbols(*model)] )
	, _state_reset(true)
    {
        for( It it = first; it != last; ++it )
        {
            // Retrieve the sequence object
            boost::shared_ptr<const typename Algo::sequence_type> seq( new typename Algo::sequence_type( fasta::get_seq(*it) ) );

            // Store the sequence
            _seq.push_back( seq );
            // Create a baum-welch algo for it
            //boost::shared_ptr<typename Algo::baumwelch_algo_type> bw(new typename Algo::baumwelch_algo_type( model, f, fb, b, bb, seq ) );
            boost::shared_ptr<typename Algo::baumwelch_algo_type> bw(new typename Algo::baumwelch_algo_type( model, seq ) );
            _baumwelch.push_back( bw );
            
        }
        _reset();
    }
    
public:
    /*
     * Baum-Welch 'a' (transitions) estimation, Multiple-sequence adapter for HMM/HSMM
     * Rabiner eqs. (40a), (79), (109)
     */
    P nexta(size_t k, size_t l)
    {
        // probabilities from/to reserved states are fixed
        //if( _model->IsReservedState(k) ||
        //	_model->IsReservedState(l)    )
        //	return _model->a(k, l);

        if( boost::core::is_same<
                typename v2::model::model_traits<typename Algo::model_type>::hidden_process
                , v2::model::semi_markov_chain>::value )
        {
            if( k==l )
                return _model->a(k, l);
        }
        
        // Rabiner Eq. (109) (p. 273)
        P sigman_xi = 0.0;
        P sigman_gamma = 0.0;
        
        _state_reset = false;

        const size_t N = _seq.size();
        for( size_t n = 0; n < N; ++n )
        {
            _baumwelch[n]->reset();
            const double Pn = exp(_baumwelch[n]->Pseq());
            //std::cout<< "P(On) = "<< Pn<< std::endl;

            const P sigmai_xi    = _baumwelch[n]->Sigmai_xi(k, l);
            const P sigmai_gamma = _baumwelch[n]->Sigmal_xi(k);
            // Both in probability space
            
            if( sigmai_gamma - sigmai_xi <= -1e-10 )
            {
                std::cout<< "nexta("<< k<< "->"<< l<< ") :: sigma_i_xi="<< sigmai_xi<< " sigma_l_xi= "<< sigmai_gamma<< std::endl;
            }
            assert( sigmai_gamma - sigmai_xi > -1e-10 ); // sigmai_gamma >= sigmai_xi

            //std::cout<< "-------------------------------- ("<< k<< "->"<< l<< ", "<< n<< " / "<< N<< ")---------------------------------"<< std::endl;
            //_baumwelch[n]->calc_all();
            //calc_all(n);
            //debug_print(n);
            //std::cout<< "------------------------------------/-----------------------------"<< std::endl;


            //std::cout<< "BW for seq "<< n<< ": sigmai_xi="<< sigmai_xi<< ", sigmal_xi="<< sigmai_gamma<< std::endl;
            
            sigman_xi += sigmai_xi / Pn ;
            sigman_gamma += sigmai_gamma / Pn ;
        }
        if( sigman_gamma < 1e-200 )
        {
            std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate 'a' for transition "<< k<< "->"<< l<<std::endl;
            return _model->a(k, l);
        }

        std::cout<< "BW total for transition "<< k<< "->"<< l<< ": sigman_xi="<< sigman_xi<< ", sigmal_xi="<< sigman_gamma<< std::endl;
	
        const P nexta = log(sigman_xi / sigman_gamma);
        assert( nexta < 0.001 );
        return (nexta < 0.0) ? nexta : 0.0;
    }

public:
    P nextb(size_t k, size_t s)
    {
        assert( !_model->is_silent(k) );
        P val = _nextB[k][s];
        if( val != _empty ) return val;
        
        _state_reset = false;
        
        if( !Algo::use_emissions_classes::value )
        {
            
            //std::cout<< "(k="<< k<< ", s="<< s<< ") ";
            _nextB[k][s] = val = _produce_nextb(k,
                                                s );
            
            if( val > -1e-6)
            {
                std::cout<< "Warning: About to return nextb==1 ("<< exp(val)<< ") for state="<< k<< ", symbol="<< s<< " (noclasses)"<< std::endl;
            }
            
        }
        else
        {
            const size_t emclass = _model->GetEmissionsClass(k);
            
            val = _produce_nextb_class(emclass,
                                       s );
            
            if( val > -1e-6)
            {
                std::cout<< "Warning: About to return nextb==1 ("<< exp(val)<< ") for state="<< k<< ", symbol="<< s<< " (classes)"<< std::endl;
            }
            
            // Store the value for all members of this class
            typename Algo::model_type::emissions_class_members_iterator mem_it, mem_end;
            boost::tie( mem_it, mem_end ) = _model->GetEmissionsClassMembersRange(emclass);
            
            for( ; mem_it != mem_end; ++mem_it )
                _nextB[*mem_it][s] = val;
        }
        
        return val;
    }
  
public:
    P nextscale( size_t k )
    {
        assert(!_model->IsReservedState(k));
        
        _state_reset = false;
        
        P sigman_scale = 0.0;
	
        const size_t N = _seq.size();
        for( size_t n = 0; n < N; ++n )
        {
            sigman_scale += _baumwelch[n]->nextscale(k);
        }
        return sigman_scale / P(N);
    }
    
public:
    P nextshape( size_t k )
    {
        assert(!_model->IsReservedState(k));
        
        _state_reset = false;
        
        P sigman_shape = 0.0;
	
        const size_t N = _seq.size();
        for( size_t n = 0; n < N; ++n )
        {
            sigman_shape += _baumwelch[n]->nextshape(k);
        }
        return sigman_shape / P(N);
    }
    
protected:
    /*
     * Baum-Welch 'e' (emissions) estimation, Multiple-sequence adapter for HMM/HSMM, no emissions classes
     * Rabiner eq. (110) (p. 273)
     */
    P _produce_nextb(size_t k, size_t s)
    {
        assert( !Algo::use_emissions_classes::value );
        assert( !_model->is_silent(k) );

        P sigman_gamma_obs_s = 0.0;
        P sigman_gamma = 0.0;
        
        _state_reset = false;

        // Sum parameter contributions over all sequences
        const size_t N = _seq.size();

        // Iterate over training sequences
        for( size_t n = 0; n < N; ++n )
        {
            _baumwelch[n]->reset();
            const P sigmai_gamma_obs_s = _baumwelch[n]->Sigmai_gamma_observe_s(k, s );
            const P sigmai_gamma = _baumwelch[n]->Sigmai_gamma(k);
            const P lnPn = _baumwelch[n]->Pseq();
            const P weight = exp(-lnPn); // lnPn = log(P) => 1/exp(lnPn) = exp(-lnPn)

            //if( s==4 )
            //    std::cout<< "n="<< n<< ", s="<< s<< "): "<< sigmai_gamma_obs_s<< " / "<< sigmai_gamma<< " (weight="<< weight<< ")"<< std::endl;
            
            if( std::isnan(sigmai_gamma_obs_s) || std::isnan(sigmai_gamma) )
            {
                std::cout<< "Warning: Got NaN during 'e' calculation for k="<< k<< ", s="<< s<< ", n="<< n<< ", sigmai_gamma_obs_s="<< sigmai_gamma_obs_s<< ", sigmai_gamma="<< sigmai_gamma<< std::endl;
                return -std::numeric_limits<P>::max();
            } 

            assert( sigmai_gamma - sigmai_gamma_obs_s > -1e-10 );
            
            if( sigmai_gamma < 1e-20 )
            {
                //std::cout<< "About to generate NaN (k="<< k<< ", s="<< s<< ", n="<< n<< ")..."<< std::endl;
                continue;  // Ignore this sequence for the purpose of calculating this value...
            }
            
            sigman_gamma_obs_s += sigmai_gamma_obs_s * weight;
            sigman_gamma       += sigmai_gamma       * weight;
        }
        
        if( sigman_gamma < 1e-200 )
        {
            std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate 'b' for state="<< k<< ", sym="<< s<<std::endl;
            return _model->e(k, s );
        }
        
        if( sigman_gamma_obs_s < 0.0 ) sigman_gamma_obs_s = 0.0;

        //if( s==4 )
        //    std::cout<< "(k="<< k<< ", s="<< s<< ") sigman_gamma="<< sigman_gamma<< " sigman_gamma_obs_s="<< sigman_gamma_obs_s<< std::endl;
        
        const P nextb = log(sigman_gamma_obs_s / sigman_gamma);
        if( nextb >= 1e-12 )
        {
            std::cout<< "Error: Sigma n gamma obs s > Sigma n gamma (Sigma n gamma obs s="<< sigman_gamma_obs_s<< ", Sigma n gamma="<< sigman_gamma<< std::endl;
            assert( nextb < 1e-12 );
        }
        if( nextb >= -1e-6)
        {
            std::cout<< "Warning: About to return nextb= "<< exp(nextb)<< " for state="<< k<< ", symbol="<< s<< " (Sigma n gamma obs s="<< sigman_gamma_obs_s<< ", sigma n gamma="<< sigman_gamma<< ")"<< std::endl;
        }
        return (nextb <= 0.0) ? nextb : 0.0;
    }
    
protected:
    /*
     * Baum-Welch 'e' (emissions) estimation, Multiple-sequence adapter for HMM/HSMM, with emissions classes
     * Rabiner eq. (110) (p. 273)
     */
    P _produce_nextb_class(size_t k, size_t s)
    {
        
        assert( Algo::use_emissions_classes::value );
        assert( !_model->is_silent(k) );
        
        P sigman_gamma_obs_s = 0.0;
        P sigman_gamma = 0.0;
        
        _state_reset = false;
        
        const bool d = false;
        //if( s==1 ) d = true;
        
        //const typename Algo::symbol_type sym = _model->get_symbol(s);
        
        typename Algo::model_type::emissions_class_members_iterator mem_it, mem_end;
        boost::tie( mem_it, mem_end ) = _model->GetEmissionsClassMembersRange(k);

        size_t members = 0;
        const size_t N = _seq.size();
        
        // Iterate over emissions class members
        for( ; mem_it != mem_end; ++mem_it )
        {
            const typename Algo::model_type::StateId mem = *mem_it;
            
            // Iterate over training sequences
            for( size_t n = 0; n < N; ++n )
            {
                const P sigmai_gamma_obs_s = _baumwelch[n]->Sigmai_gamma_observe_s( mem, s );
                const P sigmai_gamma = _baumwelch[n]->Sigmai_gamma( mem );
                
                if( std::isnan(sigmai_gamma_obs_s) || std::isnan(sigmai_gamma) )
                {
                    std::cout<< "Warning: Got NaN during 'e' calculation for k="<< k<< ", s="<< s<< ", n="<< n<< ", sigmai_gamma_obs_s="<< sigmai_gamma_obs_s<< ", sigmai_gamma="<< sigmai_gamma<< std::endl;
                    return -std::numeric_limits<P>::max();
                } 
                
                assert( sigmai_gamma - sigmai_gamma_obs_s > -1e-10 ); // sigmai_gamma >= sigmai_gamma_obs_n
                
                //std::cout<< std::endl<< "(   k="<< k<< ", sym="<< s<<", n="<< n<< ", sigmai_gamma_obs_s="<< sigmai_gamma_obs_s<< ", sigmai_gamma="<< sigmai_gamma<< std::endl;
                
                if( sigmai_gamma < 1e-20 )
                {
                    //std::cout<< "About to generate NaN (k="<< k<< ", s="<< s<< ", n="<< n<< ")..."<< std::endl;
                    continue;  // Ignore this sequence for the purpose of calculating this value...
                }
                //sigman_gamma_obs_s += sigmai_gamma_obs_s;
                //sigman_gamma += sigmai_gamma;
                sigman_gamma_obs_s += sigmai_gamma_obs_s/sigmai_gamma; //todo fix div by 0
                sigman_gamma += 1.0; // Weigh member states equally. TODO: Consider changing this
            }
            
            ++members;
        }

        assert(!std::isnan(sigman_gamma_obs_s));
        assert(!std::isnan(sigman_gamma));
        
        //sigman_gamma_obs_s /= double(sigman_count);
        //sigman_gamma /= double(sigman_count);
        
        if( sigman_gamma < 1e-200 )
        {
            std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate 'e' for state="<< k<< ", sym="<< s<< std::endl;
            return _model->e(k, s );
        }
        if(d) 
        {
            std::cout<< "k="<< k<< ", sym="<< s<< ", members="<< members<< ", sigman_gamma_obs_s="<< sigman_gamma_obs_s<< ", sigman_gamma="<< sigman_gamma<< ", N="<< N<< std::endl;
        }
        
        if( sigman_gamma_obs_s < 0.0 ) sigman_gamma_obs_s = 0.0;

        //std::cout<< "(k="<< k<< ", s="<< s<< ") sigman_gamma="<< sigman_gamma<< " sigman_gamma_obs_s="<< sigman_gamma_obs_s<< std::endl;
        
        const P nextb = log(sigman_gamma_obs_s / sigman_gamma);
        if( nextb >= 1e-12 )
        {
            std::cout<< "Error: Sigma n gamma obs s > Sigma n gamma (k="<< k<< ", s="<< s<< ", Sigma n gamma obs s="<< sigman_gamma_obs_s<< ", Sigma n gamma="<< sigman_gamma<< std::endl;
            assert( nextb < 1e-12 );
        }
        if( nextb > -1e-2 || std::isnan(nextb) )
        {
            std::cout<< "Error: Calculated invalid 'e' ("<< exp(nextb)<< ") for (k="<< k<< ", s="<< s<< ", Sigma n gamma obs s="<< sigman_gamma_obs_s<< ", Sigma n gamma="<< sigman_gamma<< std::endl;
        }
        else std::cout<< "(e="<< exp(nextb)<< ", k="<< k<< ", s="<< s<< ") ";
        return (nextb <= 0.0) ? nextb : 0.0;
    }
    
public:
    P Pseq(size_t n)
    {
        return _baumwelch[n]->Pseq();
    }

public:
    P Pseq()
    {
        // Eq. (107)
        P sigmai_Pi = 0.0;
        
        //typename _baumwelch_t::iterator it, it_end;
        //it = _baumwelch.begin();
        //it_end = _baumwelch.end();
        //for( ; it != it_end; ++it ) Pi += (*it)->Pseq();
        //return Pi;

        for( size_t i = 0; i< _seq.size(); ++i )
            sigmai_Pi += Pseq(i);

        return sigmai_Pi;
    }
	

protected:
    void _reset()
    {
        {
            typename _baumwelch_t::iterator it, it_end;
            it = _baumwelch.begin();
            it_end = _baumwelch.end();
            for( ; it != it_end; ++it ) (*it)->reset();
        }
        
        {
            const size_t m_end = _nextA.shape()[0];
            const size_t n_end = _nextA.shape()[1];
            
            for( size_t m=0; m< m_end; ++m )
                for( size_t n=0; n< n_end; ++n )
                {       _nextA[m][n] = _empty;        }
        }
        
        {
            const size_t m_end = _nextB.shape()[0];
            const size_t n_end = _nextB.shape()[1];
            
            for( size_t m=0; m< m_end; ++m )
                for( size_t n=0; n< n_end; ++n )
                {         _nextB[m][n] = _empty;        }
        }
        
        _state_reset = true;
    }
  
public:
    void reset()
    {
        // std::cout<< "MultipleBW::reset()"<< std::endl;
        if( _state_reset ) return;
        //std::cout<< "MultipleBW::reset() *"<< std::endl;
        _reset();
    }

public:
    void calc_all(size_t n)
    {
        _baumwelch[n]->calc_all();
        //_forward[n]->calc_all();
        //_forward_begin[n]->calc_all();
        //_backward[n]->calc_all();
        //_backward_begin[n]->calc_all();
    }
    
public:
    void test_28()
    {
        const size_t N = _seq.size();
        for( size_t n = 0; n < N; ++n )
        {
            tests::test_28<Algo>( *(_baumwelch[n]), *_model, *(_seq[n]) );
        }
    }
    
    void test_38()
    {
        // TODO - Restore this test
        //const size_t N = _seq.size();
        //for( size_t n = 0; n < N; ++n )
        //{
        //    tests::test_38<Algo>( *(_baumwelch[n]), *_model, *(_seq[n]) );
        //}
    }
    void test_40c()
    {
        // TODO - Restore this test
        //const size_t N = _seq.size();
        //for( size_t n = 0; n < N; ++n )
        //{
        //    tests::test_40c<Algo>( *(_baumwelch[n]), *_model );
        //}
    }
    
    void test_D59_forw_back()
    {
        // TODO - Restore this test
        //const size_t N = _seq.size();
        //for( size_t n = 0; n < N; ++n )
        //{
        //    ::tests::test_D59_forw_back<Algo>( *(_forward[n]), *(_backward_begin[n]) );
        //}
    }
    
public:
    void debug_print() const
    {
        std::cout<< "MultipleSequenceBW<Algo>"<< std::endl;
        std::cout<< " nextA: "<< std::endl;
        ::debug_print(_nextA);
        std::cout<< " nextB: "<< std::endl;
        ::debug_print(_nextB);
        
        std::cout<< "BW for sequences:"<< std::endl;
        const size_t N = _seq.size();
        for( size_t n = 0; n < N; ++n )
        {
            _baumwelch[n]->debug_print();
        }
    }

public:
    void debug_print(size_t n) const
    {
        std::cout<< "MultipleSequenceBW<Algo>"<< std::endl;
        std::cout<< " nextA: "<< std::endl;
        ::debug_print(_nextA);
        std::cout<< " nextB: "<< std::endl;
        ::debug_print(_nextB);
        
        std::cout<< "BW for sequence "<< n<< ":"<< std::endl;

        _baumwelch[n]->debug_print();
    }

public:
    P peek_forward(        size_t n, state_t i, time_t t ) { return _baumwelch[n]->peek_forward(        i,t); }
    P peek_forward_begin(  size_t n, state_t i, time_t t ) { return _baumwelch[n]->peek_forward_begin(  i,t); }
    P peek_backward(       size_t n, state_t i, time_t t ) { return _baumwelch[n]->peek_backward(       i,t); }
    P peek_backward_begin( size_t n, state_t i, time_t t ) { return _baumwelch[n]->peek_backward_begin( i,t); }
    P peek_gamma(          size_t n, state_t i, time_t t ) { return _baumwelch[n]->peek_gamma(          i,t); }
    P peek_sigma_t_xi(     size_t n, state_t i, state_t j) { return _baumwelch[n]->peek_sigma_t_xi(     i,j); }

public:
    
	
    const boost::shared_ptr<const typename Algo::model_type> _model;
    
    typedef std::vector< boost::shared_ptr<typename Algo::baumwelch_algo_type> > _baumwelch_t;
    _baumwelch_t _baumwelch;
    
    typedef std::vector< boost::shared_ptr<const typename Algo::sequence_type> > _seq_t;
    _seq_t _seq;	
    
    typedef boost::multi_array<P,2> nextA_type;
    nextA_type _nextA;
    
    typedef boost::multi_array<P,2> nextB_type;
    nextB_type _nextB;
    
    bool _state_reset;
};



template<typename Algo, typename Reestimator>
class EM
{
public:
    using state_t = typename Algo::model_type::StateId;
    using P = typename Algo::probability_type;

public:
    EM( typename Algo::model_type& model, boost::shared_ptr<Reestimator> est)	
        : _model(&model, null_deleter() )
        , _est(est)
        , _nextA(boost::extents[num_states(*_model)][num_states(*_model)] )
        , _nextB(boost::extents[num_states(*_model)][num_symbols(*_model)] )
        , _nextScale(boost::extents[num_states(*_model)] )
        , _nextShape(boost::extents[num_states(*_model)] )
    {
        _reset();
    }

public:
    // Constructor for python export
    EM( typename Algo::model_type& model, const std::string fasta_path )
        : _model(&model, null_deleter() )
        , _nextA(boost::extents[num_states(*_model)][num_states(*_model)] )
        , _nextB(boost::extents[num_states(*_model)][num_symbols(*_model)] )
        , _nextScale(boost::extents[num_states(*_model)] )
        , _nextShape(boost::extents[num_states(*_model)] )
    {
        fasta::read_fasta(fasta_path, _seqs);
        
        _est = boost::shared_ptr<Reestimator>( new Reestimator( _model, _seqs.begin(), _seqs.end() ) );
        
        _reset();
    }

public:
    // Constructor for python export
    EM( typename Algo::model_type& model, const std::string fasta_path, const boost::python::list& indexes )
        : _model(&model, null_deleter() )
        , _nextA(boost::extents[num_states(*_model)][num_states(*_model)] )
        , _nextB(boost::extents[num_states(*_model)][num_symbols(*_model)] )
        , _nextScale(boost::extents[num_states(*_model)] )
        , _nextShape(boost::extents[num_states(*_model)] )
    {
        std::vector<size_t> vec;
        python_to_stl_vector<size_t, const boost::python::list>(indexes, vec);
	
        fasta::read_fasta(fasta_path, _seqs, vec.begin(), vec.end());
        
        _est = boost::shared_ptr<Reestimator>( new Reestimator( _model, _seqs.begin(), _seqs.end() ) );
        
        _reset();
    }
    
public:
    P reestimate_a()
    {
        _reset(); // clear our buffers
        
        // Use the estimator to calculated improved values for all parameters
	
        for( size_t k=0; k< num_states(*_model); ++k )
        {
	    if( k == _model->GetTerminalState()) continue;
	    
	    for( size_t l=0; l< num_states(*_model); ++l )
	    {
                //if( l == _model->GetTerminalState() ) continue;
                _nextA[k][l] = _est->nexta(k, l);
	    }
        }
        // Update the model with the improved params (updating the model must be done at once)
        for( size_t k=0; k< num_states(*_model); ++k )
        {
            if( k == _model->GetTerminalState()) continue;
	    
            for( size_t l=0; l< num_states(*_model); ++l )
            {
		//if( l == _model->GetTerminalState() ) continue;
		_model->SetTransitionLogspace( k, l, _nextA[k][l] );
            }
        }
        return _est->Pseq();
    }
    
    P reestimate_b()
    {
        _reset(); // clear our buffers
        
        // Use the estimator to calculated improved values for all parameters
        for( size_t k=0; k< num_states(*_model); ++k )
        {
            if( _model->is_silent(k) ) continue;
            
            for( size_t sidx=0; sidx< num_symbols(*_model); ++sidx )
            {
                _nextB[k][sidx] = _est->nextb(k, sidx);
            }
        }
        // Update the model with the improved params (updating the model must be done at once)
        for( size_t k=0; k< num_states(*_model); ++k )
        {
            if( _model->is_silent(k) ) continue;
            
            for( size_t sidx=0; sidx< num_symbols(*_model); ++sidx )
            {
                const P newp = _nextB[k][sidx];
                //#ifdef EXCESSIVE_DEBUG_PRINTING
                if( !isvalidprobability( exp( newp )) )
                {
                    std::cout<< "ERROR: Got invalid probability value "<< exp(newp)<< " while calculating 'e' for state="<< k<< ", symbol="<< sidx<< std::endl;
                }
                if( newp > -1e-6 )
                {
                    std::cout<< "WARNING: Got probability of 1 ("<< exp(newp)<< ") while calculating 'e' for state="<< k<< ", symbol="<< sidx<< std::endl;
                }
                //#endif
                assert( isvalidprobability( exp( newp )) );
                
                _model->SetEmissionProbabilityLogspace( k,
                                                        _model->get_symbol(sidx),
                                                        newp );
            }
        }
        return _est->Pseq();
    }
    
    // TODO - Find out why enable_if doesn't work in this case
    P reestimate_scale()
    {
        assert(Algo::model_type::is_explicit_duration_model::value);
        _reset(); // clear our buffers
        
        // Use the estimator to calculated improved values for all parameters
        for( size_t k=0; k< num_states(*_model); ++k )
        {
            if( _model->IsReservedState(k) ) continue;
            
            _nextScale[k] = _est->nextscale(k);
        }
        // Update the model with the improved params (updating the model must be done at once)
        for( size_t k=0; k< num_states(*_model); ++k )
        {
            if( _model->IsReservedState(k) ) continue;
            
            _model->SetEta( k,
                            _nextScale[k] );
        }
        return _est->Pseq();
    }
    
    // TODO - Find out why enable_if doesn't work in this case
    P reestimate_shape()
    {
        assert(Algo::model_type::is_explicit_duration_model::value);
        _reset(); // clear our buffers
        
        // Use the estimator to calculated improved values for all parameters
        for( size_t k=0; k< num_states(*_model); ++k )
        {
            if( _model->IsReservedState(k) ) continue;
            
            _nextShape[k] = _est->nextshape(k);
        }
        // Update the model with the improved params (updating the model must be done at once)
        for( size_t k=0; k< num_states(*_model); ++k )
        {
            if( _model->IsReservedState(k) ) continue;
            
            _model->SetMu( k,
                           _nextShape[k] );
        }
        return _est->Pseq();
    }
    
    
protected:
    void _reset()
    {
        _est->reset();
        {
            const size_t m_end = _nextA.shape()[0];
            const size_t n_end = _nextA.shape()[1];
            
            for( size_t m=0; m< m_end; ++m )
                for( size_t n=0; n< n_end; ++n )
                {	  _nextA[m][n] = _empty;	}
        }
        
        {
            const size_t m_end = _nextB.shape()[0];
            const size_t n_end = _nextB.shape()[1];
            
            for( size_t m=0; m< m_end; ++m )
                for( size_t n=0; n< n_end; ++n )
                {	  _nextB[m][n] = _empty;	}
        }
        
        {
            const size_t m_end = _nextScale.shape()[0];
            
            for( size_t m=0; m< m_end; ++m )
            {	  _nextScale[m] = _empty;	}
        }
        
        {
            const size_t m_end = _nextShape.shape()[0];
            
            for( size_t m=0; m< m_end; ++m )
            {	  _nextShape[m] = _empty;	}
        }
        
        _state_reset = true;
    }
    
public:
    void reset()
    {
        //std::cout<< "EM::reset()"<< std::endl;
        if( _state_reset ) return;
        //std::cout<< "EM::reset() *"<< std::endl;
        _reset();
    }
    
    void debug_print() const
    {
        std::cout<< "EM<Algo>"<< std::endl;
        std::cout<< " nextA: ";
        ::debug_print_exp(_nextA);
        std::cout<< " nextB: ";
        ::debug_print_exp(_nextB);
        std::cout<< " nextScale: ";
        ::debug_print(_nextScale);
        std::cout<< " nextShape: ";
        ::debug_print(_nextShape);
        std::cout<< " Estimator: ";
        _est->debug_print();
    }

public:
    P peek_forward(        size_t n, state_t i, time_t t  ) { return _est->peek_forward(        n, i, t); }
    P peek_forward_begin(  size_t n, state_t i, time_t t  ) const { return _est->peek_forward_begin(  n, i, t); }
    P peek_backward(       size_t n, state_t i, time_t t  ) const { return _est->peek_backward(       n, i, t); }
    P peek_backward_begin( size_t n, state_t i, time_t t  ) const { return _est->peek_backward_begin( n, i, t); }
    P peek_gamma(          size_t n, state_t i, time_t t  ) const { return _est->peek_gamma(          n, i, t); }
    P peek_sigma_t_xi(             size_t n, state_t i, state_t j ) const { return _est->peek_sigma_t_xi( n, i, j); }
    
protected:
    boost::shared_ptr<typename Algo::model_type> _model;
    
    boost::shared_ptr<Reestimator> _est;
    
    typedef boost::multi_array<typename Algo::probability_type,2> nextA_type;
    nextA_type _nextA;
    
    typedef boost::multi_array<typename Algo::probability_type,2> nextB_type;
    nextB_type _nextB;
    
    typedef boost::multi_array<typename Algo::probability_type,1> nextScale_type;
    nextScale_type _nextScale;
    
    typedef boost::multi_array<typename Algo::probability_type,1> nextShape_type;
    nextShape_type _nextShape;
    
    bool _state_reset;
    enum {_empty=-1}; // Remove this when all memoized funcs use MemoizedFunc
    
    fasta::seq_cont_t _seqs;
};


} // namespace baumwelch_algo

// template<class Algo>
// struct FuncTraits<baumwelch_algo::detail::GammaImpl<Algo> >
// {
//     typedef typename baumwelch_algo::detail::GammaImpl<Algo>::arg_type arg_type;
//     typedef typename baumwelch_algo::detail::GammaImpl<Algo>::val_type val_type;
    
//     struct empty_val
//     {
//         val_type operator()() const { return -1; /* gamma is an expected length (in internal time units), so must be non-negative */ }
//     };
// };

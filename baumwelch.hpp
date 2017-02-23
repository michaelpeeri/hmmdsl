#pragma once
#include <vector>
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
#include "asa121.hpp"
#include "common.hpp"
#include "read_fasta.hpp"
#include "tests.hpp"
#include "python_helpers.hpp"


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
			boost::shared_ptr<typename Algo::forward_algo_type> forward, 
			boost::shared_ptr<typename Algo::forward_begin_algo_type> forward_begin, 
			boost::shared_ptr<typename Algo::backward_algo_type> backward, 
			boost::shared_ptr<typename Algo::backward_begin_algo_type> backward_begin, 
			boost::shared_ptr<const typename Algo::sequence_type> seq )
    : _model(model)
	, _forward(forward)
	, _forward_begin(forward_begin)
	, _backward(backward)
	, _backward_begin(backward_begin)
	, _seq(seq)
	, _gamma(boost::extents[num_states(*model)][length(*seq)+2])
	, _xi(boost::extents[num_states(*model)][num_states(*model)][length(*seq)+2] )
	, _nextA(boost::extents[num_states(*model)][num_states(*model)] )
	, _nextB(boost::extents[num_states(*model)][num_symbols(*model)] )
  {
	  _reset();
  }

public:
	P gamma( size_t k, size_t i)
	{
		P val =  _gamma[k][i];
		if( val != _empty ) return val;
		
		_state_reset = false;
		_gamma[k][i] = val = _produce_gamma(k, i);
		return val;
	}

public:
	double Sigmai_gamma( size_t k )
	{
		P sigmai = 0.0;
		for( size_t i = 0; i< length(*_seq)+2; ++i )
		{
			const typename Algo::probability_type gammaki = gamma( k, i );
			sigmai += exp(gammaki);
		}
		return sigmai;
	}

public:
	double Sigmai_gamma_observe_s( size_t k, size_t s )
	{
		P sigmai = 0.0;
		const typename Algo::symbol_type sym = _model->get_symbol(s);
		for( size_t i = 1; i< length(*_seq)+1; ++i ) // Only iterate over positions where symbols are emitted
		{
			if( (*_seq)[i-1] == sym )
			{
				const typename Algo::probability_type gammaki = gamma( k, i );
				sigmai += exp(gammaki);
			}
		}
		return sigmai;
	}

public:
	void calc_all()
	{
		// TODO Impl. this
		assert(false);
		throw;
	}
	
protected:
	/*
	 * The probability of being in state k at position i (given the sequence and model)
	 */
	P _produce_gamma( size_t k, size_t i )
	{
		// gamma for initial/terminal states is 0 (they are never visited on positions 0..N-1)
		//if( _model->IsReservedState(k) )
		//	return -std::numeric_limits<typename Algo::probability_type>::max();

		const typename Algo::probability_type Pseq = _forward->calc();
		//const typename Algo::forward_algo_type::arg_type args = {k,i};

		const P gamma =
			(*_forward)( typename Algo::forward_algo_type::arg_type(k, i) )
			+ (*_backward)( typename Algo::backward_algo_type::arg_type(k, i) )
			- Pseq;
		assert( gamma < 1e-10 );
		if( i == 0 )                 assert( gamma == (k == _model->GetInitialState() ) ? 0.0 : -std::numeric_limits<typename Algo::probability_type>::max() );
		if( i == length(*_seq) + 1 ) assert( gamma == (k == _model->GetTerminalState()) ? 0.0 : -std::numeric_limits<typename Algo::probability_type>::max() );
		
		return gamma;
	}
	
			

public:
	P xi( size_t k, size_t l, size_t i )
	{
		P val = _xi[k][l][i];
		if( val != _empty ) return val;
		
		_state_reset = false;
		_xi[k][l][i] = val = _produce_xi(k, l, i);
		return val;
	}
public:
	double Sigmai_xi( size_t k, size_t l)
	{
		double sigmai = 0.0;
		for( size_t i = 0; i< length(*_seq) + 2; ++i )/* TODO - fix sigma ends! */
		{
			const P xikli = xi( k, l, i );
			sigmai += exp(xikli);
		}
		return sigmai;
	}
	
public:
	double Sigmal_xi( size_t k, size_t i)
	{
		double sigmal = 0.0;
		for( size_t l = 0; l< num_states(*_model); ++l )
		{
			const P xikli = xi( k, l, i );
			sigmal += exp( xikli );
		}
		return sigmal;
	}
	
			
			
protected:
	/* Calculate xi - The frequency of k->l transitions in position i */
	P _produce_xi( size_t k, size_t l, size_t i )
	{
		// I->l  or k->T transitions do not occuer in positions [0..len)
		//if( _model->IsReservedState(k) ||
		//	_model->IsReservedState(l)    )
		//	return -std::numeric_limits<typename Algo::probability_type>::max(); 
		
		P xi_kl =
			(*_forward)( typename Algo::forward_algo_type::arg_type(k, i) )
			+ _model->a(k, l)
			+ (*_backward)( typename Algo::backward_algo_type::arg_type(l, i+1) )
			- _forward->calc();

		if( i < length(*_seq) )
			xi_kl += _model->e(l, (*_seq)[i+1-1]);

//		if( k==12 && l == 13 )
//		{
//			std::cout<<	"xi: k=12 l=13 i="<< i
//					 << "  f="<< (*_forward)( k, i)
//					 << "  a="<< _model->a(k, l) 
//					 << "  e="<< _model->e(l, (*_seq)[i+1])
//					 << "  b="<<(*_backward)(l, i+1) 
//					 << "  p=="<< _forward->calc() 
//					 << "  xi="<< xi_kl<< std::endl;
//		}
		
		return xi_kl;
	}

public:
	P nexta( size_t k, size_t l)
	{
		P val =  _nextA[k][l];
		if( val != _empty ) return val;
		
		_state_reset = false;
		_nextA[k][l] = val = _produce_nexta(k, l);
		return val;		
	}

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

public:
	struct UnsupportedOperation {};
	
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
	

public:
	P Pseq()
	{
		return _forward->calc();
	}

protected:
	void _reset()
	{
		// initialize _gamma
		{
			const size_t m_end = _gamma.shape()[0];
			const size_t n_end = _gamma.shape()[1];
			
			for( size_t m=0; m< m_end; ++m )
				for( size_t n=0; n< n_end; ++n )
				{	  _gamma[m][n] = _empty;	}
		}
		
		// initialize _xi
		{
			const size_t m_end = _xi.shape()[0];
			const size_t n_end = _xi.shape()[1];
			const size_t o_end = _xi.shape()[2];
			
			for( size_t m=0; m< m_end; ++m )
				for( size_t n=0; n< n_end; ++n )
					for( size_t o=0; o< o_end; ++o )
					{	  _xi[m][n][o] = _empty;	}
		}
		
		
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

		_forward->reset();
		_backward->reset();
		_state_reset = true;
	}
	
public:
	void reset()
	{
	  //std::cout<< "BW::reset()"<< std::endl;
		if( _state_reset ) return;
		_reset();
	}
	
			
protected:
	P _produce_nexta( size_t k, size_t l)
	{
		// probabilities from/to reserved states are fixed
		//if( _model->IsReservedState(k) ||
		//	_model->IsReservedState(l)   )
		//	return _model->a(k, l);

		const double sigmai_xi = Sigmai_xi( k, l );
		const double sigmai_gamma = Sigmai_gamma( k );
 		assert( sigmai_gamma - sigmai_xi > -1e-10 );

		if( k==0 )
		{
			std::cout<< "Debug: nexta(0, "<< l<< "): sigmai_xi="<< sigmai_xi<< ", sigmai_gamma="<< sigmai_gamma<< std::endl;
		}
		
		
		if( sigmai_gamma < 1e-200 )
		{
			std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate parameters for state "<< k<< std::endl;
			return _model->a(k, l);
		}
		
		const P nexta = log(sigmai_xi / sigmai_gamma);
		assert(nexta < 0.001 );
		return (nexta < 0.0) ? nexta : 0.0;
	}

protected:
	P _produce_nextb( size_t k, size_t s )
	{
		const double sigmai_gamma_obs_s = Sigmai_gamma_observe_s( k, s );
		const double sigmai_gamma = Sigmai_gamma( k );
		assert( sigmai_gamma - sigmai_gamma_obs_s > -1e-10  );
		assert( !Algo::use_emissions_classes::value );

		if( sigmai_gamma < 1e-200 )
		{
			std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate parameters for state "<< k<< std::endl;
			return _model->e(k, s);
		}

		const P nextb = log(sigmai_gamma_obs_s / sigmai_gamma);
		assert(nextb < 0.001);
		return (nextb < 0.0) ? nextb : 0.0;
	}

protected:
	P _produce_nextb_class( size_t k, size_t s )
	{
		assert( Algo::use_emissions_classes::value );
		
		double sigman_gamma_obs_s = 0.0;
		double sigman_gamma = 0.0;

		typename Algo::model_type::emissions_class_members_iterator mem_it, mem_end;
		boost::tie( mem_it, mem_end ) = _model->GetEmissionsClassMembersRange(k);

		for( ; mem_it != mem_end; ++mem_it )
		{
			const typename Algo::model_type::StateId mem = *mem_it;
			const double sigmai_gamma_obs_s = Sigmai_gamma_observe_s( mem, s );
			const double sigmai_gamma = Sigmai_gamma( mem );

			assert( sigmai_gamma - sigmai_gamma_obs_s > -1e-10 );
			
			sigman_gamma_obs_s += sigmai_gamma_obs_s;
			sigman_gamma += sigmai_gamma;			
		}
		//sigman_gamma_obs_s /= double(sigman_count);
		//sigman_gamma /= double(sigman_count);

		const P nextb = log(sigman_gamma_obs_s / sigman_gamma);
		assert( nextb < 0.001 );
		return (nextb <= 0.0) ? nextb : 0.0;
	}
			
protected:
  const boost::shared_ptr<const typename Algo::model_type> _model;
  const boost::shared_ptr<typename Algo::forward_algo_type> _forward;
  const boost::shared_ptr<typename Algo::forward_begin_algo_type> _forward_begin;
  const boost::shared_ptr<typename Algo::backward_algo_type> _backward;
  const boost::shared_ptr<typename Algo::backward_begin_algo_type> _backward_begin;
  const boost::shared_ptr<const typename Algo::sequence_type> _seq;

  typedef boost::multi_array<P,2> gamma_type;
  gamma_type _gamma;

  typedef boost::multi_array<P,3> xi_type;
  xi_type _xi;

  typedef boost::multi_array<P,2> nextA_type;
  nextA_type _nextA;

  typedef boost::multi_array<P,2> nextB_type;
  nextB_type _nextB;

	bool _state_reset;

public:
  void debug_print() const
  {
    std::cout<< "BaumWelchAlgorithm<Algo>"<< std::endl;
    std::cout<< " gamma: ";
    ::debug_print(_gamma);
    std::cout<< " xi: ";
    ::debug_print_exp(_xi);
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
	
template<class Algo>
class GammaImpl : public UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, GammaImpl<Algo> >
{
protected:
	typedef GammaImpl<Algo> self_type;
	typedef typename Algo::probability_type P;
	typedef UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, GammaImpl<Algo> > base_type;
	typedef gamma_environment<Algo> env_;
	
public:
	typedef typename base_type::return_type return_type;
	typedef typename base_type::arg_type arg_type;	
	typedef typename gamma_environment<Algo>::type env_type;
	
protected:
	typedef MemoizedFunc<self_type, typename gamma_environment<Algo>::type, arg_type> host_type;
	host_type& _host;
	const boost::shared_ptr<const typename Algo::sequence_type> _seq;
	const boost::shared_ptr<const typename Algo::model_type> _model;
	const boost::shared_ptr<typename Algo::forward_algo_type> _forward;
	const boost::shared_ptr<typename Algo::forward_begin_algo_type> _forward_begin;
	const boost::shared_ptr<typename Algo::backward_algo_type> _backward;
	const boost::shared_ptr<typename Algo::backward_begin_algo_type> _backward_begin;
	
protected:
	GammaImpl( host_type& host, const env_type& env)
		: _host(host) 
		, _seq(boost::get<env_::seq>(env))
		, _model(boost::get<env_::model>(env))
		, _forward(boost::get<env_::forw_end>(env))
		, _forward_begin(boost::get<env_::forw_begin>(env))
		, _backward(boost::get<env_::back_end>(env))
		, _backward_begin(boost::get<env_::back_begin>(env))
		{}
	
protected:	
	P _produce_fli( const arg_type& args )
	{
		const size_t l = boost::fusion::at_c<0>(args);
		const size_t i = boost::fusion::at_c<1>(args);

		// Currently, gamma is not supported for initial and terminal states.
		if( _model->IsReservedState(l) ) return 0.0;
		//if( i == 0 )
		//	return (l==_model->GetInitialState() ) ? 1.0 : 0.0;
		
		//const size_t lastpos = length(*_seq)+1;
		//if( i == lastpos )
		//	return (l==_model->GetTerminalState()) ? 0 : -std::numeric_limits<P>::max();

		const arg_type prev(l, i-1);
		/*
		{
			P preva = 0.0;
			if( i > 0 ) preva = _host(prev);
			
			std::cout<< "l="<< l<< "\ti="<< i;
			P fb = -1;
			P bb = -1;
			if( i>0 ) {
			  P fb = (*_forward_begin )(prev); 		std::cout<< "\tfb="<< fb;
			  P bb = (*_backward_begin)(prev); 		std::cout<< "\tbb="<< bb;
			}
			const P fe = (*_forward       )(args); 		std::cout<< "\tfe="<< fe;
			const P be = (*_backward      )(args); 		std::cout<< "\tbe="<< be;
			const P begin = exp( fb + bb ); 		std::cout<< "\tbegin="<< begin;
			const P end   = exp( fe + be ); 		std::cout<< "\tend="<< end;
			std::cout<< "\tpreva="<< preva;
			const P curra   = begin - end; 		std::cout<< "\tcurr="<< curra;
			const P sum = preva + curra; 		std::cout<< "\tsum="<< sum;
			std::cout<< std::endl;
			}*/
		
		P curr;
		if( i > 0 )
		  {
		    curr =
		      exp( (*_forward_begin)(prev) + (*_backward_begin)(prev) ) -
		      exp( (*_forward      )(prev) + (*_backward      )(prev) );
		  }
		else
		  {
		    curr = 0;
		  }
		
		if( i == 0 )
			return curr;
		else
			return curr + _host(prev);
	}	
	
public:
	inline return_type operator()( const arg_type& args ) 		{	return _produce_fli(args);	}
	
private:
	GammaImpl(const GammaImpl&);
	const GammaImpl& operator=(const GammaImpl&);
	GammaImpl& operator=(GammaImpl&);
};
	
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
	typedef typename base_type::return_type P;
	
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
		for( typename Algo::model_type::StateId k = 0; k< num_states(*_model); ++k )
		{
			fli( arg_type(k, lastpos) );
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

public:

	HSMMBaumWelchAlgorithm(boost::shared_ptr<const typename Algo::model_type> model, 
						   boost::shared_ptr<typename Algo::forward_algo_type> forward, 
						   boost::shared_ptr<typename Algo::forward_begin_algo_type> forward_begin, 
						   boost::shared_ptr<typename Algo::backward_algo_type> backward, 
						   boost::shared_ptr<typename Algo::backward_begin_algo_type> backward_begin, 
						   boost::shared_ptr<const typename Algo::sequence_type> seq )
		: _model(model)
		, _forward(forward)
		, _forward_begin(forward_begin)
		, _backward(backward)
		, _backward_begin(backward_begin)
		, _seq(seq)
		, _gamma( typename gamma_type::env_type( model, seq, forward, forward_begin, backward, backward_begin ) )
		, _nextA(boost::extents[num_states(*model)][num_states(*model)] )
		, _nextB(boost::extents[num_states(*model)][num_symbols(*model)] )
		, _nextScale(boost::extents[num_states(*model)] )
		, _nextShape(boost::extents[num_states(*model)] )
	{
		_reset();
	}
	enum {_empty=1}; // Remove this when all functions use MemoizedFunc

public:
	inline P gamma( size_t k, size_t i)
	{
		return _gamma(typename gamma_type::arg_type(k,i));
	}

public:
	double Sigmai_gamma( size_t k )
	{
		P sigmai = 0.0;
		for( size_t i = 1; i< length(*_seq)+1; ++i ) // Only iterate over positions where symbols are emitted
		{
		  const typename Algo::probability_type gammaki = gamma( k, i);
			sigmai += gammaki;
		}
		return sigmai;
	}

public:
	double Sigmai_gamma_observe_s( size_t k, size_t s )
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
	P xi( size_t k, size_t l, size_t i )
	{
		return 0.0;
	}
public:
	double Sigmai_xi( size_t k, size_t l)
	{
		return 0.0;
	}
	
public:
	double Sigmal_xi( size_t k, size_t i)
	{
		return 0.0;
	}
	
public:
	P nexta( size_t k, size_t l)
	{
		P val =  _nextA[k][l];
		if( val != _empty ) return val;
		
		_state_reset = false;
		_nextA[k][l] = val = _produce_nexta(k, l);
		return val;		
	}

public:
	P nextb( size_t k, size_t s)
	{
		assert(false);
		
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

public:
	P nextscale( size_t k )
	{
		P val =  _nextScale[k];
		if( val != _empty ) return val;
		
		_state_reset = false;
		_nextScale[k] = val = _produce_nextscale(k);
		return val;		
	}

public:
	P nextshape( size_t k )
	{
		P val =  _nextShape[k];
		if( val != _empty ) return val;
		
		_state_reset = false;
		_nextShape[k] = val = _produce_nextshape(k);
		return val;		
	}

public:
	P Pseq()
	{
		return _forward->calc();
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

		_forward->reset();
		_forward_begin->reset();
		_backward->reset();
		_backward_begin->reset();
		_gamma.reset();
		_state_reset = true;
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
		_gamma.calc_all();
	}
	
			
protected:
	P _produce_nexta( size_t k, size_t l)
	{
		// probabilities from/to reserved states are fixed
		//if( _model->IsReservedState(k) ||
		//	_model->IsReservedState(l)   )
		//	return _model->a(k, l);
		const double sigmai_xi = Sigmai_xi( k, l );
		const double sigmai_gamma = Sigmai_gamma( k );
 		assert( sigmai_gamma - sigmai_xi > -1e-10 );

		/*
		if( k==0 )
		{
			std::cout<< "Debug: nexta(0, "<< l<< "): sigmai_xi="<< sigmai_xi<< ", sigmai_gamma="<< sigmai_gamma<< std::endl;
		}
		*/
		
		if( sigmai_gamma < 1e-200 )
		{
			std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate parameters for state "<< k<< std::endl;
			return _model->a(k, l);
		}
		
		const P nexta = log(sigmai_xi / sigmai_gamma);
		assert(nexta < 0.001 );
		return (nexta < 0.0) ? nexta : 0.0;
	}

protected:
	P _produce_nextb( size_t k, size_t s )
	{
		const double sigmai_gamma_obs_s = Sigmai_gamma_observe_s( k, s );
		const double sigmai_gamma = Sigmai_gamma( k );
		assert( sigmai_gamma - sigmai_gamma_obs_s > -1e-10  );
		assert( !Algo::use_emissions_classes::value );

		if( sigmai_gamma < 1e-200 )
		{
			std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate parameters for state "<< k<< std::endl;
			return _model->e(k, s);
		}

		const P nextb = log(sigmai_gamma_obs_s / sigmai_gamma);
		assert(nextb < 0.001);
		return (nextb < 0.0) ? nextb : 0.0;
	}

protected:
	P _produce_nextb_class( size_t k, size_t s )
	{
		assert( Algo::use_emissions_classes::value );
		
		double sigman_gamma_obs_s = 0.0;
		double sigman_gamma = 0.0;


		typename Algo::model_type::emissions_class_members_iterator mem_it, mem_end;
		boost::tie( mem_it, mem_end ) = _model->GetEmissionsClassMembersRange(k);
		

		for( ; mem_it != mem_end; ++mem_it )
		{
			const typename Algo::model_type::StateId mem = *mem_it;
			const double sigmai_gamma_obs_s = Sigmai_gamma_observe_s( mem, s );
			const double sigmai_gamma = Sigmai_gamma( mem );


			assert( sigmai_gamma - sigmai_gamma_obs_s > -1e-10 );
			
			sigman_gamma_obs_s += sigmai_gamma_obs_s;
			sigman_gamma += sigmai_gamma;			
		}

		//sigman_gamma_obs_s /= double(sigman_count);
		//sigman_gamma /= double(sigman_count);

		const P nextb = log(sigman_gamma_obs_s / sigman_gamma);
		assert( nextb < 0.001 );
		return (nextb <= 0.0) ? nextb : 0.0;
	}

protected:
	P _produce_nextscale( size_t j )
	{
		// Numerator
		typedef LogspaceDouble<> logspace_t;
		//std::cout<< "nextscale j="<< j << std::endl;
		
		double sigmat_alpha_beta = 0.0;
		for( size_t t = 1; t< length(*_seq)+1; ++t ) // Only iterate over positions where symbols are emitted
		{
			const typename Algo::forward_algo_type::arg_type j_t(j, t);
			const typename Algo::forward_algo_type::arg_type j_t_1(j, t-1);
			const logspace_t f((*_forward)(j_t), true);
			const logspace_t b((*_backward)(j_t), true);

			//std::cout<< boost::format("t=%d, f=%e, b=%e, f*b=%e") % t % f % b % (f*b)<< std::endl;
			
			sigmat_alpha_beta += f * b;
		}
		const double exp1 = sigmat_alpha_beta * _model->GetMu(j);
		//std::cout<< boost::format("sigmat_alpha_beta=%e, Mu=%e, numerator=%e") % sigmat_alpha_beta % _model->GetMu(j) % exp1 << std::endl;
		const P numerator = exp1;

		// Denominator
		P sigma_t = 0.0;
		for( size_t t = 1; t< length(*_seq)+1; ++t ) // Only iterate over positions where symbols are emitted
		{
			const typename Algo::forward_algo_type::arg_type j_t(j, t);
			const typename Algo::forward_algo_type::arg_type j_t_1(j, t-1);
			const logspace_t backward( (*_backward)(j_t), true);

			P sigma_tau = 0.0;
			for( size_t tau = 1; tau <= t && tau <= Algo::max_duration; ++tau )
			{
				const logspace_t d_j_tau( _model->p(j, tau), true);

				P sigma_i = 0.0;
				for( size_t i=0; i< num_states(*_model); ++i )
				{
					if( i == j ) continue;

					const logspace_t forward( (*_forward)(typename Algo::forward_algo_type::arg_type( i, t-tau ) ), true );
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

protected:
	P _produce_nextshape( size_t j )
	{
		typedef LogspaceDouble<> logspace_t;

		// Calculate the value of the b constant (Levinson eq. (21))
		// Numerator
		P sigma_t = 0.0;
		for( size_t t = 1; t < length(*_seq)+1; ++t )
		{
			const typename Algo::forward_algo_type::arg_type j_t(j, t);
			const logspace_t backward( (*_backward)(j_t), true);

			P sigma_tau = 0.0;
			for( size_t tau = 1; tau <= t && tau <= Algo::max_duration; ++tau )
			{
				const P log_eta_tau = log(_model->GetEta(j) * tau);
				const logspace_t d_j_tau( _model->p(j, tau), true );
				
				P sigma_i = 0.0;
				for( size_t i=0; i< num_states(*_model); ++i )
				{
					if( i == j ) continue;
					const logspace_t forward( (*_forward)(typename Algo::forward_algo_type::arg_type(i, t-tau) ), true );
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
		double sigmat_alpha_beta = 0.0;
		for( size_t t = 1; t< length(*_seq)+1; ++t ) // Only iterate over positions where symbols are emitted
		{
			const typename Algo::forward_algo_type::arg_type j_t(j, t);
			const logspace_t f((*_forward)(j_t), true);
			const logspace_t b((*_backward)(j_t), true);

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
  const boost::shared_ptr<typename Algo::forward_algo_type> _forward;
  const boost::shared_ptr<typename Algo::forward_begin_algo_type> _forward_begin;
  const boost::shared_ptr<typename Algo::backward_algo_type> _backward;
  const boost::shared_ptr<typename Algo::backward_begin_algo_type> _backward_begin;
  const boost::shared_ptr<const typename Algo::sequence_type> _seq;

  typedef detail::Gamma<Algo,detail::GammaImpl<Algo> > gamma_type;
  gamma_type _gamma;

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
    std::cout<< " gamma: ";
    detail::debug_print(_gamma);
    std::cout<< " nextA: ";
    ::debug_print(_nextA);
    std::cout<< " nextB: ";
    ::debug_print(_nextB);
    std::cout<< " nextScale: ";
    ::debug_print(_nextScale);
    std::cout<< " nextShape: ";
    ::debug_print(_nextShape);
  }


public:
	typedef typename gamma_type::result_in_logspace gamma_in_logspace;

};

template<class Algo>
void debug_print( const HSMMBaumWelchAlgorithm<Algo>& obj )
{
  obj.debug_print();
}

} // namespace baumwelch_algo

namespace tests {

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
	

template<typename Algo>
void test_38(typename Algo::baumwelch_algo_type& bw, const typename Algo::model_type& model, const typename Algo::sequence_type& seq)
{
  using namespace boost::accumulators;
  //std::cout<< "test (38)"<< std::endl;

  const typename Algo::probability_type epsilon = 1e-10;

  accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;

  for( size_t i=0; i< length(seq)+1; ++i )
  {
    for( size_t k=0; k< num_states(model); ++k )
	{
		const typename Algo::probability_type gammaki = bw.gamma( k, i );
		const double sigmalxi = bw.Sigmal_xi( k, i );

		const double delta = exp(gammaki) - sigmalxi;
		if( fabs(delta) > 0.001 )
		{			std::cout<< "test_38: delta>0.001: state="<< k<< "\ti="<< i<< ":\tgammaki="<< exp(gammaki)<< "\tsigmal_xi:"<< sigmalxi<< std::endl;		}
		
		// Update statistics
		acc( delta );
	}
	
  }
  const double min_delta = min(acc);
  const double max_delta = max(acc);
  // All results are in logspace;
  assert( fabs(min_delta) < epsilon );
  assert( fabs(max_delta) < epsilon );
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
	
	const double epsilon = 1e-10;
	
	accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;
	
	for( size_t k=0; k< num_states(model); ++k )
	{
		const double sigmai_gamma = bw.Sigmai_gamma(k);
		double sigmas_gamma = 0.0;
		
		for( size_t s=0; s< num_symbols(model); ++s )
		{
			//const typename Algo::symbol_type sym = model.get_symbol(s);
			const double sigmai_gamma_obs_s = bw.Sigmai_gamma_observe_s(k, s);
			assert( ( sigmai_gamma_obs_s == 0.0 ) || (!model.is_silent(k) ) );
			sigmas_gamma += sigmai_gamma_obs_s;
		}
		if( model.is_silent(k) ) continue;
		
		const double delta = sigmai_gamma - sigmas_gamma;
		if( fabs(delta) > epsilon )
		{
			std::cout<< "test_40c: delta for symbol "<< k<< " exceeds limit (sigmai_gamma="<< sigmai_gamma<< ", sigmas_gamma="<< sigmas_gamma<< ")"<< std::endl;
		}
		
		// Update statistics
		acc( delta );
	}
	
	// All results are in logspace;
	const double min_err = min(acc);
	const double max_err = max(acc);
	assert( fabs(min_err) < epsilon );
	assert( fabs(max_err) < epsilon );
	//std::cout<< min(acc)<< " <= Error <= "<< max(acc)<< "\t(average = "<< mean(acc)<< ")"<< std::endl;
}

} // namespace tests

namespace baumwelch_algo
{
	
template<typename Algo>
class MultipleSequenceBW
{
protected:
  typedef typename Algo::probability_type P;
  enum {_empty=-1}; // Remove this when all memoized funcs use MemoizedFunc
	
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
			typename Algo::forward_algo_type::env_type env = boost::make_tuple(seq, model);
			// Create a forward algo for it
			boost::shared_ptr<typename Algo::forward_algo_type> f(new typename Algo::forward_algo_type( env ) );
			_forward.push_back( f );
			// Create a forward algo for it
			boost::shared_ptr<typename Algo::forward_begin_algo_type> fb(new typename Algo::forward_begin_algo_type( env ) );
			_forward_begin.push_back( fb );
			bind_begin_end_algorithms( f, fb );
			// Create a backward algo for it
			boost::shared_ptr<typename Algo::backward_algo_type> b(new typename Algo::backward_algo_type( env ) );
			_backward.push_back( b );
			// Create a backward algo for it
			boost::shared_ptr<typename Algo::backward_begin_algo_type> bb(new typename Algo::backward_begin_algo_type( env ) );
			_backward_begin.push_back( bb );
			bind_begin_end_algorithms( b, bb );
			// Create a baum-welch algo for it
			boost::shared_ptr<typename Algo::baumwelch_algo_type> bw(new typename Algo::baumwelch_algo_type( model, f, fb, b, bb, seq ) );
			_baumwelch.push_back( bw );

		}
		_reset();
	}

public:
	P nexta(size_t k, size_t l)
	{
		// probabilities from/to reserved states are fixed
		//if( _model->IsReservedState(k) ||
		//	_model->IsReservedState(l)    )
		//	return _model->a(k, l);

		// Rabiner Eq. (109) (p. 273)
		double sigman_xi = 0.0;
		double sigman_gamma = 0.0;

		_state_reset = false;

		const size_t N = _seq.size();
		for( size_t n = 0; n < N; ++n )
		{
			//const double Pn = exp(_forward[n]->calc());
			//const double Pn = double(N);
			const double sigmai_xi = _baumwelch[n]->Sigmai_xi(k, l);
			const double sigmai_gamma = _baumwelch[n]->Sigmai_gamma(k);

			assert( sigmai_gamma - sigmai_xi > -1e-10 );

			sigman_xi += sigmai_xi /*/ Pn*/ ;
			sigman_gamma += sigmai_gamma /*/ Pn*/ ;
		}
		if( sigman_gamma < 1e-200 )
		{
			std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate parameters for state "<< k<< std::endl;
			return _model->a(k, l);
		}
		
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
			
			_nextB[k][s] = val = _produce_nextb(k,
												s );
		}
		else
		{
			const size_t emclass = _model->GetEmissionsClass(k);
			
			val = _produce_nextb_class(emclass,
									   s );

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

		double sigman_scale = 0.0;
		
		const size_t N = _seq.size();
		for( size_t n = 0; n < N; ++n )
		{
			sigman_scale += _baumwelch[n]->nextscale(k);
		}
		return sigman_scale / double(N);
	}
	
public:
	P nextshape( size_t k )
	{
		assert(!_model->IsReservedState(k));

		_state_reset = false;

		double sigman_shape = 0.0;
		
		const size_t N = _seq.size();
		for( size_t n = 0; n < N; ++n )
		{
			sigman_shape += _baumwelch[n]->nextshape(k);
		}
		return sigman_shape / double(N);
	}

protected:	
	P _produce_nextb(size_t k, size_t s)
	{
		assert( !Algo::use_emissions_classes::value );
		assert( !_model->is_silent(k) );
		// Rabiner Eq. (110) (p. 273)
		double sigman_gamma_obs_s = 0.0;
		double sigman_gamma = 0.0;

		_state_reset = false;

		//const typename Algo::symbol_type sym = _model->get_symbol(s);

		const size_t N = _seq.size();
		for( size_t n = 0; n < N; ++n )
		{
			const double sigmai_gamma_obs_s = _baumwelch[n]->Sigmai_gamma_observe_s(k, s );
			const double sigmai_gamma = _baumwelch[n]->Sigmai_gamma(k);

			assert( sigmai_gamma - sigmai_gamma_obs_s > -1e-10 );
			
			sigman_gamma_obs_s += sigmai_gamma_obs_s ;
			sigman_gamma += sigmai_gamma ;
		}

		if( sigman_gamma < 1e-200 )
		{
			std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate parameters for state "<< k<< std::endl;
			return _model->e(k, s );
		}
		
		if( sigman_gamma_obs_s < 0.0 ) sigman_gamma_obs_s = 0.0;
		
		const P nextb = log(sigman_gamma_obs_s / sigman_gamma);
		if( nextb >= 1e-12 )
		  {
		    std::cout<< "Error: Sigma n gamma obs s > Sigma n gamma (Sigma n gamma obs s="<< sigman_gamma_obs_s<< ", Sigma n gamma="<< sigman_gamma<< std::endl;
		    assert( nextb < 1e-12 );
		  }
		return (nextb <= 0.0) ? nextb : 0.0;
	}

protected:
	P _produce_nextb_class(size_t k, size_t s)
	{
		assert( Algo::use_emissions_classes::value );
		assert( !_model->is_silent(k) );
		// Rabiner Eq. (110) (p. 273)
		double sigman_gamma_obs_s = 0.0;
		double sigman_gamma = 0.0;

		_state_reset = false;
		
		bool d = false;
		//if( s==1 ) d = true;

		//const typename Algo::symbol_type sym = _model->get_symbol(s);

		typename Algo::model_type::emissions_class_members_iterator mem_it, mem_end;
		boost::tie( mem_it, mem_end ) = _model->GetEmissionsClassMembersRange(k);

		const size_t N = _seq.size();

		size_t members = 0;
		for( ; mem_it != mem_end; ++mem_it )
		{
			const typename Algo::model_type::StateId mem = *mem_it;
			
			for( size_t n = 0; n < N; ++n )
			{
				const double sigmai_gamma_obs_s = _baumwelch[n]->Sigmai_gamma_observe_s( mem, s );
				const double sigmai_gamma = _baumwelch[n]->Sigmai_gamma( mem );
				
				assert( sigmai_gamma - sigmai_gamma_obs_s > -1e-10 ); // sigmai_gamma >= sigmai_gamma_obs_n

				//std::cout<< std::endl<< "(   k="<< k<< ", sym="<< s<<", n="<< n<< ", sigmai_gamma_obs_s="<< sigmai_gamma_obs_s<< ", sigmai_gamma="<< sigmai_gamma<< std::endl;
				
				//sigman_gamma_obs_s += sigmai_gamma_obs_s;
				//sigman_gamma += sigmai_gamma;
				sigman_gamma_obs_s += sigmai_gamma_obs_s/sigmai_gamma; //todo fix div by 0
				sigman_gamma += 1.0;
			}
			
			++members;
		}

		//sigman_gamma_obs_s /= double(sigman_count);
		//sigman_gamma /= double(sigman_count);

		if( sigman_gamma < 1e-200 )
		{
			std::cout<< "Warning: sigman_gamma==0: insufficient data to estimate parameters for state "<< k<< std::endl;
			return _model->e(k, s );
		}
		if(d) 
		{
			std::cout<< "k="<< k<< ", sym="<< s<< ", members="<< members<< ", sigman_gamma_obs_s="<< sigman_gamma_obs_s<< ", sigman_gamma="<< sigman_gamma<< ", N="<< N<< std::endl;
		}

		if( sigman_gamma_obs_s < 0.0 ) sigman_gamma_obs_s = 0.0;
			
		const P nextb = log(sigman_gamma_obs_s / sigman_gamma);
		if( nextb >= 1e-12 )
		  {
		    std::cout<< "Error: Sigma n gamma obs s > Sigma n gamma (k="<< k<< ", s="<< s<< ", Sigma n gamma obs s="<< sigman_gamma_obs_s<< ", Sigma n gamma="<< sigman_gamma<< std::endl;
		    assert( nextb < 1e-12 );
		  }
		return (nextb <= 0.0) ? nextb : 0.0;
	}

public:
	P Pseq(size_t n)
	{
		return _forward[n]->calc();
	}

public:
	P Pseq()
	{
		// Eq. (107)
		P Pi = 0.0;

		typename _forward_t::iterator it, it_end;
		it = _forward.begin();
		it_end = _forward.end();
		for( ; it != it_end; ++it ) Pi += (*it)->calc();
		return Pi;
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
				{	  _nextA[m][n] = _empty;	}
		}

		{
			const size_t m_end = _nextB.shape()[0];
			const size_t n_end = _nextB.shape()[1];
			
			for( size_t m=0; m< m_end; ++m )
				for( size_t n=0; n< n_end; ++n )
				{	  _nextB[m][n] = _empty;	}
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
		const size_t N = _seq.size();
		for( size_t n = 0; n < N; ++n )
		{
			tests::test_38<Algo>( *(_baumwelch[n]), *_model, *(_seq[n]) );
		}
	}
	void test_40c()
	{
		const size_t N = _seq.size();
		for( size_t n = 0; n < N; ++n )
		{
			tests::test_40c<Algo>( *(_baumwelch[n]), *_model );
		}
	}
	
	void test_D59_forw_back()
	{
		const size_t N = _seq.size();
		for( size_t n = 0; n < N; ++n )
		{
			::tests::test_D59_forw_back<Algo>( *(_forward[n]), *(_backward_begin[n]) );
		}
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

	
protected:
	const boost::shared_ptr<const typename Algo::model_type> _model;
	typedef std::vector< boost::shared_ptr<typename Algo::forward_algo_type> > _forward_t;
	_forward_t _forward;

	typedef std::vector< boost::shared_ptr<typename Algo::forward_begin_algo_type> > _forward_begin_t;
	_forward_begin_t _forward_begin;

	typedef std::vector< boost::shared_ptr<typename Algo::backward_algo_type> > _backward_t;
	_backward_t _backward;
	
	typedef std::vector< boost::shared_ptr<typename Algo::backward_begin_algo_type> > _backward_begin_t;
	_backward_begin_t _backward_begin;

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
	typename Algo::probability_type reestimate_a()
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

	typename Algo::probability_type reestimate_b()
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
				_model->SetEmissionProbabilityLogspace( k,
													   _model->get_symbol(sidx),
													   _nextB[k][sidx] );
			}
		}
		return _est->Pseq();
	}

	// TODO - Find out why enable_if doesn't work in this case
	typename Algo::probability_type reestimate_scale()
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
	typename Algo::probability_type reestimate_shape()
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

template<class Algo>
struct FuncTraits<baumwelch_algo::detail::GammaImpl<Algo> >
{
	enum {out_of_range_value = -1};
};

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
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/at_c.hpp> 
#include <boost/mpl/at.hpp>
#include "memoized.hpp"





namespace backward_algo
{

namespace detail
{
	
// Metafunction to return environment for BackwardAlgorithm<Algo>
template<typename Algo>
struct get_backward_env
{
	typedef boost::tuple<boost::shared_ptr<const typename Algo::sequence_type>, boost::shared_ptr<const typename Algo::model_type> > type;
	enum {seq, model};
};

template<class Impl>
struct can_calculate_sequence_probability
{
	typedef boost::mpl::bool_<false> value;
};
	

}

typedef boost::fusion::vector<size_t,size_t> backward_algo_arg_t;


template<typename Algo, typename BackwardImpl>
struct BackwardAlgorithm : MemoizedFunc<BackwardImpl, typename detail::get_backward_env<Algo>::type, backward_algo_arg_t >
{
public:
	typedef typename detail::get_backward_env<Algo>::type env_type;
	typedef boost::fusion::vector<size_t,size_t> arg_type;

protected:
	typedef MemoizedFunc<BackwardImpl,typename detail::get_backward_env<Algo>::type, arg_type> base_type;
	boost::shared_ptr<const typename Algo::sequence_type> _seq;
	boost::shared_ptr<const typename Algo::model_type> _model;

public:
	typedef typename base_type::return_type P;


public:
	using base_type::debug_print;
	using base_type::operator();
	using base_type::reset;
	using base_type::_reset;
	using base_type::fli;

public:
	BackwardAlgorithm(env_type& env)
		: base_type( env, get_standard_container(boost::fusion::vector<size_t,size_t>(
													 num_states(*(boost::get<1>(env)))
													 ,length(*(boost::get<0>(env)))+2 ) ))
		, _seq(boost::get<0>(env))
		, _model(boost::get<1>(env))
	{
		_reset();
	}

public:
	void calc_all()
	{
		const size_t firstpos =  0;
		for( typename Algo::model_type::StateId k = 0; k< num_states(*_model); ++k )
		{
			fli( arg_type(k, firstpos) );
		}
	}
	
public:
	P calc()
	{
		assert( detail::can_calculate_sequence_probability<BackwardImpl>::value::value );
		
		const size_t firstpos =  0;
		const typename Algo::model_type::StateId initstate = _model->GetInitialState();
		P p = 0.0;
		
		// TODO - Fix this (in forward algo, the same code appears to work for both types of models)
		if( Algo::model_type::is_explicit_duration_model::value )
		{
			// for HSMMs:
			P sigma_k = -std::numeric_limits<P>::max();
			for( typename Algo::model_type::StateId k = 0; k< num_states(*_model); ++k )
			{
				P next = _model->a(initstate, k) + fli( arg_type(k, firstpos) );
				
				sigma_k = logspace_sum( sigma_k,
										next );
			}
			p = sigma_k;
		}
		else
		{
			// for HMMs:
			p = fli( arg_type(initstate, firstpos) );
		}
	  
#ifndef NO_TESTS
		if ( p <= -std::numeric_limits<P>::max() )
		{
			std::cout<< "BackwardAlgorithm: P(seq|model) == 0"<< std::endl;
			debug_print();
			assert( p > -std::numeric_limits<P>::max() );
		}
#endif
		return p;
  }

private:
	BackwardAlgorithm(const BackwardAlgorithm&);
	const BackwardAlgorithm& operator=(const BackwardAlgorithm&);
	BackwardAlgorithm& operator=(BackwardAlgorithm&);
};


template<class Algo, class Impl>
void debug_print( const BackwardAlgorithm<Algo,Impl>& obj )
{
	obj.debug_print();
}


template<class Algo>
class BackwardAlgorithmHMM : public UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, BackwardAlgorithmHMM<Algo> >
{
	
protected:
	typedef BackwardAlgorithmHMM<Algo> self_type;
	typedef typename Algo::probability_type P;
	typedef UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, BackwardAlgorithmHMM<Algo> > base_type;
	typedef typename detail::get_backward_env<Algo>::type env_type;
	typedef detail::get_backward_env<Algo> env_;

public:
	typedef typename base_type::return_type return_type;
	typedef typename base_type::arg_type arg_type;	
	
protected:
	typedef MemoizedFunc<self_type, typename detail::get_backward_env<Algo>::type, arg_type> host_type;
	host_type& _host;
	const boost::shared_ptr<const typename Algo::sequence_type> _seq;
	const boost::shared_ptr<const typename Algo::model_type> _model;

public:
	void set_end_func( typename Algo::backward_algo_type& func ) {}
	void set_end_func( boost::shared_ptr<typename Algo::backward_algo_type> func ) {}
	void set_begin_func( typename Algo::backward_begin_algo_type& func ) {}
	void set_begin_func( boost::shared_ptr<typename Algo::backward_begin_algo_type> func ) {}
	
protected:
	BackwardAlgorithmHMM( host_type& host, const env_type& env)
		: _host(host) 
		, _seq(boost::get<env_::seq>(env))
		, _model(boost::get<env_::model>(env))
	{}

protected:
	P _produce_fli( const arg_type& args )
	{
		const size_t l = boost::fusion::at_c<0>(args);
		const size_t i = boost::fusion::at_c<1>(args);
		const size_t lastpos =  length(*_seq)+1;
		const typename Algo::model_type::StateId termstate = _model->GetTerminalState();
		
		if( i == lastpos )
			return (l==termstate) ? 0.0 : -std::numeric_limits<P>::max();
		
		P sum = -std::numeric_limits<P>::max();
		const typename Algo::symbol_type nextsym = (*_seq)[i+1-1];
		
		for( typename Algo::model_type::StateId k = 0; k< num_states(*_model); ++k )
		{
			P next = _model->a(l, k) + _host.fli( arg_type(k, i+1) );
			
			if( i < lastpos - 1 )
				next += _model->e(k, nextsym );
			
			sum = logspace_sum( sum,
								next );
		}
		
		return sum;
	}
	

public:
	return_type operator()( const arg_type& args ) 		{	return _produce_fli(args);	}

private:
	BackwardAlgorithmHMM(const BackwardAlgorithmHMM&);
	const BackwardAlgorithmHMM& operator=(const BackwardAlgorithmHMM&);
	BackwardAlgorithmHMM& operator=(BackwardAlgorithmHMM&);
};

namespace detail
{
template<typename Algo>
struct can_calculate_sequence_probability<BackwardAlgorithmHMM<Algo> >
{
	typedef boost::mpl::bool_<true> value;
};
} // namespace detail


template<class Algo>
class BackwardAlgorithmHSMMEnd : public UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, BackwardAlgorithmHSMMEnd<Algo> >
{
	
protected:
	typedef BackwardAlgorithmHSMMEnd<Algo> self_type;
	typedef typename Algo::probability_type P;
	typedef UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, BackwardAlgorithmHSMMEnd<Algo> > base_type;
	typedef typename detail::get_backward_env<Algo>::type env_type;
	typedef detail::get_backward_env<Algo> env_;

public:
	typedef typename base_type::return_type return_type;
	typedef typename base_type::arg_type arg_type;	
	
protected:
	typedef MemoizedFunc<self_type, typename detail::get_backward_env<Algo>::type, arg_type> host_type;
	host_type& _host;
	const boost::shared_ptr<const typename Algo::sequence_type> _seq;
	const boost::shared_ptr<const typename Algo::model_type> _model;

protected:	
	boost::shared_ptr<typename Algo::backward_begin_algo_type> _back_begin;
public:
	void set_begin_func( typename Algo::backward_begin_algo_type& func ) {	_back_begin = boost::shared_ptr<typename Algo::backward_begin_algo_type>(&func, null_deleter()); }
	void set_begin_func( boost::shared_ptr<typename Algo::backward_begin_algo_type> func ) {	_back_begin = func; }
	

protected:
	BackwardAlgorithmHSMMEnd( host_type& host, const env_type& env)
		: _host(host) 
		, _seq(boost::get<env_::seq>(env))
		, _model(boost::get<env_::model>(env))
	{}

protected:
	P _produce_fli( const arg_type& args )
	{
		const size_t l = boost::fusion::at_c<0>(args);
		const size_t i = boost::fusion::at_c<1>(args);
		const size_t lastpos =  length(*_seq)+1;
		const typename Algo::model_type::StateId termstate = _model->GetTerminalState();

		// At the terminal position, only the terminal state may end
		if( i == lastpos )
			return (l==termstate) ? 0.0 : -std::numeric_limits<P>::max();

		if( (i == 0) && (l != _model->GetInitialState()) ) return -std::numeric_limits<P>::max();
		if( l == _model->GetTerminalState() ) return -std::numeric_limits<P>::max();
		
		// Rabiner eq. (76)
		P sum = -std::numeric_limits<P>::max();
		for( typename Algo::model_type::StateId k = 0; k< num_states(*_model); ++k )
		{
			P next = _model->a(l, k) + (*_back_begin)( arg_type(k, i) );

			//std::cout<< "Back_End   (l="<< l<< ", i="<< i<< ", k="<< k<< "): a(l,k)="<< _model->a(l, k)<< ", back*(k,i)="<< (*_back_begin)( arg_type(k, i) )<< std::endl;
			
			sum = logspace_sum( sum,
								next );
		}
		
		return sum;
	}
	

public:
	inline return_type operator()( const arg_type& args ) 		{	return _produce_fli(args);	}

private:
	BackwardAlgorithmHSMMEnd(const BackwardAlgorithmHSMMEnd&);
	const BackwardAlgorithmHSMMEnd& operator=(const BackwardAlgorithmHSMMEnd&);
	BackwardAlgorithmHSMMEnd& operator=(BackwardAlgorithmHSMMEnd&);
};

template<class Algo>
class BackwardAlgorithmHSMMBegin : public UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, BackwardAlgorithmHSMMBegin<Algo> >
{
	
protected:
	typedef BackwardAlgorithmHSMMBegin<Algo> self_type;
	typedef typename Algo::probability_type P;
	typedef UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, BackwardAlgorithmHSMMBegin<Algo> > base_type;
	typedef typename detail::get_backward_env<Algo>::type env_type;
	typedef detail::get_backward_env<Algo> env_;

public:
	typedef typename base_type::return_type return_type;
	typedef typename base_type::arg_type arg_type;	
	
protected:
	typedef MemoizedFunc<self_type, typename detail::get_backward_env<Algo>::type, arg_type> host_type;
	host_type& _host;
	const boost::shared_ptr<const typename Algo::sequence_type> _seq;
	const boost::shared_ptr<const typename Algo::model_type> _model;

protected:	
	boost::shared_ptr<typename Algo::backward_algo_type> _back_end;
public:
	void set_end_func( typename Algo::backward_algo_type& func ) {	_back_end = boost::shared_ptr<typename Algo::backward_algo_type>(&func, null_deleter()); }
	void set_end_func( boost::shared_ptr<typename Algo::backward_algo_type> func ) {	_back_end = func; }

protected:
	BackwardAlgorithmHSMMBegin( host_type& host, const env_type& env)
		: _host(host) 
		, _seq(boost::get<env_::seq>(env))
		, _model(boost::get<env_::model>(env))
	{}

protected:
	P _produce_fli( const arg_type& args )
	{
		const size_t l = boost::fusion::at_c<0>(args);
		const size_t i = boost::fusion::at_c<1>(args);
		const size_t lastpos =  length(*_seq)+1;
		const typename Algo::model_type::StateId termstate = _model->GetTerminalState();

		// TODO - CHECK THIS
		if( i+1 == lastpos )
			return (l==termstate) ? 0.0 : -std::numeric_limits<P>::max();
		if( l == termstate ) return -std::numeric_limits<P>::max();
		if( l == _model->GetInitialState() ) return -std::numeric_limits<P>::max();

		// Rabiner eq. (77)
		P sum = -std::numeric_limits<P>::max();
		for( size_t d=1; d<=Algo::max_duration; ++d)
		{
			// The next state begins at position i+d+1, but this may not occur after lastpos (such paths do not contribute to the sum of probabilities)
			if( i+d+1 > lastpos ) continue; // TODO - CHECK THIS

			P next = (*_back_end)(arg_type(l, i+d))
				+ _model->p(l, d);

			//std::cout<< "Back_Begin (l="<< l<< ", i="<< i<< ", d="<< d<< ") : p(l,d)="<< _model->p(l, d)<< ", back(l,i+d)="<< (*_back_end)(arg_type(l, i+d));

			if( !_model->is_silent(l) )
			{
				 P sigmab = 0.0; //-std::numeric_limits<P>::max();
			
				for( size_t s=i+1; s<= i+d; ++s )
				{
					const P esym = _model->e( l, (*_seq)[s-1] );
					//std::cout<< " e('"<< (*_seq)[s-1]<< "')="<< esym;
					sigmab += esym;
				}
				//std::cout<< ", sigmab="<< sigmab;
				next += sigmab;
				
			}

			sum = logspace_sum( sum,
								next );
		}
		//std::cout<< std::endl;
		
		return sum;
	}
	

public:
	inline return_type operator()( const arg_type& args ) 		{	return _produce_fli(args);	}

private:
	BackwardAlgorithmHSMMBegin(const BackwardAlgorithmHSMMBegin&);
	const BackwardAlgorithmHSMMBegin& operator=(const BackwardAlgorithmHSMMBegin&);
	BackwardAlgorithmHSMMBegin& operator=(BackwardAlgorithmHSMMBegin&);
};

namespace detail
{
template<typename Algo>
struct can_calculate_sequence_probability<BackwardAlgorithmHSMMBegin<Algo> >
{
	typedef boost::mpl::bool_<true> value;
};
}

} // namespace backward_algo

template<class Algo>
struct FuncTraits<backward_algo::BackwardAlgorithmHMM<Algo> >
{
	typedef typename backward_algo::BackwardAlgorithmHMM<Algo>::arg_type arg_type;
	typedef typename backward_algo::BackwardAlgorithmHMM<Algo>::return_type val_type;

	struct empty_val
	{
		val_type operator()() const { return 1; }
	};
};
template<class Algo>
struct FuncTraits<backward_algo::BackwardAlgorithmHSMMBegin<Algo> >
{
	typedef typename backward_algo::BackwardAlgorithmHSMMBegin<Algo>::arg_type arg_type;
	typedef typename backward_algo::BackwardAlgorithmHSMMBegin<Algo>::return_type val_type;

	struct empty_val
	{
		val_type operator()() const { return 1; }
	};
};
template<class Algo>
struct FuncTraits<backward_algo::BackwardAlgorithmHSMMEnd<Algo> >
{
	typedef typename backward_algo::BackwardAlgorithmHSMMEnd<Algo>::arg_type arg_type;
	typedef typename backward_algo::BackwardAlgorithmHSMMEnd<Algo>::return_type val_type;

	struct empty_val
	{
		val_type operator()() const { return 1; }
	};
};

/*
template<class Algo>
struct FuncTraits<backward_algo::BackwardAlgorithmHMM<Algo> >
{
	enum {out_of_range_value = 1};
};
template<class Algo>
struct FuncTraits<backward_algo::BackwardAlgorithmHSMMBegin<Algo> >
{
	enum {out_of_range_value = 1};
};
template<class Algo>
struct FuncTraits<backward_algo::BackwardAlgorithmHSMMEnd<Algo> >
{
	enum {out_of_range_value = 1};
};
*/

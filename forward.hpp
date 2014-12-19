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
#include "common.hpp"
#include "memoized.hpp"
#include "v2.hpp"
#include "algo.hpp"
#include "sequential_processing.hpp"
#include "computation.hpp"


namespace forward_algo
{

namespace detail
{
        
// Metafunction to return environment for ForwardAlgorithm<Algo>
template<typename Algo>
struct get_forward_env
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

typedef boost::fusion::vector<size_t,size_t> forward_algo_arg_t;


template<typename Algo, typename ForwardImpl>
class ForwardAlgorithm : public MemoizedFunc<ForwardImpl,typename detail::get_forward_env<Algo>::type, forward_algo_arg_t >
{
public:
    typedef boost::fusion::vector<size_t,size_t> arg_type;
    typedef typename detail::get_forward_env<Algo>::type env_type;
    
protected:
    typedef MemoizedFunc<ForwardImpl,typename detail::get_forward_env<Algo>::type, arg_type> base_type;
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
    ForwardAlgorithm(env_type& env)
    : base_type( env, get_standard_container(boost::fusion::vector<size_t,size_t>(
                                                 num_states(*(boost::get<1>(env)))
                                                 , length(*(boost::get<0>(env)))+2     ) ))
    , _seq(boost::get<0>(env))
    , _model(boost::get<1>(env))
    {
        _reset();
    }
    
public:
    void calc_all()
    {
        const size_t lastpos =  length(*_seq)+1;
        fli( arg_type(_model->GetTerminalState(), lastpos) );
    }
    
public:
    P calc()
    {
        assert( detail::can_calculate_sequence_probability<ForwardImpl>::value::value );
        const size_t lastpos =  length(*_seq)+1;
        const typename Algo::model_type::StateId termstate = _model->GetTerminalState();
        
        const P sigma_k = fli( arg_type(termstate, lastpos) );
        
#ifndef NO_TESTS
        if ( sigma_k <= -std::numeric_limits<P>::max() )
        {
            //std::cout<< "ForwardAlgorithm: P(seq|model) == 0"<< std::endl;
            //debug_print();
            //::debug_print(*_model);
            //assert( sigma_k > -std::numeric_limits<P>::max() );
        }
#endif
        return sigma_k;
    }
    
private:
    ForwardAlgorithm(const ForwardAlgorithm&);
    const ForwardAlgorithm& operator=(const ForwardAlgorithm&);
    ForwardAlgorithm& operator=(ForwardAlgorithm&);
};


template<class Algo, class Impl>
void debug_print( const ForwardAlgorithm<Algo, Impl>& obj )
{
    obj.debug_print();
}       

template<class Algo>
class ForwardAlgorithmHMM : public UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, ForwardAlgorithmHMM<Algo> >
{
protected:
    typedef ForwardAlgorithmHMM<Algo> self_type;
    typedef typename Algo::probability_type P;
    typedef UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, ForwardAlgorithmHMM<Algo> > base_type;
    typedef typename detail::get_forward_env<Algo>::type env_type;
    typedef detail::get_forward_env<Algo> env_;
    
public:
    typedef typename base_type::return_type return_type;
    typedef typename base_type::arg_type arg_type;  
    
protected:
    typedef MemoizedFunc<self_type, typename detail::get_forward_env<Algo>::type, arg_type> host_type;
    host_type& _host;
    const boost::shared_ptr<const typename Algo::sequence_type> _seq;
    const boost::shared_ptr<const typename Algo::model_type> _model;
    
public:
    void set_end_func( typename Algo::forward_algo_type& func ) {}
    void set_end_func( boost::shared_ptr<typename Algo::forward_algo_type> func ) {}
    void set_begin_func( typename Algo::forward_begin_algo_type& func ) {}
    void set_begin_func( boost::shared_ptr<typename Algo::forward_begin_algo_type> func ) {}
    
protected:
    ForwardAlgorithmHMM( host_type& host, const env_type env)
        : _host(host) 
        , _seq(boost::get<env_::seq>(env))
        , _model(boost::get<env_::model>(env))
    {}

protected:      
    P _produce_fli( const arg_type& args )
    {
        const size_t l = boost::fusion::at_c<0>(args);
        const size_t i = boost::fusion::at_c<1>(args);
        
        if( i == 0 )
            return (l == _model->GetInitialState()) ? 0 : -std::numeric_limits<P>::max();
        
        const size_t lastpos =  length(*_seq)+1;
        if( (i == lastpos) && (l != _model->GetTerminalState()) ) return -std::numeric_limits<P>::max();
        
        P sum = -std::numeric_limits<P>::max();
        for( typename Algo::model_type::StateId k = 0; k< num_states(*_model); ++k )
        {
            const P prev = _host.fli( typename Algo::forward_algo_type::arg_type(k, i-1) ) + _model->a(k, l);
            
            sum = logspace_sum( sum,
                                prev );
        }
        
        if( !_model->is_silent(l) )
            sum += _model->e( l, (*_seq)[i-1] );
        
        return sum;
    }       
    
public:
    return_type operator()( const arg_type& args )          {       return _produce_fli(args);      }
    
private:
    ForwardAlgorithmHMM(const ForwardAlgorithmHMM&);
    const ForwardAlgorithmHMM& operator=(const ForwardAlgorithmHMM&);
    ForwardAlgorithmHMM& operator=(ForwardAlgorithmHMM&);
};

namespace detail
{
template<typename Algo>
struct can_calculate_sequence_probability<ForwardAlgorithmHMM<Algo> >
{
    typedef boost::mpl::bool_<true> value;
};
}


template<class Algo>
class ForwardAlgorithmHSMMEnd : public UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, ForwardAlgorithmHSMMEnd<Algo> >
{
protected:
    typedef ForwardAlgorithmHSMMEnd<Algo> self_type;
    typedef typename Algo::probability_type P;
    typedef UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, ForwardAlgorithmHMM<Algo> > base_type;
    typedef typename detail::get_forward_env<Algo>::type env_type;
    typedef detail::get_forward_env<Algo> env_;
    
public:
    typedef typename base_type::return_type return_type;
    typedef typename base_type::arg_type arg_type;  
    
protected:
    typedef MemoizedFunc<self_type, typename detail::get_forward_env<Algo>::type, arg_type> host_type;
    host_type& _host;
    const boost::shared_ptr<const typename Algo::sequence_type> _seq;
    const boost::shared_ptr<const typename Algo::model_type> _model;
    
protected:
    boost::shared_ptr<typename Algo::forward_begin_algo_type> _forw_begin;
    
public:
    void set_begin_func( typename Algo::forward_begin_algo_type& func ) {   _forw_begin = boost::shared_ptr<typename Algo::forward_begin_algo_type>(&func, null_deleter()); }
    void set_begin_func( boost::shared_ptr<typename Algo::forward_begin_algo_type> func ) { _forw_begin = func; }
    

protected:
    ForwardAlgorithmHSMMEnd( host_type& host, const env_type& env)
        : _host(host) 
        , _seq(boost::get<env_::seq>(env))
        , _model(boost::get<env_::model>(env))
    {}

protected:      
    P _produce_fli( const arg_type& args )
    {
        const size_t l = boost::fusion::at_c<0>(args);
        const size_t i = boost::fusion::at_c<1>(args);
        
        //std::cout<< "Forw_End   (l="<< l<< ", i="<< i<< ")"<< std::endl;
        
        // At i=0, ai=1.0 if i is the initial state, and 0.0 otherwise
        if( i == 0 )
        {
            return (l == _model->GetInitialState()) ? 0 : -std::numeric_limits<P>::max();
        }
        else
        {
            if( l == _model->GetInitialState() ) 
                return -std::numeric_limits<P>::max();
        }
        
        const size_t lastpos =  length(*_seq)+1;
        
        // At i=lastpos, ai=0.0 if i is not the terminal state (and calculated as usual otherwise)
        if( (i == lastpos) && (l != _model->GetTerminalState()) ) return -std::numeric_limits<P>::max();
        if( (i <  lastpos) && (l == _model->GetTerminalState()) ) return -std::numeric_limits<P>::max();
        
        // Sum over all states k
        P sum = -std::numeric_limits<P>::max();
        for( size_t d=1; d<=Algo::max_duration; ++d )
        {
            if( d>=i+1 ) continue; // TODO - CHECK THIS
            
            // Rabiner eq. (75)
            P prev = (*_forw_begin)( arg_type(l, i-d) )
                + _model->p(l, d);
            
            // Calculate sigmab, the contribution of the emissions over all d symbols
            if(!_model->is_silent(l))  // silent states have no emissions probs
            {
                P sigmab = 0.0; //-std::numeric_limits<P>::max();
                for( size_t s = i-d+1; s <= i; ++s )
                {                                       
                    const P esym = _model->e( l, (*_seq)[s-1] );
                    //std::cout<< " e('"<< (*_seq)[s-1]<< "')="<< esym;
                    sigmab += esym;
                }
                
                
                prev += sigmab;
                //std::cout<< ", sigmab="<< sigmab;
            }
            //std::cout<< std::endl;
            
            sum = logspace_sum( sum,
                                prev );
        }
        
        if( l==0 && i==1)
        {
            //std::cout<< "Forw_End   (l="<< l<< ", i="<< i<< ", sum="<< sum<< std::endl;
        }
        return sum;
    }       
    
public:
    inline return_type operator()( const arg_type& args )           {       return _produce_fli(args);      }
    
private:
    ForwardAlgorithmHSMMEnd(const ForwardAlgorithmHSMMEnd&);
    const ForwardAlgorithmHSMMEnd& operator=(const ForwardAlgorithmHSMMEnd&);
    ForwardAlgorithmHSMMEnd& operator=(ForwardAlgorithmHSMMEnd&);
};
    
    
template<class Algo>
class ForwardAlgorithmHSMMBegin : public UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, ForwardAlgorithmHSMMBegin<Algo> >
{
protected:
    typedef ForwardAlgorithmHSMMBegin<Algo> self_type;
    typedef typename Algo::probability_type P;
    typedef UnaryFunc<typename Algo::probability_type, boost::fusion::vector<size_t,size_t>, ForwardAlgorithmHMM<Algo> > base_type;
    typedef typename detail::get_forward_env<Algo>::type env_type;
    typedef detail::get_forward_env<Algo> env_;
    
public:
    typedef typename base_type::return_type return_type;
    typedef typename base_type::arg_type arg_type;  
    
protected:
    typedef MemoizedFunc<self_type, typename detail::get_forward_env<Algo>::type, arg_type> host_type;
    host_type& _host;
    const boost::shared_ptr<const typename Algo::sequence_type> _seq;
    const boost::shared_ptr<const typename Algo::model_type> _model;
    
protected:
    boost::shared_ptr<typename Algo::forward_algo_type> _forw_end;
public:
    void set_end_func( typename Algo::forward_algo_type& func ) {   _forw_end = boost::shared_ptr<typename Algo::forward_algo_type>(&func, null_deleter());         }
    void set_end_func( boost::shared_ptr<typename Algo::forward_algo_type> func ) { _forw_end = func; }
    
protected:
    ForwardAlgorithmHSMMBegin( host_type& host, const env_type& env)
        : _host(host) 
        , _seq(boost::get<env_::seq>(env))
        , _model(boost::get<env_::model>(env))
    {}
    
protected:      
    P _produce_fli( const arg_type& args )
    {
        const size_t l = boost::fusion::at_c<0>(args);
        const size_t i = boost::fusion::at_c<1>(args);
        
        // Base case - i==0
        if( i+1 == 1 )
        {
            // What's the probability of starting with state=l at position 1?
            return _model->a(_model->GetInitialState(), l);
        }
        else if( l == _model->GetInitialState()) return -std::numeric_limits<P>::max();
        
        //std::cout<< "Forw_Begin (l="<< l<< ", i="<< i<< ")"<< std::endl;
        
        const size_t lastpos =  length(*_seq)+1;
        // Constraint - ending at any state other than the terminal state is not allowed
        if( (i+1 == lastpos) && (l != _model->GetTerminalState()) ) return -std::numeric_limits<P>::max();
        // Constraint - the terminal state cannot occur at any state other than the terminal state
        if( (i+1 <  lastpos) && (l == _model->GetTerminalState()) ) return -std::numeric_limits<P>::max();
        
        
        // Rabiner eq. (74)
        P sum = -std::numeric_limits<P>::max();
        for( typename Algo::model_type::StateId k = 0; k< num_states(*_model); ++k )
        {
            const P prev = (*_forw_end)( typename Algo::forward_algo_type::arg_type(k, i) ) + _model->a(k, l);
            
            sum = logspace_sum( sum,
                                prev );
        }
        
        return sum;
    }       
    
public:
    inline return_type operator()( const arg_type& args )           {       return _produce_fli(args);      }
    
private:
    ForwardAlgorithmHSMMBegin(const ForwardAlgorithmHSMMBegin&);
    const ForwardAlgorithmHSMMBegin& operator=(const ForwardAlgorithmHSMMBegin&);
    ForwardAlgorithmHSMMBegin& operator=(ForwardAlgorithmHSMMBegin&);
};

namespace detail
{
template<typename Algo>
struct can_calculate_sequence_probability<ForwardAlgorithmHSMMEnd<Algo> >
{
    typedef boost::mpl::bool_<true> value;
};
}

} // namespace forward_algo
template<class Algo>
struct FuncTraits<forward_algo::ForwardAlgorithmHMM<Algo> >
{
    typedef typename forward_algo::ForwardAlgorithmHMM<Algo>::arg_type arg_type;
    typedef typename forward_algo::ForwardAlgorithmHMM<Algo>::return_type val_type;
    
    struct empty_val
    {
        val_type operator()() const { return 1; }
    };
};
template<class Algo>
struct
FuncTraits<forward_algo::ForwardAlgorithmHSMMBegin<Algo> >
{
    typedef typename forward_algo::ForwardAlgorithmHSMMBegin<Algo>::arg_type arg_type;
    typedef typename forward_algo::ForwardAlgorithmHSMMBegin<Algo>::return_type val_type;
    
    struct empty_val
    {
        val_type operator()() const { return 1; }
    };
};
template<class Algo>
struct
FuncTraits<forward_algo::ForwardAlgorithmHSMMEnd<Algo> >
{
    typedef typename forward_algo::ForwardAlgorithmHSMMEnd<Algo>::arg_type arg_type;
    typedef typename forward_algo::ForwardAlgorithmHSMMEnd<Algo>::return_type val_type;
    
    struct empty_val
    {
        val_type operator()() const { return 1; }
    };
};

/*

template<class Algo>
struct FuncTraits<forward_algo::ForwardAlgorithmHMM<Algo> >
{
        enum {out_of_range_value = 1};
};
template<class Algo>
struct
FuncTraits<forward_algo::ForwardAlgorithmHSMMBegin<Algo> >
{
        enum {out_of_range_value = 1};
};
template<class Algo>
struct
FuncTraits<forward_algo::ForwardAlgorithmHSMMEnd<Algo> >
{
        enum {out_of_range_value = 1};
};
*/



namespace v2
{

/*
namespace tag {
  struct forward {};
}

namespace tag {
  struct backward_begin {};
}
*/



namespace tag {
struct forward : public algorithm {};
struct forward_begin : public algorithm {};
}


namespace detail
{



template<typename Model>
class forward_hsmm :
    public Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t, time_t>
      , Model
    >
  , public memoize<
        forward_hsmm<Model>
      , boost::fusion::vector<state_t, time_t>
      , boost::fusion::vector<probability_t>
      , sequential_processing::sole_sequence_scope // TODO: Switch back to single_sequence_scope (i.e. part of the multi_sequence_scope)
    >
  , public applies_to_data<
        sequential_processing::sole_sequence_scope
    >
  , public depends_on<
        boost::mpl::set<tag::forward_begin>
      , sequential_processing::sole_sequence_scope
      , tag::forward
  >
{
    BOOST_MPL_ASSERT((boost::is_same<typename model::model_traits<Model>::hidden_process, model::semi_markov_chain>));


public:
    typedef Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t, time_t>
      , Model
    > base;

    typedef typename base::result_type result_type;
    typedef typename base::arg_type arg_type;
    typedef tag::forward tag;
    typedef probability_t P;
    
    template<typename Comp>
    result_type produce(Comp& comp, const arg_type& args)
    {

        const state_t l = boost::fusion::at_c<0>(args);
        const time_t  i = boost::fusion::at_c<1>(args);
        
        //std::cout<< "Forw_End   (l="<< l<< ", i="<< i<< ")"<< std::endl;
        
        // At i=0, ai=1.0 if i is the initial state, and 0.0 otherwise
        if( i == 0 )
        {
            return (l == this->model_c(comp).GetInitialState()) ? 0 : -std::numeric_limits<P>::max();
        }
        else
        {
            if( l == this->model_c(comp).GetInitialState() ) 
                return -std::numeric_limits<P>::max();
        }
        
        const size_t lastpos =  length(data_c(comp))+1;
        
        // At i=lastpos, ai=0.0 if i is not the terminal state (and calculated as usual otherwise)
        if( (i == lastpos) && (l != this->model_c(comp).GetTerminalState()) ) return -std::numeric_limits<P>::max();
        if( (i <  lastpos) && (l == this->model_c(comp).GetTerminalState()) ) return -std::numeric_limits<P>::max();
        
        // Sum over all states k
        P sum = -std::numeric_limits<P>::max();
        for( size_t d=1; d<=base::model_type::max_duration; ++d )
        {
            if( d>=i+1 ) continue; // TODO - CHECK THIS

            // Rabiner eq. (75)
            P prev = boost::fusion::at_c<0>( apply(comp, v2::tag::forward_begin(), typename Comp::template algo<v2::tag::forward_begin>::type::arg_type(l,i-d) ) )
                + this->model_c(comp).p(l, d);
          
            // v1 impl:
            //P prev = (*_forw_begin)( arg_type(l, i-d) )
            //    + _model->p(l, d);

            
            // Calculate sigmab, the contribution of the emissions over all d symbols
            if(!this->model_c(comp).is_silent(l))  // silent states have no emissions probs
            {
                P sigmab = 0.0; //-std::numeric_limits<P>::max();
                for( size_t s = i-d+1; s <= i; ++s )
                {                                       
                    const P esym = this->model_c(comp).e( l, data_c(comp)[s-1] ); // TODO - Check the []!
                    //std::cout<< " e('"<< data_c(comp)[s-1]<< "')="<< esym;
                    sigmab += esym;
                }
                
                
                prev += sigmab;
                //std::cout<< ", sigmab="<< sigmab;
            }
            //std::cout<< std::endl;
            
            sum = logspace_sum( sum,
                                prev );
        }
        
        //if( l==0 && i==1)
        //{
        //    //std::cout<< "Forw_End   (l="<< l<< ", i="<< i<< ", sum="<< sum<< std::endl;
        //}
        return sum;
    }


public:
  template<typename Comp>
  boost::fusion::vector<size_t,size_t> get_extents(Comp& comp)
  {
      //std::cout<< "get_extents "<< length(data(comp))+1<< "  "<< num_states( this->model(comp) )<< std::endl;
      return boost::fusion::vector<size_t,size_t>(
          num_states( this->model_c(comp) )
          , fasta::get_len(data_c(comp))+2  /* TODO - Find real length! */
          );
  }    

public:
    static result_type get_unassigned_value() { return result_type(1.0 /* Assumption: probability_t(1.0) can safely be compared exactly (without requiring a nonzero tolerance). */ ); }

public:
    const char* get_debug_id() const { return "forward[hsmm]"; }

public:
    void reset(size_t level=0) {};

}; // class forward_hsmm



template<typename Model>
class forward_begin_hsmm :
    public Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t, time_t>
      , Model
    >
  , public memoize<
        forward_begin_hsmm<Model>
      , boost::fusion::vector<state_t, time_t>
      , boost::fusion::vector<probability_t>
      , sequential_processing::sole_sequence_scope // TODO: Switch back to single_sequence_scope (i.e. part of the multi_sequence_scope)
    >
  , public applies_to_data<
        sequential_processing::sole_sequence_scope
    >
  , public depends_on<
        boost::mpl::set<tag::forward>
      , sequential_processing::sole_sequence_scope
      , tag::forward_begin
  >
{
    BOOST_MPL_ASSERT((boost::is_same<typename model::model_traits<Model>::hidden_process, model::semi_markov_chain>));


public:
    typedef Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t, time_t>
      , Model
    > base;

    typedef typename base::result_type result_type;
    typedef typename base::arg_type arg_type;
    typedef tag::forward_begin tag;
    typedef probability_t P;
    
    template<typename Comp>
    result_type produce(Comp& comp, const arg_type& args)
    {

        const state_t l = boost::fusion::at_c<0>(args);
        const time_t  i = boost::fusion::at_c<1>(args);
        
        // Base case - i==0
        if( i+1 == 1 )
        {
            // What's the probability of starting with state=l at position 1?
            return this->model_c(comp).a(this->model_c(comp).GetInitialState(), l);
        }
        else if( l == this->model_c(comp).GetInitialState()) return -std::numeric_limits<P>::max();
        
        //std::cout<< "Forw_Begin (l="<< l<< ", i="<< i<< ")"<< std::endl;
        
        const size_t lastpos =  length(data_c(comp))+1;
        // Constraint - ending at any state other than the terminal state is not allowed
        if( (i+1 == lastpos) && (l != this->model_c(comp).GetTerminalState()) ) return -std::numeric_limits<P>::max();
        // Constraint - the terminal state cannot occur at any state other than the terminal state
        if( (i+1 <  lastpos) && (l == this->model_c(comp).GetTerminalState()) ) return -std::numeric_limits<P>::max();
        
        
        // Rabiner eq. (74)
        P sum = -std::numeric_limits<P>::max();
        for( typename base::model_type::StateId k = 0; k< num_states(this->model_c(comp)); ++k )
        {

            P prev = boost::fusion::at_c<0>( apply(comp, v2::tag::forward(), typename Comp::template algo<v2::tag::forward>::type::arg_type(k,i) ) )
                + this->model_c(comp).a(k, l);

            //const P prev = (*_forw_end)( typename Algo::forward_algo_type::arg_type(k, i) ) + _model->a(k, l);
            
            sum = logspace_sum( sum,
                                prev );
        }
        
        return sum;
    }


public:
  template<typename Comp>
  boost::fusion::vector<size_t,size_t> get_extents(Comp& comp)
  {
      //std::cout<< "get_extents "<< length(data(comp))+1<< "  "<< num_states( this->model(comp) )<< std::endl;
      return boost::fusion::vector<size_t,size_t>(
          num_states( this->model_c(comp) )
          , fasta::get_len(data_c(comp))+2  /* TODO - Find real length! */
          );
  }    

public:
    static result_type get_unassigned_value() { return result_type(1.0 /* Assumption: probability_t(1.0) can safely be compared exactly (without requiring a nonzero tolerance). */ ); }

public:
    const char* get_debug_id() const { return "forward_begin[hsmm]"; }

public:
    void reset(size_t level=0) {};

}; // class forward_begin_hsmm


} // namespace detail


// Register this implementation

template <class Model>
struct implementations<
    tag::forward
  , Model
  , typename boost::enable_if<
      typename boost::is_same<
          typename model::model_traits<Model>::hidden_process
        , model::semi_markov_chain
        >::type
    >::type
  >
{
  typedef typename detail::forward_hsmm<Model> impl_type;
};


template <class Model>
struct implementations<
    tag::forward_begin
  , Model
  , typename boost::enable_if<
      typename boost::is_same<
          typename model::model_traits<Model>::hidden_process
        , model::semi_markov_chain
        >::type
    >::type
  >
{
  typedef typename detail::forward_begin_hsmm<Model> impl_type;
};


} // namespace v2

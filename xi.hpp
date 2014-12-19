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
#include "common.hpp"
#include "v2.hpp"
#include "algo.hpp"
#include "sequential_processing.hpp"
//#include "memoized.hpp"
#include "computation.hpp"
//#include "read_fasta.hpp"
#include "forward.hpp"
#include "backward.hpp"


namespace v2
{




namespace tag {
struct xi : public algorithm {};
struct sigma_i_xi : public algorithm {};
struct sigma_li_xi : public algorithm {};
}


namespace detail
{



template<typename Model>
class xi_hsmm :
    public Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t, state_t, time_t>
      , Model
    >
  , public memoize<
        xi_hsmm<Model>
      , boost::fusion::vector<state_t, state_t, time_t>
      , boost::fusion::vector<probability_t>
      , sequential_processing::sole_sequence_scope // TODO: Switch back to single_sequence_scope (i.e. part of the multi_sequence_scope)
    >
  , public applies_to_data<
        sequential_processing::sole_sequence_scope
    >
  , public depends_on<
        boost::mpl::set<tag::forward, tag::backward_begin> 
      , sequential_processing::sole_sequence_scope
      , tag::xi
  >
{
    BOOST_MPL_ASSERT((boost::is_same<typename model::model_traits<Model>::hidden_process, model::semi_markov_chain>));


public:
    typedef Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t, state_t, time_t>
      , Model
    > base;

    typedef typename base::result_type result_type;
    typedef typename base::arg_type arg_type;
    typedef tag::xi tag;
    typedef probability_t P;
    
    template<typename Comp>
    result_type produce(Comp& comp, const arg_type& args)
    {
        const state_t k = boost::fusion::at_c<0>(args);
        const state_t l = boost::fusion::at_c<1>(args);
        const time_t  i = boost::fusion::at_c<2>(args);

        //const typename Algo::forward_algo_type::arg_type k_i(k, i);
        //const typename Algo::backward_begin_algo_type::arg_type l_i_1(l, i);
        //const logspace_t f((*_forward)(j_t), true);
                
        // Rabiner eq. (79)
        const P p1 = boost::fusion::at_c<0>( apply(comp, v2::tag::forward(), typename Comp::template algo<v2::tag::forward>::type::arg_type(k,i) ) );
        const P p2 = this->model_c(comp).a(k, l);
        const P p3 = boost::fusion::at_c<0>( apply(comp, v2::tag::backward_begin(), typename Comp::template algo<v2::tag::backward_begin>::type::arg_type(l,i) ) );

        assert_range_probability_logspace(p1);
        assert_range_probability_logspace(p2);
        assert_range_probability_logspace(p3);

        //const P p1 = (*_forward)(k_i);
        //const P p2 = a_k_l;
        //const P p3 = (*_backward_begin)(l_i_1);
        
        //std::cout<< "xi("<< k<< "->"<< l<< ", i="<< i<< ") = "<< p1<< " + "<< p2<< " + "<< p3<< " = "<< p1+p2+p3<< std::endl;

        return p1 + p2 + p3; // logspace
    }


public:
  template<typename Comp>
  boost::fusion::vector<size_t,size_t,size_t> get_extents(Comp& comp)
  {
      //std::cout<< "get_extents "<< length(data(comp))+1<< "  "<< num_states( this->model(comp) )<< std::endl;
      return boost::fusion::vector<size_t,size_t,size_t>(
          num_states( this->model(comp) )
        , num_states( this->model(comp) )
        , fasta::get_len(data(comp))+1
          );
  }    

public:
    static result_type get_unassigned_value() { return result_type(1.0 /* Assumption: probability_t(1.0) can safely be compared exactly (without requiring a nonzero tolerance). */ ); }

public:
    const char* get_debug_id() const { return "xi[hsmm]"; }

public:
    void reset(size_t level=0) {};

}; // class xi_hsmm

} // namespace detail

// Register this implementation

template <class Model>
struct implementations<
    tag::xi
  , Model
  , typename boost::enable_if<
      typename boost::is_same<
          typename model::model_traits<Model>::hidden_process
        , model::semi_markov_chain
        >::type
    >::type
  >
{
  typedef typename detail::xi_hsmm<Model> impl_type;
};


namespace detail
{

template<typename Model>
class sigma_i_xi_hsmm :
    public Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t, state_t>
      , Model
    >
  , public memoize<
        sigma_i_xi_hsmm<Model>
      , boost::fusion::vector<state_t, state_t>
      , boost::fusion::vector<probability_t>
      , sequential_processing::sole_sequence_scope // TODO: Switch back to single_sequence_scope (i.e. part of the multi_sequence_scope)
    >
  , public applies_to_data<
        sequential_processing::sole_sequence_scope
    >
  , public depends_on<
        boost::mpl::set<tag::xi>
      , sequential_processing::sole_sequence_scope
      , tag::sigma_i_xi
  >
{
    BOOST_MPL_ASSERT((boost::is_same<typename model::model_traits<Model>::hidden_process, model::semi_markov_chain>));


public:
    typedef Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t, state_t>
      , Model
    > base;

    typedef typename base::result_type result_type;
    typedef typename base::arg_type arg_type;
    typedef tag::sigma_i_xi tag;
    typedef probability_t P;
    
    template<typename Comp>
    result_type produce(Comp& comp, const arg_type& args)
    {
        const state_t k = boost::fusion::at_c<0>(args);
        const state_t l = boost::fusion::at_c<1>(args);
        const time_t lastpos = length(data_c(comp))+1; // TODO - Verify boundary


        //const typename Algo::forward_algo_type::arg_type k_i(k, i);
        //const typename Algo::backward_begin_algo_type::arg_type l_i_1(l, i);
        //const logspace_t f((*_forward)(j_t), true);
                
        // Rabiner eq. (79)
        P sigma_i = 0.0;
        for( size_t i=0; i< lastpos; ++i )
        {
            //if( i > 0 )
            //{
            //const P p0 = boost::fusion::at_c<0>( apply(comp, v2::tag::sigma_i_xi(), typename Comp::template algo<v2::tag::sigma_i_xi>::type::arg_type(k,l,i-1) ) );
            //assert_range_probability_realspace(p0);
            //}

            const P p1 = boost::fusion::at_c<0>( apply(comp, v2::tag::xi(), typename Comp::template algo<v2::tag::xi>::type::arg_type(k,l,i) ) );
            assert_range_probability_logspace(p1);
            sigma_i += exp(p1);
        }
            
        //const P p1 = (*_forward)(k_i);
        //const P p2 = a_k_l;
        //const P p3 = (*_backward_begin)(l_i_1);
        
        //std::cout<< "xi("<< k<< "->"<< l<< ", i="<< i<< ") = "<< p1<< " + "<< p2<< " + "<< p3<< " = "<< p1+p2+p3<< std::endl;

        return sigma_i;
    }


public:
  template<typename Comp>
  boost::fusion::vector<size_t,size_t> get_extents(Comp& comp)
  {
      //std::cout<< "get_extents "<< length(data(comp))+1<< "  "<< num_states( this->model(comp) )<< std::endl;
      return boost::fusion::vector<size_t,size_t>(
          num_states( this->model(comp) )
        , num_states( this->model(comp) )
          );
  }    

public:
    static result_type get_unassigned_value() { return result_type(-1.0 /* Assumption: probability_t(-1.0) can safely be compared exactly (without requiring a nonzero tolerance). */ ); }

public:
    const char* get_debug_id() const { return "sigma_i_xi[hsmm]"; }

public:
    void reset(size_t level=0) {};

}; // class sigma_i_xi_hsmm


} // namespace detail



// Register this implementation

template <class Model>
struct implementations<
    tag::sigma_i_xi
  , Model
  , typename boost::enable_if<
      typename boost::is_same<
          typename model::model_traits<Model>::hidden_process
        , model::semi_markov_chain
        >::type
    >::type
  >
{
  typedef typename detail::sigma_i_xi_hsmm<Model> impl_type;
};


namespace detail
{

template<typename Model>
class sigma_li_xi_hsmm :
    public Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t>
      , Model
    >
  , public memoize<
        sigma_li_xi_hsmm<Model>
      , boost::fusion::vector<state_t>
      , boost::fusion::vector<probability_t>
      , sequential_processing::sole_sequence_scope // TODO: Switch back to single_sequence_scope (i.e. part of the multi_sequence_scope)
    >
  , public applies_to_data<
        sequential_processing::sole_sequence_scope
    >
  , public depends_on<
        boost::mpl::set<tag::sigma_i_xi>
      , sequential_processing::sole_sequence_scope
      , tag::sigma_li_xi
  >
{
    BOOST_MPL_ASSERT((boost::is_same<typename model::model_traits<Model>::hidden_process, model::semi_markov_chain>));


public:
    typedef Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t>
      , Model
    > base;

    typedef typename base::result_type result_type;
    typedef typename base::arg_type arg_type;
    typedef tag::sigma_li_xi tag;
    typedef probability_t P;
    
    template<typename Comp>
    result_type produce(Comp& comp, const arg_type& args)
    {
        const state_t k = boost::fusion::at_c<0>(args);

        //const typename Algo::forward_algo_type::arg_type k_i(k, i);
        //const typename Algo::backward_begin_algo_type::arg_type l_i_1(l, i);
        //const logspace_t f((*_forward)(j_t), true);
                
        // Rabiner eq. (79)

        P sigma_l = 0.0;
        
        for( state_t l = 0; l < num_states(this->model_c(comp)); ++l )
        {
            const P p0 = boost::fusion::at_c<0>( apply(comp, v2::tag::sigma_i_xi(), typename Comp::template algo<v2::tag::sigma_i_xi>::type::arg_type(k,l) ) );
            assert_range_probability_realspace(p0);

            sigma_l += p0;
        }

        //const P p1 = (*_forward)(k_i);
        //const P p2 = a_k_l;
        //const P p3 = (*_backward_begin)(l_i_1);
        
        //std::cout<< "xi("<< k<< "->"<< l<< ", i="<< i<< ") = "<< p1<< " + "<< p2<< " + "<< p3<< " = "<< p1+p2+p3<< std::endl;

        return sigma_l;
    }


public:
  template<typename Comp>
  boost::fusion::vector<size_t> get_extents(Comp& comp)
  {
      //std::cout<< "get_extents "<< length(data(comp))+1<< "  "<< num_states( this->model(comp) )<< std::endl;
      return boost::fusion::vector<size_t>(
          num_states( this->model(comp) )
          );
  }    

public:
    static result_type get_unassigned_value() { return result_type(-1.0 /* Assumption: probability_t(-1.0) can safely be compared exactly (without requiring a nonzero tolerance). */ ); }

public:
    const char* get_debug_id() const { return "sigma_li_xi[hsmm]"; }

public:
    void reset(size_t level=0) {};

}; // class sigma_li_xi_hsmm


} // namespace detail



// Register this implementation

template <class Model>
struct implementations<
    tag::sigma_li_xi
  , Model
  , typename boost::enable_if<
      typename boost::is_same<
          typename model::model_traits<Model>::hidden_process
        , model::semi_markov_chain
        >::type
    >::type
  >
{
  typedef typename detail::sigma_li_xi_hsmm<Model> impl_type;
};


} // namespace v2

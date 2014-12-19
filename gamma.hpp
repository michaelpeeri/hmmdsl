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
#include "computation.hpp"
#include "forward.hpp"
#include "backward.hpp"



namespace v2
{

namespace tag {
    struct gamma : public algorithm {};
}


namespace detail
{



template<typename Model>
class gamma_hsmm :
    public Algo<
        boost::fusion::vector<probability_t>
      , boost::fusion::vector<state_t, time_t>
      , Model
    >
  , public memoize<
        gamma_hsmm<Model>
      , boost::fusion::vector<state_t, time_t>
      , boost::fusion::vector<probability_t>
      , sequential_processing::sole_sequence_scope // TODO: Switch back to single_sequence_scope (i.e. part of the multi_sequence_scope)
    >
  , public applies_to_data<
        sequential_processing::sole_sequence_scope
    >
  , public depends_on<
    boost::mpl::set<tag::forward, tag::forward_begin, tag::backward, tag::backward_begin> 
      , sequential_processing::sole_sequence_scope
      , tag::gamma
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
    typedef tag::gamma tag;
    typedef probability_t P;
    
    template<typename Comp>
    result_type produce(Comp& comp, const arg_type& args) // TODO - Use recurrence to optimize!
    {
        const state_t l = boost::fusion::at_c<0>(args);
        const time_t  i = boost::fusion::at_c<1>(args);
        
        // Currently, gamma is not supported for initial and terminal states.
        if( i == 0 )
            return (l == this->model_c(comp).GetInitialState() ) ? 1.0 : 0.0;
        if( i == length(data_c(comp))+1 )
            return (l == this->model_c(comp).GetTerminalState() ) ? 1.0 : 0.0;

        // If we reached this point, occupying an initial or terminal state is impossible
        if( this->model_c(comp).IsReservedState(l) ) return 0.0;
        
        //const size_t lastpos = length(*_seq)+1;
        //if( i == lastpos )
        //	return (l==_model->GetTerminalState()) ? 0 : -std::numeric_limits<P>::max();

        //std::cout<< "gam("<< l<< ", "<< i<< ") ";
        //const arg_type prev(l, i-1);

        // TODO - memoize this!
        
        P P_l_begins_before_i = 0.0;
        for( size_t tau=0; tau< i; ++tau)
        {
            //const arg_type l_tau(l, tau);

            //const P p1 = (*_forward_begin)(l_tau);
            //if( p1 > -std::numeric_limits<P>::max() ) std::cout<< l_tau<< "[1] ";
            //const P p2 = (*_backward_begin)(l_tau);
            //if( p2 > -std::numeric_limits<P>::max() ) std::cout<< l_tau<< "[2] ";
            
            const P P_l_begins_before_tau = exp(
                boost::fusion::at_c<0>( apply(comp, v2::tag::forward_begin(),  typename Comp::template algo<v2::tag::forward_begin >::type::arg_type(l,tau) ) )
              + boost::fusion::at_c<0>( apply(comp, v2::tag::backward_begin(), typename Comp::template algo<v2::tag::backward_begin>::type::arg_type(l,tau) ) )
                );
            // v1 impl.:
            //const P P_l_begins_before_tau = exp( (*_forward_begin)(l_tau) + (*_backward_begin)(l_tau) );

            P_l_begins_before_i += P_l_begins_before_tau;
        }

        P P_l_ends_before_i = 0.0;
        for( size_t tau=0; tau< i; ++tau)
        {
            //const arg_type l_tau(l, tau);

            const P P_l_ends_before_tau = exp(
                boost::fusion::at_c<0>( apply(comp, v2::tag::forward(),  typename Comp::template algo<v2::tag::forward >::type::arg_type(l,tau) ) )
              + boost::fusion::at_c<0>( apply(comp, v2::tag::backward(), typename Comp::template algo<v2::tag::backward>::type::arg_type(l,tau) ) )
                );
            // v1 impl.:
            //const P P_l_ends_before_tau = exp( (*_forward      )(l_tau) + (*_backward      )(l_tau) );

            P_l_ends_before_i += P_l_ends_before_tau;
        }

        const P P_l_occurs_at_i = P_l_begins_before_i - P_l_ends_before_i;

        //std::cout<< "gamma(l="<< l<< ", i="<< i<< ") = "<<P_l_occurs_at_i << " (P_l_before_i="<< P_l_begins_before_i<< ", P_l_ends_before_i="<< P_l_ends_before_i<< ")"<< std::endl;
        
        if( P_l_occurs_at_i< -1e-10 )
        {
            std::cout<< "Error: Got gamma(l,i)<0 for l="<< l<< ", i="<< i<< ", P_l_before_i="<< P_l_begins_before_i<< ", P_l_ends_before_i="<< P_l_ends_before_i<< std::endl;
            return 0.0;
        }
        return P_l_occurs_at_i >= 0.0 ? P_l_occurs_at_i : 0.0;
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
    static result_type get_unassigned_value() { return result_type(-1.0 /* Assumption: probability_t(-1.0) can safely be compared exactly (without requiring a nonzero tolerance). */ ); }

public:
    const char* get_debug_id() const { return "gamma[hsmm]"; }

public:
    void reset(size_t level=0) {};

}; // class gamma_hsmm

} // namespace detail

// Register this implementation

template <class Model>
struct implementations<
    tag::gamma
  , Model
  , typename boost::enable_if<
      typename boost::is_same<
          typename model::model_traits<Model>::hidden_process
        , model::semi_markov_chain
        >::type
    >::type
  >
{
  typedef typename detail::gamma_hsmm<Model> impl_type;
};


} // namespace v2


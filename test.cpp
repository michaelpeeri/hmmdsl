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
#include <iostream>
#include <vector>
#include <boost/assign/std/vector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include "common.hpp"
#include "model.hpp"
#include "viterbi.hpp"
#include "forward.hpp"
#include "backward.hpp"
#include "baumwelch.hpp"
#include "tests.hpp"

template<typename Model>
void InitializeOccasionallyDishonestCasino( Model& model)
{
  enum { Fair=0, Loaded=1 };
  model.AddState( Fair,   "F"    );
  model.AddState( Loaded, "L"    );
  model.SetTransition( Fair,   Fair,   0.95 );
  model.SetTransition( Fair,   Loaded, 0.05 );
  model.SetTransition( Loaded, Fair,   0.10 );
  model.SetTransition( Loaded, Loaded, 0.90 );
  model.AddAlphabetSymbol( "1" );
  model.AddAlphabetSymbol( "2" );
  model.AddAlphabetSymbol( "3" );
  model.AddAlphabetSymbol( "4" );
  model.AddAlphabetSymbol( "5" );
  model.AddAlphabetSymbol( "6" );
  model.SetEmissionProbability( Fair, "1", 1./6. );
  model.SetEmissionProbability( Fair, "2", 1./6. );
  model.SetEmissionProbability( Fair, "3", 1./6. );
  model.SetEmissionProbability( Fair, "4", 1./6. );
  model.SetEmissionProbability( Fair, "5", 1./6. );
  model.SetEmissionProbability( Fair, "6", 1./6. );
  model.SetEmissionProbability( Loaded, "1", 1./10. );
  model.SetEmissionProbability( Loaded, "2", 1./10. );
  model.SetEmissionProbability( Loaded, "3", 1./10. );
  model.SetEmissionProbability( Loaded, "4", 1./10. );
  model.SetEmissionProbability( Loaded, "5", 1./10. );
  model.SetEmissionProbability( Loaded, "6", 1./2. );

  //enum { Init=0, Fair=1, Loaded=2, End=3 };
  //model.AddState( Init,   "Init" );
  //model.AddState( End,    "End"  );
  //model.AddTransition( Init,   Fair,   0.51  );
  //model.AddTransition( Init,   Loaded, 0.52  );
  //model.AddTransition( Fair,   End,    0.53  );
  //model.AddTransition( Loaded, End,    0.54  );
  //  model.SetStartState( Init );
  //  model.SetEndState( End );
}



template<typename Algo>
double get_prob( const typename Algo::sequence_type& seq, const typename Algo::model_type& model )
{
  typename Algo::forward_algo_type forward(seq, model);
  return forward.calc();
}

template<typename Algo>
void initialize_random_sequence( typename Algo::sequence_type& seq, size_t length )
{
  seq.clear();
  for( size_t n = length; n>0; --n )
    {
      switch (rand() % 6 )
	{
	case 0:	        seq.push_back("1") ; break;
	case 1:	        seq.push_back("2") ; break;
	case 2:	        seq.push_back("3") ; break;
	case 3:	        seq.push_back("4") ; break;
	case 4:	        seq.push_back("5") ; break;
	case 5:	        seq.push_back("6") ; break;
	default: assert(false);
	}
    }
}

template<typename Algo>
void multitest(const typename Algo::model_type& model)
{
  typename Algo::probability_type pbest = -std::numeric_limits<typename Algo::probability_type>::max();
  typename Algo::sequence_type seqbest;

  typename Algo::probability_type pworst = 1;
  typename Algo::sequence_type seqworst;

  bool printnow = false;

  for( size_t n =0; n< 10000000; ++n )
  {
    typename Algo::sequence_type seq;
    initialize_random_sequence<Algo>( seq, 20 );
    
    typename Algo::probability_type p = get_prob<Algo>( seq, model );

    if( p > pbest )
      {
	pbest = p;
	seqbest = seq;
	printnow = true;
      }

    if( p < pworst )
      {
	pworst = p;
	seqworst = seq;
	printnow = true;
      }


    if( ((n % 50000) == 0) && (n>0) ) printnow = true;

    if( printnow )
      {
	std::cout<< "iter "<< n<< " +: "<< stringify( seqbest,  0, 30 )<< "\t"<< pbest<< std::endl;
	std::cout<< "iter "<< n<< " -: "<< stringify( seqworst,  0, 30 )<< "\t"<< pworst<< std::endl;
	printnow = false;
      }
  }
}


template<typename Algo>
void posterior(typename Algo::forward_algo_type& f, typename Algo::backward_algo_type& b, const typename Algo::sequence_type& seq)
{
  const double Pseq = f.calc();

  std::cout<< "#i p0 p1"<< std::endl;

  for( size_t i=0; i< length(seq); ++i )
  {
    const double p0 = exp( f.fli(0/*Fair*/  , i) + b.bli(0/*Fair*/  , i) - Pseq );
    const double p1 = exp( f.fli(1/*Loaded*/, i) + b.bli(1/*Loaded*/, i) - Pseq );
    std::cout<< i<< "\t"<< p0<< "\t"<< p1<< std::endl;
  }
}

int main(int argc, char* argv[])
{
  using namespace boost::assign;
  using boost::shared_ptr;
  typedef
    HiddenMarkovAlgorithm<
	std::vector<std::string>,
		double_logscale,
		MatrixModel,
		ViterbiAlgorithm,
		ForwardAlgorithm,
		BackwardAlgorithm,
		BaumWelchAlgorithm> algo_t;

  shared_ptr<algo_t::model_type> model( new algo_t::model_type() );
  InitializeOccasionallyDishonestCasino(*model);

  shared_ptr<algo_t::sequence_type> seq( new algo_t::sequence_type() );
  //seq += "1", "2", "3", "4", "5", "6", "1", "2", "3", "5", "3", "4", "6", "6", "6", "2", "6", "6", "6", "6", "6", "3", "6", "6", "2", "4", "1", "5", "3", "3", "2", "6", "2", "4", "1", "3", "6", "6", "6", "3", "6", "6", "2", "6", "6", "6", "3", "3", "1", "6", "2", "4", "1", "1", "2", "5", "3", "1", "5", "4", "6", "6", "1", "6", "2", "6", "3", "6", "6", "5", "6", "2", "3", "6", "4", "1", "6", "1", "5", "6", "2", "6", "3", "1", "3", "3", "6", "2", "5", "6", "1", "1", "6", "2", "5", "6", "3", "1", "6", "2", "4", "2", "6", "6", "1", "6", "3", "6", "3", "3", "1", "6", "2", "5", "1";
  
  (*seq) += "3", "1", "5", "1", "1", "6", "2", "4", "6", "4", "4", "6", "6", "4", "4", "2", "4", "5", "3", "1", "1", "3", "2", "1", "6", "3", "1", "1", "6", "4", "1", "5", "2", "1", "3", "3", "6", "2", "5", "1", "4", "4", "5", "4", "3", "6", "3", "1", "6", "5", "6", "6", "2", "6", "5", "6", "6", "6", "6", "6", "6", "5", "1", "1", "6", "6", "4", "5", "3", "1", "3", "2", "6", "5", "1", "2", "4", "5", "6", "3", "6", "6", "6", "4", "6", "3", "1", "6", "3", "6", "6", "6", "3", "1", "6", "2", "3", "2", "6", "4", "5", "5", "2", "3", "6", "2", "6", "6", "6", "6", "6", "6", "2", "5", "1", "5", "1", "6", "3", "1", "2", "2", "2", "5", "5", "5", "4", "4", "1", "6", "6", "6", "5", "6", "6", "5", "6", "3", "5", "6", "4", "3", "2", "4", "3", "6", "4", "1", "3", "1", "5", "1", "3", "4", "6", "5", "1", "4", "6", "3", "5", "3", "4", "1", "1", "1", "2", "6", "4", "1", "4", "6", "2", "6", "2", "5", "3", "3", "5", "6", "3", "6", "6", "1", "6", "3", "6", "6", "6", "4", "6", "6", "2", "3", "2", "5", "3", "4", "4", "1", "3", "6", "6", "1", "6", "6", "1", "1", "6", "3", "2", "5", "2", "5", "6", "2", "4", "6", "2", "2", "5", "5", "2", "6", "5", "2", "5", "2", "2", "6", "6", "4", "3", "5", "3", "5", "3", "3", "3", "6", "2", "3", "3", "1", "2", "1", "6", "2", "5", "3", "6", "4", "4", "1", "4", "4", "3", "2", "3", "3", "5", "1", "6", "3", "2", "4", "3", "6", "3", "3", "6", "6", "5", "5", "6", "2", "4", "6", "6", "6", "6", "2", "6", "3", "2", "6", "6", "6", "6", "1", "2", "3", "5", "5", "2", "4", "5", "2", "4", "2";
//seq += "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1";
// seq += "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6";


  algo_t::viterbi_algo_type viterbi(seq, model);
 debug_print(*model);

 viterbi.calc();
 debug_print(viterbi);

 std::vector<std::string> bt;
 viterbi.get_backtrace(bt);

 print_seqs( *seq, bt );

 boost::shared_ptr<algo_t::forward_algo_type> forward( new algo_t::forward_algo_type(seq, model) );
 forward->calc();
 debug_print(*forward);
 
 std::cout<< forward->calc()<< std::endl;
 
 boost::shared_ptr<algo_t::backward_algo_type> backward( new algo_t::backward_algo_type(seq, model) );
 backward->calc();
 debug_print(*backward);
 
 std::cout<< backward->calc()<< std::endl;

 tests::test_28<algo_t>( *forward, *backward, *model, *seq );

 posterior<algo_t>(*forward, *backward, *seq);

 BaumWelchAlgorithm<algo_t> bw(model, forward, backward, seq );

 EM<algo_t, BaumWelchAlgorithm<algo_t> > em( *model, bw );
 
 // Set initial values
 model->SetTransition( 0, 0, 0.50 );
 model->SetTransition( 0, 1, 0.50 );
 model->SetTransition( 1, 0, 0.10 );
 model->SetTransition( 1, 1, 0.90 );
 bw.reset();

 for( size_t k = 0; k<200; ++k )
 {
	 std::cout<< "iter "<< k<< "a";
	 
	 //debug_print(*model);

	 tests::test_28<algo_t>( bw, *model, *seq );
	 tests::test_38<algo_t>( bw, *model, *seq );
	 tests::test_43b<algo_t, BaumWelchAlgorithm<algo_t> >( bw, *model );

	 
	 std::cout<< "\t"<< exp(model->a( 0, 0 ));
	 std::cout<< "\t"<< exp(model->a( 0, 1 ));
	 std::cout<< "\t"<< exp(model->a( 1, 0 ));
	 std::cout<< "\t"<< exp(model->a( 1, 1 ));
	 std::cout<< "\t"<< forward->calc();

	 std::cout<< std::endl;

	 em.reestimate_a();
 }
 

 for( size_t k = 0; k<200; ++k )
 {
	 std::cout<< "iter "<< k<< "e";
	 
	 //if( k%10 == 0 ) debug_print(*model);

	 tests::test_28<algo_t>( bw, *model, *seq );
	 tests::test_38<algo_t>( bw, *model, *seq );
	 tests::test_43c<algo_t,false>( *model );

	 
	 //std::cout<< "\t"<< exp(model->a( 0, 0 ));
	 //std::cout<< "\t"<< exp(model->a( 0, 1 ));
	 //std::cout<< "\t"<< exp(model->a( 1, 0 ));
	 //std::cout<< "\t"<< exp(model->a( 1, 1 ));

	 std::cout<< "\t"<< forward->calc();
	 
	 std::cout<< std::endl;

	 em.reestimate_b();
 }

 debug_print(*model);



/*
 typedef std::vector<boost::shared_ptr<algo_t::sequence_type> > seqs_t;
 seqs_t seqs;
 seqs.push_back( boost::shared_ptr<algo_t::sequence_type>( seq  ) );
 MultipleSequenceBW<algo_t> mbw( model, seqs.begin(), seqs.end() );

 tests::test_43b<algo_t, MultipleSequenceBW<algo_t> >( mbw, *model );
*/
 return 0;
}

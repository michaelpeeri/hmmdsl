#include <iostream>
#include <fstream>
#include <vector>
#include <boost/assign/std/vector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/python.hpp>
#include "test_fasta.hpp"
#include "tests.hpp"
#include "read_fasta.hpp"


template<typename Model, typename StateId>
void SetBackgroundEmissionProbs( Model& model, const StateId& state)
{
	// Source: Wikipedia
	model.SetEmissionProbability( state, 'A', 0.078 );
	model.SetEmissionProbability( state, 'C', 0.019 );
	model.SetEmissionProbability( state, 'D', 0.053 );
	model.SetEmissionProbability( state, 'E', 0.063 );
	model.SetEmissionProbability( state, 'F', 0.039 );
	model.SetEmissionProbability( state, 'G', 0.072 );
	model.SetEmissionProbability( state, 'H', 0.023 );
	model.SetEmissionProbability( state, 'I', 0.053 );
	model.SetEmissionProbability( state, 'K', 0.059 );
	model.SetEmissionProbability( state, 'L', 0.091 );
	model.SetEmissionProbability( state, 'M', 0.023 );
	model.SetEmissionProbability( state, 'N', 0.043 );
	model.SetEmissionProbability( state, 'P', 0.052 );
	model.SetEmissionProbability( state, 'Q', 0.042 );
	model.SetEmissionProbability( state, 'R', 0.051 );
	model.SetEmissionProbability( state, 'S', 0.068 );
	model.SetEmissionProbability( state, 'T', 0.059 );
	model.SetEmissionProbability( state, 'V', 0.066 );
	model.SetEmissionProbability( state, 'W', 0.014 );
	model.SetEmissionProbability( state, 'Y', 0.032 );
}

#ifdef MainProgram

template<typename Model>
void InitializeModel( Model& model)
{
	typedef typename Model::StateId StateId;
	
	enum { S0=0, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, N1, P1, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, S21  };
	model.AddState( S0,   "S"    );
	model.AddState( S1,   "S"    );
	model.AddState( S2,   "S"    );
	model.AddState( S3,   "S"    );
	model.AddState( S4,   "S"    );
	model.AddState( S5,   "S"    );
	model.AddState( S6,   "S"    );
	model.AddState( S7,   "S"    );
	model.AddState( S8,   "S"    );
	model.AddState( S9,   "S"    );
	model.AddState( S10,   "S"    );
	model.AddState( N1,    "-"    );
	model.AddState( P1,    "%"    );
	model.AddState( S11,   "S"    );
	model.AddState( S12,   "S"    );
	model.AddState( S13,   "S"    );
	model.AddState( S14,   "S"    );
	model.AddState( S15,   "S"    );
	model.AddState( S16,   "S"    );
	model.AddState( S17,   "S"    );
	model.AddState( S18,   "S"    );
	model.AddState( S19,   "S"    );
	model.AddState( S20,   "S"    );
	model.AddState( S21,   "S"    );
	model.AddAlphabetSymbol( 'E' );
	model.AddAlphabetSymbol( 'D' );
	model.AddAlphabetSymbol( 'R' );
	model.AddAlphabetSymbol( 'H' );
	model.AddAlphabetSymbol( 'K' );
	model.AddAlphabetSymbol( 'I' );
	model.AddAlphabetSymbol( 'L' );
	model.AddAlphabetSymbol( 'V' );
	model.AddAlphabetSymbol( 'S' );
	model.AddAlphabetSymbol( 'T' );
	model.AddAlphabetSymbol( 'Q' );
	model.AddAlphabetSymbol( 'W' );
	model.AddAlphabetSymbol( 'G' );
	model.AddAlphabetSymbol( 'P' );
	model.AddAlphabetSymbol( 'F' );
	model.AddAlphabetSymbol( 'Y' );
	model.AddAlphabetSymbol( 'W' );
	model.AddAlphabetSymbol( 'C' );
	model.AddAlphabetSymbol( 'M' );
	model.AddAlphabetSymbol( 'A' );
	model.AddAlphabetSymbol( 'N' );
	model.SetDefaultProbabilities();
	model.SetTransition( S0,   S1,   1.0 );
	model.SetTransition( S1,   S2,   1.0 );
	model.SetTransition( S2,   S3,   1.0 );
	model.SetTransition( S3,   S4,   1.0 );
	model.SetTransition( S4,   S5,   1.0 );
	model.SetTransition( S5,   S6,   1.0 );
	model.SetTransition( S6,   S7,   1.0 );
	model.SetTransition( S7,   S8,   1.0 );
	model.SetTransition( S8,   S9,   1.0 );
	model.SetTransition( S9,   S9,   0.5 );
	model.SetTransition( S9,   N1,   0.5 );
	model.SetTransition( N1,   N1,   0.333 );
	model.SetTransition( N1,   P1,   0.333 );
	model.SetTransition( N1,   S10,  0.334 );
	model.SetTransition( S10,  S10,  0.5 );
	model.SetTransition( S10,  P1,   0.5 );
	model.SetTransition( P1,   P1,   0.5 );
	model.SetTransition( P1,   S11,  0.5 );
	model.SetTransition( S11,  S11,  0.5 );
	model.SetTransition( S11,  S12,  0.5 );
	model.SetTransition( S12,  S13,   1.0 );
	model.SetTransition( S13,  S14,   1.0 );
	model.SetTransition( S14,  S15,   1.0 );
	model.SetTransition( S15,  S16,   1.0 );
	model.SetTransition( S16,  S17,   1.0 );
	model.SetTransition( S17,  S18,   1.0 );
	model.SetTransition( S18,  S19,   1.0 );
	model.SetTransition( S19,  S20,   1.0 );
	model.SetTransition( S20,  S21,   1.0 );
	model.SetTransition( S21,  S21,   1.0 );
	SetBackgroundEmissionProbs( model, (size_t)S0 );
	SetBackgroundEmissionProbs( model, (size_t)S1 );
	SetBackgroundEmissionProbs( model, (size_t)S2 );
	SetBackgroundEmissionProbs( model, (size_t)S3 );
	SetBackgroundEmissionProbs( model, (size_t)S4 );
	SetBackgroundEmissionProbs( model, (size_t)S5 );
	SetBackgroundEmissionProbs( model, (size_t)S6 );
	SetBackgroundEmissionProbs( model, (size_t)S7 );
	SetBackgroundEmissionProbs( model, (size_t)S8 );
	SetBackgroundEmissionProbs( model, (size_t)S9 );
	SetBackgroundEmissionProbs( model, (size_t)S10 );
	SetBackgroundEmissionProbs( model, (size_t)S11 );
	SetBackgroundEmissionProbs( model, (size_t)S12 );
	SetBackgroundEmissionProbs( model, (size_t)S13 );
	SetBackgroundEmissionProbs( model, (size_t)S14 );
	SetBackgroundEmissionProbs( model, (size_t)S15 );
	SetBackgroundEmissionProbs( model, (size_t)S16 );
	SetBackgroundEmissionProbs( model, (size_t)S17 );
	SetBackgroundEmissionProbs( model, (size_t)S18 );
	SetBackgroundEmissionProbs( model, (size_t)S19 );
	SetBackgroundEmissionProbs( model, (size_t)S20 );
	SetBackgroundEmissionProbs( model, (size_t)S21 );
	
	model.SetEmissionProbability( N1, 'D', 0.46 );
	model.SetEmissionProbability( N1, 'E', 0.54 );

	//model.SetEmissionProbability( P1, 'R', 0.38 );
	//model.SetEmissionProbability( P1, 'H', 0.17 );
	//model.SetEmissionProbability( P1, 'K', 0.45 );

	model.SetEmissionProbability( P1, 'S', 0.54 );
	model.SetEmissionProbability( P1, 'T', 0.46 );
}
#endif // MainProgram


template<typename Model>
void InitializeNullModel( Model& model, size_t len )
{
	// NOTE: Currently only supports sequences of length==len !
	typedef typename Model::StateId StateId;
	
	enum { I=0, T=1, S0=2 };

	for( size_t n = 0; n< len; ++n )
		model.AddState( S0+n,   "S"    );

	model.AddAlphabetSymbol( 'E' );
	model.AddAlphabetSymbol( 'D' );
	model.AddAlphabetSymbol( 'R' );
	model.AddAlphabetSymbol( 'H' );
	model.AddAlphabetSymbol( 'K' );
	model.AddAlphabetSymbol( 'I' );
	model.AddAlphabetSymbol( 'L' );
	model.AddAlphabetSymbol( 'V' );
	model.AddAlphabetSymbol( 'S' );
	model.AddAlphabetSymbol( 'T' );
	model.AddAlphabetSymbol( 'Q' );
	model.AddAlphabetSymbol( 'W' );
	model.AddAlphabetSymbol( 'G' );
	model.AddAlphabetSymbol( 'P' );
	model.AddAlphabetSymbol( 'F' );
	model.AddAlphabetSymbol( 'Y' );
	model.AddAlphabetSymbol( 'W' );
	model.AddAlphabetSymbol( 'C' );
	model.AddAlphabetSymbol( 'M' );
	model.AddAlphabetSymbol( 'A' );
	model.AddAlphabetSymbol( 'N' );
	model.SetDefaultProbabilities();
	model.SetTransition( I,   S0,   1.0 );

	for( size_t n = 0; n< len-1; ++n )
	{
		model.SetTransition( S0+n,  S0+n+1,   1.0 );
		SetBackgroundEmissionProbs( model, S0+n );
		//model.SetMu(  S0+n, 30.0);
		//model.SetEta( S0+n, 30.0);
	}
	
	model.SetTransition( S0+len-1,  T,   1.0 );
	//model.SetMu(  S0+len-1, 30.0);
	//model.SetEta( S0+len-1, 30.0);
	SetBackgroundEmissionProbs( model, S0+len-1 );
}

#ifdef MainProgram
int main(int argc, char* argv[])
{
	using boost::shared_ptr;

	if( argc < 2 )
	{
		std::cout<< "Usage: test_fasta <string>"<< std::endl;
		return -1;
	}

	shared_ptr<algo_t::model_type> model( new algo_t::model_type() );
	InitializeModel(*model);

	return test_model<algo_t>(argv[1], *model);
}
#endif // MainProgram


template<typename Algo>
int test_model(std::string fasta_path, typename Algo::model_type& _model)
{
	using namespace boost::assign;
	using boost::shared_ptr;


    //debug_print(*model);
	shared_ptr<typename Algo::model_type> model( &_model, null_deleter() );

	shared_ptr<typename Algo::model_type> null_model( new typename Algo::model_type );
	InitializeNullModel(*null_model);

    //debug_print(*null_model);

	//shared_ptr<algo_t::sequence_type> seq( new algo_t::sequence_type() );
  
	//(*seq) += "A", "K", "Q", "T", "H", "S", "A", "I", "Q", "E", "E", "E", "E", "R", "R", "L", "N", "L", "G", "V";
	//*seq = "AKQTHSAIQEEEERRLNLGV";

	fasta::seq_cont_t seqs;
	fasta::read_fasta(fasta_path, seqs);


	//const double Pmodel = forward->calc();
	//std::cout<< "P(O|model) = "<< Pmodel<< " ("<< (Pmodel - Pnull_model)/log(2.0) << " bits)";
	


#ifndef NO_TESTS
	{
		// Single-sequence tests

		shared_ptr<typename Algo::sequence_type> seq( new typename Algo::sequence_type( fasta::get_seq(seqs[0]) ) );

		typename Algo::forward_algo_type::env_type env = boost::make_tuple(seq, model);
		boost::shared_ptr<typename Algo::forward_algo_type> forward( new typename Algo::forward_algo_type( env ) );
		boost::shared_ptr<typename Algo::forward_begin_algo_type> forward_begin( new typename Algo::forward_begin_algo_type( env ) );
		bind_begin_end_algorithms( forward, forward_begin );
		boost::shared_ptr<typename Algo::backward_algo_type> backward( new typename Algo::backward_algo_type( env ) );
		boost::shared_ptr<typename Algo::backward_begin_algo_type> backward_begin( new typename Algo::backward_begin_algo_type( env ) );
		bind_begin_end_algorithms( backward, backward_begin );
		tests::test_28<Algo>( *forward, *backward_begin, *model, *seq );
		tests::test_D59_forw_back<Algo>( *forward, *backward_begin );

		std::cout<< "--- Forward -------" << std::endl;
		debug_print( *forward );
		std::cout<< "--- Backward -------" << std::endl;
		debug_print( *backward );
		
		typename Algo::baumwelch_algo_type bw(model, forward, forward_begin, backward, backward_begin, seq );
		tests::test_28<Algo>( bw, *model, *seq );
		tests::test_38<Algo>( bw, *model, *seq );
		tests::test_43b<Algo>( bw, *model );

		//debug_print( bw );

		std::cout<< "Finished single-sequence tests"<< std::endl;
	}
#endif 

/*
	{
		for( size_t i = 0; i<model->num_states(); ++i )
		{
			if( model->is_silent(i) ) continue;
			
			typename Algo::model_type::emissions_class_members_iterator mem_it, mem_end;
			boost::tie( mem_it, mem_end ) = model->GetEmissionsClassMembersRange(i);
			std::cout<< "In class "<< i<< ": ";

			const size_t a1 = *mem_it;
			++mem_it;
			const bool e1a = (mem_it == mem_end);
			const bool e1b = (mem_it != mem_end);

			const size_t a2 = *mem_it;
			++mem_it;
			const bool e2a = (mem_it == mem_end);
			const bool e2b = (mem_it != mem_end);

		
			for( ; mem_it != mem_end; ++mem_it )
			{
				const typename Algo::model_type::StateId mem = *mem_it;
				std::cout<< mem<< " ";
			}
			
			std::cout<< std::endl;
		}
	}
*/

	baumwelch_algo::MultipleSequenceBW<Algo> bw(model, seqs.begin(), seqs.end() );
	boost::shared_ptr<baumwelch_algo::MultipleSequenceBW<Algo> > pbw( &bw, null_deleter() );
	baumwelch_algo::EM<Algo, baumwelch_algo::MultipleSequenceBW<Algo> > em( *model, pbw );


	double p, p_prev = bw.Pseq();
	const size_t min_iters_1 = 3;
	const size_t min_iters_2 = 10;

	std::cout<< "Learning A (min. iters: "<< min_iters_1<< ")"<< std::endl;
	std::cout<< "iter 0...\tP(seqs|model) = "<< p_prev<< std::endl;

	for( size_t iter = 1; true; ++iter )
	{
		std::cout<< "iter "<< iter<< "a..."<< std::flush;
		
#ifndef NO_TESTS
		std::cout<< "."<< std::flush;
		bw.test_28();
		bw.test_38();
		tests::test_model_valid<Algo>( *model );
		bw.test_D59_forw_back();
#endif
		p = em.reestimate_a();
		//assert( (iter == 0) || (p >= p_prev) );

		std::cout<< "\tP(seqs|model) = "<< p;
		std::cout<< "\t\t("<< (p - p_prev) << ")";
		if( p < p_prev ) std::cout<< "\tWarning: re-estimation impaired model!";
		std::cout<< std::endl;

		if( iter%20==15 )     debug_print(*model);
		
		if( (iter > min_iters_1) && (p - p_prev < 0.1) ) break;

		p_prev = p;
	}


	std::cout<< "Learning E (min. iters: "<< min_iters_2<< ")"<< std::endl;
	// Model set - now train emissions

	for( size_t iter = 1; true; ++iter )
	{
		std::cout<< "iter "<< iter<< "e..."<< std::flush;
		
#ifndef NO_TESTS
		std::cout<< "."<< std::flush;
		bw.test_28();
		bw.test_38();
		bw.test_40c();
		tests::test_model_valid<Algo>( *model );
		bw.test_D59_forw_back();
#endif
		p = em.reestimate_b();
		//assert( (iter == 0) || (p >= p_prev) );

		std::cout<< "\tP(seqs|model) = "<< p;
		std::cout<< "\t\t("<< (p - p_prev) << ")";
		if( p < p_prev ) std::cout<< "\tWarning: re-estimation impaired model!";
		std::cout<< std::endl;

		if( iter%10==5 )     debug_print(*model);
		
		if( (iter > min_iters_2) && (p - p_prev < 0.1) ) break;
		
		p_prev = p;
	}


	fasta::seq_cont_t::iterator it, it_end;
	it = seqs.begin();
	it_end = seqs.end();
	for( ; it != it_end; ++it )
	{
		shared_ptr<typename Algo::sequence_type> seq( new typename Algo::sequence_type( fasta::get_seq(*it) ) );

		// Calculate probability using null model
		typename Algo::forward_algo_type::env_type null_env = boost::make_tuple(seq, null_model);
		typename Algo::forward_algo_type forward_null( null_env );
		typename Algo::forward_begin_algo_type forward_begin_null( null_env );
		bind_begin_end_algorithms( forward_null, forward_begin_null );
		
		const double Pnull_model = forward_null.calc();
		
		typename Algo::forward_algo_type::env_type env = boost::make_tuple(seq, model);
		typename Algo::forward_algo_type forward( env );
		typename Algo::forward_begin_algo_type forward_begin( env );
		bind_begin_end_algorithms( forward, forward_begin );

		const double Pmodel = forward.calc();

		std::cout<< "P(O|model') = "<< Pmodel<< " ("<< (Pmodel - Pnull_model)/log(2.0) << " bits)"<< std::endl;;
	
		typename Algo::viterbi_algo_type viterbi( seq, model);
		viterbi.calc();
		
		std::vector<std::string> bt;
		viterbi.get_backtrace(bt);
		print_seqs( *seq, bt );
	}
	
/*
	const double Pmodel2 = forward->calc();
	std::cout<< "\tP(O|model') = "<< Pmodel2<< " ("<< (Pmodel2 - Pnull_model)/log(2.0) << " bits)"<< std::endl;
	
	algo_t::viterbi_algo_type viterbi(seq, model);
	viterbi.calc();
	
	std::vector<std::string> bt;
	viterbi.get_backtrace(bt);
    print_seqs( *seq, bt );
*/
    return 0;  
}

// Explicit instantiation
template
int test_model<algo_t>(std::string fasta_path, algo_t::model_type& _model);
template
int test_model<algo_hmm_t>(std::string fasta_path, algo_hmm_t::model_type& _model);

template<typename Algo>
typename boost::enable_if<typename Algo::model_type::is_explicit_duration_model,void>::type
print_duration(const typename Algo::model_type& model)
{
	std::cout<< "Mean: ";
	double mean_len = 0.0;
	for( size_t i=0; i< num_states(model); ++i )
	{
		if( model.IsReservedState(i)) continue;
		const double mean = model.GetMu(i) / model.GetEta(i);
		std::cout<< i<< ": "<< mean<< " ";
		mean_len += mean;
	}
	
	std::cout<< "(total="<< mean_len<< ")"<< std::endl;
}

template<typename Algo>
typename boost::disable_if<typename Algo::model_type::is_explicit_duration_model,void>::type
print_duration(const typename Algo::model_type& model) {}



template<typename Algo>
double train_model(std::string fasta_path, typename Algo::model_type& _model, const double converge_threshold)
{
	using namespace boost::assign;
	using boost::shared_ptr;

    //debug_print(*model);
	shared_ptr<typename Algo::model_type> model( &_model, null_deleter() );

	shared_ptr<typename Algo::model_type> null_model( new typename Algo::model_type );
	InitializeNullModel(*null_model);


	fasta::seq_cont_t seqs;
	fasta::read_fasta(fasta_path, seqs);


#ifndef NO_TESTS
	{
		// Model tests
		tests::test_model_valid<Algo>( *model );

		// Single-sequence tests
		shared_ptr<typename Algo::sequence_type> seq( new typename Algo::sequence_type( fasta::get_seq(seqs[0]) ) );

		typename Algo::forward_algo_type::env_type env = boost::make_tuple(seq, model);
		boost::shared_ptr<typename Algo::forward_algo_type> forward( new typename Algo::forward_algo_type(env) );
		boost::shared_ptr<typename Algo::forward_begin_algo_type> forward_begin( new typename Algo::forward_begin_algo_type(env) );
		bind_begin_end_algorithms( forward, forward_begin );
		boost::shared_ptr<typename Algo::backward_algo_type> backward( new typename Algo::backward_algo_type(env) );
		boost::shared_ptr<typename Algo::backward_begin_algo_type> backward_begin( new typename Algo::backward_begin_algo_type(env) );
		bind_begin_end_algorithms( backward, backward_begin );
//		tests::test_28<Algo>( *forward, *backward_begin, *model, *seq ); // TODO - update test to support HSMMs
		tests::test_D59_forw_back<Algo>( *forward, *backward_begin );

		//std::cout<< "--- Forward -------" << std::endl;
		//debug_print( *forward );

		forward->calc_all();
		forward_begin->calc_all();
		backward->calc_all();
		backward_begin->calc_all();

		std::cout<< "Pforw_end   (O|model) = "<< forward->calc() <<std::endl;
		//std::cout<< "Pforw_begin (O|model) = "<< forward_begin->calc() <<std::endl;
		//std::cout<< "Pback_end   (O|model) = "<< backward->calc() <<std::endl;
		std::cout<< "Pback_begin (O|model) = "<< backward_begin->calc() <<std::endl;

		//std::cout<< "--- Forward End   -------" << std::endl;
		//debug_print( *forward );
		
		//std::cout<< "--- Forward Begin -------" << std::endl;
		//debug_print( *forward_begin );

		//std::cout<< "--- Backward End   -------" << std::endl;
		//debug_print( *backward );
		
		//std::cout<< "--- Backward Begin -------" << std::endl;
		//debug_print( *backward_begin );

		typename Algo::baumwelch_algo_type bw(model, forward, forward_begin, backward, backward_begin, seq );
		//tests::test_28<Algo>( bw, *model, *seq );
		//tests::test_38<Algo>( bw, *model, *seq );
		//tests::test_43b<Algo>( bw, *model );

		//debug_print( bw );

		std::cout<< "Finished single-sequence tests"<< std::endl;
	}
#endif 

	baumwelch_algo::MultipleSequenceBW<Algo> bw(model, seqs.begin(), seqs.end() );
	boost::shared_ptr<baumwelch_algo::MultipleSequenceBW<Algo> > pbw( &bw, null_deleter() );
	baumwelch_algo::EM<Algo, baumwelch_algo::MultipleSequenceBW<Algo> > em( *model, pbw );


	double p, p_prev = bw.Pseq();
	const size_t min_iters = 3;
	size_t stable_count = 0;

	// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //
	debug_print(*model,4);
	if( Algo::model_type::is_explicit_duration_model::value )
		print_duration<Algo>(*model);
	// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //

	std::cout<< "Learning (min. iters: "<< min_iters<< ")"<< std::endl;
	std::cout<< "iter 0...\tP(seqs|model) = "<< p_prev<< std::endl;
	std::ofstream trainlog("train.log", std::ios_base::out | std::ios_base::app );
	typedef enum 
	{
		estimate_e, estimate_shape, estimate_scale
	} estimate_mode;
	estimate_mode mode;

	for( size_t iter=1; true; ++iter)
	{
		if( iter<5 ) mode = estimate_e;
		else
		{
			switch((iter-5)%4)
			{
			case 0:
				mode = estimate_e;
				break;
			case 1:
			case 2:
				mode = estimate_scale;
				break;
			case 3:
				mode = estimate_shape;
				break;
			}
		}

		switch(mode)
		{
		case estimate_e:
			std::cout<< "iter "<< iter<< "e..."<< std::flush;
			break;
		case estimate_shape:
			std::cout<< "iter "<< iter<< "sh..."<< std::flush;
			break;
		case estimate_scale:
			std::cout<< "iter "<< iter<< "sc..."<< std::flush;
			break;
		default:
			assert(false);
			
		}

		
#ifndef NO_TESTS
		std::cout<< "."<< std::flush;
		tests::test_model_valid<Algo>( *model );
		//bw.test_28();
		//bw.test_38();
		//bw.test_40c();
		//bw.test_D59_forw_back();
#endif

		switch(mode)
		{
		case estimate_e:
			p = em.reestimate_b();
			break;
		case estimate_shape:
			p = em.reestimate_shape();
			break;
		case estimate_scale:
			p = em.reestimate_scale();
			break;
		default:
			assert(false);
			
		}

		if( mode == estimate_shape || mode == estimate_scale)
		{
			// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //
			if( Algo::model_type::is_explicit_duration_model::value )
				print_duration<Algo>(*model);
			// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //// DEBUG //
		}
		

		std::cout<< "\tP(seqs|model) = "<< p;
		std::cout<< "\t\t("<< (p - p_prev) << ")";
		if( p < p_prev ) std::cout<< "\tWarning: re-estimation impaired model!";
		std::cout<< std::endl;

		trainlog<< iter<< " "<< p<< std::endl;

		if( iter%10==5 )     debug_print(*model,0);
		else debug_print(*model,4);

		if(p - p_prev < converge_threshold)
			++stable_count;
		else
			stable_count = 0;

		// Break if we haven't improved in a while
		if( (iter > min_iters) && (stable_count > /*1*/3) ) break;
		// Break if we've completely divereged and P=0 (to prevent an infinite loop)
		if( p == -std::numeric_limits<typename Algo::probability_type>::max() ) break;

		p_prev = p;
		
		//em.debug_print();
	}

#ifndef NO_TESTS
		tests::test_model_valid<Algo>( *model );
	//	bw.test_28();
	//	bw.test_38();
	//	bw.test_40c();
	//	bw.test_D59_forw_back();
#endif

    return p;  
}


// Explicit instantiation
template
double train_model<algo_t>(std::string fasta_path, algo_t::model_type& _model, const double converge_threshold);
template
double train_model<algo_hmm_t>(std::string fasta_path, algo_hmm_t::model_type& _model, const double converge_threshold);


template<typename Algo>
double score_sequences(std::string fasta_path, typename Algo::model_type& _model, bool fasta_format, int filter_type, double filter_cutoff)
{
	using namespace boost::assign;
	using boost::shared_ptr;


    //debug_print(*model);
	shared_ptr<typename Algo::model_type> model( &_model, null_deleter() );


#ifndef NO_TESTS
	{
		tests::test_model_valid<Algo>( *model );
	}
#endif 


	fasta::seq_cont_t seqs;
	fasta::read_fasta(fasta_path, seqs);


	fasta::seq_cont_t::iterator it, it_end;
	it = seqs.begin();
	it_end = seqs.end();
	for( ; it != it_end; ++it )
	{
		shared_ptr<typename Algo::sequence_type> seq( new typename Algo::sequence_type( fasta::get_seq(*it) ) );

		shared_ptr<typename Algo::model_type> null_model( new typename Algo::model_type );
		InitializeNullModel(*null_model, length(*seq));

#ifndef NO_TESTS
	{
		tests::test_model_valid<Algo>( *null_model );
	}
#endif 

		// Calculate probability using null model
		typename Algo::forward_algo_type::env_type null_env = boost::make_tuple(seq, null_model);
		typename Algo::forward_algo_type forward_null( null_env );
		typename Algo::forward_begin_algo_type forward_begin_null( null_env );
		bind_begin_end_algorithms( forward_null, forward_begin_null );
		typename Algo::backward_algo_type backward_null( null_env );
		typename Algo::backward_begin_algo_type backward_begin_null( null_env );
		bind_begin_end_algorithms( backward_null, backward_begin_null );
		//null_model->debug_print();
		
		typename Algo::forward_algo_type::env_type env = boost::make_tuple(seq, model);
		typename Algo::forward_algo_type forward( env );
		typename Algo::forward_begin_algo_type forward_begin( env );
		bind_begin_end_algorithms( forward, forward_begin );
		typename Algo::backward_algo_type backward( env );
		typename Algo::backward_begin_algo_type backward_begin( env );
		bind_begin_end_algorithms( backward, backward_begin );


#ifndef NO_TESTS
	{
		tests::test_D59_forw_back<Algo>( forward, backward_begin );
		tests::test_D59_forw_back<Algo>( forward_null, backward_begin_null );
	}
#endif 
	
		const double Pnull_model = forward_null.calc();
		const double Pmodel = forward.calc();


		const double score = (Pmodel - Pnull_model)/log(2.0);
		
		if(     (filter_type == 0)
			|| ((filter_type == 1) && ( score >  filter_cutoff ))
			|| ((filter_type == 2) && ( score <= filter_cutoff )) )
		{
			if( fasta_format ) std::cout<< ">";
		
			std::cout<< fasta::get_acc(*it)<< "\t"<< Pnull_model<< "\t"<< Pmodel<< "\t"<< score << std::endl;

			if( fasta_format )  std::cout<< fasta::get_seq(*it)<< std::endl; 
	
			//typename Algo::viterbi_algo_type viterbi( seq, model);
			//viterbi.calc();
			
			//std::vector<std::string> bt;
			//viterbi.get_backtrace(bt);
			//print_seqs( *seq, bt );
		}
	}

    return 0.0;
}

// Explicit instantiation
template
double score_sequences<algo_t>(std::string fasta_path, algo_t::model_type& _model, bool fasta_format, int filter_type, double filter_cutoff);
template
double score_sequences<algo_hmm_t>(std::string fasta_path, algo_hmm_t::model_type& _model, bool fasta_format, int filter_type, double filter_cutoff);


template<typename Algo>
double score_sequences2(std::string fasta_path, typename Algo::model_type& _model, const boost::python::list& indexes)
{
	using boost::shared_ptr;


    //debug_print(*model);
	shared_ptr<typename Algo::model_type> model( &_model, null_deleter() );

#ifndef NO_TESTS
	{
		tests::test_model_valid<Algo>( *model );
	}
#endif 

	std::vector<size_t> vec;
	python_to_stl_vector<size_t, const boost::python::list>(indexes, vec);

	fasta::seq_cont_t seqs;
	fasta::read_fasta(fasta_path, seqs, vec.begin(), vec.end());

	double best_score = -1e200;

	std::cout<< "#P(O|Model)\tP(O|Null_model)\tLogOdds"<< std::endl;

	fasta::seq_cont_t::iterator it, it_end;
	it = seqs.begin();
	it_end = seqs.end();
	for( ; it != it_end; ++it )
	{
		shared_ptr<typename Algo::sequence_type> seq( new typename Algo::sequence_type( fasta::get_seq(*it) ) );

		shared_ptr<algo_hmm_t::model_type> null_model( new algo_hmm_t::model_type );
		InitializeNullModel(*null_model, length(*seq));

#ifndef NO_TESTS
		{
			tests::test_model_valid<algo_hmm_t>( *null_model );
		}
#endif 

		// Calculate probability using null model
		algo_hmm_t::forward_algo_type::env_type null_env = boost::make_tuple(seq, null_model);
		algo_hmm_t::forward_algo_type forward_null( null_env );
		algo_hmm_t::forward_begin_algo_type forward_begin_null( null_env );
		bind_begin_end_algorithms( forward_null, forward_begin_null );
		algo_hmm_t::backward_algo_type backward_null( null_env );
		algo_hmm_t::backward_begin_algo_type backward_begin_null( null_env );
		bind_begin_end_algorithms( backward_null, backward_begin_null );
		//null_model->debug_print();
		
		typename Algo::forward_algo_type::env_type env = boost::make_tuple(seq, model);
		typename Algo::forward_algo_type forward( env );
		typename Algo::forward_begin_algo_type forward_begin( env );
		bind_begin_end_algorithms( forward, forward_begin );
		typename Algo::backward_algo_type backward( env );
		typename Algo::backward_begin_algo_type backward_begin( env );
		bind_begin_end_algorithms( backward, backward_begin );


#ifndef NO_TESTS
		{
			tests::test_D59_forw_back<Algo>( forward, backward_begin );
			tests::test_D59_forw_back<algo_hmm_t>( forward_null, backward_begin_null );
		}
#endif 
	
		const double Pnull_model = forward_null.calc();
		const double Pmodel = forward.calc();


		const double score = (Pmodel - Pnull_model)/log(2.0);

		std::cout<< Pmodel<< "\t"<< Pnull_model<< "\t"<< score<< "\t"<< fasta::get_seq(*it)<< "\t"<< fasta::get_acc(*it)<< "\t"<< fasta::get_desc(*it)<< std::endl;
		//std::cout<< Pmodel<< "\t"<< Pnull_model<< "\t"<< score<< std::endl;
		

		if( score > best_score ) best_score = score;
	}

    return best_score;
}

// Explicit instantiation
template
double score_sequences2<algo_t>(std::string fasta_path, algo_t::model_type& _model, const boost::python::list& indexes);
template
double score_sequences2<algo_hmm_t>(std::string fasta_path, algo_hmm_t::model_type& _model, const boost::python::list& indexes);


/*
size_t test_vector(boost::python::list& vec)
{
	vec.append(5);
	return boost::python::len(vec);
}
*/

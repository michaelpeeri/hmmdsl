#pragma once
#include <string>
//#include <boost/python.hpp>
#include "common.hpp"
#include "model.hpp"
#include "viterbi.hpp"
#include "forward.hpp"
#include "backward.hpp"
#include "baumwelch.hpp"

// Definitions for the HMM algorithm
typedef	HiddenMarkovAlgorithm<
	std::string,
	double_logscale,
	MatrixModel,
	ViterbiAlgorithm,
	forward_algo::ForwardAlgorithm,
	forward_algo::ForwardAlgorithmHMM,
	forward_algo::ForwardAlgorithm,
	forward_algo::ForwardAlgorithmHMM,
	backward_algo::BackwardAlgorithm,
	backward_algo::BackwardAlgorithmHMM,
	backward_algo::BackwardAlgorithm,
	backward_algo::BackwardAlgorithmHMM,
	baumwelch_algo::BaumWelchAlgorithm>
algo_hmm_t;

// Definitions for the HMM algorithm
typedef	HiddenMarkovAlgorithm<
	std::string,
	double_logscale,
	HSMMMatrixModel,
	ViterbiAlgorithm,
	forward_algo::ForwardAlgorithm,
	forward_algo::ForwardAlgorithmHSMMEnd,
	forward_algo::ForwardAlgorithm,
	forward_algo::ForwardAlgorithmHSMMBegin,
	backward_algo::BackwardAlgorithm,
	backward_algo::BackwardAlgorithmHSMMEnd,
	backward_algo::BackwardAlgorithm,
	backward_algo::BackwardAlgorithmHSMMBegin,
	baumwelch_algo::HSMMBaumWelchAlgorithm,
	true,
	60>
algo_hsmm_t;

// Choose which algorithm to use:
typedef algo_hsmm_t algo_t;


///////////////////////////////////////////////////////////////////////////////////////
// Global functions to be exported through the python module:

template<typename algo_t>
int test_model(std::string fasta_path, typename algo_t::model_type& _model);


template<typename Algo>
double train_model(std::string fasta_path, typename Algo::model_type& _model, const double converge_threshold);

template<typename Algo>
double score_sequences(std::string fasta_path, typename Algo::model_type& _model, bool fasta_format, int filter_type, double filter_cutoff);

template<typename Algo>
double score_sequences2(std::string fasta_path, typename Algo::model_type& _model, const boost::python::list& indexes);


template<typename Model>
void InitializeNullModel( Model& model, size_t len = 20 );

//size_t test_vector(boost::python::list& vec);


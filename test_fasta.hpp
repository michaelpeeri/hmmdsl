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
    false, // Use emissions class?
    100>   // D (maximum state duration)
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

//template<typename Algo>
//double show_durations(typename Algo::model_type& _model );

template<typename Algo>
std::pair<std::string, double> find_stereotype(typename Algo::model_type& _model, size_t length );

template<typename Algo>
void print_stereotype_stats(typename Algo::model_type& _model, size_t length );

template<typename Model>
void InitializeNullModel( Model& model, size_t len = 20 );

//size_t test_vector(boost::python::list& vec);


template<typename Algo>
void align_viterbi(std::string fasta_path, typename Algo::model_type& _model, const boost::python::list& indexes, int format=0);

template<typename Algo>
void print_model_statistics(typename Algo::model_type& _model);

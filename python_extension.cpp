#include <boost/python.hpp>
#include "test_fasta.hpp"
#include "common.hpp"
#include "model.hpp"
#include "viterbi.hpp"
#include "forward.hpp"
#include "backward.hpp"
#include "baumwelch.hpp"
#include "emit.hpp"
#include "read_matrix.hpp"
#include "perturb_model.hpp"

using namespace boost::python;

char const* greet(size_t x)
{
	static char const* msgs[] = { "Hello", "Boost.Python", "World!" };
	if( x > 2)
		throw std::range_error("greet: index out of range");
	
	return msgs[x];
}

template<typename Algo>
void translator_SymbolNotInAlphabet(const typename MatrixModel<Algo>::SymbolNotInAlphabet& err)
{
	PyErr_SetString(PyExc_LookupError, "HMM MatrixModel: Symbol not in alphabet");
}

BOOST_PYTHON_MODULE(hmmdsl_py)
{
	def("greet", greet, "return one of 3 greetings" );

	typedef baumwelch_algo::EM<algo_t,baumwelch_algo::MultipleSequenceBW<algo_t> > em_t;
	typedef baumwelch_algo::EM<algo_hmm_t,baumwelch_algo::MultipleSequenceBW<algo_hmm_t> > em_hmm_t;
	
	algo_t::probability_type (algo_t::model_type::*pfxe)(const algo_t::model_type::StateId& state, const algo_t::symbol_type& sym) const = &algo_t::model_type::e;
	algo_hmm_t::probability_type (algo_hmm_t::model_type::*pfxe_hmm)(const algo_hmm_t::model_type::StateId& state, const algo_hmm_t::symbol_type& sym) const = &algo_hmm_t::model_type::e;
	void (algo_t::model_type::*pfx_debug_print)() const = &algo_t::model_type::debug_print;
	void (algo_hmm_t::model_type::*pfx_debug_print_hmm)() const = &algo_hmm_t::model_type::debug_print;

	class_<algo_t::model_type>("Model")
		.def("AddState", &algo_t::model_type::AddState)
		.def("SetTransition", &algo_t::model_type::SetTransition)
		.def("SetTransitionLogspace", &algo_t::model_type::SetTransitionLogspace)
		.def("AddAlphabetSymbol", &algo_t::model_type::AddAlphabetSymbol)
		.def("SetEmissionProbability", &algo_t::model_type::SetEmissionProbability)
		.def("SetEmissionProbabilityLogspace", &algo_t::model_type::SetEmissionProbabilityLogspace)
		.def("GetStateName", &algo_t::model_type::GetStateName)
		.def("num_symbols", &algo_t::model_type::num_symbols)
		.def("num_states", &algo_t::model_type::num_states)
		.def("get_symbol", &algo_t::model_type::get_symbol)
		.def("SetDefaultProbabilities", &algo_t::model_type::SetDefaultProbabilities)
		.def("a", &algo_t::model_type::a)
		.def("e", pfxe)
		.def("debug_print", pfx_debug_print )
		.def("is_silent", &algo_t::model_type::is_silent)
		.def("GetInitialState", &algo_t::model_type::GetInitialState)
		.def("GetTerminalState", &algo_t::model_type::GetTerminalState)
		.def("SetFixedEmissionsState", &algo_t::model_type::SetFixedEmissionsState)
		.def("MergeEmissionsClass", &algo_t::model_type::MergeEmissionsClass)
		.def("GetEta", &algo_t::model_type::GetEta)
		.def("GetMu", &algo_t::model_type::GetMu)
		.def("SetEta", &algo_t::model_type::SetEta)
		.def("SetMu", &algo_t::model_type::SetMu)
		.def("emissions_entropy", &algo_t::model_type::emissions_entropy)
		.def("split_state", &algo_t::model_type::split_state)
		;

	class_<algo_hmm_t::model_type>("HMMModel")
		.def("AddState", &algo_hmm_t::model_type::AddState)
		.def("SetTransition", &algo_hmm_t::model_type::SetTransition)
		.def("SetTransitionLogspace", &algo_hmm_t::model_type::SetTransitionLogspace)
		.def("AddAlphabetSymbol", &algo_hmm_t::model_type::AddAlphabetSymbol)
		.def("SetEmissionProbability", &algo_hmm_t::model_type::SetEmissionProbability)
		.def("SetEmissionProbabilityLogspace", &algo_hmm_t::model_type::SetEmissionProbabilityLogspace)
		.def("GetStateName", &algo_hmm_t::model_type::GetStateName)
		.def("num_symbols", &algo_hmm_t::model_type::num_symbols)
		.def("num_states", &algo_hmm_t::model_type::num_states)
		.def("get_symbol", &algo_hmm_t::model_type::get_symbol)
		.def("SetDefaultProbabilities", &algo_hmm_t::model_type::SetDefaultProbabilities)
		.def("a", &algo_hmm_t::model_type::a)
		.def("e", pfxe_hmm)
		.def("debug_print", pfx_debug_print_hmm )
		.def("is_silent", &algo_hmm_t::model_type::is_silent)
		.def("GetInitialState", &algo_hmm_t::model_type::GetInitialState)
		.def("GetTerminalState", &algo_hmm_t::model_type::GetTerminalState)
		.def("SetFixedEmissionsState", &algo_hmm_t::model_type::SetFixedEmissionsState)
		.def("MergeEmissionsClass", &algo_hmm_t::model_type::MergeEmissionsClass)
		.def("emissions_entropy", &algo_hmm_t::model_type::emissions_entropy)
		;

	class_<em_t>("EM", init<algo_t::model_type&, const std::string>() )
		.def(init<algo_t::model_type&, const std::string, const boost::python::list&>() )
		.def("reestimate_a", &em_t::reestimate_a)
		.def("reestimate_b", &em_t::reestimate_b)
		.def("reestimate_scale", &em_t::reestimate_scale)
		.def("reestimate_shape", &em_t::reestimate_shape)
		.def("debug_print", &em_t::debug_print)
		;

	class_<em_hmm_t>("HMMEM", init<algo_hmm_t::model_type&, const std::string>() )
		.def(init<algo_hmm_t::model_type&, const std::string, const boost::python::list&>() )
		.def("reestimate_a", &em_hmm_t::reestimate_a)
		.def("reestimate_b", &em_hmm_t::reestimate_b)
		.def("reestimate_scale", &em_hmm_t::reestimate_scale)
		.def("reestimate_shape", &em_hmm_t::reestimate_shape)
		.def("debug_print", &em_hmm_t::debug_print)
		;

	class_<random_emitter<algo_t> >("Emitter", init<const algo_t::model_type&>() )
		.def("emit", &random_emitter<algo_t>::emit)
		.def("emit_fixed_length", &random_emitter<algo_t>::emit_fixed_length)
		.def("test_random_transition", &random_emitter<algo_t>::test_random_transition)
		;

	class_<boost::multi_array<algo_t::probability_type,2> >("Matrix")
		;
	
	class_<perturb_model<algo_t> >("Perturb", init<algo_t::model_type&>() )
		.def("perturb_emissions", &perturb_model<algo_t>::perturb_emissions)
		;

	class_<perturb_model<algo_hmm_t> >("Perturb_hmm", init<algo_hmm_t::model_type&>() )
		.def("perturb_emissions", &perturb_model<algo_hmm_t>::perturb_emissions)
		;

	def("test_model", test_model<algo_t> );
	def("train_model", train_model<algo_t> );
	def("score_sequences", score_sequences<algo_t> );
	def("score_sequences2", score_sequences2<algo_t> );
	def("InitializeNullModel", InitializeNullModel<algo_t::model_type> );
//	def("test_vector", test_vector ); // Test for output to python arrays
	def("test_read_matrix_sum", read_matrix::test_read_matrix_sum<double> );
	def("read_matrix", read_matrix::read_matrix<double> );
	def("count_fasta", fasta::count_fasta );
	def("relax_emissions", relax_emissions<algo_t>);
	register_exception_translator<MatrixModel<algo_t>::SymbolNotInAlphabet>(&translator_SymbolNotInAlphabet<algo_t>);

	def("test_model_hmm", test_model<algo_hmm_t> );
	def("train_model_hmm", train_model<algo_hmm_t> );
	def("score_sequences_hmm", score_sequences<algo_hmm_t> );
	def("score_sequences2_hmm", score_sequences2<algo_hmm_t> );
	def("InitializeNullModel_hmm", InitializeNullModel<algo_hmm_t::model_type> );
//	def("test_vector", test_vector ); // Test for output to python arrays
	def("relax_emissions_hmm", relax_emissions<algo_hmm_t>);
	register_exception_translator<MatrixModel<algo_hmm_t>::SymbolNotInAlphabet>(&translator_SymbolNotInAlphabet<algo_hmm_t>);

}


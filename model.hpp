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
#include <math.h>
#include <vector>
#include <string>
#include <set>
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#undef BOOST_DISABLE_ASSERTS
#include <boost/bimap.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/mpl/logical.hpp>
#include <iostream>
#include <cassert>
#include "common.hpp"
#include "tests.hpp"
#include "v2.hpp"

template<typename Model>
typename Model::StateId begin_state(const Model& model) { return 0; }

template<typename Model>
typename Model::StateId end_state(const Model& model) { return 1; }


template<typename Algo>
struct IModelBuilder
/* Interface for initializing model at run-time */
{
public:
    //typedef typename Algo::model_type::state_id_type StateId;
    typedef size_t StateId;
    typedef typename Algo::symbol_type Symbol;
    typedef typename Algo::probability_type Probability;
    
    virtual void AddState(const StateId& state, const std::string& name) = 0;
    virtual const std::string GetStateName(const StateId& state) const = 0;
    virtual void SetTransition(const StateId& from, const StateId& to, Probability transition) = 0;
    virtual void AddAlphabetSymbol(  const Symbol& sym ) = 0;
    virtual void SetEmissionProbability( const StateId& state, const Symbol& sym, Probability emit_prob ) = 0;

    virtual typename Algo::sequence_type GetAlphabet() const = 0;
    
    /*
      virtual void SetStartState(const StateId& state) = 0;
      virtual void SetEndState(const StateId& state) = 0;
    */
    
    virtual void resize() = 0;	
};

namespace detail
{
	template<typename Iter, typename Value>
	class EmissionsClassMembersIterator
	{
		// Model of: IncrementableIterator<StateId>
	public:
		typedef Value value_type;

		// Required for: DefaultConstructible
		EmissionsClassMembersIterator()
			: _class(0) // Mark as uninitialized
			{}

		EmissionsClassMembersIterator( Iter pos, Iter end, const size_t myclass )
			: _pos(pos)
			, _end(end)
			, _class(myclass)
			{ _find_first(); }

		// Required for: CopyConstructible
		EmissionsClassMembersIterator( const EmissionsClassMembersIterator& other )
			: _pos(other._pos)
			, _end(other._end)
			, _class(other._class)
			{}
		
		
		// Required for: Assignable
		EmissionsClassMembersIterator& operator=( const EmissionsClassMembersIterator& other )
		{
			_pos = other._pos;
			_end = other._end;
			_class = other._class;
			return *this;
		}
		bool operator==( const EmissionsClassMembersIterator& other ) const 
		{
			return( (_pos == other._pos) && (_class == other._class) );
		}
		bool operator!=( const EmissionsClassMembersIterator& other ) const 
		{
			return( (_pos != other._pos) || (_class != other._class) );
		}

		// Required for: ReadableIterator
		Value operator*() const 
		{
			assert( _class != 0 ); // Verify an uninitialized iterator is not dereferenced
			return *_pos;
		}

		// Required for: IncrementableIterator
		EmissionsClassMembersIterator& operator++()
		{
			assert( _pos != _end );
			_advance();
			
			return *this;
		}

	protected:
		void _find_first()
		{
			if( ( _pos != _end ) && ( *_pos != _class ))
				_advance();
		}
		
	protected:
		void _advance()
		{
			while( true )
			{
				++_pos;
				if( (_pos == _end) || ( *_pos  == _class) )
					break;
			}
		}
		
		
		
	protected:
		Iter _pos;
		Iter _end;
		size_t _class;
	};

	template<typename Value>
	class ClassesSet : public std::set<Value>
	{
	public:
		template<typename ExtIter>
		ClassesSet( ExtIter begin, ExtIter end )
		{
			// Populate the set
			for( ; begin != end; ++begin )
				insert( *begin );
		}
	public:
		ClassesSet() {}
	};
		
	template<typename Value>
	class EmissionsClassesIterator
	{
		// Model of: IncrementableIterator<StateId>
	public:
		typedef Value value_type;

		template<typename ExtIter>
		EmissionsClassesIterator( const ClassesSet<Value>& classes )
			: _classes(classes)
		{
			_pos = _classes.begin();
		}

		// Required for: CopyConstructible
		EmissionsClassesIterator( const EmissionsClassesIterator& other )
                : _classes(other._classes)
                , _pos(other._pos)
                {}
		
		// Required for: Assignable
		EmissionsClassesIterator& operator=( const EmissionsClassesIterator& other )
		{
			_classes = other._classes;
			_pos = other._pos;
			return *this;
		}
		bool operator==( const EmissionsClassesIterator& other ) const 
		{
			return( (_pos == other._pos) && (&_classes == &other._classes) );
		}

		// Required for: ReadableIterator
		Value operator*() const 
		{
			return *_pos;
		}

		// Required for: IncrementableIterator
		EmissionsClassesIterator& operator++()
		{
			assert( _pos != _classes.end() );
			++_pos;
			return *this;
		}
		
		
		
	protected:
		const ClassesSet<Value>& _classes;
		typename ClassesSet<Value>::const_iterator _pos;
	};
	
} // namespace detail	


struct HMMMatrixModel {}; // Base type (used for template specialization)

template<typename Algo>
class MatrixModel : public IModelBuilder<Algo>, public HMMMatrixModel
{  
public:
	typedef size_t StateId;
	typedef size_t state_id_type;
	typedef typename Algo::symbol_type Symbol;
	typedef typename Algo::probability_type Probability;
	typedef MatrixModel<Algo> self_type;
  
	struct SymbolNotInAlphabet{}; // exception
	typedef boost::mpl::bool_<false> is_explicit_duration_model;
	
protected:
	typedef IModelBuilder<Algo> Base;
	typedef MatrixModel<Algo> Model;
 	typedef std::vector<size_t> emission_classes_t;

public:
	MatrixModel()
	{
		_init_fixed_parts();
	}

public:
  virtual ~MatrixModel() {}

private:
	void _init_fixed_parts()
	{
		// Add the initial state
		_states.push_back("I");
		_emission_class.push_back(0);
		
		// Add the terminal state
		_states.push_back("T");
		_emission_class.push_back(0);

		resize();
			
	}

public:	
	void AddState(const StateId& state, const std::string& name)
	{
		assert( !IsReservedState(state) );
		
		_states.push_back(name);
		_emission_class.push_back( state ); // By default, each state belongs to a class of its own
	   
		resize();
		_update_emission_classes_for_iter_impl();
	}
public:	
	void SetTransition(const StateId& from, const StateId& to, Probability transition)
	{
		assert( isvalidprobability( transition ) );
		SetTransitionLogspace( from, to, log(transition) );
	}
public:	
	void SetTransitionLogspace(const StateId& from, const StateId& to, Probability transition)
	{
		assert( transition <= 0.0);
		_a[from][to] = transition;
	}
public:	
	void AddAlphabetSymbol(  const Symbol& sym )
	{
		const size_t nextid = _symbols.size();
		_symbols.insert( typename symbols_map_t::value_type( sym, nextid ) );
		resize();
	}
public:	
	void SetEmissionProbability( const StateId& state, const Symbol& sym, Probability emit_prob )
	{
		assert( !IsReservedState(state) );
		assert( isvalidprobability( emit_prob ) );
		SetEmissionProbabilityLogspace( state, sym, log(emit_prob) );
	}
	
public:	
	void SetEmissionProbabilityLogspace( const StateId& state, const Symbol& sym, Probability emit_prob )
	{
		assert( emit_prob <= 0.0 );
		typename symbols_map_t::left_map::const_iterator it;
		it = _symbols.left.find( sym );
		assert( it != _symbols.left.end() );
		if( it == _symbols.left.end() )
		{
			assert( it != _symbols.left.end() );
			std::cout<< "Error: symbol '"<< sym<< "' not in alphabet!"<< std::endl;
			throw SymbolNotInAlphabet();
		}
		
		_e[state][it->second] = emit_prob;
	}
public:	
	void SetFixedEmissionsState( const StateId& state, bool fixed )
	{
		if( fixed )
			_emission_class[state] = 0;
		else
			_emission_class[state] = state;
		_update_emission_classes_for_iter_impl();
	}
public:	
	void MergeEmissionsClass( const StateId &s1, const StateId& s2 )
	{
		if( _emission_class[s1] < _emission_class[s2] )
			_emission_class[s2] = _emission_class[s1];
		else
			_emission_class[s1] = _emission_class[s2];
		_update_emission_classes_for_iter_impl();
	}

protected:
	void _update_emission_classes_for_iter_impl()
	{
		_emission_class_set_for_iter_impl.clear();

		emission_classes_t::const_iterator it, it_end;
		it = _emission_class.begin();
		it_end = _emission_class.end();
		for( ; it != it_end; ++it )
			_emission_class_set_for_iter_impl.insert( *it );
	}

protected:
	void resize()
	{
		_a.resize( boost::extents[_states.size()][_states.size()] );
		_e.resize( boost::extents[_states.size()][_symbols.size()] );

		//_state_learning_mode.resize( boost::extents[ _states.size()][ num_state_parameters<self_type>::value ] );


		// TODO - remove the following line?
		_state_learning_mode.resize( boost::extents[ _states.size()][ 3 /* TODO FIX THIS*/ ]);


	}

protected:
	void clear()
	{
		_symbols.clear();
		_states.clear();
		_emission_class.clear();
		_emission_class_set_for_iter_impl.clear();
		resize();

		_init_fixed_parts();
	}
	
	

public:
	size_t num_symbols() const { return _symbols.size(); }
	size_t num_states() const { return _states.size(); }
	
public:
	enum { InitialState=0, TerminalState=1 };
	StateId GetInitialState() const { return StateId(InitialState); }
	StateId GetTerminalState() const { return StateId(TerminalState); }
	static bool IsReservedState(StateId s) { return (( s == InitialState) || ( s == TerminalState )); }

	bool AreEmissionsFixed( StateId s ) const
	{
		// Class 0 is reserved for silent and fixed states
		return _emission_class[s] == 0;
	}
	size_t GetEmissionsClass( StateId s ) const
	{
		return _emission_class[s];
	}
	

protected:
	friend class detail::EmissionsClassMembersIterator<emission_classes_t::const_iterator, StateId>;
public:
	typedef detail::EmissionsClassMembersIterator<emission_classes_t::const_iterator, StateId> emissions_class_members_iterator;

	std::pair<emissions_class_members_iterator,emissions_class_members_iterator> GetEmissionsClassMembersRange( StateId state ) const
	{
		emission_classes_t::const_iterator it = _emission_class.begin();
		emission_classes_t::const_iterator it_end = _emission_class.end();
		
		return std::make_pair( emissions_class_members_iterator( it,     it_end, state), 
							   emissions_class_members_iterator( it_end, it_end, state) );
	}


protected:
	friend class detail::EmissionsClassesIterator<StateId>;
public:
	typedef detail::EmissionsClassesIterator<StateId> emissions_classes_iterator;

	std::pair<emissions_classes_iterator,emissions_classes_iterator> GetEmissionsClassesRange() const
	{
		return std::make_pair( emissions_classes_iterator( _emission_class_set_for_iter_impl, _emission_class_set_for_iter_impl.begin() ), 
							   emissions_classes_iterator( _emission_class_set_for_iter_impl, _emission_class_set_for_iter_impl.end() ) );
	}

			
public:
  typedef enum { DefaultLearning=0, PresetValuesNoLearning=1 } state_learning_mode_t;

protected:
	
  typedef typename boost::bimap<Symbol, size_t> symbols_map_t;
  symbols_map_t _symbols;
  
  boost::multi_array<Probability, 2> _a;
  boost::multi_array<Probability, 2> _e;
  
  std::vector<std::string> _states;
  emission_classes_t _emission_class;
  detail::ClassesSet<StateId> _emission_class_set_for_iter_impl;
  
  
  boost::multi_array<state_learning_mode_t, 2> _state_learning_mode;
  
  bool isvalidsymbolid( size_t val ) const { return (val < _symbols.size()); }
  bool isvalidstateid( StateId val ) const { return (val < _states.size()); }
  
public:
	void SetDefaultProbabilities()
	{
		const Probability base = -std::numeric_limits<typename Algo::probability_type>::max();
		{
			const size_t m_end = _a.shape()[0];
			const size_t n_end = _a.shape()[1];
			
			for( size_t m=0; m< m_end; ++m )
				for( size_t n=0; n< n_end; ++n )
				{	  _a[m][n] = base;	}
		}
		
		{
			const size_t m_end = _e.shape()[0];
			const size_t n_end = _e.shape()[1];
			
			for( size_t m=0; m< m_end; ++m )
				for( size_t n=0; n< n_end; ++n )
				{	  _e[m][n] = base;	}
		}
		
	}
	

public:
	Probability a(const StateId& from, const StateId& to) const
	{
		assert( isvalidstateid( from ) ); assert( isvalidstateid( to ) );
		const Probability p = _a[from][to];
		return p;
	}
	Probability e(const StateId& state, const Symbol& sym) const
	{
		assert( isvalidstateid( state ) );
		typename symbols_map_t::left_map::const_iterator it;
		it = _symbols.left.find( sym );
		if( it == _symbols.left.end() )
		{
		  assert( it != _symbols.left.end() );
                  std::cerr<< "Error: symbol '"<< sym<< "' not in alphabet! ("<< __FILE__<< ":"<< __LINE__<< ")"<< std::endl;
		  throw SymbolNotInAlphabet();
		}
		
		assert( isvalidsymbolid( it->second ) );
		
		const Probability p = _e[state][it->second];
		return p;
	}
	Probability e(const StateId& state, const size_t& s) const
	{
		assert( isvalidstateid( state ) );
		assert( !is_silent(state) );
		assert( isvalidsymbolid( s ) );
		
		const Probability p = _e[state][s];
		return p;
	}

	const Symbol get_symbol(const size_t& idx) const
	{
		assert( isvalidsymbolid( idx ) );
		typename symbols_map_t::right_map::const_iterator it;
		it = _symbols.right.find( idx );
		if( it == _symbols.right.end() )
		{
			std::cout<< "Error: index "<< idx<< " does not reference an alphabet symbol!"<< std::endl;
			throw SymbolNotInAlphabet();
		}
		return it->second;
	}

	const size_t get_symbol_idx(const Symbol& sym) const
	{
		//assert( isvalidsymbolid( idx ) );
		typename symbols_map_t::left_map::const_iterator it;
		it = _symbols.left.find( sym );
		if( it == _symbols.left.end() )
		{
			std::cout<< "Error: symbol "<< sym<< " does not reference an alphabet symbol!"<< std::endl;
			throw SymbolNotInAlphabet();
		}
		return it->second;
	}

public:	
	Probability p(const StateId& state, size_t duration) const
	{
		// TODO Impl. this
		assert(false);
		throw;
	}
	

	const std::string GetStateName(const StateId& idx) const
	{
		assert( isvalidstateid( idx ) );
		return _states.at(idx);
	}
	
	bool is_silent( const StateId& state ) const
	{
		return IsReservedState(state);
	}
	
	void debug_print() const
	{
		debug_print(0);
	}
	
	void debug_print(size_t mode) const
	{
		std::cout<< " a: ";
		::debug_print_exp(_a);
		std::cout<< " e: ";
		::debug_print_exp(_e);
	}

	double emissions_entropy(const StateId& state) const
	{
		// Return the Shannon entropy of the p.m.f. defined by the emissions params for this state
		const size_t N = num_symbols();
		const double log2 = log(2.0);

		double sigma_i = 0.0;
		for( size_t i=0; i<N; ++i)
		{
			const Probability logp = _e[state][i];
			const double p = exp(logp);
			if( p < 1e-100)
				continue;
			//std::cout<< "k="<< state<< ", c="<< i<< ", p="<< p<< ", log(p)="<< logp<< ", sigma_i="<< sigma_i<< std::endl;//
			sigma_i += p * logp / log2;
		}
		//std::cout<< "k="<< state<< ", H(e)="<< -sigma_i<< std::endl;//
		
		return -sigma_i;
	}

public:
	void normalize_emissions()
	{
		const size_t N = num_symbols();
		const size_t S = num_states();
		for( size_t s=0; s<S; ++s )
		{
			if( is_silent(s) ) continue;
			
			Probability sigma_p = 0.0;
			for( size_t n=0; n<N; ++n )
				sigma_p += exp(_e[s][n]);

			const Probability factor = log(sigma_p);

			for( size_t n=0; n<N; ++n )
				_e[s][n] -= factor;
		}
		
	}
	


public:
  size_t /*state_learning_mode_t*/ GetStateLearningMode(const StateId& state, size_t param_class) const
	{
	  assert( isvalidstateid( state ) ); assert( param_class < 2 /*num_state_parameters<self_type>::value*/ );
	  const state_learning_mode_t mode = _state_learning_mode[state][param_class];
	  return (size_t)mode;
	}
public:
  void SetStateLearningMode( const StateId& state, size_t param_class, size_t /*state_learning_mode_t*/ mode )
  {
    std::cout<< state<< ", "<< param_class<< ", "<< mode<< std::endl;
    assert( isvalidstateid( state ) ); assert( param_class < 2 /*num_state_parameters<self_type>::value*/ );
    // TODO - validate 'mode' parameterf
    _state_learning_mode[state][param_class] = (state_learning_mode_t)mode;
  }

public:
    void SetEta(const StateId& state, double eta)
    {
        assert( isvalidstateid( state ) );
    }
    
public:
    void SetMu(const StateId& state, double mu)
    {
        assert( isvalidstateid( state ) );
    }
    
public:
    double GetEta(const StateId& state) const
    {
        assert( isvalidstateid( state ) );
        // The exponential distribution with param theta is equivalent to gamma(k=1.0, theta=1/theta)
        return _e[state][state];
    }
    
public:
    double GetMu(const StateId& state) const
    {
        assert( isvalidstateid( state ) );
        // The exponential distribution with param theta is equivalent to gamma(k=1.0, theta=1/theta)
        return 1.0;
    }
    
public:
    self_type& operator=(const self_type& other)
    {
        clear();
        // TODO - impl. this
        assert(false);
        return *this;
    }
    

public:
    typename Algo::sequence_type GetAlphabet() const
    {
        typename Algo::sequence_type alphabet( _symbols.size(), ' ');

        typename symbols_map_t::right_map::const_iterator it, it_end;
        it = _symbols.right.begin();
        it_end = _symbols.right.end();
        size_t id = 0;
        for( ; it != it_end; ++it, ++id )
        {
            alphabet[id] = it->second;
        }

        return alphabet;
    }
	
};	

struct HSMMMatrixModelBase {};  // Base type (used for template specialization)

template<typename Algo>
class HSMMMatrixModel : public MatrixModel<Algo>, public HSMMMatrixModelBase
{
public:

  //typedef typename Algo::model_type::state_id_type StateId;
  typedef MatrixModel<Algo> base_type;
  typedef HSMMMatrixModel<Algo> self_type;
  
  typedef typename base_type::StateId StateId;
  typedef typename base_type::Symbol Symbol;
  typedef typename base_type::Probability Probability;
  typedef boost::mpl::bool_<true> is_explicit_duration_model;
  enum { max_duration = Algo::max_duration };

    HSMMMatrixModel() {}

public:
  size_t get_max_duration() const { return max_duration; }

protected:
	boost::multi_array<Probability, 2> _duration_params;

protected:
	enum {_eta=0, _mu=1 };

	inline double _get_eta(const StateId& state) const { return _duration_params[state][_eta]; }
	inline double _get_mu (const StateId& state) const { return _duration_params[state][_mu];  }
	using base_type::_states;

protected:
	void resize()
	{
		base_type::resize();
		_duration_params.resize( boost::extents[_states.size()][2] );
	}

protected:
	void clear()
	{
		_duration_params.clear();
		base_type::clear();
	}

public:
  virtual ~HSMMMatrixModel() {}
	

public:
	void SetEta(const StateId& state, double eta)
	{
	  assert( base_type::isvalidstateid( state ) );
		//assert( !IsReservedState(state) );
		_duration_params[state][_eta] = eta;
	}

public:
	void SetMu(const StateId& state, double mu)
	{
	  assert( base_type::isvalidstateid( state ) );
		//assert( !IsReservedState(state) );
		_duration_params[state][_mu] = mu;
	}

public:
	double GetEta(const StateId& state) const
	{
	  assert( base_type::isvalidstateid( state ) );
		//assert( !IsReservedState(state) );
		return _duration_params[state][_eta];
	}

public:
	double GetMu(const StateId& state) const
	{
	  assert( base_type::isvalidstateid( state ) );
		//assert( !IsReservedState(state) );
		return _duration_params[state][_mu];
	}


public:	
	Probability p(const StateId& state, size_t duration) const
	{		
		using namespace boost::math;
		assert(duration <= Algo::max_duration );

		if( base_type::IsReservedState(state) ) return (duration==1) ? 0.0 : -std::numeric_limits<Probability>::max();

		gamma_distribution<Probability> dist(_get_mu(state),
											 1./_get_eta(state) );

		// TODO make sure this is the most accurate discretization of the gamma function.
		//      This function must meet the following conditions:
		//      1. Sum to 1.0 (over 1>=duration>=T)
		//      2. Have the shortest distance (e.g. MSE) to the gamma distribution (since we estimate parameters as if the real gamma is used)

		double ubound = double(duration)+0.5;
		double lbound = double(duration)-0.5;
		if( duration == 1 )
			lbound = 0.0;

		return log(
			cdf(dist, ubound   ) -
			cdf(dist, lbound   )   
			);
/*
		return log(
			cdf(dist, double(duration  )   ) -
			cdf(dist, double(duration-1)   )   
			);
*/
	}


public:
  size_t /*state_learning_mode_t*/ GetStateLearningMode(const StateId& state, size_t param_class) const
	{
	  assert( base_type::isvalidstateid( state ) ); assert( param_class < 3 /* TODO FIX THIS */ /*num_state_parameters<self_type>::value*/ );
	  const typename base_type::state_learning_mode_t mode = base_type::_state_learning_mode[state][param_class];
	  return (size_t)mode;
	}
public:
  void SetStateLearningMode( const StateId& state, size_t param_class, size_t /*state_learning_mode_t*/ mode )
  {
    std::cout<< state<< ", "<< param_class<< ", "<< mode<< std::endl;
    assert( base_type::isvalidstateid( state ) ); assert( param_class < 3 /*TODO FIX THIS*/ /*num_state_parameters<self_type>::value*/ );
    // TODO - validate 'mode' parameterf
    base_type::_state_learning_mode[state][param_class] = (typename base_type::state_learning_mode_t)mode;
  }

public:
	void SetDefaultProbabilities()
	{
		base_type::SetDefaultProbabilities();

		{
			const size_t m_end = _duration_params.shape()[0];
			for( size_t m=0; m< m_end; ++m )
			{
				_duration_params[m][_mu]  =  4.00;
				_duration_params[m][_eta] =  2.00; 
			}
		}
		
	}

public:
	void debug_print() const
	{
		debug_print(0);
	}

public:
	void debug_print(size_t mode) const
	{
		if( mode == 0 )
			base_type::debug_print(mode);
		
		if( mode == 0 || mode == 4)
		{
			std::cout<< " duration (eta, mu): ";
			::debug_print(_duration_params);
		}
		
	}

public:
	void split_state(const StateId& orig_state, double balance=0.5)
	{
		const StateId new_state1 = orig_state;
		const StateId new_state2 = base_type::num_states();
		assert(balance>0.0 && balance<1.0);
		
		// Add the new state
		base_type::AddState( new_state2, "split" );
		resize();

		// Update transitions
		for(StateId i=0; i<base_type::num_states(); ++i)
		{
		  base_type::SetTransitionLogspace(new_state2,
								  i,
					      base_type::a(orig_state, i) );
		  base_type::SetTransitionLogspace(orig_state,
								  i,
								  (i==new_state2) ? 0.0 : -std::numeric_limits<Probability>::max() );
		  base_type::SetTransitionLogspace(i,
								  new_state2,
								  (i==orig_state) ? 0.0 : -std::numeric_limits<Probability>::max() );

			
		}

		// Copy emissions to the new state
		for( size_t s=0; s<base_type::num_symbols(); ++s)
		{
		  base_type::SetEmissionProbabilityLogspace( new_state2,
											base_type::get_symbol(s),
							base_type::e(orig_state, s) );
		}

		// Update duration params
		const double orig_eta = GetEta(orig_state);
		const double orig_mu = GetMu(orig_state);
		const double new_mu1 = orig_mu * balance;
		const double new_mu2 = orig_mu * (1.0 - balance);
		
		SetEta( new_state1, orig_eta );
		SetMu ( new_state1, new_mu1  );
		
		SetEta( new_state2, orig_eta );
		SetMu ( new_state2, new_mu2  );

#ifndef NO_TESTS
		{
			tests::test_model_valid<Algo>(*this);
		}
#endif //NO_TESTS

	}

public:
	self_type& operator=(const self_type& other)
	{
		clear();
		// TODO impl. this
		assert(false);
		return *this;
	}
	
	
};


template<class Algo>
void debug_print( const MatrixModel<Algo>& obj, size_t mode=0 )
{
	std::cout<< "MatrixModel<Algo>"<< std::endl;
	obj.debug_print(mode);
}

template<typename Algo>
size_t num_symbols(const MatrixModel<Algo>& model )
{
	return model.num_symbols();
}

template<typename Algo>
size_t num_states(const MatrixModel<Algo>& model )
{
	return model.num_states();
}

template<class Algo>
void debug_print( const HSMMMatrixModel<Algo>& obj, size_t mode=0 )
{
	std::cout<< "HSMMMatrixModel<Algo>"<< std::endl;
	obj.debug_print(mode);
}

template<typename Algo>
size_t num_symbols(const HSMMMatrixModel<Algo>& model )
{
	return model.num_symbols();
}

template<typename Algo>
size_t num_states(const HSMMMatrixModel<Algo>& model )
{
	return model.num_states();
}

namespace detail
{
	

typedef std::vector<size_t> order_map_t;
	
template<typename Algo>
order_map_t build_foreign_order_map(const std::string& order, const typename Algo::model_type& model)
{
	assert(order.size() == num_symbols(model) );

	order_map_t mapping(order.size());
	
	for( size_t i=0; i< order.length(); ++i )
	{
		const std::string::value_type foreign = order[i];
		const size_t local = model.get_symbol_idx(foreign);
		
		mapping[local] = i;
	}
	return mapping;
}
 
} // namespace detail

template<typename Algo>
void relax_emissions(typename Algo::model_type& model, double level, const boost::multi_array<typename Algo::probability_type,2>& substitutions, const std::string& order)
{
	assert(level <= 1.0); assert(level >= 0.0);

	typedef typename Algo::probability_type P;
	typedef LogspaceDouble<P> logspace_t;
	
	const size_t N = num_states(model);
	const size_t S = num_symbols(model);

	const detail::order_map_t local_to_foreign = detail::build_foreign_order_map<Algo>(order, model);

	for( size_t state=0; state < N; ++state)
	{
		if( model.is_silent(state) ) continue;
		
		for( size_t symbol=0; symbol<S; ++symbol )
		{
			const P orig = logspace_t(model.e(state, symbol), true);
			const size_t symbol_f = local_to_foreign[symbol];

			//std::cout<< "  --> symbol: "<< model.get_symbol(symbol)<< " orig: "<< orig<< std::endl;

			P sigma_others = 0.0;
			//std::cout<< "others: ";
			for( size_t other=0; other<S; ++other )
			{
				const size_t other_f = local_to_foreign[other];
				
				const P p = P(logspace_t(model.e(state, other), true)) * substitutions[other_f][symbol_f];
				//std::cout<< model.get_symbol(other)<< "->"<< p<< " ";
				
				sigma_others += p;
			}
			//std::cout<< std::endl;

			const P relaxed =
				orig * (1-level)
				+ sigma_others * level;

			//std::cout<< relaxed << " = ("<< orig<< " * "<< (1-level)<< ") + ("<< sigma_others<< " * "<< level<< ")"<< std::endl;

			model.SetEmissionProbability(state,
										 model.get_symbol(symbol),
										 relaxed); // accepts probability in real space
		}
	}

	model.normalize_emissions();
	
#ifndef NO_TESTS
	tests::test_emissions_valid_pdf<Algo>(model);
#endif //NO_TESTS	

}

namespace v2
{


namespace model
{

/*
 * Define model traits for HMM models
 */
template <typename Model>
struct model_traits<
    Model
  , typename boost::enable_if<
      boost::mpl::and_<
	boost::is_base_of< HMMMatrixModel, Model>                         // model derives from HMM...
      , boost::mpl::not_< boost::is_base_of< HSMMMatrixModelBase, Model> > // ...but not from HSMM
      >
    >::type
>
{
  typedef markov_chain hidden_process;
};

/*
 * Define model traits for HSMM models
 */
template <typename Model>
struct model_traits<
    Model
  , typename boost::enable_if<
      typename boost::is_base_of< HSMMMatrixModelBase, Model>
    >::type
>
{
  typedef semi_markov_chain hidden_process;
};


} // namespace model

} // namespace v2

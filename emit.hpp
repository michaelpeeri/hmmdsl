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
#include "random.hpp"
#include "model.hpp"

template<typename Algo>
class random_emitter
{
public:
	random_emitter( boost::shared_ptr<const typename Algo::model_type> model )
		: _model(model)
		, _random( 0, 1 )
	{
#ifndef NO_TESTS
		tests::test_model_valid<Algo>( *model );
#endif // NO_TESTS
	}

public:
	// Constructor for python-owned model
	random_emitter( const typename Algo::model_type& model )
		: _model( &model, null_deleter() )
		, _random( 0, 1 )
	{
#ifndef NO_TESTS
		tests::test_model_valid<Algo>( model );
#endif // NO_TESTS
	}
	
public:						  	
	typename Algo::sequence_type emit()
	{
		typename Algo::sequence_type seq;
		

		const typename Algo::model_type::StateId terminal = _model->GetTerminalState();
		typename Algo::model_type::StateId state = _model->GetInitialState();
		
		size_t len = 0;
		size_t d = 1; // default for regular HMMs
		
		while( (state != terminal) || (len == 0) )
		{
			if( Algo::model_type::is_explicit_duration_model::value )
			{
				d = random_duration(state);
			}
			
			if( !_model->is_silent( state ) )
			{
				for( size_t tau=1; tau <= d; ++tau )
					seq.append( 1, random_emission( state ) );
			}
			
			state = random_transition( state );
			++len;
		}
		
		return seq;
	}

protected:
	typename Algo::model_type::StateId random_transition( typename Algo::model_type::StateId source )
	{
		const double r = _random();
		
		double p = r;
		typename Algo::model_type::StateId dest = 0;
		while( p > 0)
		{
			p -= exp(_model->a( source, dest ));
			++dest;
		}
		return dest-1;
	}

public:
	double test_random_transition( size_t sample_size = 1000000 )
	{
		const typename Algo::model_type::StateId states_end = num_states(*_model);
		using namespace boost::accumulators;
		accumulator_set< typename Algo::probability_type, stats< tag::min, tag::max, tag::mean > > acc;
		
		for( typename Algo::model_type::StateId source = 0; source < states_end; ++source )
		{
			if( source == _model->GetTerminalState() ) continue; // No transitions for terminal state
			
			std::vector<typename Algo::model_type::StateId> counts(states_end);

			for( size_t times=0; times < sample_size; ++times )
			{
				const typename Algo::model_type::StateId dest = random_transition( source );
				assert( dest >= 0 );
				assert( dest < states_end );
				++counts[ dest ];
			}

			for( typename Algo::model_type::StateId dest = 0; dest < states_end; ++dest )
			{
				const double p = double(counts[dest]) / double(sample_size);
				const double error = p - exp(_model->a( source, dest ));
				acc( error );
			}
		}
		const double min_err = min(acc);
		const double max_err = max(acc);
		
		//std::cout<< min_err<< " <= Error <= "<< max_err<< "\t(average = "<< mean(acc)<< ")"<< std::endl;		
		return (-min_err > max_err) ? -min_err : max_err;
	}

protected:
	typename Algo::symbol_type random_emission( typename Algo::model_type::StateId source )
	{
		const double r = _random();
		
		double p = r;
		size_t s = 0;
		while( p > 0)
		{
			p -= exp(_model->e( source, s ));
			++s;
		}
		return _model->get_symbol(s-1);
	}

protected:
	typename Algo::model_type::StateId random_duration( typename Algo::model_type::StateId source )
	{
		const double r = _random();
		
		double p = r;
		size_t d = 1;
		while( p > 0)
		{
			p -= exp(_model->p( source, d ));
			++d;
		}
		return d-1;
	}
	
public:
	typename Algo::sequence_type emit_fixed_length( size_t length )
	{
		typename Algo::sequence_type seq;
		
		//typename Algo::model_type::StateId state = _model->GetInitialState();
		
		for( size_t i = 0; i < length; ++i )
		{
			break;
		}
		
		return seq;
	}
	
	
						  
protected:
	boost::shared_ptr<const typename Algo::model_type> _model;

	detail::random_gen<double> _random;
};

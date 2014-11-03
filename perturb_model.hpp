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
class perturb_model
{
public:
	perturb_model( boost::shared_ptr<typename Algo::model_type> model )
		: _model(model)
		, _random( 0, 1 )
	{
#ifndef NO_TESTS
		tests::test_model_valid<Algo>( *model );
#endif // NO_TESTS
	}

public:
	// Constructor for python-owned model
	perturb_model( typename Algo::model_type& model )
		: _model( &model, null_deleter() )
		, _random( 0, 1 )
	{
#ifndef NO_TESTS
		tests::test_model_valid<Algo>( model );
#endif // NO_TESTS
	}

public:
	void perturb_emissions( double amount=0.2 )
	{
		assert(amount <= 1.0); assert(amount >= 0.0);
		typedef typename Algo::probability_type P;
		typedef LogspaceDouble<P> logspace_t;

		const size_t N = num_states(*_model);
		const size_t S = num_symbols(*_model);

		for( size_t state=0; state < N; ++state)
		{
			if( _model->is_silent(state) ) continue;
			for( size_t symbol=0; symbol<S; ++symbol )
			{
				const P orig = logspace_t(_model->e(state, symbol), true);

				const P random = _random();

				const P perturbed =
					orig * (1-amount)
					+ random * amount;

				_model->SetEmissionProbability(state,
											   _model->get_symbol(symbol),
											   perturbed); // accepts probability in real space

			}
		}

		_model->normalize_emissions();
		
#ifndef NO_TESTS
		tests::test_emissions_valid_pdf<Algo>(*_model);
#endif //NO_TESTS	
		
	}

						  
protected:
	boost::shared_ptr<typename Algo::model_type> _model;

	detail::random_gen<double> _random;

};

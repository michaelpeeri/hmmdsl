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
#include <boost/shared_ptr.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/nondet_random.hpp>
#include "common.hpp"



namespace detail
{
	

template<typename Type, typename RandomGenerator=boost::mt19937>
class random_gen
{
public:
	random_gen( Type from, Type to)
		: _rand01( from, to )
		, _random( _rng, _rand01 )
	{
		_seed();
	}

public:
	inline Type operator()() { return _random(); }
	
protected:
	void _seed()
	{
            // Bootstrap main RNG using OS entropy pool (uses /dev/urandom by default, which seems more than good enough for seeding)
            //
            // See Rule 3 on:
            // David Jones, "Good Practice in (Pseudo) Random Number Generation for Bioinformatics Applications"
            // http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf
            // (Accessed 2014-11-19)
            //
            // For applications that strongly rely on random numbers (e.g. MC), see also Click et al., PMID PMC2992609, http://dx.doi.org/10.1002/jcc.21638 
            //
            boost::random_device seed_gen;
            _rng.seed(seed_gen);
	}

	// Primary RNG
	RandomGenerator _rng;
	typedef boost::uniform_real<> uniform_real_t;
	uniform_real_t _rand01;
	boost::variate_generator< RandomGenerator&, uniform_real_t> _random;

};


} // namespace detail

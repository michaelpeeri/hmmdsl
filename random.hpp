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

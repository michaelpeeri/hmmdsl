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
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/nondet_random.hpp>
#include "logspace.hpp"

namespace test_logspace
{
	
namespace detail
{

template<typename RandomGenerator=boost::mt19937>
class random_gen
{
public:
	random_gen()
		: _rand01( 0, 1 )
		, _random( _rng, _rand01 )
	{
		_seed();
	}


protected:
	void _seed()
	{
		boost::random_device seed_gen;
		_rng.seed(seed_gen);
	}

public:
	typedef enum
	{ plus, plusequals, times, timesequals
	} operation;
	
		

public:
	operation random_operator()
	{
		const double r = _random();
//		if( r < 0.05 )
//			return plus;
//		else if( r < 0.10 )
//			return plusequals;
//		else 
		if( r < 0.50 )
			return times;
		else 
			return timesequals;
	}

public:
	double random01()
	{
		return _random();
	}
	


						  
protected:
	// Primary RNG
	RandomGenerator _rng;
	typedef boost::uniform_real<> uniform_real_t;
	uniform_real_t _rand01;
	boost::variate_generator< RandomGenerator&, uniform_real_t> _random;

};
	
} // namespace detail;


template<typename double_t1, typename double_t2>
void test_operations(size_t N=100)
{
	detail::random_gen<> random;
	const double epsilon = 1e-12;
	
	double_t1 f1= 1.0;
	double_t2 f2= 1.0;
	
	for(size_t n=0; n<N; ++n)
	{
		detail::random_gen<>::operation op = random.random_operator();
		const double_t1 operand1 = random.random01();
		const double_t2 operand2(operand1*1.001);

		if( (f1 == double_t1(0.0)) && (f2 == double_t2(0.0)) )
		{
			std::cout<< "All values at 0 after "<< n<< " iterations"<< std::endl;
			break;
		}

		switch(op)
		{
		case detail::random_gen<>::plus:
			f1 = f1 + operand1;
			f2 = f2 + operand2;
			break;

		case detail::random_gen<>::plusequals:
			f1 += operand1;
			f2 += operand2;
			break;

		case detail::random_gen<>::times:
			f1 = f1 * operand1;
			f2 = f2 * operand2;
			break;

		case detail::random_gen<>::timesequals:
			f1 *= operand1;
			f2 *= operand2;
			break;

		default:
			assert(false);
		}

		const double result1 = f1;
		const double result2 = f2;
		if( ( fabs(result1-result2)/f1 > 0.01 ) ||
			( fabs(result1-result2)/f2 > 0.01 ) )
		{
			std::cout<< "Error: iter="<< n<< ", result1: "<< result1<< ", result2: "<< result2<< std::endl;
			//assert( fabs(result1-result2) < epsilon );
		}
		

		//if( n>0 && (n%100000 == 0))
		if( n>0 && (n%10 == 0))
		{
			std::cout<< "test_operations(iter="<< n<< "): error="<< fabs(double(f1)-double(f2))<< ", f1="<< double(f1)<< ", f2=";
			f2.debug_print();
			std::cout<< std::endl;
		}
		
	}

	const double result1 = f1;
	const double result2 = f2;

	if( fabs(result1-result2) >= epsilon )
	{
		std::cout<< "result1: "<< result1<< ", result2: "<< result2<< std::endl;
		assert( fabs(result1-result2) < epsilon );
	}
}

} // namespace test_logspace


#ifdef MainProgram

int main()
{
	std::cout<< "Logspace Test"<< std::endl;
	typedef LogspaceDouble<> logspace_t;
	
	logspace_t a = 0.5;
	logspace_t b = 0.2;
	logspace_t c = 0.1;
	
	logspace_t d = a*b + c;
	d += b;
	d *= 0.1;
	
	std::cout<< "(0.5*0.2) + 0.1 + 0.2 = "<<d<< std::endl;

	test_logspace::test_operations<double, logspace_t>(10000);

	return 0;
}

#endif //MainProgram

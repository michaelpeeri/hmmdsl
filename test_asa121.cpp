#include <iostream>
#include <cassert>
#include "asa121.hpp"

int main()
{
	int err = 0;

	for(double x= 0.000000012; x< 0.00025 ; x += 0.0000135 )
	{
		std::cout<< x<< " "<< trigam(x, &err)<< std::endl;
		assert(err == 0);
	}

	for(double x= 0.00025; x< 10.0 ; x += 0.0001525 )
	{
		std::cout<< x<< " "<< trigam(x, &err)<< std::endl;
		assert(err == 0);
	}

	for(double x= 10.0001257; x< 2000.0 ; x += 0.1931525 )
	{
		std::cout<< x<< " "<< trigam(x, &err)<< std::endl;
		assert(err == 0);
	}
	
	return 0;
}

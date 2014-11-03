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

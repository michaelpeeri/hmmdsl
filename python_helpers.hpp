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
#include <vector>
#include <boost/python.hpp>


template<typename T, typename PythonContainer>
void python_to_stl_vector(const PythonContainer& l, std::vector<T>& vec)
{
	vec.clear();
	
	const size_t N = boost::python::len(l);

	vec.reserve(N);
	for(size_t n=0; n<N ; ++n)
	{
		vec.push_back( boost::python::extract<T>(  l[n] ) );
	}
	
}

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

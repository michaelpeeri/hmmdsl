#include "common.hpp"


size_t length(const std::string& t) { return t.length(); }

template<>
std::string stringify<std::vector<std::string> >(const std::vector<std::string>& a, size_t from, size_t to )
{
	std::string out;
	if( to > a.size() ) to = a.size();
	for( size_t pos = from; pos< to; ++pos )
		out.append( a[pos] );
	return out;
}

template<>
std::string stringify<std::string>( const std::string& a, size_t from, size_t to )
{
	std::string out;
	if( to > a.size() ) to = a.size();
	for( size_t pos = from; pos< to; ++pos )
		out.append( 1, a[pos] );
	return out;
}


bool not_suspected_negative_subscript(size_t k) { return k < (std::numeric_limits<size_t>::max() - 100); }

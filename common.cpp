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

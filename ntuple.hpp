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
#include <boost/mpl/vector.hpp>
#include <boost/mpl/push_back.hpp>

namespace ntuple
{
	

template<class T, size_t N>
struct ntuple
{
	typedef typename boost::mpl::push_back<typename ntuple<T,N-1>::type, T>::type type;
};

template<class T>
struct ntuple<T,1>
{
	typedef typename boost::mpl::vector<T>::type type;
};


}

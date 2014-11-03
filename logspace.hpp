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
#include <math.h>
#include <limits>



template<typename P>
static P logspace_sum( P p, P q )
{
  // TODO - use optimization from Durbin p. 78)
	const P val = log( exp(p) + exp(q) );
	if( val == -std::numeric_limits<P>::infinity() )
		return -std::numeric_limits<P>::max();
	return val;
}


namespace detail
{

template<typename float_type>
struct ToFloat
{
	inline float_type operator()(const float_type& logspace_val) const { return exp(logspace_val); }
};

template<typename float_type>
struct FromFloat
{
	inline float_type operator()(const float_type& realspace_val) const { return log(realspace_val); }
};
	
}



template<class float_type=double, typename ToFloat=detail::ToFloat<float_type>, typename FromFloat=detail::FromFloat<float_type> >
class LogspaceDouble
{
protected:
	typedef LogspaceDouble<float_type,ToFloat,FromFloat> self_type;
	
public:
	LogspaceDouble( float_type val ) :
		_val(_fromfloat(val))
		,_floatval(0.0)
		,_floatval_ready(false)
	{}

public:
	LogspaceDouble( float_type val, bool ) :
		_val(val)
		,_floatval(0.0)
		,_floatval_ready(false)
	{}
	
public:
	LogspaceDouble( const LogspaceDouble<float_type,ToFloat,FromFloat>& other)
		: _val(other._val)
		, _floatval(other._floatval)
		, _floatval_ready(other._floatval_ready)
	{}
	
public:
	inline operator float_type() const
	{
		return _getfloat();
	}

public:
	const LogspaceDouble& operator=( float_type val)
	{
		_val = _fromfloat(val);
		_floatval_ready = false;
	}
	
public:
	self_type& operator*(const self_type& other)
	{
		_val += other._val;
		_floatval_ready = false;
		return *this;
	}

public:
	self_type& operator*=(const self_type& other)
	{
		_val += other._val;
		_floatval_ready = false;
		return *this;
	}

public:
	self_type& operator+(const self_type& other)
	{
		_val = _fromfloat( _getfloat() + other._getfloat() );
		if( _val == -std::numeric_limits<float_type>::infinity() )
			_val = -std::numeric_limits<float_type>::max();
		_floatval_ready = false;
		return *this;
	}

public:
	self_type& operator+=(const self_type& other)
	{
		_val = _fromfloat(_getfloat() + other._getfloat() );
		if( _val == -std::numeric_limits<float_type>::infinity() )
			_val = -std::numeric_limits<float_type>::max();
		_floatval_ready = false;
		return *this;
	}
	
	
protected:
	float_type _getfloat() const
	{
		if( _floatval_ready )
			return _floatval;
		
		_floatval = _tofloat(_val);
		_floatval_ready = true;
		return _floatval;
	}

public:
	void debug_print() const
	{
		std::cout<< "exp("<< _val<< ")="<< _getfloat();
	}
	
	
protected:
	float_type _val;
	mutable float_type _floatval;
	mutable bool _floatval_ready;

	ToFloat _tofloat;
	FromFloat _fromfloat;
};

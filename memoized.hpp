#pragma once
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#undef BOOST_DISABLE_ASSERTS
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/at_c.hpp> 
#include <boost/mpl/at.hpp>

template<class Func>
struct FuncTraits{};

// TODO - Add support for N!=2
template<class Func, class Env, class Args>
struct MemoizedFunc : public Func
{
public:
	typedef typename Func::return_type P;

protected:
	enum {_empty = FuncTraits<Func>::out_of_range_value	};

public:
	P fli( const Args& args )
	{
		boost::array<size_t,2> args2;
		args2[0] = boost::fusion::at_c<0>(args);
		args2[1] = boost::fusion::at_c<1>(args);

		// Negative arguments are not supported (since we're memoizing values!)
		assert(not_suspected_negative_subscript(args2[0]));
		assert(not_suspected_negative_subscript(args2[1]));
				
		P val = _v(args2);
	  
		if( val == _empty )
		{
			_state_reset = false;
			_v(args2) = val = Func::_produce_fli( args );
			assert( val != _empty);
		}
		return val;
  }

  inline P operator()( const Args& args) { return fli(args); }

protected:

	template<class Extents>
	MemoizedFunc( const Env& env, Extents ext )
		: Func(*this, env )
		,   _v(ext)
	{
		_reset();
	}
protected:
	void _reset()
	{
		// TODO Make this work
		//    for( typename v_type::iterator it = _v.begin(); it != _v.end(); ++it )
		//    *it = 0.0;
		const size_t m_end = _v.shape()[0];
		const size_t n_end = _v.shape()[1];
		
		for( size_t m=0; m< m_end; ++m )
			for( size_t n=0; n< n_end; ++n )
			{
				_v[m][n] = _empty;
			}
		
		//std::fill(_v.begin(), _v.end(), P(_empty) );
		//std::fill(_ptr.begin(), _ptr.end(), prev_ptr_t(_empty) );
		_state_reset = true;
	}
public:
	void reset()
	{
		if( _state_reset ) return;
		_reset();
	}


public:
	void debug_print() const
	{
		std::cout<< "MemoizedFunc<Algo>"<< std::endl;
		std::cout<< " v: ";
		::debug_print(_v);
	}
	

protected:
  
	typedef boost::multi_array<P,2> v_type;
	v_type _v;
	bool _state_reset;
	
};


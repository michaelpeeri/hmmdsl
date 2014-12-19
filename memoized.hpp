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
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#undef BOOST_DISABLE_ASSERTS
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/at_c.hpp> 
#include <boost/fusion/include/accumulate.hpp>
#include <boost/mpl/at.hpp>
#include "v2.hpp"

template<class Func>
struct FuncTraits{};

// TODO - Add support for N!=2
template<class Func, class Env, class Args>
struct MemoizedFunc : public Func
{
public:
    typedef typename Func::return_type P;
    typedef typename FuncTraits<Func>::val_type val_type;
    
protected:
    const val_type _empty;
    
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
            _v(args2) = val = this->_produce_fli( args );
            assert( val != _empty);
        }
        return val;
    }
    
    inline P operator()( const Args& args) { return fli(args); }
    
protected:
    
    template<class Extents>
    MemoizedFunc( const Env& env, Extents ext )
        : Func(*this, env )
        ,   _empty(typename FuncTraits<Func>::empty_val()())
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


namespace v2
{


/*
 * helper functor: set consecutive items in a vector<size_t>
 */
struct setVectorItems
{
  setVectorItems( std::vector<size_t>& vec ) : _vec(vec), _pos(0) {}
  
  template <typename T>
  void operator()(const T& val) const {     _vec[_pos++] = val;      }
private:
  std::vector<size_t>& _vec;
  mutable size_t _pos;
};

struct EInvalidValueProduced {};
struct EInvalidElementRequested {};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Note: Dependent types cannot be inferred based on the CRTP class 'Impl', because it is not fully instantiated yet (circular dependency).
// So Args is passed explicitly as an additional parameter
template <typename Impl, typename Args, typename Result, typename Level>
class memoize
{
public:

  /*
  template <typename Context>
  result_type calc(args_t const & args, Context const & context)
  {
    std::cout<< "Memoization works!"<< std::endl;
    this->get_extents(args, context);
    return this->calc2(args, context);
  }
  */

public:
  memoize() : _storage(Impl::get_unassigned_value()) {}

private:
  struct extents_helper
  {
    /* ExtentsCons: generate fusion sequences of the kind Elem^N
       e.g.: ExtentsCons<int,3> = [int,int,int]
     */
  template <typename Elem, size_t N>
  struct ExtentsCons
  {
    typedef boost::fusion::cons<Elem, typename ExtentsCons<Elem,N-1>::type> type;
  };
  template <typename Elem>
  struct ExtentsCons<Elem,0>
  {
    typedef boost::fusion::nil type;
  };

  // Testing only
  typedef typename ExtentsCons<size_t, 5>::type cons5;
  typedef typename boost::fusion::result_of::size<cons5>::type size_integral_constant5;
  BOOST_MPL_ASSERT_RELATION( size_integral_constant5::value, ==, 5 );
  // Testing only
  };

  typedef Impl impl_type;

  struct helper
  {
  typedef Args arg_type; // hide arg_type from derived classes
  BOOST_MPL_ASSERT(( boost::fusion::traits::is_sequence<arg_type> )); // To be used with memoize, the argument must be a fusion sequence.
  typedef Result result_type;
  };

private:
  //  const Result _unassigned;


public:
  typedef typename 
  boost::fusion::result_of::as_vector<
  typename extents_helper::template ExtentsCons<
        size_t
      , boost::fusion::result_of::size<typename helper::arg_type>::type::value
    >::type
  >::type dimensions_t;

  // Note: this assert is meant to make sure Dimensions is a fusion sequence.
  // TODO: replace with a concept check?
  typedef typename boost::fusion::result_of::size<dimensions_t>::type _size_of_dimensions_t;
  BOOST_MPL_ASSERT_RELATION( _size_of_dimensions_t::value, >=, 0 );
  BOOST_MPL_ASSERT(( boost::fusion::traits::is_sequence<dimensions_t > ));

  typedef typename boost::fusion::result_of::size<typename helper::arg_type>::type _size_of_arg_type;
  BOOST_MPL_ASSERT_RELATION( _size_of_arg_type::value, ==, _size_of_dimensions_t::value );


protected:
  /*
   * Hold the storage of previously computed values
   * Encapsulate implementation details
   */
  struct storage
  {
  public:
    storage( const Result& unassigned ) : _unassigned(unassigned) {}

  private:
    template <typename T, size_t Dims>
    struct container_type
    {
      typedef boost::multi_array<
        T,Dims> type;
    };

    /*
     * Store a single value; provide a compatible interface to the impl. used for storage for one or more dimensions
     */
    template <typename T>
    struct single_storage
    {
      typedef T element;
      // TODO - add concept checks for T
      //typedef boost::mpl::int_<0>::value dimensionality;
      BOOST_STATIC_CONSTANT(std::size_t, dimensionality = 0);
      typedef T* iterator;
      typedef const T* const_iterator;
    public:
      T& operator()(const boost::fusion::vector<>&)      {      return _val;      }
    public:
      const T& operator()(const boost::fusion::vector<>&) const      {  return _val;      }
    public:
      void resize(const boost::fusion::vector<>&) { /**/ }
    public:
      iterator begin() { return &_val; }
    public:
      const_iterator begin() const { return &_val; }
    public:
      iterator data() { return &_val; }
    public:
      const_iterator data() const { return &_val; }
    public:
      iterator end() { return (&_val)+1; }
    public:
      const_iterator end() const { return (&_val)+1; }
    public:
      size_t num_elements() const { return 1; }
    protected:
      T _val;
    public:
      const size_t* shape() const { static size_t _shape[1] = { 0 }; return _shape; } /* Strictly speaking a 0-length array is needed... */
    };
    
    // Specialization 1: Use a custom class for 0-dimensional storage
    template <typename T>
    struct container_type<T,0>
    {
      typedef single_storage<T> type;
    };
    /*
    // Specialization 2: Use std::vector for 1-dimensional storage of non-POD types (e.g. alignments)
    // TODO - activate this specialization if it's necessary
    template <typename T>
    enable_if< not_<is_pod<T> > > // TODO FIX THIS!
    struct container_type<T,1>
    {
      typedef std::vector<T> type;
    };
    */

    // Determine the data-type used for storage of memoized values
  private:
    typedef typename container_type<
      typename helper::result_type
      , _size_of_dimensions_t::value
      >::type container_t;
    // Trivial asserts
    BOOST_MPL_ASSERT(( boost::is_same<
                         typename helper::result_type
                       , typename container_t::element
                       > ));
    BOOST_MPL_ASSERT_RELATION( container_t::dimensionality, ==, _size_of_dimensions_t::value );

  private:
    // Instance to hold values
    container_t _values;
    const Result _unassigned;

  public:
    // Allow access to values;
    // Leave the  underlying storage impl. to define the type used to specify which item to expose.
    template <typename Address>
    typename helper::result_type& access( Address& address )
    {
      std::vector<size_t> t(_size_of_dimensions_t::value );
      boost::fusion::for_each( address, setVectorItems(t) );
      
#ifdef EXCESSIVE_DEBUG_PRINTING
      {
        std::vector<size_t>::const_iterator it, it_end;
        it = t.begin();
        it_end = t.end(); 
        std::cout<< "access([";
        for( ; it != it_end; ++it )
          {
            std::cout<< *it<< ", ";
          }
        std::cout<< "])"<< std::endl;
      }
#endif

      return _values( t ); // TODO: convert Address to something multi_array would accept (in O(1)!)
    }
    template <typename Address>
    const typename helper::result_type& access( const Address& address ) const
    {
      std::vector<size_t> t(_size_of_dimensions_t::value );
      boost::fusion::for_each( address, setVectorItems(t) );
      
#ifdef EXCESSIVE_DEBUG_PRINTING
      {
        std::vector<size_t>::const_iterator it, it_end;
        it = t.begin();
        it_end = t.end(); 
        std::cout<< "const access([";
        for( ; it != it_end; ++it )
          {
            std::cout<< *it<< ", ";
          }
        std::cout<< "])"<< std::endl;
      }
#endif

      return _values( address );
    }

  public:
    void allocate(const dimensions_t& extents)
    {
      std::vector<size_t> t(_size_of_dimensions_t::value );
      boost::fusion::for_each( extents, setVectorItems(t) );

#ifdef EXCESSIVE_DEBUG_PRINTING
      std::cout<< ">>>allocate: "<< extents<< " -> [";
      {
	std::vector<size_t>::const_iterator it, it_end;
	it = t.begin();
	it_end = t.end();
	for( ; it != it_end; ++it )
	{
	  std::cout<< *it<< ", ";
	}
      }
      std::cout<< "]"<< std::endl;
#endif

      _values.resize(t);
    }
    
  public:
    void clear()
    {
      // Reference: http://www.boost.org/doc/libs/1_53_0/libs/multi_array/example/foreach_test.cpp
      std::fill_n( _values.data(), _values.num_elements(), _unassigned );
      /*
      typename container_t::const_iterator it, it_end;
      it = _values.data();
      it_end = it + _values.num_elements();
      size_t count = 0;
      for( ; it != it_end; ++it ) ++count;
      */
#ifdef EXCESSIVE_DEBUG_PRINTING
      std::cout<< "Clear ("<< _values.num_elements()<< " items)" << std::endl;
#endif
      
      //_values.clear();
    }


  private:
    /*
     * helper functor: check consecutive items against limits in size_t[]
     */
    struct allInRange
    {
      typedef bool result_type;

      allInRange( const size_t* range ) : _range(range), _pos(0) {}
      template <typename T>
      bool operator()(const bool& state, const T& val) const { return state && (val <= _range[_pos++]); }
    private:
      const size_t* _range;
      mutable size_t _pos;
    };

  public:
    template <typename Address>
    bool check_address(Address& address) const
    {
      typedef typename boost::fusion::result_of::size<Args>::type size_of_address;
      BOOST_MPL_ASSERT_RELATION( size_of_address::value, ==, (size_t)container_t::dimensionality );

      return boost::fusion::accumulate( address, true, allInRange( _values.shape() ) );
    }

  public:
    const size_t* shape() const { return _values.shape(); }

  public:
    const Result& get_unassigned_value() const { return _unassigned; }

  } _storage;



public:
  void check_range(const Args& args ) const
  {
    /*
      if( boost::fusion::size(address) != container_t::dimensionality )
      {
      std::cout<< "Error: requrested address of wrong length ("<< boost::fusion::length(address)<< ")."<< std::endl;
      throw EInvalidElementRequested();
	}*/
    
    if( !_storage.check_address( args ) )
    {
      std::vector<size_t> t( boost::fusion::result_of::size<Args>::type::value );
      boost::fusion::for_each( args, setVectorItems(t) );
      
      std::vector<size_t>::const_iterator it, it_end;
      it = t.begin();
      it_end = t.end(); 
      std::cout<< "memoize for algo '"<< static_cast<const Impl*>(this)->get_debug_id()<<"': address outside range requested; requested address=[";
      for( ; it != it_end; ++it )
      {
	std::cout<< *it<< ", ";
      }
      std::cout<< "], extents=[";
      
      std::vector<size_t> ex( _storage.shape(), _storage.shape() + boost::fusion::result_of::size<dimensions_t>::type::value );
      it = ex.begin();
      it_end = ex.end(); 
      for( ; it != it_end; ++it )
      {
	std::cout<< *it<< ", ";
      }
      std::cout<< "]"<< std::endl;
      
      throw EInvalidElementRequested();
    }
  }
  
public:
  template<typename Comp>
  Result apply_impl(Comp& comp, const Args& args )
  {
    check_range(args);

    const Result& prev = _storage.access( args );
    //std::cout<< static_cast<Impl*>(this)->get_debug_id()<< " (prev="<< prev<< ")"<< std::endl;
    if( prev != _storage.get_unassigned_value() )
      return prev;

    //std::cout<< "memoize -> ";
    const Result& fresh = static_cast<Impl*>(this)->produce(comp, args);
    if( fresh == _storage.get_unassigned_value() ) // TODO consider changing to an assert
    {
      std::cout<< "ERROR: produce() (of algorithm '"<< static_cast<Impl*>(this)->get_debug_id() << "') returned 'unassigned' value!"<< std::endl;
      throw EInvalidValueProduced();
    }
    //std::cout<< static_cast<Impl*>(this)->get_debug_id()<< " (fresh="<< fresh<< ")"<< std::endl;
    _storage.access( args ) = fresh;
    assert( _storage.access( args ) == fresh );
    
    return fresh;
  }

public:
  template<typename Comp>
  void clear(Comp& comp)
  {
    const dimensions_t extents = static_cast<Impl*>(this)->get_extents(comp);
    
#ifdef EXCESSIVE_DEBUG_PRINTING
    {
      std::cout<< "clear got new extents: [";
      boost::fusion::for_each<>( extents, print_with_comma() );
      std::cout<< "]"<< std::endl;
    }
#endif


    _storage.allocate( extents ); // resize the storage
    _storage.clear();  // clear all values

  }


};


} // namespace v2


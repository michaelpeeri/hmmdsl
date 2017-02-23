#pragma once
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_push_back_actor.hpp>
#include <boost/spirit/include/classic_grammar_def.hpp>
#include <boost/spirit/include/classic_file_iterator.hpp> 
#include <boost/spirit/include/phoenix1_primitives.hpp>
#include <boost/spirit/include/phoenix1_functions.hpp>
#include <boost/lexical_cast.hpp>
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#undef BOOST_DISABLE_ASSERTS


namespace read_matrix {
	


typedef char                    char_t;
	
typedef boost::spirit::classic::file_iterator <char_t>  iterator_t;
	



namespace plain_matrix_grammar_dispatchers
{
	struct _disp_on_row_end
	{
		template<class T> struct result 	{			typedef void type;		};
		template<class T> void operator()(T& t) const    
		{
			t.on_row_end();
		}
	};
	phoenix::function<_disp_on_row_end> const on_row_end = _disp_on_row_end();

	struct _disp_on_value
	{
		//template<class T, class Iter1, class Iter2> struct result 	  {		  typedef void type;	  };
		template<class T, class Iter1, class Iter2> struct result 	  {		  typedef void type;	  };
		template<class T, class Iter1, class Iter2> void operator()(T& t, const Iter1& begin, const Iter2& end)  const
		{
			std::string str(begin,end);
			t.on_value(boost::lexical_cast<double>(str));
		}
	};
	phoenix::function<_disp_on_value> const on_value = _disp_on_value();

	template<class handler_t>
	struct ParsingAdaptor
	{
		ParsingAdaptor(handler_t& handler)
			: _handler(handler)
			, _col(0)
			, _row(0)
			, _num_cols(0)
			{
				// Clear the matrix before parsing
				_handler.resize( boost::extents[0][0] );
			}

		void on_row_end()
		{
			if( _row == 0 )
				_num_cols = _col;
			else
				assert(_col==_num_cols);

			//std::cout<< "*";
			
			_col = 0;
			++_row;
		}

		void on_value(double val)
		{
			if( _row == 0 )
			{
				_handler.resize( boost::extents[1][_col+1] );
				//std::cout<< "resizing to [1]["<< _col+1<< "]"<< std::endl;
			}
			else
			{
				if( _col == 0 )
				{
					_handler.resize( boost::extents[_row+1][_num_cols] );
					//std::cout<< "resizing to ["<< _row+1<< "]["<< _num_cols<< "]"<< std::endl;
				}
			}
			
			//std::cout<< "a["<< _row<< "]["<< _col<< "] <= "<< val<< std::endl;
			_handler[_row][_col] = val;

			++_col;
		}

	private:
		handler_t& _handler;
		size_t _col, _row, _num_cols;
	};
}
	

namespace spirit = boost::spirit::classic;

template <typename handler_t>
struct plain_matrix_grammar : public boost::spirit::classic::grammar<plain_matrix_grammar<handler_t> >
{
	plain_matrix_grammar(plain_matrix_grammar_dispatchers::ParsingAdaptor<handler_t>& adaptor) : _adaptor(adaptor) {}
public:
	template <typename ScannerT>
	struct definition
	{
		spirit::rule<ScannerT> row_line, comment_line, matrix_file, whitespace, elem, number;
		
		definition(plain_matrix_grammar const& self)
		{
			using namespace spirit;
			using namespace phoenix;
			using namespace plain_matrix_grammar_dispatchers;

			whitespace = +blank_p;

			number = +(digit_p | ch_p('.') | ch_p('-') | ch_p('+') | ch_p('e'));

			elem = number[on_value(var(self._adaptor), arg1, arg2)];
			
			row_line = !whitespace >> elem >> *(whitespace >> elem) >> eol_p[on_row_end(var(self._adaptor))];

			comment_line = ch_p('#') >> *(anychar_p - eol_p) >> eol_p;
			
			matrix_file
				= *(comment_line | row_line);
	  
		}
		spirit::rule<ScannerT> const& start() const {	return matrix_file;		}
	};
	
private:
	plain_matrix_grammar_dispatchers::ParsingAdaptor<handler_t>& _adaptor;
};

template<typename Scalar>
void read_matrix(const std::string& path, boost::multi_array<Scalar,2>& matrix);

template<typename Scalar>
Scalar test_read_matrix_sum(const std::string& path);
	
} // namespace read_matrix





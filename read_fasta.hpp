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


namespace fasta {
	


typedef char                    char_t;
	
typedef boost::spirit::classic::file_iterator <char_t>  iterator_t;
	



namespace fasta_grammar_dispatchers
{
	struct _disp_on_seq_begin
	{
		template<class T> struct result 	  {		  typedef void type;	  } ;
		template<class T> void operator()(T& t) const    
		{
			t.on_seq_begin();
		}
	};
	phoenix::function<_disp_on_seq_begin> const on_seq_begin = _disp_on_seq_begin();

	struct _disp_on_seq_end
	{
		template<class T> struct result 	{			typedef void type;		};
		template<class T> void operator()(T& t) const    
		{
			t.on_seq_end();
		}
		
	};
	phoenix::function<_disp_on_seq_end> const on_seq_end = _disp_on_seq_end();

	struct _disp_on_accession
	{
		template<class T, class Iter1, class Iter2> struct result 	  {		  typedef void type;	  };
		template<class T, class Iter1, class Iter2> void operator()(T& t, const Iter1& begin, const Iter2& end)  const    
		{
			t.on_accession(begin, end);
		}
	};
	phoenix::function<_disp_on_accession> const on_accession = _disp_on_accession();

	struct _disp_on_description
	{
		template<class T, class Iter1, class Iter2> struct result 		{			typedef void type;		};
		template<class T, class Iter1, class Iter2> void operator()(T& t, const Iter1& begin, const Iter2& end)  const    
		{
			t.on_description(begin, end);
		}
	};
	phoenix::function<_disp_on_description> const on_description = _disp_on_description();

	struct _disp_on_sequence_line
	{
		template<class T, class Iter1, class Iter2> struct result 	  {		  typedef void type;	  };
		template<class T, class Iter1, class Iter2> void operator()(T& t, const Iter1& begin, const Iter2& end)  const    
		{
			t.on_sequence_line(begin, end);
		}
	};
	phoenix::function<_disp_on_sequence_line> const on_sequence_line = _disp_on_sequence_line();

	template<class handler_t>
	struct ParsingAdaptor
	{
		ParsingAdaptor(handler_t& handler) : _handler(handler) 	{}
		void on_seq_begin()
		{
			_seq.clear();
			_desc.clear();
			_acc.clear();
			_handler.on_seq_begin();
		}

		void on_seq_end()
		{
			_handler.on_seq( _acc, _desc, _seq );
			_handler.on_seq_end();
		}

		template<class Iter1, class Iter2>
		void on_accession(const Iter1& begin, const Iter2& end)
		{
			_acc.append(begin, end);
			_handler.on_accession( _acc );
		}

		template<class Iter1, class Iter2>
		void on_description(const Iter1& begin, const Iter2& end)
		{
			_desc.append(begin, end);
			_handler.on_description( _desc );
		}

		template<class Iter1, class Iter2>
		void on_sequence_line(const Iter1& begin, const Iter2& end)
		{
			_seq.append(begin, end);
		}
	private:
		handler_t& _handler;
		std::string _acc, _desc, _seq;
	};
}
	

namespace spirit = boost::spirit::classic;

template <typename handler_t>
struct fasta_grammar : public boost::spirit::classic::grammar<fasta_grammar<handler_t> >
{
	fasta_grammar(fasta_grammar_dispatchers::ParsingAdaptor<handler_t>& adaptor) : _adaptor(adaptor) {}
public:
	template <typename ScannerT>
	struct definition
	{
		spirit::rule<ScannerT>header, accession, description, seq_line, fasta_file, alphabet, seq;
		
		definition(fasta_grammar const& self)
		{
			using namespace spirit;
			using namespace phoenix;
			using namespace fasta_grammar_dispatchers;
			alphabet
				=  range_p('A','Z') | range_p('a','z') | ch_p('*');
			accession
				=  *graph_p;
			description
				=  *(anychar_p - eol_p);
			header
				=  ch_p('>')
				>> epsilon_p[on_seq_begin(var(self._adaptor))]
				>> accession[on_accession(var(self._adaptor), arg1, arg2)]
				>>
				!(blank_p
				  >> description[on_description(var(self._adaptor), arg1, arg2)])
				>> eol_p;
			seq_line
				=  (+alphabet)[on_sequence_line(var(self._adaptor), arg1, arg2)]
				>> eol_p;
			seq
				=  header
				>> *seq_line
				>> epsilon_p[on_seq_end(var(self._adaptor))];

			fasta_file
				= *(seq);
	  
		}
		spirit::rule<ScannerT> const& start() const {	return fasta_file;		}
	};
	
private:
	fasta_grammar_dispatchers::ParsingAdaptor<handler_t>& _adaptor;
};
		
template<class Seq, class Acc, class Desc, class Cont=std::vector<boost::tuple<Acc,Desc,Seq> > >
struct SequenceContainer : public Cont
{
public:
	SequenceContainer() : _cummulative_seq_size(0)  {}
	void on_seq_begin() {}
	
	void on_seq_end() {}
	
	void on_accession(const Acc&) {}
	
	void on_description(const Desc&)
	{
	}
	
	void on_seq(const Acc& acc, const Desc& desc, const Seq& seq)
	{
		push_back( boost::make_tuple(acc,desc,seq) );

/*		
		std::cout<< "***"<< std::endl;
		
		typename Cont::const_iterator it=this->begin();
		typename Cont::const_iterator it_end=this->end();
		for(; it != it_end; ++it ) {
			
			std::cout<< get_seq(*it)<< std::endl;
			std::cout<< get_desc(*it)<< std::endl;
		}
		
*/		

		_cummulative_seq_size += seq.size();
	}
	size_t cummulative_seq_size() const 
	{
		return _cummulative_seq_size;
	}
private:
	size_t _cummulative_seq_size;
	
};

typedef std::string seq_t;
typedef SequenceContainer<seq_t, std::string, std::string> seq_cont_t;

const std::string& get_acc(seq_cont_t::const_reference elem);
const std::string& get_desc(seq_cont_t::const_reference elem);
const seq_t& get_seq(seq_cont_t::const_reference elem);
size_t get_len(seq_cont_t::const_reference elem);


size_t count_fasta(const std::string& path);
void read_fasta(const std::string& path, seq_cont_t& seqs);

template<class Iterator>
void read_fasta(const std::string& path, seq_cont_t& seqs, Iterator begin, Iterator end);
	
} // namespace fasta


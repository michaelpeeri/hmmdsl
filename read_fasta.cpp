#include "read_fasta.hpp"

namespace fasta
{
	
const std::string& get_acc(seq_cont_t::const_reference elem)
{
	return elem.get<0>();
}
	
const std::string& get_desc(seq_cont_t::const_reference elem)
{
	return elem.get<1>();
}

const seq_t& get_seq(seq_cont_t::const_reference elem)
{
		return elem.get<2>();
}

size_t get_len(seq_cont_t::const_reference elem)
{
	return elem.get<2>().size();
}


void read_fasta(const std::string& path, seq_cont_t& seqs)
{
	using namespace spirit;
	seqs.clear();
	fasta_grammar_dispatchers::ParsingAdaptor<seq_cont_t> trans(seqs);
	fasta_grammar<seq_cont_t> g(trans);
	
	iterator_t first(path);
	
	iterator_t last = first.make_end();
	
	parse_info<iterator_t> info = parse(first, last, g);
	
	if( (!info.full) || (info.length == 0) )
	{
		std::cerr<< "Error: Input file "<< path<< " does not appear to be a valid fasta file!"<< std::endl;
		assert(info.full);
	}
	std::cout<< "#Read "<< seqs.cummulative_seq_size()<< " letters in "<< seqs.size()<< " sequences"<< std::endl;
}

size_t count_fasta(const std::string& path)
{
	// TODO - Optimize this!
	using namespace spirit;
	seq_cont_t seqs;

	fasta_grammar_dispatchers::ParsingAdaptor<seq_cont_t> trans(seqs);
	fasta_grammar<seq_cont_t> g(trans);
	
	iterator_t first(path);
	
	iterator_t last = first.make_end();
	
	parse_info<iterator_t> info = parse(first, last, g);
	
	if( (!info.full) || (info.length == 0) )
	{
		std::cerr<< "Error: Input file "<< path<< " does not appear to be a valid fasta file!"<< std::endl;
		assert(info.full);
	}
	return seqs.size();
}

template<class Iterator>
void read_fasta(const std::string& path, seq_cont_t& seqs, Iterator begin, Iterator end)
{
	using namespace spirit;
	seq_cont_t seqs_orig;
	seqs.clear();
	fasta_grammar_dispatchers::ParsingAdaptor<seq_cont_t> trans(seqs_orig);
	fasta_grammar<seq_cont_t> g(trans);
	
	iterator_t first(path);
	
	iterator_t last = first.make_end();
	
	parse_info<iterator_t> info = parse(first, last, g);
	
	if( (!info.full) || (info.length == 0) )
	{
		std::cerr<< "Error: Input file "<< path<< " does not appear to be a valid fasta file!"<< std::endl;
		assert(info.full);
	}
	//std::cout<< "#Read "<< seqs.cummulative_seq_size()<< " letters in "<< seqs.size()<< " sequences"<< std::endl;

	for( Iterator it=begin; it != end; ++it )
	{
		seqs.push_back( seqs_orig[*it] );
	}
	// If the index-list is empty, use all sequences
	if( seqs.size() == 0)
		seqs = seqs_orig;
	
}

// Explicit instantiation
template
void read_fasta<std::vector<size_t>::iterator>(const std::string& path, seq_cont_t& seqs, std::vector<size_t>::iterator begin, std::vector<size_t>::iterator end);
	
} // namespace fasta

#include "read_matrix.hpp"

namespace read_matrix
{
	
template<typename Scalar>
void read_matrix(const std::string& path, boost::multi_array<Scalar,2>& matrix)
{
	using namespace spirit;
	typedef boost::multi_array<Scalar,2> matrix_t;

	plain_matrix_grammar_dispatchers::ParsingAdaptor<matrix_t> trans(matrix);
	plain_matrix_grammar<matrix_t> g(trans);
	
	iterator_t first(path);
	
	iterator_t last = first.make_end();
	
	parse_info<iterator_t> info = parse(first, last, g);
	
	if( (!info.full) || (info.length == 0) )
	{
		std::cerr<< "Error: Input file "<< path<< " does not appear to be a valid fasta file!"<< std::endl;
		assert(info.full);
	}
}
	

template<typename Scalar>
Scalar test_read_matrix_sum(const std::string& path)
{
	using namespace spirit;
	typedef boost::multi_array<Scalar,2> matrix_t;
	
	matrix_t matrix;
	read_matrix(path, matrix);

	Scalar sigma = 0.0;
	
	const size_t m_end = matrix.shape()[0];
	const size_t n_end = matrix.shape()[1];
	
	for( size_t m=0; m< m_end; ++m )
		for( size_t n=0; n< n_end; ++n )
		{	 sigma += matrix[m][n]; }
	

	return sigma;
}
	
// Explicit instantiation
template
double test_read_matrix_sum<double>(const std::string& path);

} // namespace read_matrix

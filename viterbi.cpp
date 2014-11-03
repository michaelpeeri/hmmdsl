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
#include "viterbi.hpp"
#include <boost/format.hpp>

namespace v2
{

  const char* AlignmentLineFormat = "%1%%|13t|%2%%|18t|%3% %4%";

std::ostream& operator<< (std::ostream& stream, const v2::decoding_t& decoding)
{
  using boost::format;
  using boost::str;
#ifdef EXCESSIVE_DEBUG_PRINTING
  stream<< "decoding_t:";
#endif

  if( decoding._uninitialized )
  {
    stream<< "decoding_t: (Null)"<< std::endl;
    return stream;
  }

  const size_t totalLength = decoding._s1.length();
  const size_t width = 50;

  for( size_t pos = 0; pos <totalLength; pos += width )
  {
    const size_t segmentLength = (pos + width <= totalLength) ? (width) : (totalLength-pos);
    stream<< "\n";
    //stream<< decoding._id1.substr(0, 8)<< " "<< pos+1<< " "<< decoding._s1.substr(pos, segmentLength)<< " "<< pos+segmentLength<< std::endl;
    stream<< str( boost::format(AlignmentLineFormat) % decoding._id1.substr(0,12) % (pos+1) % decoding._s1.substr(pos,segmentLength) % (pos+segmentLength) );
    stream<< std::endl;
    //stream<< decoding._id2.substr(0, 8)<< " "<< pos+1<< " "<< decoding._s2.substr(pos, segmentLength)<< " "<< pos+segmentLength<< std::endl;
    stream<< str( boost::format(AlignmentLineFormat) % decoding._id2.substr(0,12) % (pos+1) % decoding._s2.substr(pos,segmentLength) % (pos+segmentLength) );
    stream<< std::endl;
  }
  
  stream<< std::flush;
  return stream;
}


} // namespace v2

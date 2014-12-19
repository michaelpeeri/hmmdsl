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
#include <boost/mpl/int.hpp>

namespace v2
{


namespace sequential_processing
{
  namespace detail
  {
    struct info { typedef boost::mpl::int_<1> max_level; };
    struct info_for_1 { typedef boost::mpl::int_<0> max_level; };
  } // namespace detail

  struct sequential_processing_level {};

  struct empty_scope           : public sequential_processing_level  {};
  struct single_sequence_scope : public sequential_processing_level  {  typedef boost::mpl::int_<0> value; typedef detail::info info; typedef empty_scope prev; };
  struct multi_sequence_scope  : public sequential_processing_level  {  typedef boost::mpl::int_<1> value; typedef detail::info info; typedef single_sequence_scope prev; };
  struct everything            : public sequential_processing_level  {  typedef boost::mpl::int_<2> value; typedef detail::info info; typedef multi_sequence_scope prev; };


  struct sole_sequence_scope : public sequential_processing_level  {  typedef boost::mpl::int_<0> value; typedef detail::info_for_1 info; typedef empty_scope prev; };


template <class Level>
struct get_previous_level
{
  typedef typename Level::prev type;
};

template<>
struct get_previous_level<empty_scope>
{ 
  typedef empty_scope type;
};

} //namespace sequential_processing

} // namespace v2

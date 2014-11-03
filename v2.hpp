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
#include <string>

namespace v2
{
typedef double probability_t;
typedef size_t time_t;
typedef size_t state_t;
typedef std::string sequence_t;

namespace model
{

/*
 * Provide abstract information about a model, to allow specialization without acquaintance with the actual implementation classes
 */
template<typename Model, typename Enable=void>
struct model_traits;  // Base template (must be specialized)

/*
 * Tags for use in model_traits
 */
struct markov_chain {};       // Hidden process is a Markov chain
struct semi_markov_chain {};  // Hidden process is a semi-Markov chain


} // namespace model

struct print_with_comma
{   template <typename Item> void operator()(const Item& i) const { std::cout<< i<< ", "; } };

}  // namespace v2

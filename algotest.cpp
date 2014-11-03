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
#include <iostream>
#include <boost/tuple/tuple.hpp>

namespace hmmdsl {
/*
  Algo - accumulator
     I calculate something based on the model and data.
     I may depend on other algos.
     I cache the results of the calculation internally.
     I perform the calculation lazily.
     
  AlgoGroup - accumulator set
  Model
  Data
  
 */


template <typename T>
class Model
{
};

struct markov_dependency {};
struct semi_markov_dependency {};

template <typename T>
class HMM : public Model<T>
{
  typedef markov_dependency state_dependecy;
};

template <typename T>
class HSMM : public Model<T>
{
  typedef semi_markov_dependency state_dependecy;
};


template <typename ValueType, typename Args, typename Model>
class Algo
{
  // "Pure abstract" - Must be specialized
public:
  typedef Args args_t;
};

typedef double probability_t;
typedef size_t time_t;
typedef size_t state_t;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename ValueType, typename Args, size_t Level>
class MemoizedState
{
public:
  ValueType calc(const Args& args)
  {
  }

public:
  void reset(size_t level=0) {};
};


template<typename Model>
class delta : public Algo<boost::tuple<probability_t,size_t>, boost::tuple<time_t,state_t>,Model>, public MemoizedState< boost::tuple<probability_t,size_t>, boost::tuple<time_t,state_t>, 1 >
{
public:
  typedef Algo<probability_t,boost::tuple<time_t,state_t>,Model> base;
  
  boost::tuple<probability_t,size_t> calc(const typename base::args_t& args)
  {
    std::cout<< "BasicAlgo<HMM<T> >::calc()"<< std::endl;
    return boost::make_tuple(1.0,2);
  }
  void reset(size_t level=0) {};
};



}

int main()
{

  typedef hmmdsl::delta<hmmdsl::HMM<double> > vab1_t;
  typedef hmmdsl::delta<hmmdsl::HSMM<double> > vab2_t;  
  
  vab1_t vab1;
  vab2_t vab2;


  vab1.calc(vab1_t::args_t(1,2));
  vab2.calc(vab1_t::args_t(1,2));

  //std::cout<< vab1.calc(vab1_t::args_t(1,2)) << std::endl;
  //std::cout<< vab2.calc(vab1_t::args_t(1,2)) << std::endl;

  return 0;
}

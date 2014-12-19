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

#include "v2.hpp"
#include "computation.hpp"

namespace v2
{


// TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  //

namespace tests
{
// a -> {b}
// b -> {c, d}
// d -> {e}
//
// j -> {j}
//
// m -> {m, n}
// n -> {n, m}
//
// o -> {z}
// p -> {q, o}
// q -> {r, o}
// r -> {p, o}
//
// f1 -> {f1, f2}
// f2 -> {f2, f1}
// b1 -> {b1, b2}
// b2 -> {b2, b1}
// xi -> {f1, b2}
// sxi1 -> {sxi1, xi}
// sxi2 -> {sxi2, xi}

namespace tag {
    struct algo_a : public ::v2::tag::algorithm {};
    struct algo_b : public ::v2::tag::algorithm {};
    struct algo_c : public ::v2::tag::algorithm {};
    struct algo_d : public ::v2::tag::algorithm {};
    struct algo_e : public ::v2::tag::algorithm {};
    struct algo_j : public ::v2::tag::algorithm {};
    struct algo_m : public ::v2::tag::algorithm {};
    struct algo_n : public ::v2::tag::algorithm {};

    struct algo_o : public ::v2::tag::algorithm {};
    struct algo_p : public ::v2::tag::algorithm {};
    struct algo_q : public ::v2::tag::algorithm {};
    struct algo_r : public ::v2::tag::algorithm {};
    struct algo_z : public ::v2::tag::algorithm {};

    struct test_f1   : public ::v2::tag::algorithm {};
    struct test_f2   : public ::v2::tag::algorithm {};
    struct test_b1   : public ::v2::tag::algorithm {};
    struct test_b2   : public ::v2::tag::algorithm {};
    struct test_xi   : public ::v2::tag::algorithm {};
    struct test_sxi1 : public ::v2::tag::algorithm {};
    struct test_sxi2 : public ::v2::tag::algorithm {};
} // namespace tag

struct algo_a : 
  public depends_on<
    boost::mpl::set<tag::algo_b>  // a depends on b
    , sequential_processing::single_sequence_scope
    , tag::algo_a
    >
{};

struct algo_b : 
  public depends_on<
    boost::mpl::set<tag::algo_c, tag::algo_d> // b depends on c, d
    , sequential_processing::single_sequence_scope
    , tag::algo_b
    >
{};

struct algo_c : 
  public depends_on<
      boost::mpl::set<>
    , sequential_processing::single_sequence_scope
    , tag::algo_c
    >
{};

struct algo_d : 
  public depends_on<
      boost::mpl::set<tag::algo_e>
    , sequential_processing::single_sequence_scope
    , tag::algo_d
    >
{};

struct algo_e : 
  public depends_on<
      boost::mpl::set<>
    , sequential_processing::single_sequence_scope
    , tag::algo_e
    >
{};

struct algo_j : 
  public depends_on<
      boost::mpl::set<tag::algo_j>
    , sequential_processing::single_sequence_scope
    , tag::algo_j
    >
{};

struct algo_m : 
  public depends_on<
      boost::mpl::set<tag::algo_n>
    , sequential_processing::single_sequence_scope
    , tag::algo_m
    >
{};
struct algo_n : 
  public depends_on<
      boost::mpl::set<tag::algo_m>
    , sequential_processing::single_sequence_scope
    , tag::algo_n
    >
{};

struct algo_o : 
  public depends_on<
      boost::mpl::set<tag::algo_z>
    , sequential_processing::single_sequence_scope
    , tag::algo_o
    >
{};
struct algo_p : 
  public depends_on<
      boost::mpl::set<tag::algo_q, tag::algo_o>
    , sequential_processing::single_sequence_scope
    , tag::algo_p
    >
{};
struct algo_q : 
  public depends_on<
      boost::mpl::set<tag::algo_r, tag::algo_o>
    , sequential_processing::single_sequence_scope
    , tag::algo_q
    >
{};
struct algo_r : 
  public depends_on<
      boost::mpl::set<tag::algo_p, tag::algo_o, tag::algo_r>
    , sequential_processing::single_sequence_scope
    , tag::algo_r
    >
{};
struct algo_z : 
  public depends_on<
    boost::mpl::set<tag::algo_z>
    , sequential_processing::single_sequence_scope
    , tag::algo_z
    >
{};

struct algo_test_f1 : 
  public depends_on<
      boost::mpl::set<tag::test_f1, tag::test_f2>
    , sequential_processing::single_sequence_scope
    , tag::test_f1
    >
{};

struct algo_test_f2 : 
  public depends_on<
      boost::mpl::set<tag::test_f2, tag::test_f1>
    , sequential_processing::single_sequence_scope
    , tag::test_f2
    >
{};
struct algo_test_b1 : 
  public depends_on<
      boost::mpl::set<tag::test_b2, tag::test_b1>
    , sequential_processing::single_sequence_scope
    , tag::test_b1
    >
{};
struct algo_test_b2 : 
  public depends_on<
      boost::mpl::set<tag::test_b1, tag::test_b2>
    , sequential_processing::single_sequence_scope
    , tag::test_b2
    >
{};
struct algo_test_xi : 
  public depends_on<
      boost::mpl::set<tag::test_f1, tag::test_b2>
    , sequential_processing::single_sequence_scope
    , tag::test_xi
    >
{};
struct algo_test_sxi1 : 
  public depends_on<
      boost::mpl::set<tag::test_sxi1, tag::test_xi>
    , sequential_processing::single_sequence_scope
    , tag::test_sxi1
    >
{};
struct algo_test_sxi2 : 
  public depends_on<
      boost::mpl::set<tag::test_sxi2, tag::test_xi>
    , sequential_processing::single_sequence_scope
    , tag::test_sxi2
    >
{};


struct dummy_model {};

} // namespace tests


// Implementations must be defined in namespace v2 (not in contained namespaces)
template <>
struct implementations<tests::tag::algo_a, tests::dummy_model >
{
    typedef tests::algo_a impl_type;
};
template <>
struct implementations<tests::tag::algo_b, tests::dummy_model >
{
    typedef tests::algo_b impl_type;
};
template <>
struct implementations<tests::tag::algo_c, tests::dummy_model >
{
    typedef tests::algo_c impl_type;
};
template <>
struct implementations<tests::tag::algo_d, tests::dummy_model >
{
    typedef tests::algo_d impl_type;
};
template <>
struct implementations<tests::tag::algo_e, tests::dummy_model >
{
    typedef tests::algo_e impl_type;
};
template <>
struct implementations<tests::tag::algo_j, tests::dummy_model >
{
    typedef tests::algo_j impl_type;
};

template <>
struct implementations<tests::tag::algo_m, tests::dummy_model >
{
    typedef tests::algo_m impl_type;
};
template <>
struct implementations<tests::tag::algo_n, tests::dummy_model >
{
    typedef tests::algo_n impl_type;
};

template <>
struct implementations<tests::tag::algo_o, tests::dummy_model >
{
    typedef tests::algo_o impl_type;
};
template <>
struct implementations<tests::tag::algo_p, tests::dummy_model >
{
    typedef tests::algo_p impl_type;
};
template <>
struct implementations<tests::tag::algo_q, tests::dummy_model >
{
    typedef tests::algo_q impl_type;
};
template <>
struct implementations<tests::tag::algo_r, tests::dummy_model >
{
    typedef tests::algo_r impl_type;
};
template <>
struct implementations<tests::tag::algo_z, tests::dummy_model >
{
    typedef tests::algo_z impl_type;
};


template <>
struct implementations<tests::tag::test_f1, tests::dummy_model >
{
    typedef tests::algo_test_f1 impl_type;
};
template <>
struct implementations<tests::tag::test_f2, tests::dummy_model >
{
    typedef tests::algo_test_f2 impl_type;
};
template <>
struct implementations<tests::tag::test_b1, tests::dummy_model >
{
    typedef tests::algo_test_b1 impl_type;
};
template <>
struct implementations<tests::tag::test_b2, tests::dummy_model >
{
    typedef tests::algo_test_b2 impl_type;
};
template <>
struct implementations<tests::tag::test_xi, tests::dummy_model >
{
    typedef tests::algo_test_xi impl_type;
};
template <>
struct implementations<tests::tag::test_sxi1, tests::dummy_model >
{
    typedef tests::algo_test_sxi1 impl_type;
};
template <>
struct implementations<tests::tag::test_sxi2, tests::dummy_model >
{
    typedef tests::algo_test_sxi2 impl_type;
};


namespace tests
{

namespace tests_computation_hpp
{
/*
 * MPL unary meta-function object facade for 'dependencies_of<>'
 */
struct direct_dependencies_set
{
    template< typename Tag>
    struct apply
    {
        typedef typename detail::dependencies_of<Tag,tests::dummy_model>::type type;
    };
};
} // namespace tests_computation_hpp

// a -> {b}
// b -> {c, d}
// d -> {e}
//
// j -> {j}
//
// m -> {m, n}
// n -> {n, m}
//
// o -> {z}
// p -> {q, o}
// q -> {r, o}
// r -> {p, o}
//
// f1 -> {f1, f2}
// f2 -> {f2, f1}
// b1 -> {b1, b2}
// b2 -> {b2, b1}
// xi -> {f1, b2}
// sxi1 -> {sxi1, xi}
// sxi2 -> {sxi2, xi}

//
// New lines needed because of silly bug in BOOST_MPL_ASSERT...
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//

typedef boost::mpl::set<tag::algo_e> set771;
typedef typename detail::closure<set771, tests_computation_hpp::direct_dependencies_set>::type closure_of_set771;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set771, tag::algo_e> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set771>::value, ==, 1 );

typedef boost::mpl::set<tag::algo_d> set772;
typedef typename detail::closure<set772, tests_computation_hpp::direct_dependencies_set>::type closure_of_set772;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set772, tag::algo_d> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set772, tag::algo_e> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set772>::value, ==, 2 );

typedef boost::mpl::set<tag::algo_c> set773;
typedef typename detail::closure<set773, tests_computation_hpp::direct_dependencies_set>::type closure_of_set773;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set773, tag::algo_c> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set773>::value, ==, 1 );

typedef boost::mpl::set<tag::algo_b> set774;
typedef typename detail::closure<set774, tests_computation_hpp::direct_dependencies_set>::type closure_of_set774;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set774, tag::algo_b> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set774, tag::algo_c> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set774, tag::algo_d> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set774, tag::algo_e> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set774>::value, ==, 4 );

typedef boost::mpl::set<tag::algo_a> set775;
typedef typename detail::closure<set775, tests_computation_hpp::direct_dependencies_set>::type closure_of_set775;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set775, tag::algo_a> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set775, tag::algo_b> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set775, tag::algo_c> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set775, tag::algo_d> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set775, tag::algo_e> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set775>::value, ==, 5 );

typedef boost::mpl::set<tag::algo_j> set776;
typedef typename detail::closure<set776, tests_computation_hpp::direct_dependencies_set>::type closure_of_set776;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set776, tag::algo_j> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set776>::value, ==, 1 );


typedef boost::mpl::set<tag::algo_m> set777;
typedef typename detail::closure<set777, tests_computation_hpp::direct_dependencies_set>::type closure_of_set777;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set777, tag::algo_m> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set777, tag::algo_n> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set777>::value, ==, 2 );

typedef boost::mpl::set<tag::algo_n> set778;
typedef typename detail::closure<set778, tests_computation_hpp::direct_dependencies_set>::type closure_of_set778;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set778, tag::algo_n> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set778, tag::algo_m> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set778>::value, ==, 2 );

typedef boost::mpl::set<tag::algo_m, tag::algo_n> set779;
typedef typename detail::closure<set779, tests_computation_hpp::direct_dependencies_set>::type closure_of_set779;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set779, tag::algo_m> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set779, tag::algo_n> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set779>::value, ==, 2 );

typedef boost::mpl::set<tag::algo_p> set780;
typedef typename detail::closure<set780, tests_computation_hpp::direct_dependencies_set>::type closure_of_set780;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set780, tag::algo_o> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set780, tag::algo_p> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set780, tag::algo_q> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set780, tag::algo_r> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set780, tag::algo_z> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set780>::value, ==, 5 );

typedef boost::mpl::set<tag::test_sxi1> set790;
typedef typename detail::closure<set790, tests_computation_hpp::direct_dependencies_set>::type closure_of_set790;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set790, tag::test_sxi1> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set790, tag::test_xi> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set790, tag::test_f1> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set790, tag::test_f2> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set790, tag::test_b1> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set790, tag::test_b2> ));
BOOST_MPL_ASSERT(( boost::mpl::not_<boost::mpl::has_key<closure_of_set790, tag::algo_z>> ));
BOOST_MPL_ASSERT(( boost::mpl::not_<boost::mpl::has_key<closure_of_set790, tag::test_sxi2>> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set790>::value, ==, 6 );

typedef boost::mpl::set<tag::test_sxi2> set791;
typedef typename detail::closure<set791, tests_computation_hpp::direct_dependencies_set>::type closure_of_set791;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set791, tag::test_sxi2> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set791, tag::test_xi> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set791, tag::test_f1> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set791, tag::test_f2> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set791, tag::test_b1> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set791, tag::test_b2> ));
BOOST_MPL_ASSERT(( boost::mpl::not_<boost::mpl::has_key<closure_of_set791, tag::algo_z>> ));
BOOST_MPL_ASSERT(( boost::mpl::not_<boost::mpl::has_key<closure_of_set791, tag::test_sxi1>> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set791>::value, ==, 6 );


typedef boost::mpl::set<tag::test_xi, tag::test_sxi1> set792;
typedef typename detail::closure<set792, tests_computation_hpp::direct_dependencies_set>::type closure_of_set792;
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set792, tag::test_sxi1> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set792, tag::test_xi> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set792, tag::test_f1> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set792, tag::test_f2> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set792, tag::test_b1> ));
BOOST_MPL_ASSERT(( boost::mpl::has_key<closure_of_set792, tag::test_b2> ));
BOOST_MPL_ASSERT(( boost::mpl::not_<boost::mpl::has_key<closure_of_set792, tag::algo_z>> ));
BOOST_MPL_ASSERT(( boost::mpl::not_<boost::mpl::has_key<closure_of_set792, tag::test_sxi2>> ));
BOOST_MPL_ASSERT_RELATION( boost::mpl::size<closure_of_set792>::value, ==, 6 );


// TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  // TESTING ONLY  //
} // namespace tests

} // namespace v2

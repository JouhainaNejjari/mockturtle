/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file simulation_cec.hpp
  \brief Simulation-based CEC
  EPFL CS-472 2021 Final Project Option 2
*/

#pragma once

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include "../utils/node_map.hpp"
#include "miter.hpp"
#include "simulation.hpp"
#include <cmath>

namespace mockturtle
{

/* Statistics to be reported */
struct simulation_cec_stats
{
  /*! \brief Split variable (simulation size). */
  uint32_t split_var{ 0 };

  /*! \brief Number of simulation rounds. */
  uint32_t rounds{ 0 };
};

namespace detail
{

/*! \brief Simulates truth tables in rounds with variable limit.
 *
 * This simulator simulates truth tables.  Each primary input is assigned the
 * projection function according to the index.  The number of variables be
 * passed to the constructor of the simulator.
 */

template<class Ntk>
class simulation_cec_impl
{
public:
  using pattern_t = unordered_node_map<kitty::dynamic_truth_table, Ntk>;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit simulation_cec_impl( Ntk& ntk, simulation_cec_stats& st )
      : _ntk( ntk ),
        _st( st )
  {
  }

  bool run()
  {
    /* TODO: write your implementation here */

    // Computing split_var and rounds
    uint32_t n = _ntk.num_pis();
    uint32_t split_var = compute_splitting_var(n);
    uint32_t rounds = compute_rounds(n, split_var);

    // Reporting split_var and rounds to the user as statistics
    _st.split_var = split_var;
    _st.rounds = rounds;

    // Initializing pattern
    pattern_t patterns = init_patterns(_ntk, split_var);

    // Simulating first round
    default_simulator<kitty::dynamic_truth_table> sim(split_var);
    simulate_nodes(_ntk, patterns, sim);

    // Checking patterns
    if (!check( _ntk,patterns )){
      return(false);
    }

    // Looping over simulation rounds
    for (uint32_t round = 1; round < rounds; round++) {
      update_pattern(round, patterns,split_var);
      simulate_nodes(_ntk, patterns, sim);

      // Pattern checking
      if (!check( _ntk,patterns )){
          return(false);
          }
    }

    return true;
  }

private:
  /* you can add additional methods here */
  uint32_t compute_splitting_var(uint32_t n) {
      uint32_t split_var;
      if (n<=6)
        split_var=n;
      else {
          uint32_t v = _ntk.size();
          uint32_t m = 7;
          uint32_t cond = ( 32 + (1 << (m-3)) )*v;
          uint32_t max = 1 << 29;
          while (m<=n && cond<=max)
            m++;
          split_var = m;
      }
      return split_var;
  }

  uint32_t compute_rounds(uint32_t n, uint32_t split_var) {
      return 1 << (n-split_var);
  }

  pattern_t init_patterns( Ntk& _N, uint32_t split_var){
      pattern_t patterns(_N);

      _N.foreach_pi( [&]( auto const& n, auto p ){
       kitty::dynamic_truth_table tt (split_var);
       if (p < split_var) {
       kitty:: create_nth_var(tt , p);
       }
       patterns[n]=tt;
    } );
    return patterns;
  }

  bool check (Ntk& _N ,pattern_t& patterns){
    bool eq = true;
    _N.foreach_po( [&]( auto const& f) {
      if ( _N.is_complemented( f ) )
      {
        if ( !is_const0(~patterns[f]) ) {
          eq = false;
        }
      }
      else
      {
        if ( !is_const0(patterns[f]) ) {
        eq = false;
        }
      }
    } );
    return eq;
  }

  /*the function to update the pattern*/
  void update_pattern( uint32_t round , pattern_t& patterns, uint32_t split_var ){
    // Cleaning old patterns
    _ntk.foreach_gate( [&]( auto const& m )
    {
       patterns.erase(m);
    } );

    uint32_t r = round;
    // Updating patterns
      _ntk.foreach_pi( [&]( auto const& n, auto i )
      {
        // Splitting variables for indices after split_var
        if (i >= split_var ){
          // Case where round is odd
          if (r % 2 == 1) {
            // Updating patterns
            if ( is_const0(patterns[n]) ) patterns[n] = ~patterns[n];
          }
        // Case where round is even
        else {
          // Updating patterns
          if ( !is_const0(patterns[n]) ) patterns[n] = ~patterns[n];
        }
        r /= 2;

        }
      } );
  }

private:
  Ntk& _ntk;
  simulation_cec_stats& _st;
  /* you can add other attributes here */
};

} // namespace detail

/* Entry point for users to call */

/*! \brief Simulation-based CEC.
 *
 * This function implements a simulation-based combinational equivalence checker.
 * The implementation creates a miter network and run several rounds of simulation
 * to verify the functional equivalence. For memory and speed reasons this approach
 * is limited up to 40 input networks. It returns an optional which is `nullopt`,
 * if the network has more than 40 inputs.
 */
template<class Ntk>
std::optional<bool> simulation_cec( Ntk const& ntk1, Ntk const& ntk2, simulation_cec_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_num_pis_v<Ntk>, "Ntk does not implement the num_pis method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );

  simulation_cec_stats st;

  bool result = false;

  if ( ntk1.num_pis() > 40 )
    return std::nullopt;

  auto ntk_miter = miter<Ntk>( ntk1, ntk2 );

  if ( ntk_miter.has_value() )
  {
    detail::simulation_cec_impl p( *ntk_miter, st );
    result = p.run();
  }

  if ( pst )
    *pst = st;

  return result;
}

} // namespace mockturtle

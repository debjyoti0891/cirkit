/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2017  EPFL
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

/**
 * @file lhrs.hpp
 *
 * @brief LUT-based hierarchical reversible synthesis
 *
 * @author Mathias Soeken
 * @since  2.3
 */

#ifndef CLI_LHRS_COMMAND_HPP
#define CLI_LHRS_COMMAND_HPP

#include <memory>

#include <cli/aig_command.hpp>
#include <reversible/synthesis/lhrs/legacy/lhrs_params.hpp>

namespace cirkit
{

class lhrs_command : public aig_base_command
{
public:
  lhrs_command( const environment::ptr& env );

protected:
  rules validity_rules() const;
  void execute();

public:
  nlohmann::json log() const;

private:
  legacy::lhrs_params params;
  std::shared_ptr<legacy::lhrs_stats> stats;

  unsigned cut_size = 16u;
  unsigned lut_count = 0u;
  unsigned area_iters_init = 2u;
  unsigned flow_iters_init = 1u;

  std::string dotname_mapped;

  /* to be merged */
  std::string mapping_strategy;
  std::string assignment_strategy;
  std::string cover_method;
  std::string script;

private:
  unsigned debug_lb = 0;
};

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:

#ifndef __MX_PARMS_H__
#define __MX_PARMS_H__

//@HEADER
// ************************************************************************
// 
//               miniXyce: A simple circuit simulation benchmark code
//                 Copyright (2011) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

// Author : Karthik V Aadithya
// Mentor : Heidi K Thornquist
// Date : July 2010

#include <string>
#include <set>
#include <vector>

namespace mX_parms_utils
{
  enum parm_names{CKT_FILENAME,T_START,T_STEP,T_STOP,TOL,K,MAX_RESTART,INIT_COND,PARMS_FILE,PREV};

  void parse_command_line(int argc, std::vector<std::string> &argv, std::string &ckt_filename, double &t_start, double &t_step, double &t_stop, double &tol, int &k, int& max_restarts, 
                          std::vector<double> &init_cond, std::string &parms_file, std::set<int> &specified_parms, int p, int pid);

  std::vector<std::string> get_command_line_equivalent_from_file(std::string &filename);

  void get_parms(int argc, char* argv[], std::string &ckt_filename, double &t_start, double &t_step, double &t_stop, double &tol, int &k, int &max_restarts,
                 std::vector<double> &init_cond, bool &init_cond_specified, int p, int pid);
}

#endif

#ifndef __MX_SOLVER_DEFS_H__
#define __MX_SOLVER_DEFS_H__

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
// Date : Sept 2018

#include <vector>
#include <list>
#include <map>

namespace mX_comm_utils
{
  struct data_transfer_instruction;
}

namespace mX_vector_utils
{
  struct distributed_vector;
}
 
namespace mX_matrix_utils
{
  struct distributed_sparse_matrix;
}

namespace mX_solver_utils
{
  void gmres(mX_matrix_utils::distributed_sparse_matrix* A, mX_vector_utils::distributed_vector &b, mX_vector_utils::distributed_vector &x,
             double &tol, double &err, int k, int &iters, int &restarts);

  void rotg(double& a, double& b, double& c, double& s);
}

#endif

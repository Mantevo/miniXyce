#ifndef __MX_LINEAR_DAE_H__
#define __MX_LINEAR_DAE_H__

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

#include <map>
#include <list>
#include <vector>
#include "mX_source.h"
#include "mX_sparse_matrix.h"
#include "mX_vector.h"
#include "mX_device.h"

using namespace mX_source_utils;
using namespace mX_matrix_utils;
using namespace mX_vector_utils;
using namespace mX_device_utils;

namespace mX_DAE_utils
{
  struct mX_DAE_RHS_entry
  {
    // a single entry in the RHS of a DAE of the form "A x + B x_dot = b(t)"
      // each entry in b(t) is a list of scaled sources
        // linear combination of multiple voltage/current sources

    std::list<mX_scaled_source*> scaled_src_list;
  };

  struct mX_DAE
  {
    // a linear DAE is of the form "A x + B x_dot = b(t)"
      // A and B are distributed sparse matrices
      // b(t) is a time-varying function
        // every entry of which is a linear combination of voltage/current sources

    distributed_sparse_matrix* A;
    distributed_sparse_matrix* B;
    std::vector<mX_DAE_RHS_entry*> b;

    std::map<char, mX_device*> devices;
  };

  distributed_vector evaluate_b(double t, mX_DAE* dae);
  std::vector<double> evaluate_b_old(double t, mX_DAE* dae);
  void load_matrices(mX_DAE* dae);
  void destroy(mX_DAE* dae);
  void destroy_RHS(mX_DAE_RHS_entry* entry);
}

#endif

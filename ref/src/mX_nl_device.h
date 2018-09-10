#ifndef __MX_NL_DEVICE_DEFS_H__
#define __MX_NL_DEVICE_DEFS_H__

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
// Date : July 2018

#include "mX_vector.h"
#include "mX_sparse_matrix.h"
#include "mX_device.h"
#include <sstream>
#include <string>
#include <vector>

using namespace mX_vector_utils;
using namespace mX_matrix_utils;

namespace mX_DAE_utils
{
  struct mX_DAE;
}


namespace mX_device_utils
{

class mX_diode: public mX_device
{
        public:

        mX_diode();
        virtual ~mX_diode() {}

        int add_device( std::istringstream& input_str, int extra_nodes_ptr, mX_DAE_utils::mX_DAE* dae );

        void load_vector( distributed_vector& vec );
        void load_matrices();

        int modelCost() { return 50; }

        private:

};

std::vector<double> diode_l1(std::vector<double> &x);

}

#endif // __MX_NL_DEVICE_DEFS_H__

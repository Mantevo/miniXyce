#ifndef __MX_DEVICE_DEFS_H__
#define __MX_DEVICE_DEFS_H__

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

class mX_device
{
        // this class represents a device model, that can be time-invariant
                // or time-varying.
        public:

        mX_device( char label ) { label_ = label; }
        virtual ~mX_device() {}

        // Add a device described by this input stream, returns number of internal nodes required 
        virtual int add_device( std::istringstream& input_str, int extra_nodes_ptr, mX_DAE_utils::mX_DAE* dae ) { return 0; }

        // Return number of internal nodes introduced by this device type.
        virtual std::vector<std::pair<int,int> >& get_jacobian_stamp() { return jac_stamp_; }

        // Evaluate devices
        virtual void load_vector( distributed_vector& vec ) {}
        virtual void load_matrices() {}

        // Some indication of flops for evaluating this device model.
        virtual int modelCost() { return 0; }

        // Character (letter) indicating which device type.
        virtual char deviceLabel() { return label_; }

        protected:

        std::vector<distributed_sparse_matrix_entry* > device_entries_;

        private:

        char label_;
        std::vector<std::pair<int,int> > jac_stamp_;
};

mX_device* create_device( char device_label );

void add_device_nodes( char device_label, std::istringstream& input_str, std::vector<int>& node_list );

int num_internal_nodes( char device_label, std::istringstream& input_str );
}

#endif // __MX_DEVICE_DEFS_H__

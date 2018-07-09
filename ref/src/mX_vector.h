#ifndef __MX_VECTOR_DEFS_H__
#define __MX_VECTOR_DEFS_H__

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

#include <vector>
#include <list>
#include <map>

namespace mX_comm_utils
{
  struct data_transfer_instruction;
}

namespace mX_vector_utils
{
  struct distributed_vector
  {
    // this is a distributed vector using the same data transfer information from
      // the matrix and the row distribution from the matrix.

    std::vector<int> local_rows, overlap_rows;
    std::map<int,int> row_to_idx;
    std::vector<double> values;

    int n;
    int p;      // total number of parallel processors
    int my_pid; // this processors id

    std::list<mX_comm_utils::data_transfer_instruction*> send_instructions;
    std::list<mX_comm_utils::data_transfer_instruction*> recv_instructions;

    distributed_vector(int n_in, int p_in, int my_pid_in,
                       const std::vector<int>& local_rows_in, const std::vector<int>& overlap_rows_in,
                       const std::list<mX_comm_utils::data_transfer_instruction*>& send_inst_in,
                       const std::list<mX_comm_utils::data_transfer_instruction*>& recv_inst_in)
    : n( n_in ),
      p( p_in ),
      my_pid( my_pid_in ),
      local_rows( local_rows_in ),
      overlap_rows( overlap_rows_in ),
      send_instructions( send_inst_in ),
      recv_instructions( recv_inst_in )
    {
      int idx = 0;
      for (std::vector<int>::iterator it = local_rows.begin(); it != local_rows.end(); it++, idx++)
        row_to_idx[ *it ] = idx;
      for (std::vector<int>::iterator it = overlap_rows.begin(); it != overlap_rows.end(); it++, idx++)
        row_to_idx[ *it ] = idx;
      
      values.resize( row_to_idx.size() );
    }
  };

  void axpy(double alpha, const distributed_vector& x, double beta, distributed_vector& y);

  void scale(distributed_vector& x, double alpha);
 
  double norm2(distributed_vector& v);

  void print_vector(distributed_vector& v);
 
  void assemble_vector(distributed_vector& x);

  void init_value(distributed_vector& x, double val);

}

#endif

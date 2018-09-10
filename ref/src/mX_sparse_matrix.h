#ifndef __MX_SPARSE_MATRIX_DEFS_H__
#define __MX_SPARSE_MATRIX_DEFS_H__

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

namespace mX_DAE_utils
{
  struct mX_DAE_RHS_entry;
}

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
  struct distributed_sparse_matrix_entry
  {
    int column;    // global column index
    double value;  // value stored in the matrix
    distributed_sparse_matrix_entry* next_in_row;  // pointer to next entry in the same row

    distributed_sparse_matrix_entry()
      : column(-1), value(0.0), next_in_row(0)
    {}

    distributed_sparse_matrix_entry(int col_in, double val_in, distributed_sparse_matrix_entry* nir_in)
      : column( col_in ), value( val_in ), next_in_row( nir_in )
    {}
  };

  struct distributed_sparse_matrix
  {
    // the data structure for a distributed sparse matrix is a 1-d threaded list
      // a set of pointers called row_headers point to the first entry of each row
      // each row entry in turn points to the next entry in the same row

    // but a distributed matrix needs more data than this
      // there is a list of data transfer instructions
        // these instructions are to be followed whenever a mat-vec product is needed
    
    // each processor also stores 2 entries start_row and end_row
      // it is assumed that all processors store contiguous rows of the distributed matrix

    int start_row;
    int end_row;
    int local_nnz;
    std::vector<int> local_rows;
    std::vector<int> overlap_rows, overlap_recv, overlap_send;

    int n;
    int p;      // total number of parallel processors
    int my_pid; // this processors id

    distributed_sparse_matrix_entry ground_node;
    std::vector<distributed_sparse_matrix_entry*> row_headers;
    std::list<mX_comm_utils::data_transfer_instruction*> send_instructions;
    std::list<mX_comm_utils::data_transfer_instruction*> recv_instructions;

    distributed_sparse_matrix();

    distributed_sparse_matrix(int total_unknowns, int num_procs, int pid);
  };

  distributed_sparse_matrix_entry* distributed_sparse_matrix_insert(distributed_sparse_matrix* M,
                                                                    int row_idx, int col_idx);

  void distributed_sparse_matrix_insert(distributed_sparse_matrix* M, int row_idx, int col_idx, double val);

  void distributed_sparse_matrix_finish(distributed_sparse_matrix* A, distributed_sparse_matrix* B,
                                        const std::vector<mX_DAE_utils::mX_DAE_RHS_entry*>& b);

  void distributed_sparse_matrix_add_to(distributed_sparse_matrix* M, int row_idx, int col_idx, double val);

  void sparse_matrix_vector_product(distributed_sparse_matrix* A, mX_vector_utils::distributed_vector& x, mX_vector_utils::distributed_vector& y);

  void sparse_matrix_vector_product(distributed_sparse_matrix* A, std::vector<double> &x, std::vector<double> &y);

  void destroy_matrix(distributed_sparse_matrix* A);

  void print_matrix(distributed_sparse_matrix &A);

  void init_value(distributed_sparse_matrix& A, double val);
}


#endif

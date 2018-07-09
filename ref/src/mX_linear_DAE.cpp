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

#ifdef HAVE_MPI
#include <mpi.h> // If this routine is compiled with -DHAVE_MPI
                 // then include mpi.h
#endif
#include "mX_linear_DAE.h"
#include <algorithm>
#include <iostream>

using namespace mX_source_utils;
using namespace mX_matrix_utils;
using namespace mX_linear_DAE_utils;

distributed_vector mX_linear_DAE_utils::evaluate_b(double t, mX_linear_DAE* dae)
{
  distributed_vector global_result( dae->A->n, dae->A->p, dae->A->my_pid,
                                    dae->A->local_rows, dae->A->overlap_rows,
                                    dae->A->send_instructions, dae->A->recv_instructions );

  std::vector<mX_linear_DAE_RHS_entry*>::iterator it2;
  std::map<int,int>::iterator row_it;
  for (int i = 0; i < dae->b.size(); i++)
  {
    if ( dae->b[i] )
    {
      double sum = 0.0;
      std::list<mX_scaled_source*>::iterator it3 = dae->b[i]->scaled_src_list.begin();

      for( ; it3 != dae->b[i]->scaled_src_list.end(); it3++)
      {
        mX_source* src = (*it3)->src;
        sum += src->output(t * (*it3)->scale);
      }

      row_it = global_result.row_to_idx.find( i );
      if (row_it != global_result.row_to_idx.end())
      {
        global_result.values[row_it->second] = sum;
        std::cout << "b[" << i << "] = " << sum << " is a scaled source" << std::endl;
      }
    }
  }

  mX_vector_utils::assemble_vector( global_result );

  return global_result;
}

std::vector<double> mX_linear_DAE_utils::evaluate_b_old(double t, mX_linear_DAE* dae)
{
  // given a linear DAE "A x + B x_dot = b(t)"
  // and a particular time point t
  // this function computes and returns the vector b at that time point t

  std::vector<mX_linear_DAE_RHS_entry*>::iterator it2;
  std::vector<double> local_result, global_result;

  for (it2 = dae->b.begin(); it2 != dae->b.end(); it2++)
  {
    if (!(*it2))
    {
      local_result.push_back(0.0);
    }
    else
    {
      double sum = (double)(0);
      std::list<mX_scaled_source*>::iterator it3;

      for(it3 = (*it2)->scaled_src_list.begin(); it3 != (*it2)->scaled_src_list.end(); it3++)
      {
	mX_source* src = (*it3)->src;
	sum += src->output((double)(t)) * ((*it3)->scale);
      }

      local_result.push_back(sum);
    }
  }

#ifdef HAVE_MPI
  global_result.resize( dae->b.size(), 0.0 );
  MPI_Allreduce(&local_result[0],&global_result[0],dae->b.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return global_result;
#else
  return local_result;
#endif

}

void mX_linear_DAE_utils::destroy_RHS(mX_linear_DAE_RHS_entry* entry)
{
  if (entry)
  {
    // delete send_instructions
    while (!entry->scaled_src_list.empty())
    {
      mX_scaled_source* curr = entry->scaled_src_list.front();
      if (curr) delete curr->src; delete curr; curr=0;
      entry->scaled_src_list.pop_front();
    }
    delete entry; entry=0;
  }
}

void mX_linear_DAE_utils::destroy(mX_linear_DAE* dae)
{
  // Destroy A
  destroy_matrix( dae->A );

  // Destroy B
  destroy_matrix( dae->B );

  // Destroy b
  for (int i=0; i<dae->b.size(); ++i)
  {
    destroy_RHS( dae->b[i] );
  }
  dae->b.resize(0);

  delete dae; dae=0;
}

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
#include "mX_comm.h"
#include "mX_vector.h"
#include <iostream>
#include <cmath>

using namespace mX_comm_utils;

void mX_vector_utils::print_vector(distributed_vector& x)
{
  int n = 1, my_pid = 0;

#ifdef HAVE_MPI
  /* Find this processor number */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);

  /* Find the number of processors */
  MPI_Comm_size(MPI_COMM_WORLD, &n);
#endif

  for (int i=0; i<n; ++i)
  {
    if (i == my_pid)
    {
      if (my_pid == 0)
        std::cout << "Proc\tGID\tValue" << std::endl;

      int j=0;
      for ( ; j<x.local_rows.size(); ++j)
      {
        std::cout << my_pid << "\t" << x.local_rows[j] << "\t" << x.values[j] << std::endl;
      }

      for ( ; j<x.row_to_idx.size(); ++j)
      {
        int jj = j-x.local_rows.size();
        std::cout << my_pid << "\t" << x.overlap_rows[jj] << "\t" << x.values[j] << " [O]" << std::endl;
      }

    }
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
}

double mX_vector_utils::norm2(distributed_vector& x)
{
  double global_norm;
  double local_norm = 0.0;

  for (int i = 0; i < x.local_rows.size(); i++)
  {
    local_norm += x.values[i]*x.values[i];
  }
#ifdef HAVE_MPI
  MPI_Allreduce(&local_norm,&global_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
  global_norm = local_norm;
#endif

  return std::sqrt(global_norm);
}

double mX_vector_utils::dot(distributed_vector& x, distributed_vector& y)
{
  double global_dot;
  double local_dot = 0.0;

  for (int i = 0; i < x.local_rows.size(); i++)
  {
    local_dot += x.values[i]*y.values[i];
  }
#ifdef HAVE_MPI
  MPI_Allreduce(&local_dot,&global_dot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
  global_dot = local_dot;
#endif

  return global_dot;
}

void mX_vector_utils::assemble_vector(distributed_vector& x)
{
#ifdef HAVE_MPI
  // now's the time to send and receive the results from other processors

  std::list<data_transfer_instruction*>::iterator it1;
  std::list<int>::iterator it2;
  std::map<int,int>::iterator it3;

  for (it1 = x.send_instructions.begin(); it1 != x.send_instructions.end(); it1++)
  {

    for (it2 = (*it1)->indices.begin(); it2 != (*it1)->indices.end(); it2++)
    {
      int x_idx = 0;
      double x_vec_entry = 0.0;
      it3 = x.row_to_idx.find( *it2 );
      if (it3 != x.row_to_idx.end())
      {
        x_vec_entry = x.values[it3->second];
        //std::cout << "PID " << x.my_pid << ", SENDING row_idx = " << *it2 << ", x_idx = " << it3->second << ", x_values = " << x.values[it3->second] << std::endl;
      }

      MPI_Send(&x_vec_entry,1,MPI_DOUBLE,(*it1)->pid,*it2,MPI_COMM_WORLD);
    }
  }
  for (it1 = x.recv_instructions.begin(); it1 != x.recv_instructions.end(); it1++)
  {
    for (it2 = (*it1)->indices.begin(); it2 != (*it1)->indices.end(); it2++)
    {
      int x_idx = 0;
      double x_vec_entry = 0.0;
      MPI_Status status;
      MPI_Recv(&x_vec_entry,1,MPI_DOUBLE,(*it1)->pid,*it2,MPI_COMM_WORLD,&status);

      it3 = x.row_to_idx.find( *it2 );
      if (it3 != x.row_to_idx.end())
      {
        x.values[it3->second] += x_vec_entry;
        //std::cout << "PID " << x.my_pid << ", UPDATED row_idx = " << *it2 << ", x_idx = " << it3->second << ", x_values = " << x.values[it3->second] << std::endl;
      }
   }
 }

#endif
}

void mX_vector_utils::init_value(distributed_vector& x, double val)
{
  for (std::vector<double>::iterator it = x.values.begin(); it != x.values.end(); it++)
  {
    (*it) = val;
  }
} 

void mX_vector_utils::axpy(double alpha, const distributed_vector& x, double beta, distributed_vector& y)
{
  std::vector<double>::const_iterator itx;
  std::vector<double>::iterator ity;
  for (itx = x.values.begin(), ity = y.values.begin(); itx != x.values.end(); itx++, ity++)
  {
    (*ity) = beta*(*ity) + alpha*(*itx);
  }
}

void mX_vector_utils::scale(distributed_vector& x, double val)
{
  for (std::vector<double>::iterator it = x.values.begin(); it != x.values.end(); it++)
  {
    (*it) *= val;
  }
}


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
#include "mX_sparse_matrix.h"
#include "mX_vector.h"
#include "mX_linear_DAE.h"
#include <map>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>

using namespace mX_comm_utils;
using namespace mX_matrix_utils;
using namespace mX_vector_utils;
using namespace mX_linear_DAE_utils;

distributed_sparse_matrix::distributed_sparse_matrix()
  : start_row(0), end_row(0), local_nnz(0)
{}

void mX_matrix_utils::distributed_sparse_matrix_insert(distributed_sparse_matrix* M, int row_idx, int col_idx, double val, int proc)
{
  // Inserts row_idx,col_idx,val into M.
  if ((row_idx < 0) || (col_idx < 0))
  {
    // check for negative indices
      // useful in cases where you want to ignore the ground
        // because the ground is not really a node

    return;
  }

  if (proc == M->my_pid)
  { 
    bool inserted = false;
    bool found = false;

    distributed_sparse_matrix_entry* prev = 0;
    distributed_sparse_matrix_entry* curr = M->row_headers[row_idx];

    while (curr && !found)
    {
      if (curr->column < col_idx)
      {
        prev = curr;
        curr = curr->next_in_row;
      }
      else
      {
        if (curr->column == col_idx)
        {
          curr->value = curr->value + val;
          found = true;
        }
        else
        {
          // The entry should go here, but it's not less than or equal to, so the col_idx must be >.
          break;
        }
      }
    }

    if (!found)
    {
      distributed_sparse_matrix_entry* entry_ptr_1 = new distributed_sparse_matrix_entry(col_idx, val, curr);

      if (prev)
      {
        prev->next_in_row = entry_ptr_1;
      }
      else
      {
        M->row_headers[row_idx] = entry_ptr_1;
      }

      inserted = true;
    }

    // The number of local nonzeros has gone up only if this is a new entry.
    if (inserted)
      M->local_nnz++;
  }
} 


void mX_matrix_utils::distributed_sparse_matrix_finish(distributed_sparse_matrix* A, distributed_sparse_matrix* B,
                                                       const std::vector<mX_linear_DAE_RHS_entry*>& b)
{
  // MPI_MIN is used to assign processors to shared nodes, so put the 
  // fill value as the number of processors for nodes not required by this processor.
  // 
  // The sends and recvs are set up using the nonzero pattern of both A, B, and b.
  // This ensures we can form the circuit Jacobian and perform matrix-vector products.
  //
  int fill_value = A->p; 

  std::vector<int> local_nnz_rows(A->n,fill_value), global_nnz_rows(A->n,0);
  std::vector<int> local_nnz_cols;

  for (int i=0; i<A->n; i++)
  {
    distributed_sparse_matrix_entry* curr = A->row_headers[i];
    if (curr)
      local_nnz_rows[i] = A->my_pid;
    while (curr)
    {
      local_nnz_cols.push_back(curr->column);
      curr = curr->next_in_row;
    }

    curr = B->row_headers[i];
    if (curr)
      local_nnz_rows[i] = B->my_pid;
    while (curr)
    {
      local_nnz_cols.push_back(curr->column);
      curr = curr->next_in_row;
    } 

    mX_linear_DAE_RHS_entry* curr_b = b[i];
    if (curr_b)
    {
      local_nnz_rows[i] = A->my_pid; 
      local_nnz_cols.push_back(i);
    }
  } 

  // Order the columns needed for matrix-vector products.
  std::sort( local_nnz_cols.begin(), local_nnz_cols.end() );
  local_nnz_cols.erase( std::unique( local_nnz_cols.begin(), local_nnz_cols.end() ), local_nnz_cols.end() ); 
 
#ifdef HAVE_MPI
  MPI_Allreduce(&local_nnz_rows[0],&global_nnz_rows[0],A->n,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
#endif

  for (int i=0; i<A->n; i++)
  {
    // Determine who owns the nodes that are shared
    if (global_nnz_rows[i]==A->my_pid)
    {
      A->local_rows.push_back(i);
      std::cout << "Processor " << A->my_pid << " owns row " << i << std::endl;

      std::vector<int>::iterator it = std::find( local_nnz_cols.begin(), local_nnz_cols.end(), i );
      if ( it != local_nnz_cols.end() )
      {
        local_nnz_cols.erase( it );
      }
    } 
    else if (local_nnz_rows[i]==A->my_pid)
    {
      // This processor will need to send the result of it's matrix-vector product to global_nnz_rows[i].
      // Since there is already a send_instructions for the column, we'll call this a recv_instructions.
      // For linear circuits, matrices are symmetric, so recv_instructions will be the transpose of send.
      bool recv_instruction_posted = false;
      int pid_to_recv_info = global_nnz_rows[i];

      std::list<data_transfer_instruction*>::iterator it1 = A->recv_instructions.begin();

      for (; it1 != A->recv_instructions.end(); it1++)
      {
        if ((*it1)->pid == pid_to_recv_info)
        {
          std::list<int>::iterator it2 = std::find( (*it1)->indices.begin(), (*it1)->indices.end(), i );
          if (it2 == (*it1)->indices.end())
          {
            (*it1)->indices.push_back(i);
          }
          recv_instruction_posted = true;
        }
      }

      if (!recv_instruction_posted)
      {
        data_transfer_instruction* dti_ptr_1 = new data_transfer_instruction(pid_to_recv_info);
        dti_ptr_1->indices.push_back(i);

        A->recv_instructions.push_back(dti_ptr_1);
      }
    }
  }

  for (std::vector<int>::iterator it = local_nnz_cols.begin(); it != local_nnz_cols.end(); it++)
  {
    std::cout << "Processor " << A->my_pid << " has overlap row " << *it << std::endl;
  }

  // Set the local and ghost rows with both A and B.
  A->overlap_rows = local_nnz_cols;
  B->local_rows = A->local_rows;
  B->overlap_rows = local_nnz_cols;

  std::vector<int> local_overlap( A->p ), global_overlap( A->p );
  local_overlap[ A->my_pid ] = local_nnz_cols.size();

#ifdef HAVE_MPI
  MPI_Allreduce(&local_overlap[0],&global_overlap[0],A->p,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  int total_overlap = std::accumulate( global_overlap.begin(), global_overlap.end(), 0 );

  int ptr = 0;
  std::vector<int> local_recv_rows( 2*total_overlap ), global_recv_rows( 2*total_overlap );
  for (int proc = 0; proc < A->p; proc++)
  {
    if (proc == A->my_pid)
    {
      for (std::vector<int>::const_iterator it = local_nnz_cols.begin(); it != local_nnz_cols.end(); it++)
      {
        local_recv_rows[ ptr ] = proc;
        local_recv_rows[ ptr+1 ] = *it;
        ptr += 2;
      }
      break;
    }
    ptr += 2*global_overlap[proc]; 
  }
  
#ifdef HAVE_MPI
  MPI_Allreduce(&local_recv_rows[0],&global_recv_rows[0],2*total_overlap,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  for (int i=0; i<total_overlap; i++)
  {
    int pid_to_send_info = global_recv_rows[2*i];
    int row = global_recv_rows[2*i+1];

    if (pid_to_send_info != A->my_pid)
    {
      bool send_instruction_posted = false;
      std::vector<int>::const_iterator it = std::find( A->local_rows.begin(), A->local_rows.end(), row );
      if (it != A->local_rows.end())
      {
        std::cout << "Processor " << A->my_pid << " needs to send row " << row << " to processor " << pid_to_send_info << std::endl;

        std::list<data_transfer_instruction*>::iterator it1 = A->send_instructions.begin();

        for (; it1 != A->send_instructions.end(); it1++)
        {
          if ((*it1)->pid == pid_to_send_info)
          {
            std::list<int>::iterator it2 = std::find( (*it1)->indices.begin(), (*it1)->indices.end(), row );
            if (it2 == (*it1)->indices.end())
            {
              (*it1)->indices.push_back(row);
            }
            send_instruction_posted = true;
            break;
          }
        }

        if (!send_instruction_posted)
        {
          data_transfer_instruction* dti_ptr_1 = new data_transfer_instruction(pid_to_send_info);
          dti_ptr_1->indices.push_back(row);

          A->send_instructions.push_back(dti_ptr_1);
        }
      }
    }
  }
}


void mX_matrix_utils::distributed_sparse_matrix_add_to(distributed_sparse_matrix* M, int row_idx, int col_idx, double val)
{
  // implements M[row_idx][col_idx] += val
    // man, in the distributed matrix world, this simple thing takes such a lot of effort!

  if ((row_idx < 0) || (col_idx < 0))
  {
    // check for negative indices
      // useful in cases where you want to ignore the ground
        // because the ground is not really a node

    return;
  }

  bool inserted = false;
  bool found = false;

  distributed_sparse_matrix_entry* prev = 0;
  distributed_sparse_matrix_entry* curr = M->row_headers[row_idx];
  
  while (curr && !found)
  {
    if (curr->column < col_idx)
    {
      prev = curr;
      curr = curr->next_in_row;
    }
    else
    {
      if (curr->column == col_idx)
      {
        curr->value = curr->value + val;
        found = true;
      }
      else
      {
        // The entry should go here, but it's not less than or equal to, so the col_idx must be >.
        break;
      }
    }
  }

  if (!found)
  {
    distributed_sparse_matrix_entry* entry_ptr_1 = new distributed_sparse_matrix_entry(col_idx, val, curr);

    if (prev)
    {
      prev->next_in_row = entry_ptr_1;
    }
    else
    {
      M->row_headers[row_idx] = entry_ptr_1;
    }

    inserted = true;
  }

  // The number of local nonzeros has gone up only if this is a new entry.
  if (inserted)
    M->local_nnz++;

  return;
}


void mX_matrix_utils::sparse_matrix_vector_product(distributed_sparse_matrix* A, distributed_vector &x, distributed_vector &y)
{
  // compute the matrix vector product A*x and return it in y

  // Zero out values in y.
  init_value( y, 0.0 );

  for (int i=0; i<A->n; i++)
  {
    distributed_sparse_matrix_entry* curr = A->row_headers[i];

    while (curr)
    {
      int x_idx, y_idx;

      int col_idx = curr->column;
      double x_value = 0.0;
      std::map<int,int>::iterator it2 = x.row_to_idx.find( col_idx ); 
      if (it2 != x.row_to_idx.end())
      {
        x_value = x.values[it2->second];
        std::cout << "PID " << A->my_pid << ", col_idx = " << col_idx << ", x_idx = " << it2->second << ", x_values = " << x_value << std::endl;
      } 

      std::map<int,int>::iterator it3 = y.row_to_idx.find( i );
      if (it3 != y.row_to_idx.end())
      {
        y.values[it3->second] += (curr->value)*x_value;
        std::cout << "PID " << A->my_pid << ", row_idx = " << i << ", y_idx = " << it3->second  << ", y_values = " << y.values[it3->second] << std::endl;
      }

      // Update pointer
      curr = curr->next_in_row;
    }
  }

  // Now assemble vector in parallel before returning.
  mX_vector_utils::assemble_vector( y );

}


void mX_matrix_utils::sparse_matrix_vector_product(distributed_sparse_matrix* A, std::vector<double> &x, std::vector<double> &y)
{
  // compute the matrix vector product A*x and return it in y

#ifdef HAVE_MPI
  // ok, now's the time to follow the send instructions that each pid has been maintaining

  std::list<data_transfer_instruction*>::iterator it1;
  
  for (it1 = A->send_instructions.begin(); it1 != A->send_instructions.end(); it1++)
  {
    std::list<int>::iterator it2;

    for (it2 = (*it1)->indices.begin(); it2 != (*it1)->indices.end(); it2++)
    {
      MPI_Send(&x[*it2],1,MPI_DOUBLE,(*it1)->pid,*it2,MPI_COMM_WORLD);
    }
  }
#endif

  // and everytime a processor receives an x_vec entry
    // it stores the entry in a temporary map

  std::map<int,double> x_vec_entries;
  std::map<int,double>::iterator it3;

  for (int i = 0; i < A->n; i++)
  {
    // compute the mat_vec product for the i'th row
    if (y.size() > i)
    {
      y[i] = (double)(0);
    }
    else
    {
      y.push_back((double)(0));
    }

    distributed_sparse_matrix_entry* curr = A->row_headers[i];

    while(curr)
    {
      int col_idx = curr->column;
      
      if (std::find( A->local_rows.begin(), A->local_rows.end(), col_idx )!= A->local_rows.end())
      {
        // aha, this processor has the correct x_vec entry locally

        y[i] += (curr->value)*x[col_idx];
      }
#ifdef HAVE_MPI
      else
      {
        // this processor does not have the x_vec entry locally
          // but some other processor might have sent it to this guy
            // in which case the entry would have been stored in the local x_vec_entries map
            // so check if the entry is in the map

        it3 = x_vec_entries.find(col_idx);

        if (it3 != x_vec_entries.end())
        {
          y[i] += (double)(it3->second)*(curr->value);
        }
        else
        {
          // no, the entry is not in the map either
            // so this processor waits until someone sends the entry
            // and once it gets the entry, it does two things
              // puts the entry in the map for future reference
              // continues with the matrix vector multiplication

          double x_vec_entry;
          MPI_Status status;
          MPI_Recv(&x_vec_entry,1,MPI_DOUBLE,MPI_ANY_SOURCE,col_idx,MPI_COMM_WORLD,&status);

          x_vec_entries[col_idx] = x_vec_entry;
          y[i] += x_vec_entry*(curr->value);
        }
      }
#endif
      curr = curr->next_in_row;
    }
  }
}

void mX_matrix_utils::gmres(distributed_sparse_matrix* A, std::vector<double> &b, std::vector<double> &x0, double &tol, double &err, int k, std::vector<double> &x, int &iters, int &restarts)
{}

/*
void mX_matrix_utils::gmres(distributed_sparse_matrix* A, std::vector<double> &b, std::vector<double> &x0, double &tol, double &err, int k, std::vector<double> &x, int &iters, int &restarts)
{
  // here's the star of the show, the guest of honor, none other than Mr.GMRES
  
  // first Mr.GMRES will compute the error in the initial guess
    // if it's already smaller than tol, he calls it a day
    // otherwise he settles down to work in mysterious ways his wonders to perform

  x = x0;

  std::vector<double> temp1;
  sparse_matrix_vector_product(A,x,temp1);

  for (int i = 0; i < temp1.size(); i++)
  {
    temp1[i] -= b[i];
  }

  err = norm(temp1);
  restarts = -1;
  iters = 0;

  while (err > tol)
  {
    // at the start of every re-start
      // the initial guess is already stored in x

    restarts++;
    
    std::vector< std::vector<double> > V;
    sparse_matrix_vector_product(A,x,temp1);

    for (int i = 0; i <= A->n; i++)
    {
      temp1[i] -= b[i];
      temp1[i] *= (double)(-1);

      V.push_back(temp1);
    }

    double beta = norm(temp1);

    for (int i = 0; i < A->n; i++)
    {
      V[i][0] /= beta;
    }

    err = beta;
    iters = 0;

    std::vector<double> cosines;
    std::vector<double> sines;
    std::vector<double> g;
    std::vector< std::vector<double> > R;
    
    g.push_back(beta);

    // ok, Mr.GMRES has determined the initial values for
      // V,R,g,sines,cosines,err and iters
      // he updates these at every iteration until
        // either err becomes less than tol
        // or a new restart is required

    // note that Mr.GMRES does not think it necessary to store the upper Hessenberg matrix at each iteration
      // he computes R at each iteration from the new Hessenberg matrix column
      // and he's clever enough to determine err without having to solve for x at each iteration

    while ((err > tol) && (iters < k))
    {
      iters++;

      // Mr.GMRES is now going to update the V matrix
        // for which he will require a matrix vector multiplication

      std::vector<double> temp1;
      std::vector<double> temp2;

      for (int i = 0; i < A->n; i++)
      {
        temp1.push_back(V[i][iters-1]);  
      }
      sparse_matrix_vector_product(A,temp1,temp2);

      // Right, Mr.GMRES now has the matrix vector product
        // now he will orthogonalize this vector with the previous ones 
          // with some help from Messrs Gram and Schmidt
      
      std::vector<double> new_col_H;

      for (int i = 0; i < iters; i++)
      {
        double local_dot = 0.0;
        double global_dot;

        for (int j = 0; j < A->n; j++)
        {
          local_dot += temp2[j]*V[j][i];
        }
#ifdef HAVE_MPI
        MPI_Allreduce(&local_dot,&global_dot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
                                global_dot = local_dot;
#endif        
        for (int j = 0; j < A->n; j++)
        {
          temp2[j] -= global_dot*V[j][i];
        }

        new_col_H.push_back(global_dot);
      }

      new_col_H.push_back(norm(temp2));

      for (int i = 0; i < A->n; i++)
      {
        temp2[i]/= new_col_H.back();
        V[i].push_back(temp2[i]);
      }

      // Right, Mr.GMRES has successfully updated V
        // on the side, he has also been computing the new column of the Hessenberg matrix
      // now he's going to get the new column of R using the current sines and cosines
        // and he will also add a new sine and a new cosine for future use

      for (int i = 0; i < iters-1; i++)
      {
        double old_i = new_col_H[i];
        double old_i_plus_one = new_col_H[i+1];

        new_col_H[i] = cosines[i]*old_i + sines[i]*old_i_plus_one;
        new_col_H[i+1] = -sines[i]*old_i + cosines[i]*old_i_plus_one;
      }

      double r = std::sqrt(new_col_H[iters-1]*new_col_H[iters-1] + new_col_H[iters]*new_col_H[iters]);
      cosines.push_back(new_col_H[iters-1]/r);
      sines.push_back(new_col_H[iters]/r);
      
      double old_i = new_col_H[iters-1];
      double old_i_plus_one = new_col_H[iters];

      new_col_H[iters-1] = cosines.back()*old_i + sines.back()*old_i_plus_one;
      new_col_H.pop_back();

      R.push_back(new_col_H);

      // Right, the new column of R is ready
      // the only thing left to do is to update g
        // which will also tell Mr.GMRES what the new error is

      double old_g = g[iters-1];
      g[iters-1] = old_g*cosines.back();
      g.push_back(-old_g*sines.back());

      err = std::abs(g.back());
    }

    // ok, so either Mr.GMRES has a solution
      // or he's being forced to restart
    // either way, he needs to compute x
      // now he needs to solve Ry = g
      // after which he will say x += (V without its last column)*y

    std::vector<double> y;

    for (int i = iters-1; i >= 0; i--)
    {
      double sum = (double)(0);

      for (int j = iters-1; j > i; j--)
      {
        sum += R[j][i]*y[iters-1-j];
      }

      y.push_back((g[i] - sum)/R[i][i]);
    }

    // ok, so y is ready (although it's stored upside down)

    for (int i = 0; i < A->n; i++)
    {
      double sum = (double)(0);

      for (int j = iters-1; j >= 0; j--)
      {
        sum += y[iters-1-j]*V[i][j];
      }

      x[i] += sum;
    }

    // the new x is also ready
      // either return it or use it as an initial guess for the next restart
  }
  
  // if Mr.GMRES ever reaches here, it means he's solved the problem

  if (restarts < 0)
  {
    restarts = 0;
  }
}
*/

void mX_matrix_utils::destroy_matrix(distributed_sparse_matrix* A)
{
  if (A)
  {
      // delete row_headers
      for (int j=0, cnt=0; j<A->n; ++j, ++cnt)
      {
        distributed_sparse_matrix_entry* curr = (*A).row_headers[cnt], *next = 0;
        
        while (curr)
        {
          next = curr->next_in_row;
          delete curr; curr=0;
          curr = next;
        }
      }
      A->row_headers.resize(0);  

      // delete send_instructions
      while (!A->send_instructions.empty())
      {
        data_transfer_instruction* curr = A->send_instructions.front();
        delete curr; curr=0;
        A->send_instructions.pop_front();
      }

      delete A; A=0;
  }
}

void mX_matrix_utils::print_matrix(distributed_sparse_matrix &A)
{
  if (A.my_pid == 0)
    std::cout << "Matrix\nProc\tRow\tColumn\tValue" << std::endl;

  for (int i=0; i<A.p; ++i)
  {
    if (i == A.my_pid)
    {
      for (int j=0; j<A.n; ++j)
      {
        distributed_sparse_matrix_entry* curr = A.row_headers[j];

        while (curr)
        {
          std::cout << A.my_pid << "\t" << j << "\t" << curr->column << "\t" << curr->value << std::endl;
          curr = curr->next_in_row;
        }
      }
    }
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
}



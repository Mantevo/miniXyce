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

#ifdef HAVE_MPI
#include <mpi.h> // If this routine is compiled with -DHAVE_MPI
                 // then include mpi.h
#endif
#include "mX_comm.h"
#include "mX_sparse_matrix.h"
#include "mX_vector.h"
#include "mX_solver.h"

#include <map>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>

using namespace mX_comm_utils;
using namespace mX_matrix_utils;
using namespace mX_vector_utils;
using namespace mX_solver_utils;
using namespace mX_linear_DAE_utils;

void mX_solver_utils::gmres(distributed_sparse_matrix* A, distributed_vector &b, distributed_vector &x, double &tol, double &err, int k, int &iters, int &restarts)
{
  // x is expected to contain the initial guess, x0, for the linear system upon entry and will contain the approximate solution upon exit.

  // Define a local k based on the size of the linear system.
  int local_k = std::min( b.n, k );

  std::vector< distributed_vector > V(local_k+1, distributed_vector( b.n, b.p, b.my_pid, b.local_rows, b.overlap_rows, b.send_instructions, b.recv_instructions ) );
  sparse_matrix_vector_product(A,x,V[0]);

  // r0 = temp1 = b - A*x
  axpy( 1.0, b, -1.0, V[0] );

  err = norm2(V[0]);

  int max_restarts = restarts;
  restarts = -1;
  iters = 0;

  std::vector<double> cosines(local_k);
  std::vector<double> sines(local_k);
  std::vector<double> g(local_k+1);

  while ((err > tol) && (restarts < max_restarts))
  {
    // at the start of every re-start
      // the initial guess is already stored in x

    restarts++;
    iters = 0;
    
    // Normalize the first vector, kernel, of the Krylov subspace.
    scale( V[0], 1.0/err );

    std::vector< std::vector<double> > R;
    g[0] = err;

    // ok, Mr.GMRES has determined the initial values for
      // V,R,g,sines,cosines,err and iters
      // he updates these at every iteration until
        // either err becomes less than tol
        // or a new restart is required

    // note that Mr.GMRES does not think it necessary to store the upper Hessenberg matrix at each iteration
      // he computes R at each iteration from the new Hessenberg matrix column
      // and he's clever enough to determine err without having to solve for x at each iteration

    while ((err > tol) && (iters < local_k))
    {  
      iters++;

      // Compute the next Krylov-vector which will require a matrix vector multiplication

      sparse_matrix_vector_product( A, V[iters-1], V[iters] );

      // Perform orthogonalization, modified Gram-Schmidt
      std::vector<double> new_col_H(iters+1);

      for (int i = 0; i < iters; i++)
      {
        double new_dot = dot( V[iters], V[i] );
        axpy( -1.0*new_dot, V[i], 1.0, V[iters] );      
        new_col_H[i] = new_dot;
      }

      double new_norm = norm2( V[iters] );
      new_col_H[iters] = new_norm;
      scale( V[iters], 1.0/new_norm );
     
      // Update the least squares system to get the new column of R using the current sines and cosines

      for (int i = 0; i < iters-1; i++)
      {
        double old_i = new_col_H[i];
        double old_i_plus_one = new_col_H[i+1];

        new_col_H[i] = cosines[i]*old_i + sines[i]*old_i_plus_one;
        new_col_H[i+1] = -sines[i]*old_i + cosines[i]*old_i_plus_one;
      }

      // Compute next Given's rotation
      rotg( new_col_H[iters-1], new_col_H[iters], cosines[iters-1], sines[iters-1] );

      new_col_H.pop_back();
      R.push_back(new_col_H);

      // Right, the new column of R is ready
      // the only thing left to do is to update g
        // which will also tell Mr.GMRES what the new error is

      double old_g = g[iters-1];
      g[iters-1] = old_g*cosines[iters-1];
      g[iters] = -old_g*sines[iters-1];

      err = std::abs(g[iters]);
    }

    // Solve the least-squares system Ry = g
    // after which compute the next solution x += (V without its last column)*y

    std::vector<double> y(iters);

    for (int i = iters-1; i >= 0; i--)
    {
      double sum = (double)(0);

      for (int j = iters-1; j > i; j--)
      {
        sum += R[j][i]*y[j];
      }

      y[i] = (g[i] - sum)/R[i][i];

      // Update the solution.
      axpy( y[i], V[i], 1.0, x );
    }

    // Compute the actual residual 
    sparse_matrix_vector_product(A,x,V[0]);

    // r0 = temp1 = b - A*x
    axpy( 1.0, b, -1.0, V[0] );

    // r0 / || r_0 ||
    double beta = norm2(V[0]);
    err = beta;

    // the new x is also ready
      // either return it or use it as an initial guess for the next restart
  }
  
  // if Mr.GMRES ever reaches here, it means he's solved the problem

  if (restarts < 0)
  {
    restarts = 0;
  }
}

void mX_solver_utils::rotg(double& a, double& b, double& c, double& s)
{
  double roe = b;
  double abs_a = std::abs(a);
  double abs_b = std::abs(b);
  double abs_roe = abs_b;

  if (abs_a > abs_b)
  {
    roe = a;
    abs_roe = abs_a;
  }
  double scale = abs_a + abs_b;
  
  if (scale == 0.0)
  {
    c = 1.0;
    s = 0.0;
  }
  else
  {
    double roe_sign = roe / abs_roe;
    double r = roe_sign*scale*std::sqrt( (a/scale)*(a/scale) + (b/scale)*(b/scale) );
    c = a/r;
    s = b/r;
    double z = 1.0;
    if (abs_a > abs_b)
      z = s;
    if ((abs_b >= abs_a) && (c != 0.0))
      z = 1.0/c;
    a = r;
    b = z;
  }
} 


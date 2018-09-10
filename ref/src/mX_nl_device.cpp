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

#include "mX_nl_device.h"
#include "mX_DAE.h"
#include <cmath>
#include <vector>

using namespace mX_DAE_utils;

namespace mX_device_utils
{

mX_diode::mX_diode()
: mX_device('D')
{}

int mX_diode::add_device( std::istringstream& input_str, int extra_nodes_ptr, mX_DAE_utils::mX_DAE* dae )
{ 
  int node1, node2;
  input_str >> node1;
  input_str >> node2;
  
  // Don't process any more values, there is an independent function to evaluate.             
  //double dvalue;
  //input_str >> dvalue;
    
  // Insert A entries first. 
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node1-1,node1-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node1-1,extra_nodes_ptr) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node2-1,node2-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node2-1,extra_nodes_ptr) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,extra_nodes_ptr,node1-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,extra_nodes_ptr,node2-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,extra_nodes_ptr,extra_nodes_ptr) );
  
  // Then insert B entres.
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->B,node2-1,node2-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->B,node2-1,extra_nodes_ptr) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->B,extra_nodes_ptr,node2-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->B,extra_nodes_ptr,extra_nodes_ptr) );
  
  return 1;
}

void mX_diode::load_vector( distributed_vector& vec )
{
}

void mX_diode::load_matrices()
{ 
  std::vector<distributed_sparse_matrix_entry* >::iterator d_it =  device_entries_.begin();
  std::vector<distributed_sparse_matrix_entry* >::iterator d_end =  device_entries_.end();
  
  for ( ; d_it != d_end; )
  { 
    // Evaluate model
    std::vector<double> input_vals( 3, std::rand() );
    std::vector<double> model_vals = diode_l1( input_vals );
    
    double Gspr = model_vals[0];
    double Gd = model_vals[1];
    double Cd = model_vals[2];
    
    // A-matrix
    (*d_it)->value += Gspr; d_it++;
    (*d_it)->value -= Gspr; d_it++;
    (*d_it)->value += Gd; d_it++;
    (*d_it)->value -= Gd; d_it++;
    (*d_it)->value -= Gspr; d_it++;
    (*d_it)->value -= Gd; d_it++;
    (*d_it)->value += Gspr + Gd; d_it++;
    
    // B-matrix
    (*d_it)->value += Cd; d_it++;
    (*d_it)->value -= Cd; d_it++;
    (*d_it)->value -= Cd; d_it++;
    (*d_it)->value += Cd; d_it++;
  }
}

std::vector<double> diode_l1(std::vector<double> &x)
{
    double z0 = x[0];
    double z1 = x[1];
    double z2 = x[2];

    double z3 = +z0;
    double z4 = z2 * z1;
    double z5 = std::abs(z4);
    double z6 = std::cos(z5);
    double z7 = std::sin(z3);
    double z8 = z6 + z7;
    double z9 = std::cbrt(z8);
    double z10 = std::tanh(z9);
    double z11 = z10 * z8;
    double z12 = std::tanh(z11);
    double z13 = std::cbrt(z12);
    double z14 = z13 * z10;
    double z15 = z14 - z0;
    double z16 = z15 + z7;
    double z17 = z16 + z14;
    double z18 = std::pow(z17, 2.0);
    double z19 = std::atan(z18);
    double z20 = -z19;
    double z21 = std::sin(z20);
    double z22 = std::cbrt(z21);
    double z23 = std::sin(z22);
    double z24 = z23 + z14;
    double z25 = std::cos(z24);
    double z26 = z25 - z11;
    double z27 = std::atan(z26);
    double z28 = std::pow(z27, 3.0);
    double z29 = z28 + z28;
    double z30 = z29 * z17;
    double z31 = std::atan(z30);
    double z32 = std::tanh(z31);
    double z33 = z32 + z2;
    double z34 = -z33;
    double z35 = z34 - z12;
    double z36 = std::atan(z35);
    double z37 = std::abs(z36);
    double z38 = std::abs(z37);
    double z39 = z38 + z23;
    double z40 = std::cos(z39);
    double z41 = std::cos(z40);
    double z42 = +z41;
    double z43 = z42 - z6;
    double z44 = std::sin(z43);
    double z45 = +z44;
    double z46 = std::pow(z45, 3.0);
    double z47 = z46 * z28;
    double z48 = std::cbrt(z47);
    double z49 = z48 + z30;
    double z50 = std::cos(z49);
    double z51 = std::atan(z50);
    double z52 = std::tanh(z51);

    double z53 = z52;
    double z54 = z52;
    double z55 = z11;

    std::vector<double> y;
    y.push_back(z53);
    y.push_back(z54);
    y.push_back(z55);

    return y;
}

}

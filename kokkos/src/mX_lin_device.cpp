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

#include <iostream>

#include "mX_source.h"
#include "mX_lin_device.h"
#include "mX_DAE.h"

using namespace mX_source_utils;
using namespace mX_DAE_utils;

namespace mX_device_utils
{

//
// class mX_resistor 
//
mX_resistor::mX_resistor()
: mX_device('R')
{}

int mX_resistor::add_device( std::istringstream& input_str, int extra_nodes_ptr, mX_DAE_utils::mX_DAE* dae )
{
  int node1, node2;
  input_str >> node1;
  input_str >> node2;

  double rvalue;
  input_str >> rvalue;

  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node1-1,node1-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node2-1,node2-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node1-1,node2-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node2-1,node1-1) );

  R_.push_back( rvalue ); 

  return 0;
}

void mX_resistor::load_vector( distributed_vector& vec )
{
}

void mX_resistor::load_matrices()
{
  std::vector<double>::iterator r_it = R_.begin();
  std::vector<double>::iterator r_end = R_.end();
  std::vector<distributed_sparse_matrix_entry* >::iterator d_it =  device_entries_.begin();
  std::vector<distributed_sparse_matrix_entry* >::iterator d_end =  device_entries_.end();
   
  for ( ; r_it != r_end && d_it != d_end; r_it++ )
  {
    double rvalue = 1.0/(*r_it);
 
    (*d_it)->value += rvalue; d_it++;
    (*d_it)->value += rvalue; d_it++;
    (*d_it)->value -= rvalue; d_it++;
    (*d_it)->value -= rvalue; d_it++;
  }
}

//
// class mX_inductor
//
mX_inductor::mX_inductor()
: mX_device('L')
{}

int mX_inductor::add_device( std::istringstream& input_str, int extra_nodes_ptr, mX_DAE_utils::mX_DAE* dae )
{
  int node1, node2;
  input_str >> node1;
  input_str >> node2;
  
  double lvalue;
  input_str >> lvalue;
  
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,extra_nodes_ptr,node1-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,extra_nodes_ptr,node2-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->B,extra_nodes_ptr,extra_nodes_ptr) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node1-1,extra_nodes_ptr) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node2-1,extra_nodes_ptr) );

  L_.push_back( lvalue );

  return 1;
}

void mX_inductor::load_vector( distributed_vector& vec ) 
{
}

void mX_inductor::load_matrices() 
{
  std::vector<double>::iterator l_it = L_.begin();
  std::vector<double>::iterator l_end = L_.end();
  std::vector<distributed_sparse_matrix_entry* >::iterator d_it =  device_entries_.begin();
  std::vector<distributed_sparse_matrix_entry* >::iterator d_end =  device_entries_.end();

  for ( ; l_it != l_end && d_it != d_end; l_it++ )
  {
    (*d_it)->value += 1.0; d_it++;
    (*d_it)->value -= 1.0; d_it++;
    (*d_it)->value -= *l_it; d_it++;
    (*d_it)->value += 1.0; d_it++;
    (*d_it)->value -= 1.0; d_it++;
  }
}

//
// class mX_capacitor 
//
mX_capacitor::mX_capacitor()
: mX_device('C')
{}

int mX_capacitor::add_device( std::istringstream& input_str, int extra_nodes_ptr, mX_DAE_utils::mX_DAE* dae )
{
  int node1, node2;
  input_str >> node1;
  input_str >> node2;

  double cvalue;
  input_str >> cvalue;

  device_entries_.push_back( distributed_sparse_matrix_insert(dae->B,node1-1,node1-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->B,node2-1,node2-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->B,node1-1,node2-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->B,node2-1,node1-1) );
  
  C_.push_back( cvalue );

  return 0;
}

void mX_capacitor::load_vector( distributed_vector& vec ) 
{
}

void mX_capacitor::load_matrices()
{
  std::vector<double>::iterator c_it = C_.begin();
  std::vector<double>::iterator c_end = C_.end();
  std::vector<distributed_sparse_matrix_entry* >::iterator d_it =  device_entries_.begin();
  std::vector<distributed_sparse_matrix_entry* >::iterator d_end =  device_entries_.end();
   
  for ( ; c_it != c_end && d_it != d_end; c_it++ )
  { 
    (*d_it)->value += *c_it; d_it++;
    (*d_it)->value += *c_it; d_it++;
    (*d_it)->value -= *c_it; d_it++;
    (*d_it)->value -= *c_it; d_it++;
  }
}

//
// class mX_vsrc
//
mX_vsrc::mX_vsrc()
: mX_device('V')
{}

int mX_vsrc::add_device( std::istringstream& input_str, int extra_nodes_ptr, mX_DAE_utils::mX_DAE* dae )
{ 
  int node1, node2;
  input_str >> node1;
  input_str >> node2;

  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,extra_nodes_ptr,node1-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,extra_nodes_ptr,node2-1) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node1-1,extra_nodes_ptr) );
  device_entries_.push_back( distributed_sparse_matrix_insert(dae->A,node2-1,extra_nodes_ptr) );

  mX_source* src = parse_source(input_str);

  mX_scaled_source* scaled_src = new mX_scaled_source();
  scaled_src->src = src;
  scaled_src->scale = (double)(1);

  if(!(dae->b[extra_nodes_ptr]))
  {
    dae->b[extra_nodes_ptr] = new mX_DAE_RHS_entry();
  }

  (dae->b[extra_nodes_ptr])->scaled_src_list.push_back(scaled_src);

  return 1;
}

void mX_vsrc::load_vector( distributed_vector& vec )
{
}

void mX_vsrc::load_matrices()
{
  std::vector<distributed_sparse_matrix_entry* >::iterator d_it =  device_entries_.begin();
  std::vector<distributed_sparse_matrix_entry* >::iterator d_end =  device_entries_.end();

  for ( ; d_it != d_end; )
  {
    (*d_it)->value += 1.0; d_it++;
    (*d_it)->value -= 1.0; d_it++;
    (*d_it)->value -= 1.0; d_it++;
    (*d_it)->value += 1.0; d_it++;
  }
}

//
// class mX_isrc
//
mX_isrc::mX_isrc()
: mX_device('I')
{}

int mX_isrc::add_device( std::istringstream& input_str, int extra_nodes_ptr, mX_DAE_utils::mX_DAE* dae )
{ 
  int node1, node2;
  input_str >> node1;
  input_str >> node2;

  if (node1)
  {
    mX_source* src = parse_source(input_str);

    mX_scaled_source* scaled_src = new mX_scaled_source();
    scaled_src->src = src;
    scaled_src->scale = (double)(1);

    if(!(dae->b[node1-1]))
    {
      dae->b[node1-1] = new mX_DAE_RHS_entry();
    }

    (dae->b[node1-1])->scaled_src_list.push_back(scaled_src);
  }
  if (node2)
  {
    mX_source* src = parse_source(input_str);

    mX_scaled_source* scaled_src = new mX_scaled_source();
    scaled_src->src = src;
    scaled_src->scale = (double)(-1);

    if(!(dae->b[node2-1]))
    {
      dae->b[node2-1] = new mX_DAE_RHS_entry();
    }

    (dae->b[node2-1])->scaled_src_list.push_back(scaled_src);
  }

  return 0;
}

void mX_isrc::load_vector( distributed_vector& vec )
{
}

void mX_isrc::load_matrices()
{
}


}



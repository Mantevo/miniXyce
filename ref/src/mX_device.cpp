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

#include "mX_device.h"

namespace mX_device_utils
{

//
// class mX_resistor 
//
mX_resistor::mX_resistor()
: mX_device('R')
{}

int mX_resistor::add_device( std::istringstream& input_str )
{
  return 0;
}

void mX_resistor::loadVector( distributed_vector& vec )
{
}

void mX_resistor::loadMatrices( distributed_sparse_matrix& dfdx, distributed_sparse_matrix& dqdx ) 
{
}

//
// class mX_inductor
//
mX_inductor::mX_inductor()
: mX_device('L')
{}

int mX_inductor::add_device( std::istringstream& input_str )
{
  return 1;
}

void mX_inductor::loadVector( distributed_vector& vec ) 
{
}

void mX_inductor::loadMatrices( distributed_sparse_matrix& dfdx, distributed_sparse_matrix& dqdx ) 
{
}

//
// class mX_capacitor 
//
mX_capacitor::mX_capacitor()
: mX_device('C')
{}

int mX_capacitor::add_device( std::istringstream& input_str )
{
  return 0;
}

void mX_capacitor::loadVector( distributed_vector& vec ) 
{
}

void mX_capacitor::loadMatrices( distributed_sparse_matrix& dfdx, distributed_sparse_matrix& dqdx ) 
{
}

//
// Create device container
//
mX_device* create_device( char device_label )
{
  mX_device* new_device = 0;

  switch( device_label )
  {
    case 'R':
    {    
      new_device = new mX_resistor();
      break;
    }
    case 'C':
    {
      new_device = new mX_capacitor();
      break;
    }
    case 'L':
    {
      new_device = new mX_inductor();
      break;
    }
  }

  return new_device;
} 

void add_device_nodes( char device_label, std::istringstream& input_str, std::vector<int>& node_list )
{
  // Each device has at least 2 nodes, some have more.
  // As we add device models, those that have more will be inserted here.
  int node1, node2;
  input_str >> node1;
  input_str >> node2;
  node_list.push_back( node1 );
  node_list.push_back( node2 );
}


int num_internal_nodes( char device_label, std::istringstream& input_str )
{
  int num_nodes = 0;

  // Add internal nodes for KCL based on certain device models.
  switch( device_label )
  {
    case 'V':
    {
      num_nodes = 1;
    }
    break;

    case 'L':
    {
      num_nodes = 1;
    }
    break;
  }

  return num_nodes;
}

}


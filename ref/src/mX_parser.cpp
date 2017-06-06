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

#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include "mX_source.h"
#include "mX_linear_DAE.h"
#include "mX_sparse_matrix.h"
#include "mX_parser.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

using namespace mX_source_utils;
using namespace mX_matrix_utils;
using namespace mX_linear_DAE_utils;
  
mX_linear_DAE* mX_parse_utils::parse_netlist(std::string filename, int p, int pid, int &total_devices, int &total_unknowns, int &num_internal_nodes, std::map<std::string, int>& device_count)
{
  // here's the netlist parser I was boasting about

  // given a linear netlist description stored in a file
        // this function will construct a DAE out of it
  // and return that part of the DAE to be stored with processor pid

  // I assume that all processors will be parsing the netlist simultaneously
  // and each processor will only take the entries that are rightfully his own
  // thou shall not covet your fellow processor's distributed matrix entries

  std::ifstream infile;
  infile.open(filename.data());
  
  // Variables for the device count
  int num_voltage_sources = 0;
  int num_current_sources = 0;
  int num_resistors = 0;
  int num_capacitors = 0;
  int num_inductors = 0;

  // Make first pass over the netlist and count the number of devices
  // and collect the integer node names.
  std::vector<int> node_list;

  // Make first pass and perform device count
  while (!infile.eof())
  {
    std::string curr_line;
    getline(infile,curr_line);
   
    if ((curr_line.length() == 0) || (curr_line[0] == '%'))
    {
      continue;  // comments begin with %
    }

    std::istringstream input_str(curr_line);
    std::string first_word;
    input_str >> first_word;

    switch(first_word[0])
    {
      case 'r':
      case 'R':
      {
        // seen a resistor
        num_resistors++;
      }
      break;

      case 'c':
      case 'C':
      {
        // seen a capacitor
        num_capacitors++;
      }
      break;
           
      case 'l':
      case 'L':
      {
        // seen an inductor
        num_inductors++;
      }
      break;

      case 'v':
      case 'V':
      {
        // seen a voltage source
        num_voltage_sources++;
      } 
      break;

      case 'i':
      case 'I':
      {
        // seen a current source
        num_current_sources++;
      }
      break;
    }

    // Save the node names, which are numbers
    int node1, node2;
    input_str >> node1;
    input_str >> node2;
    node_list.push_back( node1 );
    node_list.push_back( node2 );
  }

  // Sort the node list and count how many nonzero nodes there are.
  std::sort( node_list.begin(), node_list.end() );
  node_list.erase( std::unique( node_list.begin(), node_list.end() ), node_list.end() );
  num_internal_nodes = node_list.size();
  if (node_list[0] == 0)
    num_internal_nodes--;  // Get rid of ground node, which is '0', from the node count.

  total_unknowns = num_internal_nodes + num_voltage_sources + num_inductors;
  total_devices = num_voltage_sources + num_current_sources + num_resistors + num_capacitors + num_inductors;

  if (num_voltage_sources)
    device_count["V"] = num_voltage_sources;
  if (num_current_sources)
    device_count["I"] = num_current_sources;
  if (num_resistors)
    device_count["R"] = num_resistors;
  if (num_capacitors)
    device_count["C"] = num_capacitors;
  if (num_inductors)
    device_count["L"] = num_inductors;

  // Closing and reopening the file, should be able to rewind the file, but it didn't work
  infile.close();
  infile.open(filename.data());
  //infile.seekg( 0, infile.beg );

  // Variables for the device count
  int voltage_src_number = 0;
  int inductor_number = 0; 

  int start_row = (total_unknowns/p)*(pid) + ((pid < total_unknowns%p) ? pid : total_unknowns%p);
  int end_row = start_row + (total_unknowns/p) - 1 + ((pid < total_unknowns%p) ? 1 : 0);

  mX_linear_DAE* dae = new mX_linear_DAE();
  dae->A = new distributed_sparse_matrix();
  dae->B = new distributed_sparse_matrix();

  distributed_sparse_matrix* A = dae->A;
  distributed_sparse_matrix* B = dae->B;

  // initialise the A, B and b parts that need to be stored in this processor

  A->start_row = start_row;
  A->end_row = end_row;
  B->start_row = start_row;
  B->end_row = end_row;

  for (int i = start_row; i <= end_row; i++)
  {
    distributed_sparse_matrix_entry* null_ptr_1 = 0;
    A->row_headers.push_back(null_ptr_1);
    B->row_headers.push_back(null_ptr_1);

    mX_linear_DAE_RHS_entry* null_ptr_2 = 0;
    (dae->b).push_back(null_ptr_2);
  }

  // initialisation done
  // now read each line from the input file and do what is expected

  while (!infile.eof())
  {
    std::string curr_line;
    getline(infile,curr_line);

    if ((curr_line.length() == 0) || (curr_line[0] == '%'))
    {
      continue;  // comments begin with %
    }

    std::istringstream input_str(curr_line);
    std::string first_word;
    input_str >> first_word;

    int node1, node2;
    input_str >> node1;
    input_str >> node2;

    switch(first_word[0])
    {
      case 'r':
      case 'R':

      {
        double rvalue;
        input_str >> rvalue;
        
        distributed_sparse_matrix_add_to(A,node1-1,node1-1,(double)(1)/rvalue,total_unknowns,p);
        distributed_sparse_matrix_add_to(A,node2-1,node2-1,(double)(1)/rvalue,total_unknowns,p);
        distributed_sparse_matrix_add_to(A,node1-1,node2-1,(double)(-1)/rvalue,total_unknowns,p);
        distributed_sparse_matrix_add_to(A,node2-1,node1-1,(double)(-1)/rvalue,total_unknowns,p);
      }

      break;

      case 'c':
      case 'C':

      {
        double cvalue;
        input_str >> cvalue;
        
        distributed_sparse_matrix_add_to(B,node1-1,node1-1,cvalue,total_unknowns,p);
        distributed_sparse_matrix_add_to(B,node2-1,node2-1,cvalue,total_unknowns,p);
        distributed_sparse_matrix_add_to(B,node1-1,node2-1,-cvalue,total_unknowns,p);
        distributed_sparse_matrix_add_to(B,node2-1,node1-1,-cvalue,total_unknowns,p);
      }

      break;

      case 'l':
      case 'L':

      {
        // seen an inductor
        
        inductor_number++;
        
        int k = num_internal_nodes + num_voltage_sources + inductor_number;

        double lvalue;
        input_str >> lvalue;
        
        distributed_sparse_matrix_add_to(A,k-1,node1-1,(double)(1),total_unknowns,p);
        distributed_sparse_matrix_add_to(A,k-1,node2-1,(double)(-1),total_unknowns,p);
        distributed_sparse_matrix_add_to(B,k-1,k-1,-lvalue,total_unknowns,p);
        distributed_sparse_matrix_add_to(A,node1-1,k-1,(double)(1),total_unknowns,p);
        distributed_sparse_matrix_add_to(A,node2-1,k-1,(double)(-1),total_unknowns,p);
      }

      break;

      case 'v':
      case 'V':

      {
        // seen a voltage source 
        
        voltage_src_number++;
        
        int k = num_internal_nodes + voltage_src_number;

        distributed_sparse_matrix_add_to(A,k-1,node1-1,(double)(1),total_unknowns,p);
        distributed_sparse_matrix_add_to(A,k-1,node2-1,(double)(-1),total_unknowns,p);
        distributed_sparse_matrix_add_to(A,node1-1,k-1,(double)(-1),total_unknowns,p);
        distributed_sparse_matrix_add_to(A,node2-1,k-1,(double)(1),total_unknowns,p);
        
        if ((k-1 >= start_row) && (k-1 <= end_row))
        {
          mX_source* src = parse_source(input_str);
        
          mX_scaled_source* scaled_src = new mX_scaled_source();
          scaled_src->src = src;
          scaled_src->scale = (double)(1);

          if(!(dae->b[k-1-start_row]))
          {
            dae->b[k-1-start_row] = new mX_linear_DAE_RHS_entry();
          }

          (dae->b[k-1-start_row])->scaled_src_list.push_back(scaled_src);

        }

      }

      break;

      case 'i':
      case 'I':

      {
        if ((node1-1 >= start_row) && (node1-1 <= end_row))
        {
          mX_source* src = parse_source(input_str);
        
          mX_scaled_source* scaled_src = new mX_scaled_source();
          scaled_src->src = src;
          scaled_src->scale = (double)(1);

          if(!(dae->b[node1-1-start_row]))
          {
            dae->b[node1-1-start_row] = new mX_linear_DAE_RHS_entry();
          }

          (dae->b[node1-1-start_row])->scaled_src_list.push_back(scaled_src);

          if ((node2-1 >= start_row) && (node2-1 <= end_row))
          {
            mX_scaled_source* scaled_src_2 = new mX_scaled_source();
            scaled_src_2->src = src;
            scaled_src_2->scale = (double)(-1);

            if(!(dae->b[node2-1-start_row]))
            {
              dae->b[node2-1-start_row] = new mX_linear_DAE_RHS_entry();
            }

            (dae->b[node2-1-start_row])->scaled_src_list.push_back(scaled_src_2);

          }
        }

        else
        {
          if ((node2-1 >= start_row) && (node2-1 <= end_row))
          {
            mX_source* src = parse_source(input_str);
          
            mX_scaled_source* scaled_src = new mX_scaled_source();
            scaled_src->src = src;
            scaled_src->scale = (double)(-1);

            if(!(dae->b[node2-1-start_row]))
            {
              dae->b[node2-1-start_row] = new mX_linear_DAE_RHS_entry();
            }

            (dae->b[node2-1-start_row])->scaled_src_list.push_back(scaled_src);
            
          }
        }
      }

      break;
    }
  }

  // whew! all the lines have been read
    // and each processor hopefully has his correct share of the DAE

  infile.close();

  return dae;
}

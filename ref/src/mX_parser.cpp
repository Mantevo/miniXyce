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

#include <cctype>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include "mX_source.h"
#include "mX_linear_DAE.h"
#include "mX_sparse_matrix.h"
#include "mX_parser.h"
#include "mX_device.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

using namespace mX_source_utils;
using namespace mX_matrix_utils;
using namespace mX_linear_DAE_utils;
using namespace mX_device_utils;
  
mX_linear_DAE* mX_parse_utils::parse_netlist(std::string filename, int p, int pid, int &total_devices, int &total_unknowns, int &num_external_nodes, std::map<char, int>& device_count)
{
  std::ifstream infile;
  infile.open(filename.data());
  
  // Make first pass over the netlist and count the number of devices in device_count
  // and collect the integer node names in node_list.
  total_devices = 0;
  total_unknowns = 0;
  device_count.clear();
  std::vector<int> node_list;

  // Make first pass and perform device count
  while (!infile.eof())
  {
    std::string curr_line;
    getline(infile,curr_line);
   
    if ((curr_line.length() == 0) || (curr_line[0] == '*') || (curr_line[0] == '.'))
    {
      continue;  // comments begin with *, ignore '.' lines 
    }

    std::istringstream input_str(curr_line);
    std::string first_word;
    input_str >> first_word;

    // Convert to upper-case.
    std::transform(first_word.begin(), first_word.end(), first_word.begin(),
                   [](unsigned char c) -> unsigned char { return std::toupper(c); });

    // Increment device count for this device type.
    device_count[first_word[0]]++;
    total_devices++;

    // Save the node names, which are numbers
    mX_device_utils::add_device_nodes( first_word[0], input_str, node_list );

    // Get the number of internal nodes for this device to add to the number of unknowns.
    total_unknowns += mX_device_utils::num_internal_nodes( first_word[0], input_str ); 
  }

  // Sort the node list and count how many nonzero nodes there are.
  std::sort( node_list.begin(), node_list.end() );
  node_list.erase( std::unique( node_list.begin(), node_list.end() ), node_list.end() );
  num_external_nodes = node_list.size();
  if (node_list[0] == 0)
    num_external_nodes--;  // Get rid of ground node, which is '0', from the node count.

  total_unknowns += num_external_nodes;

  // Closing and reopening the file, should be able to rewind the file, but it didn't work
  infile.close();
  infile.open(filename.data());
  //infile.seekg( 0, infile.beg );

  // Variables for the device count
  int voltage_src_number = 0;
  int inductor_number = 0; 

  // Determine how many devices each processor should get.
  int num_my_devices = total_devices/p;
  if ( p > 1 && pid < total_devices%p )
    num_my_devices++;

  std::vector<int> local_dev_count( p, 0 ), global_dev_count( p, 0 );
  local_dev_count[ pid ] = num_my_devices;
  
#ifdef HAVE_MPI
  MPI_Allreduce(&local_dev_count[0],&global_dev_count[0],p,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif

  mX_linear_DAE* dae = new mX_linear_DAE();
  dae->A = new distributed_sparse_matrix();
  dae->B = new distributed_sparse_matrix();

  distributed_sparse_matrix* A = dae->A;
  distributed_sparse_matrix* B = dae->B;

  // initialise the A, B and b parts that need to be stored in this processor
  A->n = total_unknowns;
  A->p = p;
  A->my_pid = pid;

  B->n = total_unknowns; 
  B->p = p;
  B->my_pid = pid;

  for (int i = 0; i < total_unknowns; i++)
  {
    distributed_sparse_matrix_entry* null_ptr_1 = 0;
    A->row_headers.push_back(null_ptr_1);
    B->row_headers.push_back(null_ptr_1);

    mX_linear_DAE_RHS_entry* null_ptr_2 = 0;
    (dae->b).push_back(null_ptr_2);
  }

  // In parallel, processor 0 communicates devices to all other processors
  // before taking its share.
  int current_pid = 0;
  int extra_nodes_ptr = num_external_nodes;

  if (p > 1)
    current_pid = 1;   

  // initialisation done
  // now read each line from the input file and do what is expected
  while (!infile.eof())
  {
    std::string curr_line;
    getline(infile,curr_line);

    if ((curr_line.length() == 0) || (curr_line[0] == '*') || (curr_line[0] == '.'))
    {
      continue;  // comments begin with * 
    }

    std::istringstream input_str(curr_line);
    std::string first_word;
    input_str >> first_word;

    std::transform(first_word.begin(), first_word.end(), first_word.begin(),
                   [](unsigned char c) -> unsigned char { return std::toupper(c); });

    if ( first_word[0] == 'V' )
        voltage_src_number++;
    if ( first_word[0] == 'L' )
        inductor_number++;

    if (pid == current_pid)
    {
      // Get the container for this device
      std::map<char, mX_device*>::iterator itr = dae->devices.find(first_word[0]);
      mX_device* currDevice = 0;
      if (itr == dae->devices.end())
      {
        std::cout << "Processor " << current_pid << " didn't find container for device " << first_word[0] << std::endl;
        currDevice = create_device( first_word[0] );
        dae->devices[first_word[0]] = currDevice;
      }
      else
      {
        currDevice = itr->second;
      }

      // If the device generates internal nodes, increment the extra_nodes_ptr.
      extra_nodes_ptr += currDevice->add_device( input_str, extra_nodes_ptr, dae );
    }

    // Update device count and determine if next processor is to receive new device.
    global_dev_count[ current_pid ]--;
    if ( !global_dev_count[ current_pid ] )
    {
      current_pid++;
      // If this is the last processor, change current_pid to processor 0
      if (current_pid == p)
        current_pid = 0;

#ifdef HAVE_MPI
      std::cout << "Processor " << current_pid << " extra_nodes_ptr = " << extra_nodes_ptr << std::endl;
      int my_extra_nodes_ptr = extra_nodes_ptr;
      MPI_Allreduce(&my_extra_nodes_ptr,&extra_nodes_ptr,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);      
#endif
    }
  }

  // Have each device type load their entries now.
  load_matrices( dae );

  // whew! all the lines have been read
    // and each processor hopefully has his correct share of the DAE
  infile.close();

  // Now tell the matrices to assemble themselves and determine off-processor communication.
  distributed_sparse_matrix_finish( A, B, dae->b );

  return dae;
}

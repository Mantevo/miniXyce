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

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include "mX_comm.h"
#include "mX_parser.h"
#include "mX_source.h"
#include "mX_sparse_matrix.h"
#include "mX_solver.h"
#include "mX_DAE.h"
#include "mX_parms.h"
#include "mX_timer.h"
#include "YAML_Element.hpp"
#include "YAML_Doc.hpp"

#include <Kokkos_Core.hpp>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

using namespace mX_comm_utils;
using namespace mX_parse_utils;
using namespace mX_source_utils;
using namespace mX_DAE_utils;
using namespace mX_parms_utils;
using namespace mX_solver_utils;

int main(int argc, char* argv[])
{
  // this is of course, the actual transient simulator
  int p=1, pid=0;
#ifdef HAVE_MPI  
  MPI_Init(&argc,&argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#endif
  double sim_start = mX_timer();

  Kokkos::initialize (argc, argv);

  printf ("Hello World on Kokkos execution space %s\n",
          typeid (Kokkos::DefaultExecutionSpace).name ());

  // initialize YAML doc
  YAML_Doc doc("miniXyce","1.1");

  // initialize the simulation parameters

  std::string ckt_netlist_filename;
  double t_start, t_step, t_stop, tol, res;
  int k=0, iters=0, restarts=0;

  std::vector<double> x;
  bool init_cond_specified;

  double tstart = mX_timer();
  get_parms(argc,argv,ckt_netlist_filename,t_start,t_step,t_stop,tol,k,restarts,x,init_cond_specified,p,pid);
  double tend = mX_timer() - tstart;
  doc.add("Parameter_parsing_time",tend);

  // build the DAE from the circuit netlist

  tstart = mX_timer();
  int total_devices, total_unknowns, num_external_nodes;
  std::map<char, int> device_count;

  mX_DAE* dae = parse_netlist(ckt_netlist_filename,p,pid,total_devices,total_unknowns,num_external_nodes,device_count);

  tend = mX_timer() - tstart;
  doc.add("Netlist_parsing_time",tend);

  // document circuit and matrix attributes
 
  doc.add("Netlist_file",ckt_netlist_filename.c_str());

  doc.add("Circuit_attributes","");
  doc.get("Circuit_attributes")->add("Number_of_devices",total_devices);
  if (device_count.find('R') != device_count.end())
  {
    int num_resistors = device_count['R'];
    doc.get("Circuit_attributes")->add("Resistors_(R)",num_resistors);
  }
  if (device_count.find('L') != device_count.end())
  {
    int num_inductors = device_count['L'];
    doc.get("Circuit_attributes")->add("Inductors_(L)",num_inductors);
  }
  if (device_count.find('C') != device_count.end())
  {
    int num_capacitors = device_count['C'];
    doc.get("Circuit_attributes")->add("Capacitors_(C)",num_capacitors);
  }
  if (device_count.find('D') != device_count.end())
  {
    int num_diodes = device_count['D'];
    doc.get("Circuit_attributes")->add("Diodes_(D)",num_diodes);
  if (device_count.find('V') != device_count.end())
  {
    int num_voltage_sources = device_count['V'];
    doc.get("Circuit_attributes")->add("Voltage_sources_(V)",num_voltage_sources);
  }
  if (device_count.find('I') != device_count.end())
  {
    int num_current_sources = device_count['I'];
    doc.get("Circuit_attributes")->add("Current_sources_(I)",num_current_sources);
  }
  }

  int num_my_nnz = dae->A->local_nnz, sum_nnz = dae->A->local_nnz;
  int min_nnz = num_my_nnz, max_nnz = num_my_nnz;
  int num_my_rows = dae->A->local_rows.size();
  int min_rows = num_my_rows, max_rows = num_my_rows, sum_rows = num_my_rows;
  int num_my_overlap = dae->A->overlap_rows.size();
  int min_overlap = num_my_overlap, max_overlap = num_my_overlap, sum_overlap = num_my_overlap;

#ifdef HAVE_MPI
  MPI_Allreduce(&num_my_nnz,&sum_nnz,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&num_my_nnz,&min_nnz,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&num_my_nnz,&max_nnz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&num_my_rows,&sum_rows,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&num_my_rows,&min_rows,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&num_my_rows,&max_rows,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&num_my_overlap,&sum_overlap,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&num_my_overlap,&min_overlap,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&num_my_overlap,&max_overlap,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  doc.add("Matrix_attributes","");
  doc.get("Matrix_attributes")->add("Global_rows",total_unknowns);
  doc.get("Matrix_attributes")->add("Rows_per_proc_MIN",min_rows);
  doc.get("Matrix_attributes")->add("Rows_per_proc_MAX",max_rows);
  doc.get("Matrix_attributes")->add("Rows_per_proc_AVG",(double)sum_rows/p);
  doc.get("Matrix_attributes")->add("Global_overlap_rows",sum_overlap);
  doc.get("Matrix_attributes")->add("Rows_overlap_per_proc_MIN",min_overlap);
  doc.get("Matrix_attributes")->add("Rows_overlap_per_proc_MAX",max_overlap);
  doc.get("Matrix_attributes")->add("Rows_overlap_per_proc_AVG",(double)sum_overlap/p);
  doc.get("Matrix_attributes")->add("Global_NNZ",sum_nnz);
  doc.get("Matrix_attributes")->add("NNZ_per_proc_MIN",min_nnz);
  doc.get("Matrix_attributes")->add("NNZ_per_proc_MAX",max_nnz);
  doc.get("Matrix_attributes")->add("NNZ_per_proc_AVG",(double)sum_nnz/p);
        
  distributed_sparse_matrix* A = dae->A;
  distributed_sparse_matrix* B = dae->B;

  /*
  for (int proc = 0; proc < p; proc++)
  {
    if (pid == proc)
    {
      int num_my_sends = A->send_instructions.size();
      std::cout << "Processor : " << pid << " has " << num_my_sends << " sends: ";
      std::list<data_transfer_instruction*>::iterator it1;
      for (it1 = A->send_instructions.begin(); it1 != A->send_instructions.end(); it1++)
      {
        std::cout << " ( " << (*it1)->pid << " ) " << (*it1)->indices.size() << " ";
      }
      std::cout << std::endl;
      int num_my_recvs = A->recv_instructions.size();
      std::cout << "Processor : " << pid << " has " << num_my_recvs << " recvs: ";
      for (it1 = A->recv_instructions.begin(); it1 != A->recv_instructions.end(); it1++)
      {
        std::cout << " ( " << (*it1)->pid << " ) " << (*it1)->indices.size() << " ";
      }
      std::cout << std::endl;
    }
  }

  if (pid==0)
    std::cout << "A: " << std::endl;
  print_matrix( *A );

  if (pid==0)
    std::cout << "B: " << std::endl;
  print_matrix( *B );
  */

  distributed_vector init_RHS2 = evaluate_b(t_start,dae);

  // compute the initial condition if not specified by user

  tstart = mX_timer();
        
  distributed_vector x2( total_unknowns, p, pid, A->local_rows, A->overlap_rows, A->send_instructions, A->recv_instructions );

  if (pid==0)
    std::cout << "Initial condition specified: " << init_cond_specified << std::endl;

  if (!init_cond_specified)
  {
    mX_vector_utils::init_value( x2 , 0.0 );

    int curr_restarts = restarts; 
    gmres(A,init_RHS2,x2,tol,res,k,iters,curr_restarts);

    doc.add("DCOP Calculation","");
    doc.get("DCOP Calculation")->add("Init_cond_specified", false);
    doc.get("DCOP Calculation")->add("GMRES_tolerance",tol);
    doc.get("DCOP Calculation")->add("GMRES_subspace_dim",k);
    doc.get("DCOP Calculation")->add("GMRES_iterations",iters);
    doc.get("DCOP Calculation")->add("GMRES_restarts",curr_restarts);
    doc.get("DCOP Calculation")->add("GMRES_native_residual",res);
  }
  else 
  {
    doc.add("DCOP Calculation","");
    doc.get("DCOP Calculation")->add("Init_cond_specified", true);
  }
  tend = mX_timer() - tstart;
  doc.get("DCOP Calculation")->add("DCOP_calculation_time",tend);

  // from now you won't be needing any more Ax = b solves
    // but you will be needing many (A + B/t_step)x = b solves
    // compute the nonzero entry pattern for (A + B/t_step) right now
      // so you won't have to compute it at each time step
 
  tstart = mX_timer();
  distributed_sparse_matrix ApB( total_unknowns, p, pid );

  std::vector<distributed_sparse_matrix_entry*>::iterator it1;
  int row_idx = -1;

  // First insert A.
  for (it1 = A->row_headers.begin(); it1 != A->row_headers.end(); it1++)
  {
    row_idx++;
    distributed_sparse_matrix_entry* curr = *it1;

    while (curr)
    {
      int col_idx = curr->column;
      double value = (curr->value);
      distributed_sparse_matrix_add_to(&ApB,row_idx,col_idx,value);

      curr = curr->next_in_row;
    }
  }

  // Then insert B.
  row_idx = -1;
  for (it1 = B->row_headers.begin(); it1 != B->row_headers.end(); it1++)
  {
    row_idx++;
    distributed_sparse_matrix_entry* curr = *it1;

    while (curr)
    {
      int col_idx = curr->column;
      double value = (curr->value);
      distributed_sparse_matrix_add_to(&ApB,row_idx,col_idx,value);

      curr = curr->next_in_row;
    }
  }

  double matrix_setup_tend = mX_timer() - tstart;

  // even though you use a constant timestep, in practice that time step
    // is variable, so the combination of A + B/t_step will be computed
    // for each timestep to get the correct performance numbers. 

  tstart = mX_timer();
  double t = t_start + t_step;
  double total_res_load = 0.0;
  double total_jac_load = 0.0;
  double total_jac_assembly = 0.0;
  double inner_tstart = 0.0;

  int trans_steps = 0;

  while (t < t_stop)
  {
    trans_steps++;

    // --------------------------------------------------
    // new time point t => new value for b(t)
    // --------------------------------------------------

    inner_tstart = mX_timer();

    distributed_vector b = evaluate_b(t,dae);

    total_res_load += (mX_timer() - inner_tstart);

    // --------------------------------------------------
    // load matrices to assemble the Jacobian
    // NOTE:  this cannot be loaded directly into ApB
    //        because the devices have direct access to
    //        the memory in A and B to load the values.
    // --------------------------------------------------

    inner_tstart = mX_timer();

    load_matrices( dae );

    total_jac_load += (mX_timer() - inner_tstart);

    // --------------------------------------------------
    // assemble jacobian
    // --------------------------------------------------

    inner_tstart = mX_timer();
    init_value(ApB, 0.0);

    // first insert A.
    for (it1 = A->row_headers.begin(), row_idx=0; it1 != A->row_headers.end(); it1++, row_idx++)
    {
      distributed_sparse_matrix_entry* curr = *it1;

      while (curr)
      {
        int col_idx = curr->column;
        double value = (curr->value);
        distributed_sparse_matrix_add_to(&ApB,row_idx,col_idx,value);

        curr = curr->next_in_row;
      }
    }

    // then insert B.
    for (it1 = B->row_headers.begin(), row_idx=0; it1 != B->row_headers.end(); it1++, row_idx++)
    {
      distributed_sparse_matrix_entry* curr = *it1;

      while (curr)
      {
        int col_idx = curr->column;
        double value = (curr->value)/t_step;
        distributed_sparse_matrix_add_to(&ApB,row_idx,col_idx,value);
 
        curr = curr->next_in_row;
      }
    }
   
    total_jac_assembly += (mX_timer() - inner_tstart);

    // increment t
    t += t_step;
  }

  // Document the accrued time for the transient device loads.
  // No solves are being performed in transient.
  doc.add("Transient Calculation","");
  doc.get("Transient Calculation")->add("Number_of_time_steps", trans_steps);
  doc.get("Transient Calculation")->add("Time_start", t_start);
  doc.get("Transient Calculation")->add("Time_end", t_stop);
  doc.get("Transient Calculation")->add("Time_step", t_step);
  doc.get("Transient Calculation")->add("Matrix_setup_time",matrix_setup_tend);
  doc.get("Transient Calculation")->add("Residual_load_time",total_res_load);
  doc.get("Transient Calculation")->add("Jacobian_load_time",total_jac_load);
  doc.get("Transient Calculation")->add("Jacobian_assembly_time",total_jac_assembly);

  if (pid==0) { // Only PE 0 needs to compute and report timing results

    std::string yaml = doc.generateYAML();
    std::cout << yaml;
  }

  // Clean up
  mX_DAE_utils::destroy( dae );

  // You must call finalize() after you are done using Kokkos.
  Kokkos::finalize ();
 
#ifdef HAVE_MPI  
  MPI_Finalize();
#endif

  return 0;
}

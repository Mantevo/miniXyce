# miniXyce Linear Circuit Simulator

At this time, miniXyce is a simple linear circuit simulator with a
basic parser that performs transient analysis on any circuit with
resistors (R), inductors (L), capacitors (C), and voltage/current
sources.  The parser incorporated into this version of miniXyce is a
single pass parser, where the netlist is expected to be flat
(no hierarchy via subcircuits is allowed). Simulating the system of
DAEs generates a nonsymmetric linear problem, which is solved using
un-preconditioned GMRES. The time integration method used in miniXyce
is backward Euler with a constant time-step.  The simulator outputs
all the solution variables at each time step in a 'prn' file.

The development of the first version of miniXyce resulted in something
closer to a compact application than a miniapp since more focus was put
on the simulator returning the correct answer, than modeling performance
characteristics of interest. Further analysis of Xyce has called out
particular performance issues in each of the three phases of circuit
simulation:  parsing, device evaluation/loading, and the solution of linear
equations.  These issues will inspire enhancements to, and a second version
of, miniXyce.

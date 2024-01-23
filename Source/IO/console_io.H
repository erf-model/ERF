#ifndef CONSOLE_IO_H
#define CONSOLE_IO_H

#include <iostream>
#include "AMReX.H"

//namespace ERF::io {

void print_banner(MPI_Comm /*comm*/, std::ostream& /*out*/);

void print_usage(MPI_Comm /*comm*/, std::ostream& /*out*/);

void print_error(MPI_Comm /*comm*/, const std::string& msg);

void print_summary(std::ostream&);

void print_tpls(std::ostream& /*out*/);

//} // namespace ERF::io

#endif /* CONSOLE_IO_H */

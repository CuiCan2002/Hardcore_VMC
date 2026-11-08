#include "vmc.hpp"
#include <exception>
#include <iostream>
#ifdef VMC_MPI
#include <mpi.h>
#endif
int main(int argc,char** argv) {
#ifdef VMC_MPI
    MPI_Init(&argc,&argv);
#endif
    int status=0;
    if(argc!=2) {std::cerr<<"Usage: hardcore_vmc model.inp\n";status=2;}
    else {
    try {vmc::run(vmc::Config::read(argv[1]));}
    catch(const std::exception& e) {std::cerr<<"Error: "<<e.what()<<'\n';status=1;}
    }
#ifdef VMC_MPI
    // Fail the entire communicator rather than strand peers in a collective.
    if(status)MPI_Abort(MPI_COMM_WORLD,status);
    MPI_Finalize();
#endif
    return status;
}

#include "vmc.hpp"
#include <exception>
#include <iostream>
int main(int argc,char** argv) {
    if(argc!=2) {std::cerr<<"Usage: hardcore_vmc model.inp\n";return 2;}
    try {vmc::run(vmc::Config::read(argv[1]));}
    catch(const std::exception& e) {std::cerr<<"Error: "<<e.what()<<'\n';return 1;}
}

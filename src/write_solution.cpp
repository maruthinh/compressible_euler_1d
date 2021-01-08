#include "global_declarations.h"
#include "basic_functions.h"

void Write_Solution(int ib, double *&x, double **&cv, double *&p){

    std::string oup;
    oup = "../output/";
	std::ofstream other(oup+"Res_"+test_case.c_str());
	other.flags(std::ios::dec | std::ios::scientific);
	other.precision(16);

	
	if(!other){
		std::cerr<<"File couldn't be opened to write the solution"<<std::endl;
		exit(1);
	}

	//other << "TITLE = Flow" << endl << "VARIABLES = xc, rho, u, p" << endl;
	//other << "Zone T = Omega I = " << N <<  endl ;


	for (int i=2;i<=ib;i++){
		other<<x[i]<<"\t"<<cv[0][i]<<"\t"<<cv[1][i]/cv[0][i]<<"\t"<<p[i]<<"\t"<<(p[i]/(cv[0][i]*(Gamma-1.0)))<<std::endl;

	}

	other.close();
}

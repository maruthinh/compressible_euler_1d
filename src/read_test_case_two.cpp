#include "global_declarations.h"

void Read_Test_Case_Two(){

	std::string skip;
	std::string inp;
	inp = "../input/";
	std::ifstream infile(inp+test_case.c_str());

	if(infile.fail()){
		std::cerr<<"test case file couldn't be opened to read the input parameters"<<std::endl;
		exit(1);
	}

/****to read grid information****/
	getline(infile, skip);
	getline(infile, skip);
	infile>>xll;
	infile.ignore(80,'\n');
	infile>>xul;
	infile.ignore(80,'\n');
	infile>>x01;
	infile.ignore(80,'\n');
	infile>>x02;
	infile.ignore(80,'\n');
	getline(infile, skip);	
	
/****to read initial conditions****/
	infile>>Prim_var_L[0];
	infile.ignore(80,'\n');
	infile>>Prim_var_L[1];
	infile.ignore(80,'\n');
	infile>>Prim_var_L[2];
	infile.ignore(80,'\n');
	infile>>Prim_var_M[0];
	infile.ignore(80,'\n');
	infile>>Prim_var_M[1];
	infile.ignore(80,'\n');
	infile>>Prim_var_M[2];
	infile.ignore(80,'\n');
	infile>>Prim_var_R[0];
	infile.ignore(80,'\n');
	infile>>Prim_var_R[1];
	infile.ignore(80,'\n');
	infile>>Prim_var_R[2];
	infile.ignore(80,'\n');

/****to read final time at which computations will be stopped****/
	infile>>t_fin;
	infile.ignore(80,'\n');

	std::cout<<"-------------------------------------------------------------"<<std::endl   
		 <<"the test case details..............."<<std::endl
		 <<"-------------------------------------------------------------"<<std::endl  
		 <<"lower and upper limit of grid="<<xll<<"\t"<<xul<<std::endl
		 <<"position of the initial discontinuity="<<x0<<std::endl 
		 <<"left and right states of den, vel, pres="<<Prim_var_L[0]<<"\t"<<Prim_var_L[1]<<"\t"<<Prim_var_L[2]<<"\t"
          	 <<Prim_var_M[0]<<"\t"<<Prim_var_M[1]<<"\t"<<Prim_var_M[2]<<"\t"
	         <<Prim_var_R[0]<<"\t"<<Prim_var_R[1]<<"\t"<<Prim_var_R[2]<<std::endl
		 <<"final time at which computations will be stopped="<<t_fin<<std::endl;
}


#include "global_declarations.h"

void Read_Input(){

	std::string skip;
	std::ifstream infile;

	infile.open("../input/input_file.txt");

	if(infile.fail()){
		std::cerr<<"input file couldn't be opened to read the input parameters"<<std::endl;
		exit(1);
	}

/****read number of grid points required****/
	getline(infile, skip);
	infile>>N;
	infile.ignore(80,'\n');

/****read cfl number****/
	getline(infile, skip);
	infile>>cfl;
	infile.ignore(80,'\n');

/****spatial accuracy****/
	getline(infile, skip);
	infile>>space_accuracy;
	infile.ignore(80,'\n');	
	
/****type of flux used for computations****/
	getline(infile, skip);
	infile>>flux_type;
	infile.ignore(80,'\n');

/****time accuracy****/
	getline(infile, skip);
	infile>>time_accuracy;
	infile.ignore(80,'\n');

/****test case to be solved****/
	getline(infile, skip);
	infile>>test_case;

	ib=N+2; id1=N+3, id2=N+4, imax=N+5;
	
	std::cout<<"-------------------------------------------------------------"<<std::endl   
		 <<"the input parameters read from the input file as follows"<<std::endl
		 <<"-------------------------------------------------------------"<<std::endl   
		 <<"Number of grid pointd required, N="<<N<<std::endl
              	 <<"cfl="<<cfl<<std::endl
		 <<"spatial accuracy desired="<<space_accuracy<<std::endl
		 <<"flux type to be used for computations="<<flux_type<<std::endl
		 <<"method of time accuracy desired="<<time_accuracy<<std::endl
		 <<"test case to be solved="<<test_case<<std::endl;	

}

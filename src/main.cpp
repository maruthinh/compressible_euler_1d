#include "global_declarations.h"
#include "basic_functions.h"

int main(){

	Read_Input();
	
	Declaration_1d_array();
	Declaration_2d_array();

	if(test_case=="test_case10.txt"){


		Read_Test_Case_Two();
	}
	else{

		Read_Test_Case();

	}
	
	x=Grid(xll, xul, N, ib);
	
	if(test_case=="test_case10.txt"){

		Initial_Conditions_Two(x01, x02, ib, imax, x, Prim_var_L, Prim_var_M, Prim_var_R, cv, p);
	}
	else if(test_case=="test_case11.txt" or test_case=="test_case12.txt"){
		Initial_Conditions_Oscil(x0, ib, imax, x, Prim_var_L, Prim_var_R, cv, p);
	}
	else{
		Initial_Conditions(x0, ib, imax, x, Prim_var_L, Prim_var_R, cv, p);
	}

	Boundary_Condition(2, ib, cv, p);

	double t=0;

	while(t<t_fin){

		if(time_accuracy=="explicit"){
			Solver_Explicit(ib, id1, id2, x, cv, cv_old, p, L_state, R_state, residue, dt);
		}
		else if(time_accuracy=="ssprk3"){
            Solver(ib, id1, id2, x, cv, cv_old, p, L_state, R_state, residue, dt);
		}	
		else{
			std::cout<<"wrong input on time accuracy, pls check input file"<<std::endl;
		}

		t=t+dt;
	
		std::cout<<"current time step"<<"\t"<<dt<<"\t"<<"total time"<<"\t"<<t<<std::endl;
	}
	
	std::cout<<"----------------------------------------------------------------"<<std::endl
	         <<"iterations are completed"<<std::endl
             <<"Writing solutions to the file........."<<std::endl
             <<"----------------------------------------------------------------"<<std::endl;

	Write_Solution(ib, x, cv, p);
	
	return 0;
}

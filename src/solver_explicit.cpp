#include "global_declarations.h"
#include "basic_functions.h"

/****time update is done here****/

void Solver_Explicit(int ib, int id1, int id2, double *&x, double **&cv, double **&cv_old, double *&p, double **&L_state, double **&R_state, double **&residue, double &dt){


	for(int i=1; i<=imax; i++){

		cv_old[0][i] = cv[0][i];
		cv_old[1][i] = cv[1][i];
		cv_old[2][i] = cv[2][i];
	}

	dt=Time_Step(ib, imax, cfl, x, cv, p);
	//std::cout<<"ib and imax are"<<ib<<""<<imax<<std::endl;
	//exit(0);
	//LR_State(ib, id, imax, x, cv, p, L_state, R_state);
	
	if(flux_type=="llf"){

		Flux_LLF(ib, imax, L_state, R_state, residue);
	}
	else if(flux_type=="movers"){

		Flux_Movers(ib, imax, L_state, R_state, residue);
	}
	else if(flux_type=="movers_h"){

		Flux_Movers_H(ib, imax, L_state, R_state, residue);
	}
	else if(flux_type=="movers_le"){

		Flux_Movers_LE(ib, imax, L_state, R_state, residue);
	}
	else if(flux_type=="movers_cp"){

		Flux_Movers_CP(ib, imax, L_state, R_state, residue);
	}
	else if(flux_type=="kfds"){

		Flux_KFDS(ib, id1, id2, L_state, R_state, residue);
	}
	else if(flux_type=="kfds_movers"){

		Flux_KFDS_Movers(ib, imax, L_state, R_state, residue);
	}
    else if(flux_type=="kfds2nd"){

        Flux_KFDS2nd_Order(ib, imax, cv, p, t, residue);
    }
	else if(flux_type=="kfds2nd_movers"){
		Flux_Con_Var_Defns(ib, imax, cv, p, Dcv, Acv, Dflux, Aflux);
        //std::cout<<"in solver explicit ib and imax are"<<ib<<" "<<id1<<" "<<imax<<std::endl;
        //exit(0);
        Flux_KFDS2nd_Movers(ib, id1, id2, cv, Dcv, Acv, p, t, flux, Dflux, Aflux, residue);
	}
	
	for(int i=2; i<=ib; i++){
		cv[0][i] = cv_old[0][i]-dt*residue[0][i];
		cv[1][i] = cv_old[1][i]-dt*residue[1][i];
		cv[2][i] = cv_old[2][i]-dt*residue[2][i];
						
		//update pressure
	
		p[i]=(Gamma-1.0)*(cv[2][i]-0.5*cv[1][i]*cv[1][i]/cv[0][i]);

		if (isnan(p[i])) {
			std::cout<<"\t"<<p[i]<<std::endl;
			exit(0);
		
	    }
	Boundary_Condition(2, ib, cv, p);
	}
}

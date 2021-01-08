#include "global_declarations.h"

void Initial_Conditions(double x0, int ib, int imax, double *&x, double *Prim_var_L, double *Prim_var_R, double **&cv, double *&p){

	double rho, u, p_i, E;
	
	for(int i=2;i<=ib;i++){
		if(x[i]<x0){
			rho=Prim_var_L[0];
			u=Prim_var_L[1];
			p_i=Prim_var_L[2];
			E=p_i/(rho*(Gamma-1.0)) + 0.5*u*u;

			cv[0][i] = rho;
			cv[1][i] = rho*u;
			cv[2][i] = rho*E;
			p[i]     = p_i;
		}
		else{
			rho=Prim_var_R[0];
			u=Prim_var_R[1];
			p_i=Prim_var_R[2];
			E=p_i/(rho*(Gamma-1.0)) + 0.5*u*u;

			cv[0][i] = rho;
			cv[1][i] = rho*u;
			cv[2][i] = rho*E;
			p[i]     = p_i;
		}
	}

	for (int k=0;k<3;k++){
	    cv[k][1] = cv[k][2];
	    cv[k][0] = cv[k][1];

        cv[k][id1] = cv[k][ib];
        cv[k][id2] = cv[k][id1];
	}
        p[1]   = p[2];
        p[0]   = p[1];
        p[id1] = p[ib];
	    p[id2] = p[id1];

	//debug
    /*for (int i=0;i<=imax;i++){
        std::cout<<"cv[0]["<<i<<"]="<<cv[0][i]<<"\t"<<"cv[1]["<<i<<"]="<<cv[1][i]<<"\t"<<"cv[2]["<<i<<"]="<<cv[2][i]
                <<"\t"<<"p["<<i<<"]="<<p[i]<<std::endl;
    }
    exit(0);*/
}





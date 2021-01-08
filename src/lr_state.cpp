#include "global_declarations.h"
#include "basic_functions.h"

/****to compute left and right state****/

template <typename T>
T Venki_Limiter(T del_plus, T del_minus, T dx);

template <typename T>
T VanAlbada_Limiter(T del_plus, T del_minus);

void LR_State(int ib, int id1, int id2, double *&x, double **&cv, double *&p, double **&L_state, double **&R_state){

	double **du, dx;

	du=Allocate_2D(du, ndvar, imax);

	if(space_accuracy=="1"){
		for(int i=2; i<=ib;i++){
	
			R_state[0][i] = cv[0][i+1]; //rho
			R_state[1][i] = cv[1][i+1]/cv[0][i+1]; //u
			R_state[2][i] = p[i+1];

			L_state[0][i] = cv[0][i]; //rho
			L_state[1][i] = cv[1][i]/cv[0][i]; //u
			L_state[2][i] = p[i];
		}

		for(int k=0;k<3;k++){
			R_state[k][1]    = R_state[k][2];
			R_state[k][0]    = R_state[k][1];
			R_state[k][id1]   = R_state[k][ib];
			R_state[k][id2] = R_state[k][id1];

			L_state[k][1]    = L_state[k][2];
			L_state[k][0]    = L_state[k][1];
			L_state[k][id1]   = L_state[k][ib];
			L_state[k][id2] = L_state[k][id1];
		}

		/*for (int i=0;i<=imax;i++){
			std::cout<<"L_state[0]["<<i<<"]="<<L_state[0][i]<<"\t"<<"L_state[1]["<<i<<"]="<<L_state[1][i]<<"\t"
				     <<"L_state[2]["<<i<<"]="<<L_state[2][i]<<std::endl;
		}
		exit(0);*/
	}

	else if	(space_accuracy=="2"){
		for(int i=2; i<=ib;i++){
			du[0][i] = cv[0][i+1]-cv[0][i]; //rho
			du[1][i] = cv[1][i+1]/cv[0][i]-cv[1][i]/cv[0][i]; //u
			du[2][i] = p[i+1]-p[i]; //p
		}

		for(int k=0;k<3;k++){
			du[k][1]   = du[k][2];
			du[k][0]   = du[k][1];
			du[k][id1] = du[k][ib];
			du[k][id2] = du[k][id1];
		}

		for (int i=0;i<=imax;i++){
			std::cout<<"du[0]["<<i<<"]="<<du[0][i]<<"\t"<<"du[1]["<<i<<"]="<<du[1][i]<<"\t"<<"du[2]["<<i<<"]="<<du[2][i]
				 <<p[i]<<std::endl;
		}

		for(int i=2; i<=ib;i++){
		
			dx = sqrt((x[i]-x[i-1])*(x[i]-x[i-1]));

			L_state[0][i] = cv[0][i]	      + 0.5*Venki_Limiter(du[0][i], du[0][i-1], dx); //rho
			L_state[1][i] = cv[1][i]/cv[0][i] + 0.5*Venki_Limiter(du[1][i], du[1][i-1], dx); //u
			L_state[2][i] = p[i]	          + 0.5*Venki_Limiter(du[2][i], du[2][i-1], dx); //p

			dx = sqrt((x[i+1]-x[i])*(x[i+1]-x[i]));

			R_state[0][i] = cv[0][i+1]	          - 0.5*Venki_Limiter(du[0][i+1], du[0][i], dx); //rho
			R_state[1][i] = cv[1][i+1]/cv[0][i+1] - 0.5*Venki_Limiter(du[1][i+1], du[1][i], dx); //u
			R_state[2][i] = p[i+1]		          - 0.5*Venki_Limiter(du[2][i+1], du[2][i], dx); //p

			/*L_state[0][i] = cv[0][i]	  + 0.5*VanAlbada_Limiter(du[0][i], du[0][i-1])*du[0][i-1]; //rho
			L_state[1][i] = cv[1][i]/cv[0][i] + 0.5*VanAlbada_Limiter(du[1][i], du[1][i-1])*du[1][i-1]; //u
			L_state[2][i] = p[i]	          + 0.5*VanAlbada_Limiter(du[2][i], du[2][i-1])*du[2][i-1]; //p

			
			R_state[0][i] = cv[0][i+1]	      - 0.5*VanAlbada_Limiter(du[0][i+1], du[0][i])*du[0][i+1]; //rho
			R_state[1][i] = cv[1][i+1]/cv[0][i+1] - 0.5*VanAlbada_Limiter(du[1][i+1], du[1][i])*du[1][i+1]; //u
			R_state[2][i] = p[i+1]		      - 0.5*VanAlbada_Limiter(du[2][i+1], du[2][i])*du[2][i+1]; //p*/
			
		}

		for(int k=0;k<3;k++){
			R_state[k][1]    = R_state[k][2];
			R_state[k][0]    = R_state[k][1];
			R_state[k][id1]   = R_state[k][ib];
			R_state[k][id2] = R_state[k][id1];

			L_state[k][1]   = L_state[k][2];
			L_state[k][0]   = L_state[k][1];
			L_state[k][id1] = L_state[k][ib];
			L_state[k][id2] = L_state[k][id1];
		}

		/*for (int i=0;i<=imax;i++){
    		std::cout<<"L_state[0]["<<i<<"]="<<L_state[0][i]<<"\t"<<"L_state[1]["<<i<<"]="<<L_state[1][i]<<"\t"
            		 <<"L_state[2]["<<i<<"]="<<L_state[2][i]<<std::endl;
		}
		exit(0);*/
	}
	Delete_2D(du, ndvar);
}


template <typename T>
T Venki_Limiter(T del_plus, T del_minus, T dx){

	double const k=0.5;
	double eps;
	
	eps=pow(k*dx, 3);

	return ((del_plus*del_plus+eps)*del_minus+(del_minus*del_minus+eps)*del_plus)/(del_plus*del_plus+del_minus*del_minus+2.0*eps+1e-3);
}

template <typename T>
T VanAlbada_Limiter(T del_plus, T del_minus){

	const double	omega=1e-12;
	double r;

	del_minus = Sign(del_minus)*(fabs(del_minus)+omega);

	r=del_plus/(del_minus+1e-5);

	return ((r*r+r)/(r*r+1.0));
}


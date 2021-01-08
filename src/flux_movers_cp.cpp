#include "global_declarations.h"
#include "basic_functions.h"

/****computes the interface flux using MOVERS dissipation****/


template <typename T>
T Movers_Con(T Fr, T Fl, T Ur,T Ul,T ur, T ul);

template <typename T>
T Movers_Pres(T Fr, T Fl, T Ur,T Ul,T ar, T al);


void Flux_Movers_CP(int ib, int imax, double **L_state, double **R_state, double **residue){

	double rhor, ur, pr, rhoul, El, hr, rhol, ul, pl, rhour, Er, hl;
	double fconv1, fconv2, fconv3, fdiss1, fdiss2, fdiss3, dx;
	double al, ar;
	double **iflux;
	double con_diss1, con_diss2, con_diss3;
	double pre_diss1, pre_diss2, pre_diss3;

	iflux = Allocate_2D(iflux, ndvar, imax);

	for (int i=1;i<=ib;i++){

		rhol = L_state[0][i];
		  ul = L_state[1][i];
		  pl = L_state[2][i];
		  El = (pl/(rhol*(Gamma-1.0))) + 0.5*ul*ul;
		  hl = pl*Gamma/(rhol*(Gamma-1.0)) + 0.5*ul*ul;
	       rhoul = rhol*ul;	
		  al = sqrt(Gamma*pl/rhol);

		rhor = R_state[0][i];
		  ur = R_state[1][i];
		  pr = R_state[2][i];
		  Er = (pr/(rhor*(Gamma-1.0))) + 0.5*ur*ur;
		  hr = pr*Gamma/(rhor*(Gamma-1.0)) + 0.5*ur*ur;
	       rhour = rhor*ur;	
		  ar = sqrt(Gamma*pr/rhor);
		
/****convective flux****/

		fconv1 = rhoul + rhour;
		fconv2 = pl + rhoul*ul + pr + rhour*ur;
		fconv3 = rhoul*hl + rhour*hr;

/****dissipative flux****/

		con_diss1 = Movers_Con(rhour, rhoul, rhor, rhol, ur, ul);
		con_diss2 = Movers_Con(rhour*ur, rhoul*ul, rhour, rhoul, ur, ul);
		con_diss3 = Movers_Con(rhour*Er, rhoul*El, rhor*Er, rhol*El, ur, ul);

		pre_diss1 = 0.0;
		pre_diss2 = Movers_Pres(pr, pl, rhour, rhoul, ar, al);
		pre_diss3 = Movers_Pres(pr*ur, pl*ul, rhor*Er, rhol*El, ar, al);
 
		fdiss1 = (con_diss1+pre_diss1)*(rhor-rhol);		
		fdiss2 = (con_diss2+pre_diss2)*(rhour-rhoul);
		fdiss3 = (con_diss3+pre_diss3)*(rhor*Er-rhol*El);

		iflux[0][i] = 0.5*(fconv1-fdiss1);
		iflux[1][i] = 0.5*(fconv2-fdiss2);
		iflux[2][i] = 0.5*(fconv3-fdiss3);

	}

/**** residue or RHS of the update formula****/

	for(int i=2;i<=ib;i++){

		dx=0.5*fabs(x[i+1]-x[i-1]);						

		residue[0][i] = (1.0/dx)*(iflux[0][i]-iflux[0][i-1]);
		residue[1][i] = (1.0/dx)*(iflux[1][i]-iflux[1][i-1]);
		residue[2][i] = (1.0/dx)*(iflux[2][i]-iflux[2][i-1]);
	}

	Delete_2D(iflux, ndvar);
}



template <typename T>
T Movers_Con(T Fr, T Fl, T Ur,T Ul,T ur, T ul){

	const double epsilon=1e-6;
	double S;

	if (fabs(Fr-Fl)<epsilon) return 0.0;
	else if (fabs(Ur-Ul)<epsilon) return std::max(ur,ul);
	else if(fabs(Ur-Ul)>epsilon and fabs(Fr-Fl)>epsilon) S= fabs(((Fr-Fl)/(Ur-Ul)));

	return S;
}

template <typename T>
T Movers_Pres(T Fr, T Fl, T Ur,T Ul,T ar, T al){

	double max_eig, S;
	const double epsilon=1e-6;

	max_eig = std::max(fabs(((Gamma-1.0)/Gamma))*ar, fabs(((Gamma-1.0)/Gamma))*al);

	if (fabs(Fr-Fl)<epsilon) return 0.0;
	else if (fabs(Ur-Ul)<epsilon) return max_eig;
	else if(fabs(Ur-Ul)>epsilon and fabs(Fr-Fl)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));	

	if (S<0.0) return 0.0;
	else if(S>max_eig) return max_eig;	
	
	return S;
}

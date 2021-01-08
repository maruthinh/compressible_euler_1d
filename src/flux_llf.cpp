#include "global_declarations.h"
#include "basic_functions.h"

/****computes the interface flux using MOVERS dissipation****/

template <typename T>
T MaxEigVal(T ur, T ul, T ar, T al);

void Flux_LLF(int ib, int imax, double **L_state, double **R_state, double **residue){

	double rhor, ur, pr, rhoul, El, hr, rhol, ul, pl, rhour, Er, hl;
	double fconv1, fconv2, fconv3, fdiss1, fdiss2, fdiss3, dx;
	double al, ar, L_max;
	double **iflux;

	iflux = Allocate_2D(iflux, ndvar, imax);

	for (int i=1;i<=ib;i++){

		rhol = L_state[0][i];
		  ul = L_state[1][i];
		  pl = L_state[2][i];
		  El = (pl/(rhol*(Gamma-1.0))) + 0.5*ul*ul;
		  hl = (pl*Gamma)/(rhol*(Gamma-1.0)) + 0.5*ul*ul;
	       rhoul = rhol*ul;	
		  al = sqrt(Gamma*pl/rhol);

		rhor = R_state[0][i];
		  ur = R_state[1][i];
		  pr = R_state[2][i];
		  Er = (pr/(rhor*(Gamma-1.0))) + 0.5*ur*ur;
		  hr = (pr*Gamma)/(rhor*(Gamma-1.0)) + 0.5*ur*ur;
	       rhour = rhor*ur;	
		  ar = sqrt(Gamma*pr/rhor);
		
/****convective flux****/

		fconv1 = rhoul + rhour;
		fconv2 = pl + rhoul*ul + pr + rhour*ur;
		fconv3 = rhoul*hl + rhour*hr;

/****dissipative flux****/

		L_max = MaxEigVal(ur, ul, ar, al);

 		fdiss1 = L_max*(rhor-rhol);		
		fdiss2 = L_max*(rhour-rhoul);
		fdiss3 = L_max*(rhor*Er-rhol*El);

		iflux[0][i] = 0.5*(fconv1-fdiss1);
		iflux[1][i] = 0.5*(fconv2-fdiss2);
		iflux[2][i] = 0.5*(fconv3-fdiss3);

	}

/**** residue or RHS of the update formula****/

	for(int i=2;i<=ib;i++){

		dx=0.5*(x[i+1]-x[i-1]);						

		residue[0][i] = (1.0/dx)*(iflux[0][i]-iflux[0][i-1]);
		residue[1][i] = (1.0/dx)*(iflux[1][i]-iflux[1][i-1]);
		residue[2][i] = (1.0/dx)*(iflux[2][i]-iflux[2][i-1]);
	}

	Delete_2D(iflux, ndvar);
}

template <typename T>
T MaxEigVal(T ur, T ul, T ar, T al){

	double L1r, L2r, L3r, L1l, L2l, L3l;

		L1r = fabs(ur+ar);  L1l = fabs(ul+al);
		L2r = fabs(ur);     L2l = fabs(ul);
		L3r = fabs(ur-ar);  L3l = fabs(ul-al);

	//return Max3(max(L1l,L1r), max(L2l,L2r), max(L3l,L3r));

	return Max2(Max3(L1r, L2r, L3r), Max3(L1l, L2l, L3l));
}



#include "global_declarations.h"
#include "basic_functions.h"

/****computes the interface flux using MOVERS dissipation****/

template <typename T>
T MaxEigVal(T ur, T ul, T ar, T al);

template <typename T>
T MinEigVal(T ur, T ul, T ar, T al);

template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min);


void Flux_Movers(int ib, int imax, double **L_state, double **R_state, double **residue){

	double rhor, ur, pr, rhoul, El, hr, rhol, ul, pl, rhour, Er, hl;
	double fconv1, fconv2, fconv3, fdiss1, fdiss2, fdiss3, dx;
	double al, ar, L_max, L_min;
	double **iflux;

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

		L_max = MaxEigVal(ur, ul, ar, al);
		L_min = MinEigVal(ur, ul, ar, al);
 
		fdiss1 = Movers(rhour, rhoul, rhor, rhol, L_max, L_min)*(rhor-rhol);		
		fdiss2 = Movers(pr + rhour*ur, pl + rhoul*ul, rhour, rhoul, L_max, L_min)*(rhour-rhoul);
		fdiss3 = Movers(rhour*hr, rhoul*hl, rhor*Er, rhol*El, L_max, L_min)*(rhor*Er-rhol*El);

/*		fdiss1 = Movers(rhour*hr, rhoul*hl, rhor*Er, rhol*El, L_max, L_min)*(rhor-rhol);		
		fdiss2 = Movers(rhour*hr, rhoul*hl, rhor*Er, rhol*El, L_max, L_min)*(rhour-rhoul);
		fdiss3 = Movers(rhour*hr, rhoul*hl, rhor*Er, rhol*El, L_max, L_min)*(rhor*Er-rhol*El);*/

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
T MaxEigVal(T ur, T ul, T ar, T al){

	double L1r, L2r, L3r, L1l, L2l, L3l;

		L1r = fabs(ur+ar);  L1l = fabs(ul+al);
		L2r = fabs(ur);     L2l = fabs(ul);
		L3r = fabs(ur-ar);  L3l = fabs(ul-al);

	//return Max3(max(L1l,L1r), max(L2l,L2r), max(L3l,L3r));

	return Max2(Max3(L1r, L2r, L3r), Max3(L1l, L2l, L3l));
}

template <typename T>
T MinEigVal(T ur, T ul, T ar, T al){

	double L1r, L2r, L3r, L1l, L2l, L3l;

		L1r = fabs(ur+ar);  L1l = fabs(ul+al);
		L2r = fabs(ur);     L2l = fabs(ul);
		L3r = fabs(ur-ar);  L3l = fabs(ul-al);

	//return Max3(min(L1l,L1r), min(L2l,L2r), min(L3l,L3r));

	return Max2(Min3(L1r, L2r, L3r), Min3(L1l, L2l, L3l));
}

template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S, alpha;
	const double k=0.5, epsilon=1e-6;

	if (fabs(Fr-Fl)<epsilon) return 0.0;
	else if (fabs(Ur-Ul)<epsilon) return L_min;
	else if (fabs(Ur-Ul)>epsilon and fabs(Fr-Fl)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));
	else S=L_min;

	if (S<epsilon)	return 0;
	else if ((S)>=L_max) return (L_max);
	else if((S)<=L_min) return (L_min);

	/*alpha=k*L_max;

	if(S<alpha){
		return (S*S+alpha*alpha)/(2.0*alpha);		
	}
	else{
		return S;
	}*/
}

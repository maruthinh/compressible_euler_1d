#include "global_declarations.h"
#include "basic_functions.h"

/****computes the interface flux using MOVERS dissipation****/

template <typename T>
T MaxEigVal(T ur, T ul, T ar, T al);

template <typename T>
T MinEigVal(T ur, T ul, T ar, T al);

template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min);


void Flux_Movers_LE(int ib, int imax, double **L_state, double **R_state, double **residue){

	double rhor, ur, pr, rhoul, El, hr, rhol, ul, pl, rhour, Er, hl;
	double fconv1, fconv2, fconv3, fdiss1, fdiss2, fdiss3, dx;
	double al, ar, L_max, L_min;
	double **iflux, **diss, **du;
	double mdiss1, mdiss2, mdiss3;
	double diss1, diss2, diss3;

	
	du    = Allocate_2D(du, ndvar, imax+1);
	diss  = Allocate_2D(diss, ndvar, imax);
	iflux = Allocate_2D(iflux, ndvar, imax);

	for(int i=1; i<=ib;i++){
		du[0][i] = cv[0][i+1]-cv[0][i]; //rho
		du[1][i] = cv[1][i+1]-cv[1][i]; //rho*u
		du[2][i] = cv[2][i+1]-cv[2][i]; //rho*E
	}

	du[0][0] = du[0][1];
	du[1][0] = du[1][1];
	du[2][0] = du[2][1];

	du[0][imax] = du[0][ib];
	du[1][imax] = du[1][ib];
	du[2][imax] = du[2][ib];

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

		mdiss1 = Movers(rhour, rhoul, rhor, rhol, L_max, L_min);
		mdiss2 = Movers(pr + rhour*ur, pl + rhoul*ul, rhour, rhoul, L_max, L_min);
		mdiss3 = Movers(rhour*hr, rhoul*hl, rhor*Er, rhol*El, L_max, L_min);

		diss1 = mdiss1+MinMod(du[0][i+1], du[0][i], du[0][i-1])*(fabs(std::max(ur, ul))-mdiss1);
		diss2 = mdiss2+MinMod(du[1][i+1], du[1][i], du[1][i-1])*(fabs(std::max(ur, ul))-mdiss2);
		diss3 = mdiss3+MinMod(du[2][i+1], du[2][i], du[2][i-1])*(fabs(std::max(ur, ul))-mdiss3);
 
		fdiss1 = diss1*(rhor-rhol);		
		fdiss2 = diss2*(rhour-rhoul);
		fdiss3 = diss3*(rhor*Er-rhol*El);

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
	Delete_2D(diss, ndvar);
	Delete_2D(du, ndvar);
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

	double S;
	const double epsilon=1e-6;

	if (fabs(Fr-Fl)<epsilon) return 0.0;
	else if (fabs(Ur-Ul)<epsilon) return L_min;
	else if (fabs(Ur-Ul)>epsilon and fabs(Fr-Fl)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));
	else S=L_min;
	
	if (S<epsilon)	return 0;
	else if ((S)>=L_max) return (L_max);
	else if((S)<=L_min) return (L_min);
	else return (S);
	return S;
}

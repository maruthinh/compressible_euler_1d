//
// Created by Maruthi on 08-05-2018.
//

#include "global_declarations.h"
#include "basic_functions.h"

/****computes the interface flux using MOVERS dissipation****/

void Flux_KFDS(int ib, int id1, int id2, double **L_state, double **R_state, double **residue){

    double rhor, ur, pr, rhoul, El, hr, rhol, ul, pl, rhour, Er, hl;
    double fconv1, fconv2, fconv3, fdiss1, fdiss2, fdiss3, dx;
    double al, ar, L_max;
    double **iflux;

    iflux = Allocate_2D(iflux, ndvar, imax);
    L_max = GlobMaxEigVal(ib, imax, cv, p);

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

double GlobMaxEigVal(int ib, int imax, double **&cv, double *&p){

    double rho, u, a, *EigVal;

    EigVal = new double [imax+1];

    for(int i=2; i<=ib;i++){
        rho = cv[0][i];
        u = cv[1][i]/cv[0][i];
        a = sqrt(Gamma*p[i]/rho);
        EigVal[i] = fabs(u)+a;
    }

    EigVal[1] = EigVal[2];
    EigVal[0] = EigVal[1];
    EigVal[id1] = EigVal[ib];
    EigVal[id2] = EigVal[id1];

    GlobMaxEigenVal=10e-32;

    for(int i=1; i<=imax; i++){
        if(EigVal[i]>GlobMaxEigenVal) GlobMaxEigenVal=EigVal[i];
    }
    delete [] EigVal;
    return GlobMaxEigenVal;
}

double GlobMinEigVal(int ib, int imax, double **&cv, double *&p){

    double rho, u, a, *EigVal;

    EigVal = new double [imax+1];

    for(int i=2; i<=ib;i++){
        rho = cv[0][i];
        u = cv[1][i]/cv[0][i];
        a = sqrt(Gamma*p[i]/rho);
        EigVal[i] = fabs(u)+a;
    }

    EigVal[1] = EigVal[2];
    EigVal[0] = EigVal[1];
    EigVal[id1] = EigVal[ib];
    EigVal[id2] = EigVal[id1];

    GlobMinEigenVal=10e32;

    for(int i=1; i<=imax; i++){
        if(EigVal[i]<GlobMinEigenVal) GlobMinEigenVal=EigVal[i];
    }
    delete [] EigVal;
    return GlobMinEigenVal;
}


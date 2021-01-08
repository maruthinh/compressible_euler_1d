//
// Created by Maruthi on 22-05-2018.
//
#include "global_declarations.h"
#include "basic_functions.h"

/****computes the interface flux using MOVERS dissipation****/

void Flux_Con_Var_Defns(int ib, int imax, double **&cv, double *&p, double **&Dcv, double **&Acv, double **&Dflux,
                        double **&Aflux){

    double rho, u, P, rhou, E, h;

    for (int i=2;i<=ib;i++){
        rho = cv[0][i];
        u  = cv[1][i]/rho;
        P  = p[i];
        E  = (P/(rho*(Gamma-1.0))) + 0.5*u*u;
        h  = E + P/rho;
        rhou = rho*u;

/****flux****/
        flux[0][i] = rhou;
        flux[1][i] = P + rhou*u;
        flux[2][i] = rhou*h;

        //std::cout<<"flux="<<flux[0][i]<<"\t"<<flux[1][i]<<"\t"<<flux[2][i]<<std::endl;
    }

    for (int k=0;k<3;k++){
        flux[k][1] = flux[k][2];
        flux[k][0] = flux[k][1];

        flux[k][id1] = flux[k][ib];
        flux[k][id2] = flux[k][id1];
    }

    //debug
    /*for (int i=0;i<=id2;i++){
        std::cout<<"flux[0]["<<i<<"]="<<flux[0][i]<<"\t"<<"flux[1]["<<i<<"]="<<flux[1][i]<<"\t"<<"flux[2]["<<i<<"]="
                 <<flux[2][i]<<std::endl;
    }
    exit(0);*/

/**flux differences**/
    for(int i=2; i<=ib;i++){
        Dflux[0][i] = flux[0][i+1]-flux[0][i]; //rho*u
        Dflux[1][i] = flux[1][i+1]-flux[1][i]; //p+rho*u^2
        Dflux[2][i] = flux[2][i+1]-flux[2][i]; //p*u + rho*u*E
    }

    for (int k=0;k<3;k++){
        Dflux[k][1] = Dflux[k][2];
        Dflux[k][0] = Dflux[k][1];

        Dflux[k][id1]   = Dflux[k][ib];
        Dflux[k][id2] = Dflux[k][id1];
    }

/**flux 2*averages**/
    for(int i=2; i<=ib;i++){
        Aflux[0][i] = flux[0][i+1]+flux[0][i]; //rho*u
        Aflux[1][i] = flux[1][i+1]+flux[1][i]; //p+rho*u^2
        Aflux[2][i] = flux[2][i+1]+flux[2][i]; //p*u + rho*u*E
    }

    for (int k=0;k<3;k++){
        Aflux[k][1] = Aflux[k][2];
        Aflux[k][0] = Aflux[k][1];

        Aflux[k][id1] = Aflux[k][ib];
        Aflux[k][id2] = Aflux[k][id1];
    }
/**Con var differences**/
    for(int i=2; i<=ib;i++){
        Dcv[0][i] = cv[0][i+1]-cv[0][i]; //rho
        Dcv[1][i] = cv[1][i+1]-cv[1][i]; //rho*u
        Dcv[2][i] = cv[2][i+1]-cv[2][i]; //rho*u*E
    }

    for (int k=0;k<3;k++){
        Dcv[k][1] = Dcv[k][2];
        Dcv[k][0] = Dcv[k][1];

        Dcv[k][id1]   =  Dcv[k][ib];
        Dcv[k][id2] = Dcv[k][id1];
    }

/**Con var 2*averages**/
    for(int i=2; i<=ib;i++){
        Acv[0][i] = cv[0][i+1]+cv[0][i]; //rho
        Acv[1][i] = cv[1][i+1]+cv[1][i]; //rho*u
        Acv[2][i] = cv[2][i+1]+cv[2][i]; //rho*u*E
    }

    for (int k=0;k<3;k++){
        Acv[k][1] = Acv[k][2];
        Acv[k][0] = Acv[k][1];

        Acv[k][id1] = Acv[k][ib];
        Acv[k][id2] = Acv[k][id1];
    }

    //debug
    /*for (int i=0;i<=id2;i++){
        std::cout<<"Acv[0]["<<i<<"]="<<Acv[0][i]<<"\t"<<"Acv[1]["<<i<<"]="<<Acv[1][i]<<"\t"<<"Acv[2]["<<i<<"]="
                 <<Acv[2][i]<<std::endl;
    }
    exit(0);*/
}
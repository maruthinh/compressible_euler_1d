//
// Created by Maruthi on 20-05-2018.
//

#include "global_declarations.h"
#include "basic_functions.h"

template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min);

template <typename T>
T MoversWEF(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min);

/****computes the interface flux using MOVERS dissipation****/

void Flux_KFDS2nd_Movers(int ib, int imax, double **&cv, double *&p, double *&t, double **&residue){

    double rho, u, P, rhou, E, h, a;
    double fconv1, fconv2, fconv3, fdiss1, fdiss2, fdiss3, dx;
    double fterm1, fterm2, fterm3,fterm4, fterm5;
    double al, ar, L_max, L_min, dt_dx;
    double **iflux;

    iflux = Allocate_2D(iflux, ndvar, imax);
    L_max = GlobMaxEigVal(ib, imax, cv, p);
    L_min = GlobMinEigVal(ib, imax, cv, p);

    for (int i=1;i<=imax;i++){

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
    }

/**** residue or RHS of the update formula****/
    for(int i=2;i<ib;i++){

        dx=0.5*(x[i+1]-x[i-1]);
        dt_dx = t[i]/dx;

        fterm1 = 0.5*(flux[0][i]+flux[0][i+1])-0.5*L_max*(cv[0][i+1]-cv[0][i]);
        fterm2 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[0][i+1]-cv[0][i])-0.5*(flux[0][i+1]-flux[0][i]);
        fterm3 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[0][i]-cv[0][i-1])-0.5*(flux[0][i]-flux[0][i-1]);

        fterm4 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[0][i+1]-cv[0][i])-0.5*(flux[0][i+1]-flux[0][i]);
        fterm5 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[0][i+2]-cv[0][i+1])-0.5*(flux[0][i+2]-flux[0][i+1]);

        iflux[0][i] = fterm1;//-0.5*MinModLim2(fterm2,fterm3)-0.5*MinModLim2(fterm4,fterm5);

        fterm1 = 0.5*(flux[1][i]+flux[1][i+1])-0.5*L_max*(cv[1][i+1]-cv[1][i]);
        fterm2 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[1][i+1]-cv[1][i])-0.5*(flux[1][i+1]-flux[1][i]);
        fterm3 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[1][i]-cv[1][i-1])-0.5*(flux[1][i]-flux[1][i-1]);

        fterm4 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[1][i+1]-cv[1][i])-0.5*(flux[1][i+1]-flux[1][i]);
        fterm5 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[1][i+2]-cv[1][i+1])-0.5*(flux[1][i+2]-flux[1][i+1]);

        iflux[1][i] = fterm1;//-0.5*MinModLim2(fterm2,fterm3)-0.5*MinModLim2(fterm4,fterm5);

        fterm1 = 0.5*(flux[2][i]+flux[2][i+1])-0.5*L_max*(cv[2][i+1]-cv[2][i]);
        fterm2 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[2][i+1]-cv[2][i])-0.5*(flux[2][i+1]-flux[2][i]);
        fterm3 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[2][i]-cv[2][i-1])-0.5*(flux[2][i]-flux[2][i-1]);

        fterm4 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[2][i+1]-cv[2][i])-0.5*(flux[2][i+1]-flux[2][i]);
        fterm5 = (0.5*L_max*L_max*dt_dx-0.5*L_max)*(cv[2][i+2]-cv[2][i+1])-0.5*(flux[2][i+2]-flux[2][i+1]);

        iflux[2][i] = fterm1;//-0.5*MinModLim2(fterm2,fterm3)-0.5*MinModLim2(fterm4,fterm5);
    }

    for(int i=2;i<=ib;i++){

        dx=0.5*(x[i+1]-x[i-1]);

        residue[0][i] = (1.0/dx)*(iflux[0][i]-iflux[0][i-1]);
        residue[1][i] = (1.0/dx)*(iflux[1][i]-iflux[1][i-1]);
        residue[2][i] = (1.0/dx)*(iflux[2][i]-iflux[2][i-1]);
    }

    Delete_2D(iflux, ndvar);
}

template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){
    double S;
    const double k=0.9, epsilon=1e-16;

    if (fabs(Fr-Fl)<epsilon) return 0.0;
    else if (fabs(Ur-Ul)<epsilon) return L_min;
    else if (fabs(Ur-Ul)>epsilon and fabs(Fr-Fl)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));
    else S=L_min;

    if (S<epsilon)	return 0;
    else if (S>=L_max) return L_max;
    else if(S<=L_min) return L_min;
    else return S;
}

template <typename T>
T MoversWEF(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

    double S, alpha;
    const double k=0.5, epsilon=1e-16;

    if (fabs(Fr-Fl)<epsilon) return 0.0;
    else if (fabs(Ur-Ul)<epsilon) return L_min;
    else if (fabs(Ur-Ul)>epsilon and fabs(Fr-Fl)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));
    else S=L_min;

    if (S<epsilon)	S = 0;
    else if (S>=L_max) S = L_max;
    else if (S<=L_min) S = L_min;
    alpha=k*L_max;

    if(S<alpha){
        return (S*S+alpha*alpha)/(2.0*alpha);
    }
    else{
        return S;
    }
}

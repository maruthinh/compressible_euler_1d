//
// Created by Maruthi on 22-05-2018.
//

#include "global_declarations.h"
#include "basic_functions.h"

template <typename T>
T MaxEigVal(T ur, T ul, T ar, T al);

template <typename T>
T MinEigVal(T ur, T ul, T ar, T al);

template <typename T>
T Movers(T deltaF, T deltaU, T L_max, T L_min);

/****computes the interface flux using MOVERS dissipation****/

void Flux_KFDS2nd_Movers(int ib, int id1, int id2, double **&cv, double **&Dcv, double **&Acv, double *&p, double *&t,
                        double **&flux, double **&Dflux, double **&Aflux, double **&residue){

    double dx;
    double fterm1, fterm2, fterm3,fterm4, fterm5;
    double al, ar, dt_dx;
    double **iflux, *max_eig_val, *min_eig_val, *lambda0, *lambda1, *lambda2, global_max_eig;
    double rhol, ul , pl, rhor, ur, pr;

    iflux = Allocate_2D(iflux, ndvar, imax);
    lambda1 = new double[3]; lambda2 = new double[3]; lambda0 = new double[3]; max_eig_val = new double[3];
    min_eig_val = new double[3];

    global_max_eig = GlobMaxEigVal(ib, imax, cv, p);
    //std::cout<<"global max eigen value="<<global_max_eig<<std::endl;
    for (int i=2;i<=ib;i++){
        rhol = cv[0][i-1];
        ul   = cv[1][i-1]/rhol;
        pl   = p[i-1];
        al   = sqrt(Gamma*pl/rhol);

        rhor = cv[0][i];
        ur   = cv[1][i]/rhor;
        pr   = p[i];
        ar   = sqrt(Gamma*pr/rhor);

/****dissipative flux****/
        max_eig_val[0] = MaxEigVal(ur, ul, ar, al);
        min_eig_val[0] = MinEigVal(ur, ul, ar, al);

        for(int k=0;k<3;k++){
            lambda0[k] = Movers(Dflux[k][i-1], Dcv[k][i-1], max_eig_val[0], min_eig_val[0]);
        }

        /*********************/
        rhol = cv[0][i];
        ul   = cv[1][i]/rhol;
        pl   = p[i];
        al   = sqrt(Gamma*pl/rhol);

        rhor = cv[0][i+1];
        ur   = cv[1][i+1]/rhor;
        pr   = p[i+1];
        ar   = sqrt(Gamma*pr/rhor);

/****dissipative flux****/
        max_eig_val[1] = MaxEigVal(ur, ul, ar, al);
        min_eig_val[1] = MinEigVal(ur, ul, ar, al);

        for(int k=0;k<3;k++){
            lambda1[k] = Movers(Dflux[k][i], Dcv[k][i], max_eig_val[1], min_eig_val[1]);
        }

        //if(i<ib){
            rhol = cv[0][i+1];
            ul = cv[1][i+1]/rhol;
            pl = p[i+1];
            al = sqrt(Gamma*pl/rhol);

            rhor = cv[0][i+2];
            ur = cv[1][i+2]/rhor;
            pr = p[i+2];
            ar = sqrt(Gamma*pr/rhor);

            max_eig_val[2] = MaxEigVal(ur, ul, ar, al);
            min_eig_val[2] = MinEigVal(ur, ul, ar, al);

            for(int k=0;k<3;k++){
                lambda2[k] = Movers(Dflux[k][i+1], Dcv[k][i+1], max_eig_val[2], min_eig_val[2]);
            }

            dx=0.5*fabs(x[i+1]-x[i-1]);
            //dt_dx = t[i]/dx;
            dt_dx = dt/dx;

        //llf for both first and second order
        /*for(int k=0;k<3;k++){
            fterm1 = 0.5*(Aflux[k][i])-0.5*lambda1[k]*Dcv[k][i];
            fterm2 = (0.5*lambda1[k]*lambda1[k]*dt_dx-0.5*lambda1[k])*(Dcv[k][i])-0.5*Dflux[k][i];
            fterm3 = (0.5*lambda0[k]*lambda0[k]*dt_dx-0.5*lambda0[k])*(Dcv[k][i-1])-0.5*Dflux[k][i-1];

            fterm4 = (0.5*lambda1[k]*lambda1[k]*dt_dx-0.5*lambda1[k])*(Dcv[k][i])+0.5*Dflux[k][i];
            fterm5 = (0.5*lambda2[k]*lambda2[k]*dt_dx-0.5*lambda2[k])*(Dcv[k][i+1])+0.5*(Dflux[k][i+1]);

            iflux[k][i] = fterm1-0.5*MinModLim2(fterm2,fterm3)-0.5*MinModLim2(fterm4,fterm5);
        }*/

        //movers for first order and llf for second order
        for(int k=0;k<3;k++){
            fterm1 = 0.5*Aflux[k][i]-0.5*lambda1[k]*Dcv[k][i];
            fterm2 = (0.5*max_eig_val[1]*max_eig_val[1]*dt_dx-0.5*lambda1[k])*Dcv[k][i]-0.5*Dflux[k][i];
            fterm3 = (0.5*max_eig_val[0]*max_eig_val[0]*dt_dx-0.5*lambda0[k])*Dcv[k][i-1]-0.5*Dflux[k][i-1];

            fterm4 = (0.5*max_eig_val[1]*max_eig_val[1]*dt_dx-0.5*lambda1[k])*Dcv[k][i]+0.5*Dflux[k][i];
            fterm5 = (0.5*max_eig_val[2]*max_eig_val[2]*dt_dx-0.5*lambda2[k])*Dcv[k][i+1]+0.5*Dflux[k][i+1];

            iflux[k][i] = fterm1-0.5*MinModLim2(fterm2, fterm3)-0.5*MinModLim2(fterm4, fterm5);
        }

        //movers for both first and second order
        /*for(int k=0;k<3;k++){
            fterm1 = 0.5*Aflux[k][i]-0.5*lambda1[k]*Dcv[k][i];
            fterm2 = (0.5*lambda1[k]*lambda1[k]*dt_dx-0.5*lambda1[k])*Dcv[k][i]-0.5*Dflux[k][i];
            fterm3 = (0.5*lambda0[k]*lambda0[k]*dt_dx-0.5*lambda0[k])*Dcv[k][i-1]-0.5*Dflux[k][i-1];

            fterm4 = (0.5*lambda1[k]*lambda1[k]*dt_dx-0.5*lambda1[k])*Dcv[k][i]+0.5*Dflux[k][i];
            fterm5 = (0.5*lambda2[k]*lambda2[k]*dt_dx-0.5*lambda2[k])*Dcv[k][i+1]+0.5*Dflux[k][i+1];

            iflux[k][i] = fterm1-0.5*MinModLim2(fterm2, fterm3)-0.5*MinModLim2(fterm4, fterm5);
        }*/

        //Movers-1
        /*for(int k=0;k<3;k++){
            fterm1 = 0.5*Aflux[k][i]-0.5*lambda1[2]*Dcv[k][i];
            fterm2 = (0.5*lambda1[2]*lambda1[2]*dt_dx-0.5*lambda1[2])*Dcv[k][i]-0.5*Dflux[k][i];
            fterm3 = (0.5*lambda0[2]*lambda0[2]*dt_dx-0.5*lambda0[2])*Dcv[k][i-1]-0.5*Dflux[k][i-1];

            fterm4 = (0.5*lambda1[2]*lambda1[2]*dt_dx-0.5*lambda1[2])*Dcv[k][i]+0.5*Dflux[k][i];
            fterm5 = (0.5*lambda2[2]*lambda2[2]*dt_dx-0.5*lambda2[2])*Dcv[k][i+1]+0.5*Dflux[k][i+1];

            iflux[k][i] = fterm1-0.5*MinModLim2(fterm2, fterm3)-0.5*MinModLim2(fterm4, fterm5);
        }*/

        //Movers and global max eigen value
        /*for(int k=0;k<3;k++){
            fterm1 = 0.5*Aflux[k][i]-0.5*lambda1[k]*Dcv[k][i];
            fterm2 = (0.5*global_max_eig*global_max_eig*dt_dx-0.5*lambda1[k])*Dcv[k][i]-0.5*Dflux[k][i];
            fterm3 = (0.5*global_max_eig*global_max_eig*dt_dx-0.5*lambda0[k])*Dcv[k][i-1]-0.5*Dflux[k][i-1];

            fterm4 = (0.5*global_max_eig*global_max_eig*dt_dx-0.5*lambda1[k])*Dcv[k][i]+0.5*Dflux[k][i];
            fterm5 = (0.5*global_max_eig*global_max_eig*dt_dx-0.5*lambda2[k])*Dcv[k][i+1]+0.5*Dflux[k][i+1];

            iflux[k][i] = fterm1-0.5*MinModLim2(fterm2, fterm3)-0.5*MinModLim2(fterm4, fterm5);
        }*/

        //std::cout<<iflux[0][i]<<std::endl;
    }

    for(int i=3;i<ib;i++){

        dx=0.5*(x[i+1]-x[i-1]);

        residue[0][i] = (1.0/dx)*(iflux[0][i]-iflux[0][i-1]);
        residue[1][i] = (1.0/dx)*(iflux[1][i]-iflux[1][i-1]);
        residue[2][i] = (1.0/dx)*(iflux[2][i]-iflux[2][i-1]);
    }

    Delete_2D(iflux, ndvar);
    delete[] lambda0;
    delete[] lambda1;
    delete[] lambda2;
    delete[] max_eig_val;
    delete[] min_eig_val;
    //iflux = NULL; lambda0 = NULL; lambda1 = NULL; lambda2 = NULL;
    //max_eig_val = NULL; min_eig_val = NULL;
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

/*
template <typename T>
T Movers(T deltaF, T deltaU, T L_max, T L_min){

    double S, alpha;
    const double k=0.5, epsilon=1e-6;

    if (fabs(deltaF)<epsilon) return 0.0;
    else if (fabs(deltaU)<epsilon) return L_min;
    else if (fabs(deltaU)>epsilon and fabs(deltaF)>epsilon) S=fabs(((deltaF)/(deltaU)));
    else S=L_min;

    if (S<epsilon)	return 0;
    else if ((S)>=L_max) return (L_max);
    else if((S)<=L_min) return (L_min);

    alpha=k*L_max;

    if(S<alpha){
        return (S*S+alpha*alpha)/(2.0*alpha);
    }
    else{
        return S;
    }
}*/

template <typename T>
T Movers(T deltaF, T deltaU, T L_max, T L_min){

    double S;
    const double epsilon=1e-6;

    if (fabs(deltaF)<epsilon) return 0.0;
    else if (fabs(deltaU)<epsilon) return L_min;
    else if (fabs(deltaU)>epsilon and fabs(deltaF)>epsilon) S=fabs(deltaF/deltaU);
    else S=L_min;

    if (S<epsilon)	return 0;
    else if (S>=L_max) return L_max;
    else if (S<=L_min) return L_min;
    else return S;
}
#include "global_declarations.h"
#include "basic_functions.h"

/****to compute global time step****/

double Time_Step(int ib, int imax, double cfl, double *&x, double **&cv, double *&p){

	double spectrad, dx, rho, u, a, *dtime;
	
	for(int i=2; i<=ib;i++){
        rho = cv[0][i];
        u = cv[1][i]/cv[0][i];
        a = sqrt(Gamma*p[i]/rho);
        dx = 0.5*fabs(x[i+1]-x[i-1]);
        spectrad = fabs(u)+a;
        t[i] = cfl*dx/spectrad;
	}

    t[1]   = t[2];
	t[0]   = t[1];
    t[id1] = t[ib];
	t[id2] = t[id1];

	dt=10e32;
	for(int i=2; i<=ib; i++){
		if(t[i]<dt) dt=t[i];
	}
    /*for (int i=0;i<=imax;i++){
        std::cout<<"t["<<i<<"]="<<t[i]<<"\t"<<"dt="<<dt<<std::endl;
    }
    exit(0);*/
	return dt;
}

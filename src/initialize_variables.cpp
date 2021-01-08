#include "global_declarations.h"
#include "basic_functions.h"


/****declaration of the global variables used in the code****/
int N, ib, id1, id2, imax;
std::string space_accuracy, flux_type, time_accuracy, test_case;
double cfl;
double xll, xul, x0, x01, x02, t_fin;
double *Prim_var_L, *Prim_var_M, *Prim_var_R;
double *x, *t, *xc, **cv, **cv_old, *p, **residue, **flux, **Dflux, **Aflux, **Dcv, **Acv;
double dx, dt, GlobMaxEigenVal, GlobMinEigenVal;
double **L_state, **R_state;



/**************  to allocate 1d array***********/
void Declaration_1d_array(){
	x=new double[imax]; 	t=new double[imax], xc=new double[imax];	p=new double[imax]; Prim_var_L=new double[ndvar];
	Prim_var_M=new double[ndvar];	Prim_var_R=new double[ndvar];
}

/**************  to allocate 2d array***********/
void Declaration_2d_array(){
	cv=Allocate_2D(cv, ndvar, imax); cv_old=Allocate_2D(cv_old, ndvar, imax); L_state=Allocate_2D(L_state, ndvar, imax); 
	R_state=Allocate_2D(R_state, ndvar, imax); residue=Allocate_2D(residue, ndvar, imax); flux=Allocate_2D(flux, 3, imax);
	Dflux=Allocate_2D(Dflux, 3, imax); Aflux=Allocate_2D(Aflux, 3, imax); Dcv=Allocate_2D(Dcv, 3, imax);
	Acv=Allocate_2D(Acv, 3, imax);
}

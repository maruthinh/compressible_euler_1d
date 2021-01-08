#ifndef _Header_H
#define _Header_H

#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <string>	
#include <iosfwd>

/*****constants used throught the program******/
const double Gamma=1.4, cp=1005.0, pi=4.0*atan(1.0);
const int ndvar=3;
/*******SSPRK coefficients************/
static const double ConstRK1[]={0.0, 3.0/4.0, 1.0/3.0}, ConstRK2[]={1.0, 1.0/4.0, 2.0/3.0};
static const double ConstSSPRK2a[]={0.0, 1.0/2.0}, ConstSSPRK2b[]={1.0, 1.0/2.0};

/****global variables declared here, which will be used throught the function*****/

extern int N, ib, id1, id2, imax;
extern std::string space_accuracy, flux_type, time_accuracy, test_case;
extern double cfl;
extern double xll, xul, x0, x01, x02, t_fin;
extern double *Prim_var_L, *Prim_var_M, *Prim_var_R;
extern double *x, *t, *xc, **cv, **flux, **cv_old, *p, **residue, **Dflux, **Aflux, **Dcv, **Acv;
extern double dx, dt, GlobMaxEigenVal, GlobMinEigenVal;
extern double **L_state, **R_state;

/*********************forward declaration of functions called in the main program **********************/
void Read_Input();
void Read_Test_Case();
void Read_Test_Case_Two();

void Declaration_1d_array();
void Declaration_2d_array();


double *Grid(double xll, double xul, int N, int ib);
void Initial_Conditions(double x0, int ib, int imax, double *&x, double *Prim_var_L, double *Prim_var_R, double **&cv, double *&p);
void Initial_Conditions_Two(double x01, double x02, int ib, int imax, double *&x, double *Prim_var_L, double *Prim_var_M, double *Prim_var_R, double **&cv, double *&p);
void Initial_Conditions_Oscil(double x0, int ib, int imax, double *&x, double *Prim_var_L, double *Prim_var_R, double **&cv, double *&p);
void Boundary_Condition(int imax, int ib, double **&cv, double *&p);
void LR_State(int ib, int id, int imax, double *&x, double **&cv, double *&p, double **&L_state, double **&R_state);
void Flux_Con_Var_Defns(int ib, int imax, double **&cv, double *&p, double **&Dcv, double **&Acv, double **&Dflux,
                        double **&Aflux);
void Flux_LLF(int ib, int imax, double **L_state, double **R_state, double **residue);
double GlobMaxEigVal(int ib, int imax, double **&cv, double *&p);
double GlobMinEigVal(int ib, int imax, double **&cv, double *&p);
void Flux_KFDS(int ib, int id1, int id2, double **L_state, double **R_state, double **residue);
void Flux_KFDS_Movers(int ib, int imax, double **L_state, double **R_state, double **residue);
void Flux_Movers(int ib, int imax, double **L_state, double **R_state, double **residue);
void Flux_Movers_H(int ib, int imax, double **L_state, double **R_state, double **residue);
void Flux_Movers_LE(int ib, int imax, double **L_state, double **R_state, double **residue);
void Flux_Movers_CP(int ib, int imax, double **L_state, double **R_state, double **residue);
void Flux_KFDS2nd_Order(int ib, int imax, double **&cv, double *&p, double *&t, double **&residue);
void Flux_KFDS2nd_Movers(int ib, int id1, int id2, double **&cv, double **&Dcv, double **&Acv, double *&p, double *&t,
                         double **&flux, double **&Dflux, double **&Aflux, double **&residue);
double Time_Step(int ib, int imax, double cfl, double *&x, double **&cv, double *&p);
void Solver(int ib, int id1, int id2, double *&x, double **&cv, double **&cv_old, double *&p, double **&L_state,
            double **&R_state, double **&residue, double &dt);
void Solver_Explicit(int ib, int id1, int id2, double *&x, double **&cv, double **&cv_old, double *&p, double **&L_state, double **&R_state, double **&residue, double &dt);
void Write_Solution(int ib, double *&x, double **&cv, double *&p);

#endif  //#ifndef _Header_H

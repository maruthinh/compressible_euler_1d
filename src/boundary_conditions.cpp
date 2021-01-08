#include "global_declarations.h"
 void Boundary_Condition(int left_bound, int right_bound, double **&cv, double *&p){
	
	double u, ub;
	int left_inter, right_inter;

	left_inter  = left_bound  + 1;
	right_inter = right_bound - 1;

	if (test_case=="test_case10.txt") {

/****inlet boundary condition****/

        u = cv[1][left_inter] / cv[0][left_inter];
        ub = -u;

        cv[0][left_bound] = cv[0][left_inter];
        cv[1][left_bound] = cv[0][left_bound] * ub;
        cv[2][left_bound] = cv[2][left_inter];
            p[left_bound] = p[left_inter];

        for (int k = 0; k < 3; k++) {
            cv[k][left_bound-1] = cv[k][left_bound];
        }
        p[left_bound-1] = p[left_bound];

/****outlet boundary condition****/

        u = cv[1][right_inter] / cv[0][right_inter];
        ub = -u;

        cv[0][right_bound] = cv[0][right_inter];
        cv[1][right_bound] = cv[0][right_bound] * ub;
        cv[2][right_bound] = cv[2][right_inter];
            p[right_bound] = p[right_inter];

        for (int k = 0; k < 3; k++) {
            cv[k][right_bound + 1] = cv[k][right_bound];
        }
        p[right_bound + 1] = p[right_bound];
    }
	else {
/****inlet boundary condition****/
        for (int k = 0; k < 3; k++) {
            cv[k][left_bound]   = cv[k][left_inter];
            cv[k][left_bound-1] = cv[k][left_bound];
        }

        p[left_bound]   = p[left_inter];
        p[left_bound-1] = p[left_bound];

/****outlet boundary condition****/
        for (int k = 0; k < 3; k++) {
            cv[k][right_bound]     = cv[k][right_inter];
            cv[k][right_bound + 1] = cv[k][right_bound];
        }
        p[right_bound] = p[right_inter];
        p[right_bound + 1] = p[right_bound];
    }
}

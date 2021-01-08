#include "global_declarations.h"


double *Grid(double xll, double xul, int N, int ib){

	dx=(xul-xll)/N;

    for (int i=0;i<=id2;i++){
        x[i]=0.0;
    }

	for (int i=2;i<=ib;i++){
		x[i] = xll+(i-2)*dx;
    }

	x[1]   = x[2];
    x[0]   = x[1];
	x[id1] = x[ib];
	x[id2] = x[id1];

	/*for (int i=0;i<=id2;i++){
        std::cout<<"x["<<i<<"]="<<x[i]<<std::endl;
    }
    exit(0);*/
	return x;
}

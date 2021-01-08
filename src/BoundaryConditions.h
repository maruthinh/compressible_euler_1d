void BcProb1(){

	rho[1]=1.0;
	u[1]=0.75;
	p[1]=1.0;
	a[1]=sqrt((1.4*fabs(p[1]))/rho[1]);
	H[1]=(p[1]/rho[1])*(Gamma/(Gamma-1.0))+0.5*(u[1]*u[1]);
	E[1]=(p[1]/(rho[1]*(Gamma-1.0)))+0.5*u[1]*u[1];

/*	rho[N+2]=rho[N+1];
	u[N+2]=u[N+1];
	p[N+2]=p[N+1];
	E[N+2]=E[N+1];
	a[N+2]=a[N+1];
	H[N+2]=H[N+1];*/

	rho[N+2]=0.125;
	u[N+2]=0.0;
	p[N+2]=0.1;
	a[N+2]=sqrt((1.4*fabs(p[N+2]))/rho[N+2]);
	H[N+2]=(p[N+2]/rho[N+2])*(Gamma/(Gamma-1.0))+0.5*(u[N+2]*u[N+2]);
	E[N+2]=(p[N+2]/(rho[N+2]*(Gamma-1.0)))+0.5*u[N+2]*u[N+2];
/*
	rho[2]=rho[3];
	u[2]=u[3];
	p[2]=p[3];
	E[2]=E[3];
	a[2]=a[3];
	H[2]=H[3]; 

	rho[N+1]=rho[N];
	u[N+1]=u[N];
	p[N+1]=p[N];
	E[N+1]=E[N];
	a[N+1]=a[N];
	H[N+1]=H[N]; */
	
}

void BcProb2(){

	rho[1]=1.0;
	u[1]=-2.0;
	p[1]=0.4;
	a[1]=sqrt((1.4*fabs(p[1]))/rho[1]);
	H[1]=(p[1]/rho[1])*(Gamma/(Gamma-1.0))+0.5*(u[1]*u[1]);
	E[1]=(p[1]/(rho[1]*(Gamma-1.0)))+0.5*u[1]*u[1];

/*	rho[N+2]=rho[N+1];
	u[N+2]=u[N+1];
	p[N+2]=p[N+1];
	E[N+2]=E[N+1];
	a[N+2]=a[N+1];
	H[N+2]=H[N+1];*/

	rho[N+2]=1.0;
	u[N+2]=2.0;
	p[N+2]=0.4;
	a[N+2]=sqrt((1.4*fabs(p[N+2]))/rho[N+2]);
	H[N+2]=(p[N+2]/rho[N+2])*(Gamma/(Gamma-1.0))+0.5*(u[N+2]*u[N+2]);
	E[N+2]=(p[N+2]/(rho[N+2]*(Gamma-1.0)))+0.5*u[N+2]*u[N+2];
/*
	rho[2]=rho[3];
	u[2]=u[3];
	p[2]=p[3];
	E[2]=E[3];
	a[2]=a[3];
	H[2]=H[3]; 

	rho[N+1]=rho[N];
	u[N+1]=u[N];
	p[N+1]=p[N];
	E[N+1]=E[N];
	a[N+1]=a[N];
	H[N+1]=H[N]; */
	
}

void BcProb3(){

	rho[1]=1.0;
	u[1]=0.0;
	p[1]=1000.0;
	a[1]=sqrt((1.4*fabs(p[1]))/rho[1]);
	H[1]=(p[1]/rho[1])*(Gamma/(Gamma-1.0))+0.5*(u[1]*u[1]);
	E[1]=(p[1]/(rho[1]*(Gamma-1.0)))+0.5*u[1]*u[1];

/*	rho[N+2]=rho[N+1];
	u[N+2]=u[N+1];
	p[N+2]=p[N+1];
	E[N+2]=E[N+1];
	a[N+2]=a[N+1];
	H[N+2]=H[N+1];*/

	rho[N+2]=1.0;
	u[N+2]=0.0;
	p[N+2]=0.01;
	a[N+2]=sqrt((1.4*fabs(p[N+2]))/rho[N+2]);
	H[N+2]=(p[N+2]/rho[N+2])*(Gamma/(Gamma-1.0))+0.5*(u[N+2]*u[N+2]);
	E[N+2]=(p[N+2]/(rho[N+2]*(Gamma-1.0)))+0.5*u[N+2]*u[N+2];
/*
	rho[2]=rho[3];
	u[2]=u[3];
	p[2]=p[3];
	E[2]=E[3];
	a[2]=a[3];
	H[2]=H[3]; 

	rho[N+1]=rho[N];
	u[N+1]=u[N];
	p[N+1]=p[N];
	E[N+1]=E[N];
	a[N+1]=a[N];
	H[N+1]=H[N]; */
	
}

void BcProb4(){

	rho[1]=5.99924;
	u[1]=19.5975;
	p[1]=460.894;
	a[1]=sqrt((1.4*fabs(p[1]))/rho[1]);
	H[1]=(p[1]/rho[1])*(Gamma/(Gamma-1.0))+0.5*(u[1]*u[1]);
	E[1]=(p[1]/(rho[1]*(Gamma-1.0)))+0.5*u[1]*u[1];

/*	rho[N+2]=rho[N+1];
	u[N+2]=u[N+1];
	p[N+2]=p[N+1];
	E[N+2]=E[N+1];
	a[N+2]=a[N+1];
	H[N+2]=H[N+1];*/

	rho[N+2]=5.99924;
	u[N+2]=-6.19633;
	p[N+2]=40.0950;
	a[N+2]=sqrt((1.4*fabs(p[N+2]))/rho[N+2]);
	H[N+2]=(p[N+2]/rho[N+2])*(Gamma/(Gamma-1.0))+0.5*(u[N+2]*u[N+2]);
	E[N+2]=(p[N+2]/(rho[N+2]*(Gamma-1.0)))+0.5*u[N+2]*u[N+2];
/*
	rho[2]=rho[3];
	u[2]=u[3];
	p[2]=p[3];
	E[2]=E[3];
	a[2]=a[3];
	H[2]=H[3]; 

	rho[N+1]=rho[N];
	u[N+1]=u[N];
	p[N+1]=p[N];
	E[N+1]=E[N];
	a[N+1]=a[N];
	H[N+1]=H[N]; */
	
}

void BcProb5(){

	rho[1]=1.0;
	u[1]=-19.59745;
	p[1]=1000.0;
	a[1]=sqrt((1.4*fabs(p[1]))/rho[1]);
	H[1]=(p[1]/rho[1])*(Gamma/(Gamma-1.0))+0.5*(u[1]*u[1]);
	E[1]=(p[1]/(rho[1]*(Gamma-1.0)))+0.5*u[1]*u[1];

/*	rho[N+2]=rho[N+1];
	u[N+2]=u[N+1];
	p[N+2]=p[N+1];
	E[N+2]=E[N+1];
	a[N+2]=a[N+1];
	H[N+2]=H[N+1];*/

	rho[N+2]=1.0;
	u[N+2]=-19.59745;
	p[N+2]=0.01;
	a[N+2]=sqrt((1.4*fabs(p[N+2]))/rho[N+2]);
	H[N+2]=(p[N+2]/rho[N+2])*(Gamma/(Gamma-1.0))+0.5*(u[N+2]*u[N+2]);
	E[N+2]=(p[N+2]/(rho[N+2]*(Gamma-1.0)))+0.5*u[N+2]*u[N+2];
/*
	rho[2]=rho[3];
	u[2]=u[3];
	p[2]=p[3];
	E[2]=E[3];
	a[2]=a[3];
	H[2]=H[3]; 

	rho[N+1]=rho[N];
	u[N+1]=u[N];
	p[N+1]=p[N];
	E[N+1]=E[N];
	a[N+1]=a[N];
	H[N+1]=H[N]; */
	
}

void BcProb6(){

	rho[1]=1.4;
	u[1]=0.0;
	p[1]=1.0;
	a[1]=sqrt((1.4*fabs(p[1]))/rho[1]);
	H[1]=(p[1]/rho[1])*(Gamma/(Gamma-1.0))+0.5*(u[1]*u[1]);
	E[1]=(p[1]/(rho[1]*(Gamma-1.0)))+0.5*u[1]*u[1];

/*	rho[N+2]=rho[N+1];
	u[N+2]=u[N+1];
	p[N+2]=p[N+1];
	E[N+2]=E[N+1];
	a[N+2]=a[N+1];
	H[N+2]=H[N+1];*/

	rho[N+2]=1.0;
	u[N+2]=0.0;
	p[N+2]=1.0;
	a[N+2]=sqrt((1.4*fabs(p[N+2]))/rho[N+2]);
	H[N+2]=(p[N+2]/rho[N+2])*(Gamma/(Gamma-1.0))+0.5*(u[N+2]*u[N+2]);
	E[N+2]=(p[N+2]/(rho[N+2]*(Gamma-1.0)))+0.5*u[N+2]*u[N+2];
/*
	rho[2]=rho[3];
	u[2]=u[3];
	p[2]=p[3];
	E[2]=E[3];
	a[2]=a[3];
	H[2]=H[3]; 

	rho[N+1]=rho[N];
	u[N+1]=u[N];
	p[N+1]=p[N];
	E[N+1]=E[N];
	a[N+1]=a[N];
	H[N+1]=H[N]; */
	
}

void BcProb8(){
	double M1=2.0;

	rho[1]=1.0;
	u[1]=1.0;
	p[1]=1.0/(Gamma*M1*M1);
	a[1]=sqrt((1.4*fabs(p[1]))/rho[1]);
	H[1]=(p[1]/rho[1])*(Gamma/(Gamma-1.0))+0.5*(u[1]*u[1]);
	E[1]=(p[1]/(rho[1]*(Gamma-1.0)))+0.5*u[1]*u[1];

/*	rho[N+2]=rho[N+1];
	u[N+2]=u[N+1];
	p[N+2]=p[N+1];
	E[N+2]=E[N+1];
	a[N+2]=a[N+1];
	H[N+2]=H[N+1];*/

	rho[N+2]=1/((2.0/((Gamma+1.0)*M1*M1))+((Gamma-1.0)/(Gamma+1.0)));
	u[N+2]=1/rho[N+2];
	p[N+2]=(1.0/(Gamma*M1*M1))*(((2.0*Gamma*M1*M1)/(Gamma+1.0))-((Gamma-1.0)/(Gamma+1.0)));
	a[N+2]=sqrt((1.4*fabs(p[N+2]))/rho[N+2]);
	H[N+2]=(p[N+2]/rho[N+2])*(Gamma/(Gamma-1.0))+0.5*(u[N+2]*u[N+2]);
	E[N+2]=(p[N+2]/(rho[N+2]*(Gamma-1.0)))+0.5*u[N+2]*u[N+2];

/*	rho[2]=rho[3];
	u[2]=u[3];
	p[2]=p[3];
	E[2]=E[3];
	a[2]=a[3];
	H[2]=H[3]; 

	rho[N+1]=rho[N];
	u[N+1]=u[N];
	p[N+1]=p[N];
	E[N+1]=E[N];
	a[N+1]=a[N];
	H[N+1]=H[N]; */
	
}

void BcProb9(){
	double M1=2.0;

	rho[1]=3.86;
	u[1]=-3.1266/rho[1];
	p[1]=(Gamma-1.0)*(27.0913-0.5*rho[1]*u[1]*u[1]);
	a[1]=sqrt((1.4*fabs(p[1]))/rho[1]);
	H[1]=(p[1]/rho[1])*(Gamma/(Gamma-1.0))+0.5*(u[1]*u[1]);
	E[1]=(p[1]/(rho[1]*(Gamma-1.0)))+0.5*u[1]*u[1];

/*	rho[N+2]=rho[N+1];
	u[N+2]=u[N+1];
	p[N+2]=p[N+1];
	E[N+2]=E[N+1];
	a[N+2]=a[N+1];
	H[N+2]=H[N+1];*/

	rho[N+2]=1.4;
	u[N+2]=-3.44/rho[N+2];
	p[N+2]=(Gamma-1.0)*(8.4168-0.5*rho[N+2]*u[N+2]*u[N+2]);
	a[N+2]=sqrt((1.4*fabs(p[N+2]))/rho[N+2]);
	H[N+2]=(p[N+2]/rho[N+2])*(Gamma/(Gamma-1.0))+0.5*(u[N+2]*u[N+2]);
	E[N+2]=(p[N+2]/(rho[N+2]*(Gamma-1.0)))+0.5*u[N+2]*u[N+2];

/*	rho[2]=rho[3];
	u[2]=u[3];
	p[2]=p[3];
	E[2]=E[3];
	a[2]=a[3];
	H[2]=H[3]; 

	rho[N+1]=rho[N];
	u[N+1]=u[N];
	p[N+1]=p[N];
	E[N+1]=E[N];
	a[N+1]=a[N];
	H[N+1]=H[N]; */
	
}

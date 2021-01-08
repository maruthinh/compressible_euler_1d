template <typename T> //t=0.25
void InitialCondition0(T *&xc){
	for(int i=2;i<=N+1;i++){
		if(xc[i]<0.3){
			rho[i]=1.0;
			u[i]=0.0;
			p[i]=1.0;
			
		}
		else{
			rho[i]=0.125;
			u[i]=0.0;
			p[i]=0.1;
		}
	}
	
	for (int i=0;i<=1;i++){
		rho[i]=1.0;
		u[i]=0.0;
		p[i]=1.0;
	}
	
	for (int i=N+2;i<=N+3;i++){
		rho[i]=0.125;
		u[i]=0.0;
		p[i]=0.1;
	}

	for(int i=0;i<=N+3;i++){
		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];

		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
	}
		
}


template <typename T> //t=0.2
void InitialCondition1(T *&xc){
	for(int i=2;i<=N+1;i++){
		if(xc[i]<0.3){
			rho[i]=1.0;
			u[i]=0.75;
			p[i]=1.0;
			
		}
		else{
			rho[i]=0.125;
			u[i]=0.0;
			p[i]=0.1;
		}

	}

	for (int i=0;i<=1;i++){
		rho[i]=1.0;
		u[i]=0.75;
		p[i]=1.0;
	}
	
	for (int i=N+2;i<=N+3;i++){
		rho[i]=0.125;
		u[i]=0.0;
		p[i]=0.1;
	}

	for(int i=0;i<=N+3;i++){
		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];

		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
	}
}

template <typename T> //t=0.15
void InitialCondition2(T *&xc){
	for(int i=2;i<=N+1;i++){
		if(xc[i]<0.5){
			rho[i]=1.0;
			u[i]=-2.0;
			p[i]=0.4;
			
		}
		else{
			rho[i]=1.0;
			u[i]=2.0;
			p[i]=0.4;
		}

	}

	for (int i=0;i<=1;i++){
		rho[i]=1.0;
		u[i]=-2.0;
		p[i]=0.4;
	}
	
	for (int i=N+2;i<=N+3;i++){
		rho[i]=1.0;
		u[i]=2.0;
		p[i]=0.4;
	}

	for(int i=0;i<=N+3;i++){
		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];

		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
	}
}

template <typename T> //t=0.012
void InitialCondition3(T *&xc){
	for(int i=2;i<=N+1;i++){
		if(xc[i]<0.5){
			rho[i]=1.0;
			u[i]=0.0;
			p[i]=1000.0;
			
		}
		else{
			rho[i]=1.0;
			u[i]=0.0;
			p[i]=0.01;
		}

	}

	for (int i=0;i<=1;i++){
		rho[i]=1.0;
		u[i]=0.0;
		p[i]=1000.0;
	}
	
	for (int i=N+2;i<=N+3;i++){
		rho[i]=1.0;
		u[i]=0.0;
		p[i]=0.01;
	}

	for(int i=0;i<=N+3;i++){
		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];

		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
	}

}

template <typename T> //t=0.035
void InitialCondition4(T *&xc){
	for(int i=2;i<=N+1;i++){
		if(xc[i]<0.4){
			rho[i]=5.99924;
			u[i]=19.5975;
			p[i]=460.894;
			
		}
		else{
			rho[i]=5.99924;
			u[i]=-6.19633;
			p[i]=40.0950;
		}

	}

	for (int i=0;i<=1;i++){
		rho[i]=5.99924;
		u[i]=19.5975;
		p[i]=460.894;
	}
	
	for (int i=N+2;i<=N+3;i++){
		rho[i]=5.99924;
		u[i]=-6.19633;
		p[i]=40.0950;
	}

	for(int i=0;i<=N+3;i++){
		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];

		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
	}

}


template <typename T> //t=0.012
void InitialCondition5(T *&xc){
	for(int i=2;i<=N+1;i++){
		if(xc[i]<0.8){
			rho[i]=1.0;
			u[i]=-19.59745;
			p[i]=1000.0;
			
		}
		else{
			rho[i]=1.0;
			u[i]=-19.59745;
			p[i]=0.01;
		}
	}

	for (int i=0;i<=1;i++){
		rho[i]=1.0;
		u[i]=-19.59745;
		p[i]=1000.0;
	}
	
	for (int i=N+2;i<=N+3;i++){
		rho[i]=1.0;
		u[i]=-19.59745;
		p[i]=0.01;
	}

	for(int i=0;i<=N+3;i++){
		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];

		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
	}
}

template <typename T> //t=2.0
void InitialCondition6(T *&xc){

	for(int i=0;i<=N+3;i++){
		if(i<((N+3)/2)+1){
			rho[i]=1.4;
			u[i]=0.0;
			p[i]=1.0;
			
		}
		else{
			rho[i]=1.0;
			u[i]=0.0;
			p[i]=1.0;
		}

		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];

		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
		
	}

}

template <typename T> //t=0.17
void InitialCondition7(T *&xc){
	for(int i=2;i<=N+1;i++){
		if(xc[i]<0.5){
			rho[i]=1.0;
			u[i]=0.0;
			p[i]=1.0;
			
		}
		else{
			rho[i]=0.125;
			u[i]=0.0;
			p[i]=0.1;
		}

		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];

		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
		
	}
}

template <typename T> //t=0.012
void InitialCondition8(T *&xc){
	double M1=2.0;

	for(int i=0;i<=(N+3)/2;i++){
		rho[i]=1.0;
		u[i]=1.0;
		p[i]=1.0/(Gamma*M1*M1); //0.17857142857142858
		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];
		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
			
	}

	for(int i=((N+3)/2)+1;i<=N+3;i++){
		rho[i]=1/((2.0/((Gamma+1.0)*M1*M1))+((Gamma-1.0)/(Gamma+1.0)));//2.6666666666666665
		u[i]=1/rho[i]; //0.375
		p[i]=(1.0/(Gamma*M1*M1))*(((2.0*Gamma*M1*M1)/(Gamma+1.0))-((Gamma-1.0)/(Gamma+1.0))); //0.8035714285714286
		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];

		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
	}


		
	
}

template <typename T> //t=0.012
void InitialCondition9(T *&xc){

	for(int i=2;i<=N+1;i++){
		if(xc[i]<0.5){
			rho[i]=3.86;
			u[i]=-3.1266/rho[i];
			p[i]=(Gamma-1.0)*(27.0913-0.5*rho[i]*u[i]*u[i]);
			
		}
		else{
			rho[i]=1.4;
			u[i]=-3.44/rho[i];
			p[i]=(Gamma-1.0)*(8.4168-0.5*rho[i]*u[i]*u[i]);
		}

		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];

		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
		
	}
}


template <typename T> //t=2.0 //slowly moving contact page no 344 from Toro book
void InitialCondition10(T *&xc){
	for(int i=2;i<=N+1;i++){
		if(xc[i]<0.5){
			rho[i]=1.4;
			u[i]=0.1;
			p[i]=1.0;
			
		}
		else{
			rho[i]=1.0;
			u[i]=0.1;
			p[i]=1.0;
		}

		a[i]=sqrt(Gamma*p[i]/rho[i]);
		E[i]=(p[i]/(rho[i]*(Gamma-1.0)))+0.5*u[i]*u[i];
		H[i]=(p[i]/rho[i])+E[i];
		
		U1[i]=rho[i];
		U2[i]=rho[i]*u[i];
		U3[i]=rho[i]*E[i];
	}
}

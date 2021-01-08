template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S;
	const double epsilon=1e-6;

	if (fabs(Fr-Fl)<epsilon) return 0.0;
	else if (fabs(Ur-Ul)<epsilon) return L_min;
	else if (fabs(Ur-Ul)>epsilon and fabs(Fr-Fl)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));
	else S=L_min;
	
	if (S<epsilon)	return 0;
	else if ((S)>=L_max) return (L_max);
	else if((S)<=L_min) return (L_min);
	else return (S);
	return S;
}

/*
template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S;

	if (Ur!=Ul) S=fabs(((Fr-Fl)/(Ur-Ul)));
	else S=fabs(L_min);
	
	if (fabs(S)>=L_max) return fabs(L_max);
	else if (fabs(S)<=L_min) return fabs(L_min);
	else return fabs(S);
}

*/
/*
template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S, S1, delta;
	
	delta=0.5*L_max;

	if (Ur!=Ul) S=fabs(((Fr-Fl)/(Ur-Ul)));
	else S=fabs(L_min);
	
	if (fabs(S)>=L_max) S1=fabs(L_max);
	else if (fabs(S)<=L_min) S1=abs(L_min);
	else S1=fabs(S);
	
	if (fabs(S1)<delta) return (S1*S1+delta*delta)/(2.0*delta);
	else return S1;
}
*/
/*
template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S;
	const double epsilon=1e-6;

	if (fabs(Fr-Fl)<epsilon) return 0.0;
	else if (fabs(Ur-Ul)<epsilon) return (L_min+L_max)/2;
	else if (fabs(Ur-Ul)>epsilon and fabs(Fr-Fl)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));
	else S=(L_min+L_max)/2;
	
	if (S<0.00001)	return S;
	else if (fabs(S)>=L_max) return fabs(L_min+L_max)/2;
	else if(fabs(S)<=L_min) return fabs(L_min+L_max)/2;
	else return fabs(S);
	return S;
}


template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S;
	if ((fabs(Ur)-fabs(Ul))<0.000001 && (fabs(Fr)-fabs(Fl))<0.000001) return fabs(L_min+L_max)/2.0;
	else if ((Ur!=Ul) && (fabs(Fr)-fabs(Fl))<0.000001) return 0;
	else if ((Ur!=Ul) && (Fr!=Fl)) S=fabs(((Fr-Fl)/(Ur-Ul)));
	else if (Ur==Ul) return fabs(L_min+L_max)/2.0;
	
	if (fabs(S)>=L_max) return fabs(L_max);
	else if(fabs(S)<=L_min) return fabs(L_min+L_max)/2.0;
	else return fabs(S);
	return S;
}
*/

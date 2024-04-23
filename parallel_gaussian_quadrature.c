#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>

//funzione da integrare
float function(float x){
	float y;	
	y = cos(x) - x * exp(x);
	return y;
}

//funzione che confronta termini per ordinare array
int cmpfunc (const void * a, const void * b)
{
  if (*(double*)a > *(double*)b)
    return 1;
  else if (*(double*)a < *(double*)b)
    return -1;
  else
    return 0;  
}

/*Nth Legendre Polynomial Pn(x)*/
double Pn( int n ,double x )
{
	double r,s,t ;
	int m;
	r=0; s=1;
	
	for (m=0; m<n; m++)
	{
		t=r; 
		r=s;
		s=(2*m+1 )*x*r-m*t;
		s/=(m+1);
	}
	return s;
} 

/*Lagrange terms*/
long double Li(int n, double x[n+1], int i, double X){
    int j;
    long double prod=1;
    for(j=0;j<=n;j++){
        if (j!=i){
            prod=prod*(X-x[j])/(x[i]-x[j]);     
        }
    }
    return prod;
}

/*Funzione per definire l'integrazione con Simpson's 1/3rd Rule */
long double Ci(int i, int n, double x[n], double a, double b, int N){
  long double h,integral,X,sum=0;
  int j,k;
  h=(b-a)/N;
  for(j=1;j<N;j++){
    X=a+j*h;
    if(j%2==0){
      sum=sum+2*Li(n-1,x,i,X);
    }
    else{
      sum=sum+4*Li(n-1,x,i,X);;
    }
  }
    long double Fa=Li(n-1,x,i,a);;
    long double Fb=Li(n-1,x,i,b);

  integral=(h/3.0)*(Fa+Fb+sum);
  return integral;
}

/*Funzione del metodo della bisezione[Return la radice quando la trova o 999 altrimenti]*/
double bisection(int n,double f(int n,double x),double a, double b, double eps, int maxSteps){
  double c;
  if(f(n,a)*f(n,b)<=0){  
    int iter=1;
	
    do{
      c=(a+b)/2;
      if(f(n,a)*f(n,c)>0){
	  a=c;
	}
	else if(f(n,a)*f(n,c)<0){
	  b=c;
	}
	else if(f(n,c)==0){
		return c;
	}
      iter++;
	      
    }while(fabs(a-b)>=eps&&iter<=maxSteps);
    return c;
  }
  else{
    return 999;
  }
} 
int main (int argc, char** argv)
{
	MPI_Status status;
	MPI_Request req;
	//utile per calcolare il tempo d'esecuzione
	double time;
	//MPI size e rank 
	int size, rank;
	//buffer per send e recv
	double buff_send,buff_recv,buff_recv1;
	//estremi di integrazione per il polinomio di legendre
	int la = -1, lb = 1;
	//estremi di integrazione della function ricercata
	double a, b, temp;
	//variabili per salvare (b-a)/2 e (b+a)/2
   	 double d, e;
	//Utili per differenziare in base al rank l'intervallo di lavoro per processso
	double my_a, my_b;
	int i=0,n;
	//per salvare risultati parziali e finali dei singoli processi
	double my_res = 0,result;
	//variabile utile per cicli for
	int y = 0 ;
	//window(Step-size) per bisection method
	double h=0.00001;
	//variable utile per bisection method
	double x;
	//variabile dove viene ritornata la radice dopo il bisection loop 
	double root;
	
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	
	if (argc != 4 )
	{
	if (rank == 0)
	{
	
		printf("Specificare i data point <n> e gli estremi di integrazione <a> <b>\n");
		return 1;
	}	
	}
	
	
	n = strtof(argv[1], NULL);	
	a = strtof(argv[2], NULL);
	b = strtof(argv[3], NULL);
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	time = -MPI_Wtime();
	
	//se il primo estremo di integrazione è minore del secondo, scambio le variabili
	if (a>b) 
	{
		temp = a;
		a = b;
		b = temp;
	}
	
	d=((b-a)/2);
	e=((b+a)/2);

	
	//Array per salvare le radici del Legendre polynomials
	double xi[n];
	
	if (rank == 0)
	{
		printf("\n\nLe Radici (xi's) sono:\n\n\n");
	}
			
	double w = (lb - la)/h;

	my_a = (rank) * (w/(size));
	my_b = my_a + (w/(size));
		
	double my_a_n =my_a*h -1;
	double my_b_n =my_b*h -1;
	
	for (x = my_a_n; x < my_b_n; x += h) 
	{
       	root=bisection(n,Pn,x,x+h,0.0000001,1000000);
       	//printf("Trovata\n");
      		if(root!=999)    	      	
      		{        	
      			if (rank == 0)
           	  	 {
				MPI_Irecv(&buff_recv, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0 , MPI_COMM_WORLD,&req);
			 }

      		        buff_send = root;		
           	        MPI_Send(&buff_send, 1, MPI_DOUBLE, 0 ,0, MPI_COMM_WORLD);
           		
 		 	if (rank == 0)
     			 {			
				MPI_Wait(&req, MPI_STATUS_IGNORE);	
    				xi[i]=buff_recv;
    			 }
              
               	 i++;
        	 }
	}
	
	//ricevo le n-i restanti radici	
	if (rank == 0)
	{		
    	    for(y = i ; y < n; y++ )
    	    {    				
    			  MPI_Recv(&buff_recv1, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0 , MPI_COMM_WORLD,&status);
    			  xi[y]=buff_recv1;   
	    }
     	//ordino array		
     	   qsort (xi, n, sizeof(double), cmpfunc);
           printf("\nOrdino l'array\n\n");
   	   for( int y = 0 ; y < n; y++ )
    	   {
    		    printf(" x[%d] = %f\n",y+1, xi[y]);
   	   }
		printf("\n__________________________________________\n");
 		printf("\n\nI pesi (ci's) sono:\n\n");
       }	
	 
	//Invio a tutti i processi l'array ordinato		
      	MPI_Bcast(xi, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				
	my_a = (rank) * (n/(size));
	my_b = my_a + (n/(size)); 
		
	for(i=my_a;i<my_b;i++)
	{
        	double my_root = xi[i];
		double my_peso = Ci(i,n,xi,-1,1,100000);
		printf("%d) Peso c[%d] = %7.6lf\n",rank,i+1,my_peso);
		double t = my_peso*function((d*my_root) + e);
		my_res += t;
		
   	}
	//residuo
	int rankn = size-1;
	
	if(rank==rankn) 
   	{
        	for(i=my_b;i<n;i++)
        	{
        			double my_root = xi[i];
				double my_peso = Ci(i,n,xi,-1,1,100000);
				printf("%d) Peso c[%d] = %7.6lf\n",rank,i+1,my_peso);
				double t = my_peso*function((d*my_root) + e);
				
				my_res += t;	
		}
   	}	
   		//MPI_reduce(MPI_SUM) per sommare tutti i singoli contributi
   		MPI_Reduce(&my_res, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0){
	 	result = d * result;
	//	Arrotondo Risultato 	
	//ceilif restituisce un valore a virgola mobile che rappresenta il numero intero più piccolo maggiore o uguale a x
	 	float rounded_up = ceilf(result * 100) / 100; 
	 	printf("\n__________________________________________\n");
		printf("\nRisultato Integrale: %.2f\n", rounded_up);
		printf("__________________________________________\n\n");
	}
		
 	//Tempo
 	MPI_Barrier(MPI_COMM_WORLD);
	time += MPI_Wtime();
	if (rank == 0)
	{
		printf("Wall clock time %f\n", time);
	}
	MPI_Finalize();
}

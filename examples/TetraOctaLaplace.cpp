/*======================================================================================
Rich Brower Mon Dec 26 19:22:00 EST 2016
Laplacian on a "Tetrahedral-octahedral honeycomb" for an L*L*L torus.
https://en.wikipedia.org/wiki/Tetrahedral-octahedral_honeycomb

Two constants c1 and c2:
c1 for 12 nearest neighbor  even to even  (0,+/-,+/-) 4*3 = 12
c2 for 6  nearest neighbor  odd to even.  (0,0,+/-1)  3*2 = 6

Even-Odd lattice a= 1  even-even d = \sqrt[2}
Vol_Tetra = d^3/(6 Sqrt[2]) = 1/3
Vol_Oct_Tetra= 1/6 
Sum  on unit cube = a^3 =  1/3 + 4*(1/6) = 1

DEC weights c1 =     c2  = 
FEM weigths c1 =     c2  =

======================================================================================*/

#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

#define L 48
#define maxIter 1000
#define resEpsilon  0.0001

// index x = i%L , y = (i/L)%L, z = (i/(L*L))%L
// i = x + y*L + z*L*L

inline int even(int i)
{
  return ( i%L + (i/L)%L + (i/(L*L))%L) % 2;
}

// tetra is even from even: pairs 0:3, 1:2 & 4:7, 5:6 & 8:11, 9:10
 
inline int tetraOld(int i, int mu)
{int imu = -1;
switch(mu)
  {
case 0: imu = (i+1)%L + ((i/L +1)%L)*L + ((i/L/L)%L)*L*L;  break;
case 1: imu = (i-1+L)%L + ((i/L +1)%L)*L + ((i/L/L)%L)*L*L;  break;
case 2: imu = (i+1)%L + ((i/L -1 +L)%L)*L + ((i/L/L)%L)*L*L;  break;
case 3: imu = (i-1+L)%L + ((i/L -1 +L)%L)*L + ((i/L/L)%L)*L*L;  break;
      
case 4: imu = (i+1)%L +  ((i/L)%L)*L + ((i/L/L +1)%L)*L*L;  break;
case 5: imu = (i-1+L)%L +  ((i/L)%L)*L + ((i/L/L +1)%L)*L*L;  break;
case 6: imu = (i+1)%L +   ((i/L)%L)*L + ((i/L/L -1 +L)%L)*L*L;  break;
case 7: imu = (i-1+L)%L +  ((i/L)%L)*L + ((i/L/L -1 +L)%L)*L*L;  break;

case 8: imu = i%L + ((i/L+1)%L)*L + ((i/L/L +1)%L)*L*L;  break;
case 9: imu =i%L +  ((i/L-1+L)%L)*L + ((i/L/L +1)%L)*L*L;  break;
case 10: imu = i%L +  ((i/L+1)%L)*L + ((i/L/L -1 +L)%L)*L*L;  break;
case 11: imu = i%L +  ((i/L-1+L)%L)*L + ((i/L/L -1 +L)%L)*L*L;  break;

default: printf("Tetra index out of range\n");  break;
}
return imu;
}


//0:3, 1:2 & 4:7, 5:6 & 8:11, 9:10

inline int tetra(int i, int mu)
{
int imu = -1;
switch(mu)
  {
case 0: imu = (i+1)%L + ((i/L +1)%L)*L + ((i/L/L)%L)*L*L;  break;      //0
case 1: imu = (i-1+L)%L + ((i/L +1)%L)*L + ((i/L/L)%L)*L*L;  break;    //1
case 2: imu = (i+1)%L +  ((i/L)%L)*L + ((i/L/L +1)%L)*L*L;  break;     //4
case 3: imu = (i-1+L)%L +  ((i/L)%L)*L + ((i/L/L +1)%L)*L*L;  break;   //5
case 4: imu = i%L + ((i/L+1)%L)*L + ((i/L/L +1)%L)*L*L;  break;        //8
case 5: imu =i%L +  ((i/L-1+L)%L)*L + ((i/L/L +1)%L)*L*L;  break;      //9

case 6: imu = (i-1+L)%L + ((i/L -1 +L)%L)*L + ((i/L/L)%L)*L*L;  break;  //3
case 7: imu = (i+1)%L + ((i/L -1 +L)%L)*L + ((i/L/L)%L)*L*L;  break;    //2
case 8: imu = (i-1+L)%L +  ((i/L)%L)*L + ((i/L/L -1 +L)%L)*L*L;  break; //7
case 9: imu = (i+1)%L +   ((i/L)%L)*L + ((i/L/L -1 +L)%L)*L*L;  break;  //6
case 10: imu = i%L +  ((i/L-1+L)%L)*L + ((i/L/L -1 +L)%L)*L*L;  break;  //11
case 11: imu = i%L +  ((i/L+1)%L)*L + ((i/L/L -1 +L)%L)*L*L;  break;    //10

default: printf("Tetra index out of range\n");  break;
}
return imu;
}

// tetra/octa center: even from odd and odd form even
inline int octa(int i, int mu)
{int imu = -1;
  switch(mu)
    {
    case 0: imu = (i+1)%L + ((i/L)%L)*L + ((i/L/L)%L)*L*L ;  break;
    case 1: imu = (i-1+L)%L + ((i/L)%L)*L + ((i/L/L)%L)*L*L;  break;
      
    case 2: imu = i%L +  ((i/L + 1)%L)*L + ((i/L/L)%L)*L*L;  break;
    case 3: imu =  i%L +  ((i/L -1 +L)%L)*L + ((i/L/L)%L)*L*L;  break;
      
    case 4: imu = i%L +   ((i/L)%L)*L + ((i/(L*L) +1)%L)*L*L;  break;
    case 5: imu =  i%L +  ((i/L)%L)*L + ((i/(L*L) -1 +L)%L)*L*L;  break;
      
    default: printf("Octa index out of range\n");  break;
    }
  return imu;
}


inline int mod(int x, int n)
{
  return  (n + (x % n))%n;
}

struct parameters
{
 double K[6];
  double rmin;
  double rmax;
}p;
  
void Aphi(double *phi2, const double *phi1);
double Ainv_phi(double  *phi, const double  *b, const double  *phi0);

int main()
{
  double truersq = 0.0;
  double *phi  = (double *) malloc(L*L*L*sizeof(double));
  double *phi0 = (double *) malloc(L*L*L*sizeof(double));
  double *b = (double *) malloc(L*L*L*sizeof(double));
  p.rmin = 2.0;
  p.rmax = L/4.0;
  /**** 
  p.K[0] =  1.0/6.0;
  p.K[1] =  1.0/6.0;
  p.K[2] =  1.0/6.0;
  p.K[3] =  0.75/6.0;
  p.K[4] =  0.75/6.0;
  p.K[5] =  0.75/6.0;
  *****/

  p.K[0] =  1.0/6.0;
  p.K[1] =  1.0/6.0;
  p.K[2] =  1.0/6.0;
  p.K[3] =  1.0/6.0;
  p.K[4] =  1.0/6.0;
  p.K[5] =  1.0/6.0;
  
  ofstream Corr_out, Cut_out;
  Corr_out.open("../data/Corr.dat");
   Cut_out.open("../data/Cut.dat");
   
  for(int i = 0; i<L*L*L; i++) {
    phi[i] = 0.0; phi0[i] = 0.0; b[i] = 0.0;
  }
  b[L/2 + L*(L/2) + L*L*(L/2)] = 10.0;
  
  
  int x,y,z;
  // Testing Aphi
#if 0  
  Aphi(phi0, b);
  for(int i = 0; i<L*L*L; i++) {
    x = i%L ; y = (i/L)%L; z = (i/(L*L))%L;
    cout << endl <<" phi0[i] = " <<  phi0[i] <<  "  at (x,y,z) = (" << x <<","  << y <<","   << z <<")" << endl;
  }
  for(int i = 0; i<L*L*L; i++) {
    phi0[i] = 0.0;
  }
  
#endif
  
  truersq = Ainv_phi(phi, b, phi0);
  cout <<"  truersq =  " <<  truersq  << endl;
  // Headers 
  // Cut_out << "#   rmin  = " << p.rmin << " rmax = "<< p.rmax << endl;
  
  for(int i = 0; i<L*L*L; i++) {
    x = i%L ; y = (i/L)%L; z = (i/(L*L))%L;
    if( (x + y + z)%2 ==0)
      {
	int  rsqr =  (x-L/2)*(x-L/2)  +  (y-L/2)* (y-L/2) +   (z-L/2)* (z-L/2);
	
	if(rsqr > 0)
	  { Corr_out  <<"     "   << sqrt(rsqr) <<   "   "  << phi[i] << "   " << abs(x-L/2) <<  "   " << abs(y-L/2) << "   " << abs(z-L/2) << endl;}
       	
	if(rsqr > p.rmin*p.rmin && rsqr < p.rmax*p.rmax)
	  {  Cut_out  <<"     "   << sqrt(rsqr) <<   "   "  <<phi[i] << "   " << abs(x-L/2) <<  "   " << abs(y-L/2) << "   " << abs(z-L/2) << endl;}
      }
  }


  
   for(int i = 0; i<L*L*L; i++)
     cout << i << "  " << phi[i] << endl;
  
   //  for(int i = 0; i<L*L*L; i++) cout << " i = " << i << "even = "<< even(i) << endl;
	   
   free(phi);
   free(phi0);
   free(b);

   Corr_out.close();
    Cut_out.close();
  return 0;
}

void Aphi(double *phi2, const double *phi1)
{
  // double c1 = 1.0/6.0;  //even to even
  //  double c2 = 2.0/3.0; //even to odd
  double c2 = 0;
   //double c2 = 1.0/3.0;
  double m = 1.0/(double)L;
  
  for(int i = 0; i < L*L*L; i++)
    {
      phi2[i] = 0;
       if(even(i)==0)
	{
	  for(int mu =0; mu < 6; mu++) //even to even
	    phi2[i] +=  p.K[mu]*( phi1[i] - phi1[tetra(i,mu)] );
	  
	  for(int mu =6; mu < 12; mu++) //even to even
	    phi2[i] +=  p.K[mu - 6]*( phi1[i] - phi1[tetra(i,mu)] );
	}
      
      phi2[i] += m*m*phi1[i];
      if(c2 != 0){
           for(int mu =0; mu < 6; mu++)  //even to odd and odd to even
      	phi2[i] += c2*(phi1[i] - phi1[octa(i,mu)] );
      }
    }
}


double Ainv_phi(double  *phi, const double  *b, const double  *phi0)
{
// CG solutions to Aphi = b 
//  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  int latVol = L*L*L;
  double  res[latVol], resNew[latVol],  pvec[latVol], Apvec[latVol];
  double alpha, beta, denom ;
  double rsq = 0, rsqNew = 0, bsqrt = 0, truersq=0.0;
  int k, i;


  for(i = 0; i<latVol; i++) {
    res[i] = b[i];  // for now assume start with phi0 = 0
    resNew[i] = 0.0;
    pvec[i] = res[i];
    Apvec[i] = 0.0;
    bsqrt += b[i]*b[i];
  }
  bsqrt = sqrt(bsqrt);

  // iterate till convergence
  for(k = 0; k< maxIter; k++) {
    rsq = 0;
    for (i = 0; i < latVol; i++) rsq += res[i]*res[i];
     cout << endl << " rsq = " << rsq << endl;

    Aphi(Apvec, pvec);  
    denom = 0;
    for(i=0; i< latVol; i++) denom += pvec[i]*Apvec[i];
    alpha = rsq/denom;


    for(i=0; i < latVol; i++)  phi[i] +=  alpha * pvec[i];
    for(i=0; i < latVol; i++) resNew[i] = res[i]- alpha*Apvec[i];

    // Exit if new residual is small enough
    rsqNew = 0;
    for (i = 0; i < latVol; i++) rsqNew += resNew[i]*resNew[i];
   

    if (sqrt(rsqNew) < resEpsilon*bsqrt) {
     	printf("Final rsq = %g\n", rsqNew);
      break;
    }

    // Update vec using new residual
    beta = rsqNew / rsq;
    for (i = 0; i < latVol; i++) {
      pvec[i] = resNew[i] + beta * pvec[i];
      res[i] = resNew[i];
    }
  }
 
  if(k == maxIter) {
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq); 
    return 0; /* Failed convergence */
  }
 
  Aphi(Apvec,phi);
  for(i=0; i < latVol; i++) truersq += (Apvec[i] - b[i])*(Apvec[i] - b[i]);
  
  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return truersq; /* Convergence */
}

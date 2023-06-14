/*======================================================================================

 Rich Brower Sun Apr  2 19:28:40 EDT 2023 

 Develop as prototype code for Peter's Quintessence Grid Two Sphere implementation

        y (mu = 1)
        ^               
         \                  
           *    *    *    *    *    *   
             *    *    *    *    *    *  
               *    *    *    *    *   *  --> x  (mu = 0)
              /  
             /
           z (mu = 3)

x = 0,...,Lx   y = 0,..,Ly -1
mu =  x,y,z, -x, -y, -z = 0,1,2,3,4,5

======================================================================================*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <climits>
#include <iomanip>      // std::setprecision
#include <stack>

using namespace std;

#define TWOPI  6.283185307179586
#define Lx 32
#define Ly 32
#define Three 3
#define N  Lx*Ly   //  Lx * Ly
#define Debug 0

inline int mod(int x, int n)
{
  return  (n + (x % n))%n;
}

void printSpins(int s[Ly][Lx]);
int  getMag( int  s[Ly][Lx]);  // return magntization
void HeatBath(int s[Ly][Lx],double beta);
int  FlipCluster(int s[Ly][Lx], double beta);// sign flips spins
void testMap(int s[Ly][Lx]);

inline int neighbor(int n,int mu)
{
  int x,y;
  int xp, xm, yp, ym;
  int nn;
  
  x = n%Lx; y = n/Lx;
  xp = mod(x+1,Lx); xm =  mod(x-1,Lx); yp  = mod(y+1,Ly); ym =  mod(y-1,Ly);
  
  switch(mu){
  case 0:  nn =  xp + Lx*y;break; // positive
  case 1:  nn =  x + Lx*yp;  break;
  case 2:  nn =  xp + Lx*yp; break;
    
  case 3:  nn =  xm + Lx*y; break;  // negative
  case 4:  nn =  x + Lx*ym; break;
  case 5: nn =  xm + Lx*ym; break;
    
  default: nn = -1;
  }
  return nn;
}

int main()
{ 
  int s[Ly][Lx];
  double beta; //  Triagular critical point β_crit = ln(3)/4 =  0.274653
  beta = log(3.0)/4.0;
  ofstream out_Heat;
  out_Heat.open("Heat.dat");
    ofstream out_Wolff;
  out_Wolff.open("Wolff.dat");
  
   cout<<" Ising on Lx by Ly  triagular lattice: "<< Lx <<"  " << Ly <<" critical point " << log(3.0)/4.0 << endl;
 
  srand(137);
  // Hot start
    for(int y = 0; y< Ly; y++)
      for(int x = 0; x < Lx; x++)
	s[y][x] = 2 * ( rand() % 2) - 1;

 #if Debug
    cout << " Magnetization is " << getMag(s) << endl;
   printSpins(s);
   
    testMap(s);
    cout << "After TestMap Magnetization is " << getMag(s) << endl;
    printSpins(s);
   
#endif
   
    beta = beta + 0.02;

    double AveMag =0.0;
    int heatMag;

   for(int iter = 0;iter < 100000; iter++)
    {
       HeatBath(s,beta);
        heatMag =  getMag(s);
	AveMag += (double) abs(heatMag);
       
    if(iter%1000  == 0 )
	{
	  cout << iter;
	  cout <<" In Heat Bath   " << heatMag <<"   " << endl;
	  out_Heat << AveMag/1000.0 << endl;
	  AveMag = 0.0;
	}     
    }

   int WolffMag  =  getMag(s);

   double AveAbsMag = 0.0;
   
   for(int iter = 0;iter < 100000 ; iter++)
     {
       int  signed_flip =  FlipCluster(s,beta);
       WolffMag += - 2*signed_flip;
   
      AveMag += (double) WolffMag;
       AveAbsMag += abs( (double) WolffMag);
  
       if(iter%1000  == 0 )
       {
	 cout << iter;
	 cout <<" In Wolff signed_flip =   " << signed_flip << "Mag and  Update mag "<< getMag(s) << " " << WolffMag << endl;
	 out_Wolff << AveAbsMag/1000.0 <<"  " << AveMag/1000.0 <<endl;
	   AveMag = 0.0;
	   AveAbsMag = 0.0;   
       }
     }
   return 0;
}
  
/*******************************************************
Ising exp[beta sum_<i,j> s_i s_j]  and h_i = sum_{j nn of i} s_j
Defining  prob =  exp[2*beta*h_i]
prob to aligning with h_i
is  exp[beta*h]/(exp[beta*h] + (exp[-beta*h])   = pup/(1 + pup)
*******************************************************/


void  HeatBath(int s[Ly][Lx],double beta)
{
  int h = 0;
  int xp,xm,yp, ym;
  double xran, prob;
  
      for(int y = 0; y<Ly; y++)
	for(int x= 0 ; x<Lx; x++)
	  { xp = mod(x+1,Lx); xm =  mod(x-1,Lx); yp  = mod(y+1,Ly); ym =  mod(y-1,Ly);
	    h = s[y][xp] +  s[yp][x] +  s[yp][xp] +  s[y][xm] +  s[ym][x] +  s[ym][xm];
	    xran = (double)rand()/(double)RAND_MAX;
	    prob = exp(2.0*beta* h); 
	    s[y][x] = xran < prob/(1.0 + prob) ? +1 : -1;
	  } 
}


int  FlipCluster( int s[Ly][Lx], double beta){
  int *spt;
  spt = *s;  // external field
  bool is_cluster[N];
  int  wolff_cluster_size = 0;
  int mag_sign;
  double K[2 *Three];
  for(int mu =0;mu < 2 *Three; mu++) K[mu] = beta;

   // index to site n =  x * Lx*y so that value is spt[n] == s[y][x];
  for(int n= 0; n<N; n++) is_cluster[n] = false;  // clear is_cluster
  
  // create the stack
  std::stack<int> stack;
  // choose a random site and add it to the cluster
  int site =    rand() % (N); // careful with defines
   mag_sign = spt[site];
   is_cluster[site] = true;
   stack.push(site);
   //  cout << " pushed in stack site = " << site << endl;
   
  
    //     cout << "Seed Site = " << site << endl;

 
  while (stack.size() != 0) {
    int  n = stack.top();
    stack.pop();
      wolff_cluster_size += 1;

    // flip the spin
      spt[n] = - spt[n];

    // try to add neighbors
    for (int mu = 0; mu < 2*Three; mu++) {
      int nn = neighbor(n,mu);
      // skip if the site not allinged or  already pushed into clustered
       if ( spt[nn] == spt[n] || is_cluster[nn]  )    continue;
       double rate = -2.0 * K[mu]; // link coupling
       double xran = (double)rand()/(double)RAND_MAX;	    
      if (rate >= 0.0 || xran < exp(rate)) continue;
      is_cluster[nn] = true;
      stack.push(nn);
      //     cout << " pushed in stack nn = " << nn << endl;
    }
  }

   int  signed_flip =  mag_sign*wolff_cluster_size;
    
  return signed_flip; // signed cluster
}


int  getMag( int s[Ly][Lx])
{
 int  mag = 0;
  for(int y = 0; y<Ly; y++)
    for(int x= 0 ; x<Lx; x++)
      mag += s[y][x];
  return mag;
}

void printSpins(int s[Ly][Lx])
{
	cout<<"\n-------------------------------------------- \n";

	for(int y = 0; y<Ly; y++)
	  {
	    for(int x= 0 ; x<Lx; x++)
	      cout<< s[y][x] <<"   ";
	    cout << endl; 
	  }					 
} 
  
 void testMap(int s[Ly][Lx])
{
  int *spt;
  spt = *s;
  cout<<"\n Insidet testMap \n";
	for(int y = 0; y<Ly; y++)
	  {
	    for(int x= 0 ; x<Lx; x++)
	      cout<< spt[ x + y*Lx] <<"   ";
	    cout << endl; 
	  }

	cout<<"\n Remap insdie estMap \n";
      	for(int y = 0; y<Ly; y++)
	  {
	    for(int x= 0 ; x<Lx; x++){
	      spt[x + y*Lx] = x + y*Lx;
	       cout << spt[x + y*Lx] << "   ";
	    }
	     cout << endl;
	  }
		
}
 
 

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

#include <filesystem> // c++17

using namespace std;

#define TWOPI  6.283185307179586
#define Lx 32
#define Ly 32
#define Three 3
#define N  Lx*Ly   //  Lx * Ly
#define Debug 0
#define NInit 100
#define NIter 10000
#define NInter 10

const string dir_obs = "obs/";
const string f_mag = dir_obs+"mag";
const string f_susc = dir_obs+"susc";
const string f_ss = dir_obs+"ss";
const string f_exex = dir_obs+"exex";
const string f_exey = dir_obs+"exey";
const string f_exez = dir_obs+"exez";
const string f_eyey = dir_obs+"eyey";
const string f_eyez = dir_obs+"eyez";
const string f_ezez = dir_obs+"ezez";

inline int mod(int x, int n)
{
  return  (n + (x % n))%n;
}

void printSpins(int s[Ly][Lx]);
int  getMag( int  s[Ly][Lx]);  // return magntization
void HeatBath(int s[Ly][Lx], double beta, double K[2*Three]);
int  FlipCluster(int s[Ly][Lx], double beta, double K[2*Three]);// sign flips spins

void getCorrSS( double corr[Ly][Lx], const int s[Ly][Lx] );
void getCorrEEs( double exex[Ly][Lx],
                 double exey[Ly][Lx],
                 double exez[Ly][Lx],
                 double eyey[Ly][Lx],
                 double eyez[Ly][Lx],
                 double ezez[Ly][Lx],
                 const int s[Ly][Lx] );

void writeCorr(ofstream& fout, const double corr[Ly][Lx]);

void update_mean(double& mean, const double obs, const int counter);
void update_mean(double mean[Ly][Lx], const double obs[Ly][Lx], const int counter);

void write2file(const string& filename_header, const double mean, const string& memo);
void write2file(const string& filename_header, const double mean[Ly][Lx], const string& memo);


inline int neighbor(int n,int mu)
{
  int x,y;
  int xp, xm, yp, ym;
  int nn;

  x = n%Lx; y = n/Lx;
  xp = mod(x+1,Lx); xm = mod(x-1,Lx); yp = mod(y+1,Ly); ym = mod(y-1,Ly);

  switch(mu){
  case 0:  nn =  xp + Lx*y;  break; // positive
  case 1:  nn =  x + Lx*yp;  break;
  case 2:  nn =  xp + Lx*yp; break;

  case 3:  nn =  xm + Lx*y;  break;  // negative
  case 4:  nn =  x + Lx*ym;  break;
  case 5:  nn =  xm + Lx*ym; break;

  default: nn = -1;
  }
  return nn;
}



int main()
{
  filesystem::create_directory(dir_obs);

  double mag_mean = 0.0;
  double susc_mean = 0.0;
  double ss_mean[Ly][Lx] = {{0.0}};
  double exex_mean[Ly][Lx] = {{0.0}};
  double exey_mean[Ly][Lx] = {{0.0}};
  double exez_mean[Ly][Lx] = {{0.0}};
  double eyey_mean[Ly][Lx] = {{0.0}};
  double eyez_mean[Ly][Lx] = {{0.0}};
  double ezez_mean[Ly][Lx] = {{0.0}};

  int s[Ly][Lx];
  double beta = log(3.0)/4.0; //  Triagular critical beta = ln(3)/4 = 0.274653
  const double Kx = 1.0, Ky=1.0, Kz=1.0;

  double K[2*Three];
  K[0] = Kx; K[1] = Ky; K[2] = Kz;
  K[3] = Kx; K[4] = Ky; K[5] = Kz;

  cout << " Ising on Lx by Ly triangular lattice: " << Lx << "  " << Ly
       << " critical point " << log(3.0)/4.0
       << " Kx " << Kx << " Ky " << Ky << " Kz " << Kz
       << " beta " << beta
       << endl;

  srand(137);
  // Hot start
  for(int y=0; y<Ly; y++)
    for(int x=0; x<Lx; x++)
      s[y][x] = 2 * ( rand() % 2) - 1;

  int mag = getMag(s);
  for(int iter=0; iter<NInit; iter++){
    const int signed_flip = FlipCluster(s,beta,K);
    mag += -2*signed_flip;
  }

  int counter=0;
  for(int iter=0; iter<NIter; iter++){
    counter++;
    const int signed_flip = FlipCluster(s,beta,K);
    mag += -2*signed_flip;
    // HeatBath(s,beta,K);
    // mag = getMag(s);

    // -----
    // obs part
    if(iter%NInter==0){

      // mag
      update_mean(mag_mean, mag, counter);
      write2file(f_mag, mag_mean, "count="+std::to_string(counter));

      // susc
      update_mean(susc_mean, 1.0*mag*mag/(Ly*Lx), counter);
      write2file(f_susc, susc_mean, "count="+std::to_string(counter));

      // corr_ss
      double corr[Ly][Lx];
      getCorrSS(corr, s);
      update_mean(ss_mean, corr, counter);
      write2file(f_ss, ss_mean, "count="+std::to_string(counter));

      // corr_ee
      double exex[Ly][Lx];
      double exey[Ly][Lx];
      double exez[Ly][Lx];
      double eyey[Ly][Lx];
      double eyez[Ly][Lx];
      double ezez[Ly][Lx];
      getCorrEEs(exex,exey,exez,eyey,eyez,ezez, s);
      update_mean(exex_mean, exex, counter);
      write2file(f_exex, exex_mean, "count="+std::to_string(counter));
      update_mean(exey_mean, exey, counter);
      write2file(f_exey, exey_mean, "count="+std::to_string(counter));
      update_mean(exez_mean, exez, counter);
      write2file(f_exez, exez_mean, "count="+std::to_string(counter));
      update_mean(eyey_mean, eyey, counter);
      write2file(f_eyey, eyey_mean, "count="+std::to_string(counter));
      update_mean(eyez_mean, eyez, counter);
      write2file(f_eyez, eyez_mean, "count="+std::to_string(counter));
      update_mean(ezez_mean, ezez, counter);
      write2file(f_ezez, ezez_mean, "count="+std::to_string(counter));
    }
    // -----

  }
  return 0;
}





/*******************************************************
Ising exp[beta sum_<i,j> s_i s_j]  and h_i = sum_{j nn of i} s_j
Defining  prob =  exp[2*beta*h_i]
prob to aligning with h_i
is  exp[beta*h]/(exp[beta*h] + (exp[-beta*h])   = pup/(1 + pup)
*******************************************************/


void  HeatBath(int s[Ly][Lx], double beta, double K[2*Three])
{
  int h = 0;
  int xp,xm,yp, ym;
  double xran, prob;

  for(int y=0; y<Ly; y++)
    for(int x=0; x<Lx; x++)
      {
        xp = mod(x+1,Lx); xm = mod(x-1,Lx); yp = mod(y+1,Ly); ym = mod(y-1,Ly);
        h = K[0]*s[y][xp] + K[1]*s[yp][x] + K[2]*s[yp][xp] + K[0]*s[y][xm] + K[1]*s[ym][x] + K[2]*s[ym][xm]; // K
        xran = (double)rand()/(double)RAND_MAX;
        prob = exp(2.0*beta*h);
        s[y][x] = xran < prob/(1.0 + prob) ? +1 : -1;
      }
}


int  FlipCluster( int s[Ly][Lx], double beta, double K[2*Three]){
  int *spt; spt = *s;
  bool is_cluster[N];
  int wolff_cluster_size = 0;
  int mag_sign;

  // index to site n =  x * Lx*y so that value is spt[n] == s[y][x];
  for(int n=0; n<N; n++) is_cluster[n] = false;  // clear is_cluster

  // create the stack
  std::stack<int> stack;
  // choose a random site and add it to the cluster
  int site = rand() % (N); // careful with defines
  mag_sign = spt[site];
  is_cluster[site] = true;
  stack.push(site);

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
      if ( spt[nn] == spt[n] || is_cluster[nn] ) continue;
      double rate = -2.0 * beta * K[mu]; // link coupling
      double xran = (double)rand()/(double)RAND_MAX;
      if (rate >= 0.0 || xran < exp(rate)) continue;
      is_cluster[nn] = true;
      stack.push(nn);
    }
  }
  int  signed_flip = mag_sign*wolff_cluster_size;
  return signed_flip; // signed cluster
}


int  getMag( int s[Ly][Lx])
{
  int  mag=0;
  for(int y=0; y<Ly; y++)
    for(int x=0 ; x<Lx; x++)
      mag += s[y][x];
  return mag;
}

void printSpins(int s[Ly][Lx])
{
  cout << "\n-------------------------------------------- \n";

  for(int y=0; y<Ly; y++)
    {
      for(int x=0; x<Lx; x++)
        cout<< s[y][x] <<"   ";
      cout << endl;
    }
}


void getCorrSS( double corr[Ly][Lx], const int s[Ly][Lx] ){
  for(int y1=0; y1<Ly; y1++){
    for(int x1=0; x1<Lx; x1++){
      corr[y1][x1] = 0.0;

      for(int y2=0; y2<Ly; y2++){
        for(int x2=0; x2<Lx; x2++){
          corr[y1][x1] += s[y2][x2]*s[(y2+y1)%Ly][(x2+x1)%Lx];
        }}

      corr[y1][x1] /= Ly*Lx;
    }}
}


void getCorrEEs( double exex[Ly][Lx],
                 double exey[Ly][Lx],
                 double exez[Ly][Lx],
                 double eyey[Ly][Lx],
                 double eyez[Ly][Lx],
                 double ezez[Ly][Lx],
                 const int s[Ly][Lx] ){
  for(int y1=0; y1<Ly; y1++){
    for(int x1=0; x1<Lx; x1++){
      exex[y1][x1] = 0.0;
      exey[y1][x1] = 0.0;
      exez[y1][x1] = 0.0;
      eyey[y1][x1] = 0.0;
      eyez[y1][x1] = 0.0;
      ezez[y1][x1] = 0.0;

      for(int y2=0; y2<Ly; y2++){
        for(int x2=0; x2<Lx; x2++){
          const int y2p1 = (y2+1)%Ly;
          const int x2p1 = (x2+1)%Lx;
          const int ex = s[y2][x2]*s[y2][x2p1];
          const int ey = s[y2][x2]*s[y2p1][x2];
          const int ez = s[y2][x2]*s[y2p1][x2p1];
          //
          const int y2py1 = (y2+y1)%Ly;
          const int x2px1 = (x2+x1)%Lx;
          const int y2py1p1 = (y2+y1+1)%Ly;
          const int x2px1p1 = (x2+x1+1)%Lx;
          const int ex1 = s[y2py1][x2px1]*s[y2py1][x2px1p1];
          const int ey1 = s[y2py1][x2px1]*s[y2py1p1][x2px1];
          const int ez1 = s[y2py1][x2px1]*s[y2py1p1][x2px1p1];

          exex[y1][x1] += ex*ex1;
          exey[y1][x1] += ex*ey1;
          exez[y1][x1] += ex*ez1;
          eyey[y1][x1] += ey*ey1;
          eyez[y1][x1] += ey*ez1;
          ezez[y1][x1] += ez*ez1;
        }}

      exex[y1][x1] /= Ly*Lx;
      exey[y1][x1] /= Ly*Lx;
      exez[y1][x1] /= Ly*Lx;
      eyey[y1][x1] /= Ly*Lx;
      eyez[y1][x1] /= Ly*Lx;
      ezez[y1][x1] /= Ly*Lx;
    }}
}

void writeCorr(ofstream& fout, const double corr[Ly][Lx]){
  fout << std::right << std::scientific << std::setprecision(15);
  for(int y=0; y<Ly; y++){
    for(int x=0; x<Lx; x++){
      fout << std::setw(25) << corr[y][x] << " ";
    }
    fout << endl;
  }
}


void update_mean(double& mean, const double obs, const int counter){
  mean = ( (counter-1.0)*mean + obs)/counter;
}

void update_mean(double mean[Ly][Lx], const double obs[Ly][Lx], const int counter){
  for(int y=0; y<Ly; y++) for(int x=0; x<Lx; x++){
      update_mean(mean[y][x], obs[y][x], counter);
    }
}

void write2file(const string& filename_header, const double v, const string& memo){
  const string f_tmp = filename_header+"_tmp.dat";
  const string f_fin = filename_header+".dat";
  ofstream fout(f_tmp, ios::trunc);
  fout << std::right << std::scientific << std::setprecision(15);
  fout << "#"+memo << std::endl;
  fout << std::setw(25) << v;
  fout.close();
  filesystem::copy(f_tmp, f_fin, filesystem::copy_options::update_existing);
  filesystem::remove(f_tmp);
}

void write2file(const string& filename_header, const double v[Ly][Lx], const string& memo){
  const string f_tmp = filename_header+"_tmp.dat";
  const string f_fin = filename_header+".dat";
  ofstream fout(f_tmp, ios::trunc);
  fout << "#"+memo << std::endl;
  writeCorr(fout, v);
  fout.close();
  filesystem::copy(f_tmp, f_fin, filesystem::copy_options::update_existing);
  filesystem::remove(f_tmp);
}


// void testMap(int s[Ly][Lx])
// {
//   int *spt;
//   spt = *s;
//   cout<<"\n Insidet testMap \n";
//   for(int y=0; y<Ly; y++)
//     {
//       for(int x=0 ; x<Lx; x++)
//         cout << spt[x+y*Lx] <<"   ";
//       cout << endl;
//     }

//   cout<<"\n Remap insdie estMap \n";
//   for(int y=0; y<Ly; y++)
//     {
//       for(int x=0 ; x<Lx; x++){
//         spt[x+y*Lx] = x+y*Lx;
//         cout << spt[x+y*Lx] << "   ";
//       }
//       cout << endl;
//     }

// }


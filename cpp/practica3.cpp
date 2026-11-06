#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif
using namespace std;
using ll = int64_t;
using Flat = vector<ll>;
using Mat = vector<vector<ll>>;
uint32_t rng=734251;
uint32_t azar(){rng^=rng<<13;rng^=rng>>17;rng^=rng<<5;return rng;}
double memoria(){
#ifdef _WIN32
 PROCESS_MEMORY_COUNTERS pm{};
 return GetProcessMemoryInfo(GetCurrentProcess(),&pm,sizeof(pm))?double(pm.PeakWorkingSetSize)/1048576.0:-1;
#else
 struct rusage r{};getrusage(RUSAGE_SELF,&r);
#ifdef __APPLE__
 return r.ru_maxrss/1048576.0;
#else
 return r.ru_maxrss/1024.0;
#endif
#endif
}
template<class T> vector<T> datos(size_t n){
 vector<T>a(n);for(size_t i=0;i<n;++i)a[i]=T(1+(i*7919ULL+1234567ULL)%10000000ULL);
 rng=734251;for(size_t i=n;i>1;--i)swap(a[i-1],a[azar()%i]);return a;
}
template<class T>T quickselect(vector<T>&a,size_t k){
 size_t lo=0,hi=a.size();--k;
 while(true){if(hi-lo==1)return a[lo];
  T p=a[lo+(hi-lo)/2];size_t lt=lo,i=lo,gt=hi;
  while(i<gt){if(a[i]<p)swap(a[lt++],a[i++]);else if(a[i]>p)swap(a[i],a[--gt]);else ++i;}
  if(k<lt)hi=lt;else if(k>=gt)lo=gt;else return a[k];
 }
}
Flat suma(const Flat&a,const Flat&b,int s){Flat c(s*s);for(int i=0;i<s*s;i++)c[i]=a[i]+b[i];return c;}
Flat resta(const Flat&a,const Flat&b,int s){Flat c(s*s);for(int i=0;i<s*s;i++)c[i]=a[i]-b[i];return c;}
Flat cortar(const Flat&a,int n,int fi,int fj){int h=n/2;Flat c(h*h);for(int i=0;i<h;i++)copy_n(a.begin()+(fi*h+i)*n+fj*h,h,c.begin()+i*h);return c;}
Flat tradicional_flat(const Flat&a,const Flat&b,int n){Flat c(n*n,0);for(int i=0;i<n;i++)for(int k=0;k<n;k++){ll v=a[i*n+k];for(int j=0;j<n;j++)c[i*n+j]+=v*b[k*n+j];}return c;}
Flat strassen_flat(const Flat&a,const Flat&b,int n){
 if(n<=64)return tradicional_flat(a,b,n);
 int h=n/2;
 Flat a11=cortar(a,n,0,0),a12=cortar(a,n,0,1),a21=cortar(a,n,1,0),a22=cortar(a,n,1,1);
 Flat b11=cortar(b,n,0,0),b12=cortar(b,n,0,1),b21=cortar(b,n,1,0),b22=cortar(b,n,1,1);
 Flat m1=strassen_flat(suma(a11,a22,h),suma(b11,b22,h),h);
 Flat m2=strassen_flat(suma(a21,a22,h),b11,h);
 Flat m3=strassen_flat(a11,resta(b12,b22,h),h);
 Flat m4=strassen_flat(a22,resta(b21,b11,h),h);
 Flat m5=strassen_flat(suma(a11,a12,h),b22,h);
 Flat m6=strassen_flat(resta(a21,a11,h),suma(b11,b12,h),h);
 Flat m7=strassen_flat(resta(a12,a22,h),suma(b21,b22,h),h);
 Flat c11=suma(resta(suma(m1,m4,h),m5,h),m7,h),c12=suma(m3,m5,h),c21=suma(m2,m4,h),c22=suma(resta(suma(m1,m3,h),m2,h),m6,h);
 Flat c(n*n);for(int i=0;i<h;i++)for(int j=0;j<h;j++){c[i*n+j]=c11[i*h+j];c[i*n+j+h]=c12[i*h+j];c[(i+h)*n+j]=c21[i*h+j];c[(i+h)*n+j+h]=c22[i*h+j];}return c;
}
Mat suma(const Mat&a,const Mat&b){int n=a.size();Mat c(n,vector<ll>(n));for(int i=0;i<n;i++)for(int j=0;j<n;j++)c[i][j]=a[i][j]+b[i][j];return c;}
Mat resta(const Mat&a,const Mat&b){int n=a.size();Mat c(n,vector<ll>(n));for(int i=0;i<n;i++)for(int j=0;j<n;j++)c[i][j]=a[i][j]-b[i][j];return c;}
Mat cortar(const Mat&a,int fi,int fj){int h=a.size()/2;Mat c(h,vector<ll>(h));for(int i=0;i<h;i++)copy_n(a[i+fi*h].begin()+fj*h,h,c[i].begin());return c;}
Mat tradicional_2d(const Mat&a,const Mat&b){int n=a.size();Mat c(n,vector<ll>(n,0));for(int i=0;i<n;i++)for(int k=0;k<n;k++){ll v=a[i][k];for(int j=0;j<n;j++)c[i][j]+=v*b[k][j];}return c;}
Mat strassen_2d(const Mat&a,const Mat&b){int n=a.size();if(n<=64)return tradicional_2d(a,b);int h=n/2;
 Mat a11=cortar(a,0,0),a12=cortar(a,0,1),a21=cortar(a,1,0),a22=cortar(a,1,1);
 Mat b11=cortar(b,0,0),b12=cortar(b,0,1),b21=cortar(b,1,0),b22=cortar(b,1,1);
 Mat m1=strassen_2d(suma(a11,a22),suma(b11,b22)),m2=strassen_2d(suma(a21,a22),b11);
 Mat m3=strassen_2d(a11,resta(b12,b22)),m4=strassen_2d(a22,resta(b21,b11));
 Mat m5=strassen_2d(suma(a11,a12),b22),m6=strassen_2d(resta(a21,a11),suma(b11,b12));
 Mat m7=strassen_2d(resta(a12,a22),suma(b21,b22));
 Mat c11=suma(resta(suma(m1,m4),m5),m7),c12=suma(m3,m5),c21=suma(m2,m4),c22=suma(resta(suma(m1,m3),m2),m6);
 Mat c(n,vector<ll>(n));for(int i=0;i<h;i++)for(int j=0;j<h;j++){c[i][j]=c11[i][j];c[i][j+h]=c12[i][j];c[i+h][j]=c21[i][j];c[i+h][j+h]=c22[i][j];}return c;
}
int main(int argc,char**argv){if(argc<3){cerr<<"Uso: practica3 select bits n [k] | matrix tradicional|strassen flat|2d n\n";return 2;}
 try{string modo=argv[1];auto start=chrono::steady_clock::now();double secs=0;string ok="si",extra="";
 if(modo=="select"){
  if(argc<4)throw invalid_argument("Faltan bits y n");int bits=stoi(argv[2]);size_t n=stoull(argv[3]);if(!n||n>10000000)throw invalid_argument("n fuera del dominio");size_t k=argc>4?stoull(argv[4]):(n+1)/2;if(k<1||k>n)throw invalid_argument("k fuera de rango");
  if(bits==32){auto a=datos<int32_t>(n);start=chrono::steady_clock::now();auto v=quickselect(a,k);secs=chrono::duration<double>(chrono::steady_clock::now()-start).count();size_t less=0,eq=0;for(auto x:a){less+=x<v;eq+=x==v;}ok=(less<k&&k<=less+eq)?"si":"no";extra=to_string(v);}
  else if(bits==64){auto a=datos<int64_t>(n);start=chrono::steady_clock::now();auto v=quickselect(a,k);secs=chrono::duration<double>(chrono::steady_clock::now()-start).count();size_t less=0,eq=0;for(auto x:a){less+=x<v;eq+=x==v;}ok=(less<k&&k<=less+eq)?"si":"no";extra=to_string(v);}
  else throw invalid_argument("Solo 32 o 64 bits");
  cout<<fixed<<setprecision(6)<<"cpp,quickselect,"<<bits<<",- ,"<<n<<","<<k<<","<<secs<<","<<memoria()<<","<<ok<<","<<extra<<"\n";
 }else if(modo=="matrix"){
  if(argc<5)throw invalid_argument("Faltan algoritmo, representacion y n");string alg=argv[2],rep=argv[3];int n=stoi(argv[4]);if(n<16||n>512||(n&(n-1)))throw invalid_argument("n debe ser potencia de dos de 16 a 512");
  if(rep=="flat"){
   Flat a(n*n),b(n*n);for(int i=0;i<n*n;i++){a[i]=(i*17+3)%11-5;b[i]=(i*37+9)%11-5;}
   start=chrono::steady_clock::now();auto c=alg=="strassen"?strassen_flat(a,b,n):tradicional_flat(a,b,n);secs=chrono::duration<double>(chrono::steady_clock::now()-start).count();
   if(alg!="strassen"&&alg!="tradicional")throw invalid_argument("Algoritmo desconocido");
   ll suma=0;for(ll x:c)suma+=x;extra=to_string(suma);
   for(int i: {0,n/2,n-1})for(int j:{0,n/2,n-1}){ll v=0;for(int k=0;k<n;k++)v+=a[i*n+k]*b[k*n+j];if(v!=c[i*n+j])ok="no";}
  }else if(rep=="2d"){
   Mat a(n,vector<ll>(n)),b(n,vector<ll>(n));for(int i=0;i<n;i++)for(int j=0;j<n;j++){int z=i*n+j;a[i][j]=(z*17+3)%11-5;b[i][j]=(z*37+9)%11-5;}
   start=chrono::steady_clock::now();auto c=alg=="strassen"?strassen_2d(a,b):tradicional_2d(a,b);secs=chrono::duration<double>(chrono::steady_clock::now()-start).count();
   if(alg!="strassen"&&alg!="tradicional")throw invalid_argument("Algoritmo desconocido");
   ll suma=0;for(auto&r:c)for(ll x:r)suma+=x;extra=to_string(suma);
   for(int i: {0,n/2,n-1})for(int j:{0,n/2,n-1}){ll v=0;for(int k=0;k<n;k++)v+=a[i][k]*b[k][j];if(v!=c[i][j])ok="no";}
  }else throw invalid_argument("Representacion desconocida");
  cout<<fixed<<setprecision(6)<<"cpp,"<<alg<<",64,"<<rep<<","<<n<<",-,"<<secs<<","<<memoria()<<","<<ok<<","<<extra<<"\n";
 }else throw invalid_argument("Modo desconocido");return ok=="si"?0:4;
 }catch(const exception&e){cerr<<e.what()<<"\n";return 1;}
}

// http://cap-lore.com/MathPhys/Simplex/dilog/
#include <stdio.h>
#include <complex.h>
#include <math.h>
typedef double _Complex C;
typedef double R;

static R  sq(R x){return x*x;}
static C csq(C x){return x*x;}
static R m2(C x){return sq(creal(x)) + sq(cimag(x));}
static int t(C x){return fabs(creal(x)) < 1.e-16 && fabs(cimag(x)) < 1.e-16;}
static R const pi = 3.141592653589793238, p2 = 3.141592653589793238*3.141592653589793238;

C dilog(C z){R r = m2(z);
  if(r > 1) {C r = 1/z; return - dilog(r) - p2/6 - 0.5*csq(clog(-r));}
  if(creal(z) > 0.6) return - dilog(1-z) + p2/6 - (z==1?0:clog(1-z)*clog(z));
  if(r > 0.4) if(creal(z) < 0.15) return -dilog(z/(z-1))-0.5*csq(clog(1-z));
              else return 0.5*dilog(z*z) - dilog(-z);
  if(m2(z + 0.268) < 0.0049){C s = 0; C w = z; R a = 1, b = 3;
     while(1) {C T = w/a; s += T; if(t(T)) break; w *= z; a += b; b += 2;}
     return s;}
 {C z2=z*z, s = z + 1.4375*z2 + 0.75*(1-z2)*clog(1-z), f = z*z2; R J=6, K=18, L=18;
    while(1){C a = f/(J*J); s+=a; if(t(a)) break; f *= z; J+=K; K+=L; L+=6;}
    return 4*s/(1+4*z+z2);}}

C dilogD(C z){long double _Complex s = 0; C w = z; R a = 1, b = 3;
  if(m2(z)>1) {C r = 1/z; return - dilogD(r) - p2/6 - 0.5*csq(clog(-r));}
  while(1) {C T = w/a; s += T; if(t(T)) break; w *= z; a += b; b += 2;}
  return s;}

int main(void)
{
    printf("%g\n", creal(dilog(0.5+0*I)));
    return 0;
}

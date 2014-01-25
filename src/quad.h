#ifndef __QUAD_H
#define __QUAD_H

#if defined(__ICC) || defined(__INTEL_COMPILER)
    #define __float128 _Quad

    #define isinfq(x) (x/10 == x)
    #define isnanq(x) (x != x)

    _Quad __logq(_Quad);
    _Quad __expq(_Quad);
    _Quad __sqrtq(_Quad);
    _Quad __log1pq(_Quad);

    #define logq __logq
    #define expq __logq
#elif defined(__GNUC__) || defined(__GNUG__)
    #include <quadmath.h>
#endif

#endif

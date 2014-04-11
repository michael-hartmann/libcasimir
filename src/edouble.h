#ifndef __EDOUBLE_H
#define __EDOUBLE_H

#ifdef EXTENDED_DOUBLE
    #define CASIMIR_ARITHMETICS "long double"
    #define edouble long double

    #define logq      logl
    #define expq      expl
    #define sqrtq     sqrtl
    #define log1pq    log1pl
    #define fabsq     fabsl
    #define sinq      sinl
    #define cosq      cosl
    #define copysignq copysignl
    #define isinfq(x) (x/10 == x)
    #define isnanq(x) (x != x)
#else
    #if defined(__ICC) || defined(__INTEL_COMPILER)
        #define CASIMIR_ARITHMETICS "icc _Quad"
        #define edouble _Quad
    
        /* define isinf and isnan */
        #define isinfq(x) (x/10 == x)
        #define isnanq(x) (x != x)
    
        /* define prototypes. without these prototypes icc will return nan. */
        _Quad __logq(_Quad);
        #define logq __logq
    
        _Quad __cosq(_Quad);
        #define cosq __cosq
    
        _Quad __sinq(_Quad);
        #define sinq __sinq
    
        _Quad __expq(_Quad);
        #define expq __expq
    
        _Quad __sqrtq(_Quad);
        #define sqrtq __sqrtq
    
        _Quad __log1pq(_Quad);
        #define log1pq __log1pq
    
        _Quad __fabsq(_Quad);
        #define fabsq __fabsq
    
        _Quad __copysignq(_Quad, _Quad);
        #define copysignq __copysignq
    
    #elif defined(__GNUC__) || defined(__GNUG__)
        #define CASIMIR_ARITHMETICS "gcc __float128"

        #include <quadmath.h>
        #define edouble __float128
    #else
        #error "I'm sorry, but quad precision is only supported with gcc or icc at the moment."
    #endif
#endif

#if defined(__ICC) || defined(__INTEL_COMPILER)
    #define COMPILER "icc"
#elif defined(__GNUC__) || defined(__GNUG__)
    #define COMPILER "gcc"
#else
    #define COMPILER "unknown"
#endif

#ifdef MATRIX_QUAD
    #define CASIMIR_MATRIX "plain"
#else
    #define CASIMIR_MATRIX "log"
#endif

#ifdef INTEGRATION_QUAD
    #define CASIMIR_INTEGRATION "plain"
#else
    #define CASIMIR_INTEGRATION "log"
#endif

#endif

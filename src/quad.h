#ifndef __QUAD_H
#define __QUAD_H

#if defined(__ICC) || defined(__INTEL_COMPILER)
    #define quad_t _Quad

    _Quad __logq(_Quad);
    _Quad __expq(_Quad);

    #define logq __logq
    #define expq __logq
#elif defined(__GNUC__) || defined(__GNUG__)
    #include <quadmath.h>
#endif

#endif

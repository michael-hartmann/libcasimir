#include <stdio.h>

#include "libcasimir.h"
#include "integration.h"

#include "unittest.h"

int test_Lambda(void);
int test_integration(void);


int test_Lambda(void)
{
    unittest_t test;
    unittest_init(&test, "Lambda", "Test Lambda function for various parameters");

    AssertAlmostEqual(&test, casimir_lnLambda(50,50,0),  -3.921875301871158);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,1),  -11.76572394002363);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,50), -367.6612508574347);
    AssertAlmostEqual(&test, casimir_lnLambda(50,20,10), -72.40153130583653);
    AssertAlmostEqual(&test, casimir_lnLambda(5,7,3),    -12.77235621226475);
    AssertAlmostEqual(&test, casimir_lnLambda(16,6,4),   -20.7139882421892);
    AssertAlmostEqual(&test, casimir_lnLambda(100,6,4),  -28.88322376001037);

    return test_results(&test, stderr);
}

int test_integration(void)
{
    casimir_integrals_t cint;
    unittest_t test;
    unittest_init(&test, "Integration", "Test integration for various parameters");

    casimir_integrate(&cint, 4, 4, 1, 0.01, 1);
    AssertAlmostEqual(&test, cint.A, +2.4806179125126554e17);
    AssertAlmostEqual(&test, cint.B, -2.2226323455151368e24);
    AssertAlmostEqual(&test, cint.C, -6.9457269656680333e20);
    AssertAlmostEqual(&test, cint.D, +6.9457269656680333e20);

    casimir_integrate(&cint, 4, 4, 1, 0.01, 10);
    AssertAlmostEqual(&test, cint.A, +2.4806179125126554e18);
    AssertAlmostEqual(&test, cint.B, -2.2226323455151368e25);
    AssertAlmostEqual(&test, cint.C, -6.9457269656680333e21);
    AssertAlmostEqual(&test, cint.D, +6.9457269656680333e21);

    casimir_integrate(&cint, 40, 40, 1, 0.5, 1);
    AssertAlmostEqual(&test, cint.A, +1.5754477603435539e159);
    AssertAlmostEqual(&test, cint.B, -6.3723632215476122e166);
    AssertAlmostEqual(&test, cint.C, -9.9568222699306801e162);
    AssertAlmostEqual(&test, cint.D, +9.9568222699306801e162);

    casimir_integrate(&cint, 40, 40, 40, 2, 1);
    AssertAlmostEqual(&test, cint.A, +6.4140686579381969e91);
    AssertAlmostEqual(&test, cint.B, -1.0147301906459434e95);
    AssertAlmostEqual(&test, cint.C, -2.5352219594503741e93);
    AssertAlmostEqual(&test, cint.D, +2.5352219594503736e93);

    casimir_integrate(&cint, 7, 4, 3, 17, 1);
    AssertAlmostEqual(&test, cint.A, +4.8180365200137397e-9);
    AssertAlmostEqual(&test, cint.B, -1.3731640166794149e-8);
    AssertAlmostEqual(&test, cint.C, -6.7659079909128738e-9);
    AssertAlmostEqual(&test, cint.D, +9.44463292099617e-9);

    casimir_integrate(&cint, 40, 40, 0, 5, 1);
    AssertAlmostEqual(&test, cint.B, -6.0455421304871757e85);

    casimir_integrate(&cint, 100, 41, 0, 5, 1);
    AssertAlmostEqual(&test, cint.B, 8.8689390374540308e185);

    return test_results(&test, stderr);
}

int main(int argc, char *argv[])
{
    test_Lambda();
    test_integration();
    
    return 0;
}

#include <math.h>
#include <stdio.h>

#include "libcasimir.h"
#include "integration.h"
#include "sfunc.h"
#include "givens.h"
#include "unittest.h"

#include "tests.h"

int test_logdet()
{
    unittest_t test;
    casimir_t casimir;
    casimir_mie_cache_t cache;
    const double RbyScriptL = 0.97;
    const double T = 0.1;
    const int lmax = 200;
    double logdet;

    unittest_init(&test, "logdet M", "calculate logdet");

    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_lmax(&casimir, lmax);

    logdet = casimir_logdetD(&casimir, 0, 0, &cache);
    AssertAlmostEqual(&test, logdet, -3.45236396285874);

    logdet = casimir_logdetD(&casimir, 0, 1, &cache);
    AssertAlmostEqual(&test, logdet, -2.63586999367158);

    logdet = casimir_logdetD(&casimir, 0, 10, &cache);
    AssertAlmostEqual(&test, logdet, -0.0276563864490425);

    casimir_mie_cache_init(&cache, 1*T*RbyScriptL);
    casimir_mie_cache_alloc(&casimir, &cache, lmax);

    logdet = casimir_logdetD(&casimir, 1, 1, &cache);
    AssertAlmostEqual(&test, logdet, -2.63900987016801);

    casimir_mie_cache_free(&cache);

    casimir_free(&casimir);

    return test_results(&test, stderr);
}


int test_logadd()
{
    int sign;
    double ret;
    unittest_t test;
    unittest_init(&test, "logadd", "add log");

    AssertAlmostEqual(&test, log(10+20), logadd(log(10),log(20)));
    AssertAlmostEqual(&test, log(900), logadd(log(600),log(300)));
    AssertAlmostEqual(&test, log(10), logadd(log(10),log(0)));

    // + +
    ret = logadd_s(log(2), 1, log(4), 1, &sign);
    AssertAlmostEqual(&test, log(6), ret);
    AssertEqual(&test, sign, +1);

    // - -
    ret = logadd_s(log(2), -1, log(4), -1, &sign);
    AssertAlmostEqual(&test, log(6), ret);
    AssertEqual(&test, sign, -1);

    // + -
    ret = logadd_s(log(2), +1, log(4), -1, &sign);
    AssertAlmostEqual(&test, log(2), ret);
    AssertEqual(&test, sign, -1);

    // - +
    ret = logadd_s(log(2), -1, log(4), +1, &sign);
    AssertAlmostEqual(&test, log(2), ret);
    AssertEqual(&test, sign, +1);

    // - +
    ret = logadd_s(log(2), -1, -INFINITY, +1, &sign);
    AssertAlmostEqual(&test, log(2), ret);

    // a == b
    ret = logadd_s(log(4), -1, log(4), +1, &sign);
    Assert(&test, isinf(ret));
    Assert(&test, ret < 0);

    // large values
    ret = logadd_s(1000, +1, 1004, +1, &sign);
    AssertAlmostEqual(&test, ret, 1004.018149927917809740354983318);
    AssertEqual(&test, sign, +1);

    return test_results(&test, stderr);
}

int test_givens()
{
    matrix_t *M;
    unittest_t test;
    unittest_init(&test, "QR decomposition", "Test QR decomposition using givens rotation");

    {
        M = matrix_alloc(2);
        matrix_set(M, 0,0, log(20));
        matrix_set(M, 0,1, log(2)+3000*log(10));
        matrix_set(M, 1,0, log(1)-3001*log(10));
        matrix_set(M, 1,1, log(1));

        matrix_log_balance(M);

        matrix_exp(M);

        AssertAlmostEqual(&test, matrix_logdet(M), log(19.8));

        matrix_free(M);
    }

    return test_results(&test, stderr);
}

int test_besselI(void)
{
    unittest_t test;
    unittest_init(&test, "Bessel function I_nu", "Test modified Bessel function I_nu");

    AssertAlmostEqual(&test, bessel_lnInu(0,1e-6), -7.133546631626697);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e-5), -5.982254085113174);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e-4), -4.830961536966152);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e-3), -3.679668825469134);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e-2), -2.528359779027661);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e-1), -1.375417787678169);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e0),  -0.064351991073531);
    AssertAlmostEqual(&test, bessel_lnInu(0,5e0),   3.276297109617906);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e1),   7.929768918237150);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e2),   96.77847637380128);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e3),   995.6271838273042);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e4),   9994.475891280807);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e5),   99993.32459873431);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e6),   999992.1733061878);

    AssertAlmostEqual(&test, bessel_lnInu(1,1e-6), -22.04766947825915);
    AssertAlmostEqual(&test, bessel_lnInu(1,1),    -1.225791352644727);
    AssertAlmostEqual(&test, bessel_lnInu(1,3),     1.131235470744604);
    AssertAlmostEqual(&test, bessel_lnInu(1,1e6),   999992.1733051878);

    AssertAlmostEqual(&test, bessel_lnInu(2,1e-6), -37.47261794865755);
    AssertAlmostEqual(&test, bessel_lnInu(2,1),    -2.862970265776753);
    AssertAlmostEqual(&test, bessel_lnInu(2,5),     2.622265862896675);
    AssertAlmostEqual(&test, bessel_lnInu(2,1e6),   999992.1733031878);

    AssertAlmostEqual(&test, bessel_lnInu(3,1), -4.824473578629219);

    AssertAlmostEqual(&test, bessel_lnInu(23,1e-6), -394.1439513814884);
    AssertAlmostEqual(&test, bessel_lnInu(23,5),    -31.40382021014728);
    AssertAlmostEqual(&test, bessel_lnInu(23,1e6),   999992.1730301876);

    AssertAlmostEqual(&test, bessel_lnInu(119,1e-6), -2189.202200199878);
    AssertAlmostEqual(&test, bessel_lnInu(119,0.5),  -621.0792579289692);
    AssertAlmostEqual(&test, bessel_lnInu(119,3),    -406.9458492626251);
    AssertAlmostEqual(&test, bessel_lnInu(119,30),   -129.9524456900199);
    AssertAlmostEqual(&test, bessel_lnInu(119,300),   272.6929318295042);
    AssertAlmostEqual(&test, bessel_lnInu(119,1e6),   999992.1661661842);

    AssertAlmostEqual(&test, bessel_lnInu(702,1e-6),  -14098.666835519577122094);
    AssertAlmostEqual(&test, bessel_lnInu(702,1e-4),  -10863.534779862939382744);
    AssertAlmostEqual(&test, bessel_lnInu(702,3),     -3621.4923374733442116413);
    AssertAlmostEqual(&test, bessel_lnInu(702,1234),   1034.4300403851143436433);
    AssertAlmostEqual(&test, bessel_lnInu(702,12345),  12319.387046237228462572);

    AssertAlmostEqual(&test, bessel_lnInu(1000,1e-6), -20431.4944983961827997);
    AssertAlmostEqual(&test, bessel_lnInu(1000,1e-3), -13520.2853417743050538);
    AssertAlmostEqual(&test, bessel_lnInu(1000,1e-2), -11216.5489562090494164);
    AssertAlmostEqual(&test, bessel_lnInu(1000,1e-1), -8912.81256819721365224);
    AssertAlmostEqual(&test, bessel_lnInu(1000,1),    -6609.07593552739598009);
    AssertAlmostEqual(&test, bessel_lnInu(1000,3),    -5509.91234371294526732);
    AssertAlmostEqual(&test, bessel_lnInu(1000,1e3),   527.852986878681152219);

    AssertAlmostEqual(&test, bessel_lnInu(2000,7), -10704.166550337010374);

    return test_results(&test, stderr);
}

int test_besselK(void)
{
    unittest_t test;
    unittest_init(&test, "Bessel function K_nu", "Test modified Bessel function K_nu");

    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-6), 7.133545631626864);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-5), 5.9822440851298415);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-4), 4.830861538632819);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-3), 3.6786689921357962);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-2), 2.5183764456387734);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-1), 1.2770838991417504);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e0), -0.7742086473552725);
    AssertAlmostEqual(&test, bessel_lnKnu(0,5e0), -5.5789276035723227);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e1), -10.925501193852295);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e2), -102.07679374034932);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e3), -1003.2280862868463);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e4), -10004.379378833343);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e5), -100005.53067137984);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e6), -1.0000066819639263e6);

    AssertAlmostEqual(&test, bessel_lnKnu(1,1e-6), 20.949057189590638);
    AssertAlmostEqual(&test, bessel_lnKnu(1,1), -0.08106146679532716);
    AssertAlmostEqual(&test, bessel_lnKnu(1,3),   -3.0358327192375464858953059975);
    AssertAlmostEqual(&test, bessel_lnKnu(1,1e6), -1.00000668196292633e6);

    AssertAlmostEqual(&test, bessel_lnKnu(2,1e-6), 35.863180036223355);
    AssertAlmostEqual(&test, bessel_lnKnu(2,1), 1.17170150170004);
    AssertAlmostEqual(&test, bessel_lnKnu(2,5), -5.036603312746961080665958204772);
    AssertAlmostEqual(&test, bessel_lnKnu(2,1e6), -1.0000066819609263389096196908e6);

    AssertAlmostEqual(&test, bessel_lnKnu(3,1), 2.8367092652889516);

    AssertAlmostEqual(&test, bessel_lnKnu(4,1e15), -1.00000000000001704359684e15);

    AssertAlmostEqual(&test, bessel_lnKnu(23,1e-6), 390.29380377977833);
    AssertAlmostEqual(&test, bessel_lnKnu(23,5), 27.5314997887589672718741222750056);
    AssertAlmostEqual(&test, bessel_lnKnu(23,1e6), -1.00000668168792647542217765299e6);

    AssertAlmostEqual(&test, bessel_lnKnu(119,1e-6), 2183.7257366479472175742539693253862993069);
    AssertAlmostEqual(&test, bessel_lnKnu(119,0.5), 615.6027856231534);
    AssertAlmostEqual(&test, bessel_lnKnu(119,3), 401.4690706673959);
    AssertAlmostEqual(&test, bessel_lnKnu(119,30), 124.44542144141829);
    AssertAlmostEqual(&test, bessel_lnKnu(119,300), -279.16349731660983);

    AssertAlmostEqual(&test, bessel_lnKnu(702,1e-4), 10856.28698728117152647293);
    AssertAlmostEqual(&test, bessel_lnKnu(702,3), 3614.24453577321548255948381274);
    AssertAlmostEqual(&test, bessel_lnKnu(702,1234), -1042.3815655681729711090061175312483747);
    AssertAlmostEqual(&test, bessel_lnKnu(702,12345), -12329.50281632819683895772331427);

    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e-6),  20423.89309606159906635627);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e-3),  13512.68393943972082096386);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e-2),  11208.94755387441573291358);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e-1),  8905.211165857634910126568);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,3),     5502.310936882873879713131);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e3),  -535.8007129753599475405978);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e6),  -1.0000061814642183370632e6);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e10), -1.000000001128708406232e10);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e15), -1.000000000000017043596e15);

    return test_results(&test, stderr);
}

int test_Lambda(void)
{
    unittest_t test;
    unittest_init(&test, "Lambda", "Test Lambda function for various parameters");

    AssertAlmostEqual(&test, casimir_lnLambda(50,50,0),   -3.921875301871158);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,1),   -11.76572394002363);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,50),  -367.6612508574347);
    AssertAlmostEqual(&test, casimir_lnLambda(50,20,10),  -72.40153130583653);
    AssertAlmostEqual(&test, casimir_lnLambda(5,7,3),     -12.77235621226475);
    AssertAlmostEqual(&test, casimir_lnLambda(16,6,4),    -20.7139882421892);
    AssertAlmostEqual(&test, casimir_lnLambda(100,6,4),   -28.88322376001037);
    AssertAlmostEqual(&test, casimir_lnLambda(100,100,0), -4.61013297533022);
    AssertAlmostEqual(&test, casimir_lnLambda(100,100,50), -461.1524718729809);
    AssertAlmostEqual(&test, casimir_lnLambda(100,201,10), -104.0968227550132);
    AssertAlmostEqual(&test, casimir_lnLambda(200,200,0), -5.300808027860489);
    AssertAlmostEqual(&test, casimir_lnLambda(200,100,70), -690.4926643211061);
    AssertAlmostEqual(&test, casimir_lnLambda(500,500,0), -6.215606600751781);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,0), -6.908254904273569);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,1), -20.72476496257093);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,2), -34.54127302286429);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,3), -48.35777708713769);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,10), -145.0730258071027);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,20), -283.2357984094714);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,50), -697.6921182236243);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,100), -1388.225291090218);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,499), -6856.454154873699);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,500), -6869.983981361874);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,501), -6883.51247629193);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,999), -13205.83170295854);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,1000), -13213.43260541808);

    return test_results(&test, stderr);
}

int test_Xi(void)
{
    int sign;
    unittest_t test;
    unittest_init(&test, "Xi", "Test Xi function for various parameters");

    AssertAlmostEqual(&test, casimir_lnXi(4,3,2,&sign), -0.1101206735572);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, casimir_lnXi(4,2,2,&sign), -2.394730234408415);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, casimir_lnXi(100,100,100,&sign), 587.0039751538028);
    AssertAlmostEqual(&test, casimir_lnXi(100,100,50,&sign),  696.7380895450116);
    AssertAlmostEqual(&test, casimir_lnXi(17,14,10,&sign),    45.8135805997528);

    return test_results(&test, stderr);
}

int test_mie(void)
{
    int sign;
    unittest_t test;
    unittest_init(&test, "Mie", "Test Mie functions al,bl for various parameters");

    /* b_l */
    AssertAlmostEqual(&test, casimir_lnb(5,3,&sign), -3.206110089012862);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, casimir_lnb(6,3,&sign), -6.093433624873396);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, casimir_lnb(50,1,&sign),   -365.8137152817732);
    AssertAlmostEqual(&test, casimir_lnb(50,100,&sign),  174.3104974165916);
    AssertAlmostEqual(&test, casimir_lnb(100,2,&sign),  -726.3166073149845);
    AssertAlmostEqual(&test, casimir_lnb(100,100,&sign), 104.9919945452843);
    AssertAlmostEqual(&test, casimir_lnb(100,200,&sign), 349.7964954441692);
    AssertAlmostEqual(&test, casimir_lnb(100,300,&sign), 565.9447715085943);
    AssertAlmostEqual(&test, casimir_lnb(40,0.01,&sign), -648.6664276814638);
    AssertAlmostEqual(&test, casimir_lnb(4,0.01,&sign),  -52.95166526324419);

    /* a_l */
    AssertAlmostEqual(&test, casimir_lna(3,3,&sign), 1.692450306201961);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, casimir_lna(4,3,&sign), -0.50863950281017);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, casimir_lna(70,1,&sign),     -557.4493819729695);
    AssertAlmostEqual(&test, casimir_lna(70,2,&sign),     -459.6943641467319);
    AssertAlmostEqual(&test, casimir_lna(70,70,&sign),     73.02602649528605);
    AssertAlmostEqual(&test, casimir_lna(70,100,&sign),    151.4135590544529);
    AssertAlmostEqual(&test, casimir_lna(7,0.2,&sign),    -50.34157726932342);
    AssertAlmostEqual(&test, casimir_lna(20,0.1,&sign),   -206.3146872637107);
    AssertAlmostEqual(&test, casimir_lna(20,0.01,&sign),  -300.7209163862779);
    AssertAlmostEqual(&test, casimir_lna(30,0.01,&sign),  -471.3445070668955);
    AssertAlmostEqual(&test, casimir_lna(30,0.001,&sign), -611.8021993589887);

    return test_results(&test, stderr);
}

int test_integration(void)
{
    casimir_integrals_t cint;
    unittest_t test;
    unittest_init(&test, "Integration", "Test integration for various parameters");

    casimir_integrate(&cint, 4, 4, 0, 0.01);
    AssertAlmostEqual(&test, cint.logB, 56.28387814539346);

    casimir_integrate(&cint, 4, 4, 1, 0.01);
    AssertAlmostEqual(&test, cint.signA*exp(cint.logA), +2.4806179125126554e17);
    AssertAlmostEqual(&test, cint.signB*exp(cint.logB), -2.2226323455151368e24);
    AssertAlmostEqual(&test, cint.signC*exp(cint.logC), -6.9457269656680333e20);
    AssertAlmostEqual(&test, cint.signD*exp(cint.logD), +6.9457269656680333e20);

    casimir_integrate(&cint, 40, 40, 1, 0.5);
    AssertAlmostEqual(&test, cint.signA*exp(cint.logA), +1.5754477603435539e159);
    AssertAlmostEqual(&test, cint.signB*exp(cint.logB), -6.3723632215476122e166);
    AssertAlmostEqual(&test, cint.signC*exp(cint.logC), -9.9568222699306801e162);
    AssertAlmostEqual(&test, cint.signD*exp(cint.logD), +9.9568222699306801e162);

    casimir_integrate(&cint, 40, 40, 40, 2);
    AssertAlmostEqual(&test, cint.signA*exp(cint.logA), +6.4140686579381969e91);
    AssertAlmostEqual(&test, cint.signB*exp(cint.logB), -1.0147301906459434e95);
    AssertAlmostEqual(&test, cint.signC*exp(cint.logC), -2.5352219594503741e93);
    AssertAlmostEqual(&test, cint.signD*exp(cint.logD), +2.5352219594503736e93);

    casimir_integrate(&cint, 7, 4, 3, 17);
    AssertAlmostEqual(&test, cint.signA*exp(cint.logA), +4.8180365200137397e-9);
    AssertAlmostEqual(&test, cint.signB*exp(cint.logB), -1.3731640166794149e-8);
    AssertAlmostEqual(&test, cint.signC*exp(cint.logC), -6.7659079909128738e-9);
    AssertAlmostEqual(&test, cint.signD*exp(cint.logD), +9.44463292099617e-9);

    casimir_integrate(&cint, 40, 40, 0, 5);
    AssertAlmostEqual(&test, cint.signB*exp(cint.logB), -6.0455421304871757e85);

    casimir_integrate(&cint, 100, 41, 0, 5);
    AssertAlmostEqual(&test, cint.signB*exp(cint.logB), 8.8689390374540308e185);

    return test_results(&test, stderr);
}

int main(int argc, char *argv[])
{
    test_Lambda();
    test_Xi();
    test_integration();
    test_mie();
    test_besselI();
    test_besselK();
    test_givens();
    test_logadd();
    test_logdet();
    
    return 0;
}

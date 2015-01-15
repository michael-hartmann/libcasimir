#include <math.h>
#include <stdio.h>

#include "libcasimir.h"
#include "integration.h"
#include "sfunc.h"
#include "matrix.h"
#include "unittest.h"

#include "tests.h"

int test_casimirF()
{
    unittest_t test;
    casimir_t casimir;
    double F;

    unittest_init(&test, "casimirF", "Compare free energies");

    casimir_init(&casimir, 0.85, 2.7);
    casimir_set_precision(&casimir, 1e-14);
    casimir_set_lmax(&casimir, 30);
    F = casimir_F(&casimir, NULL);
    AssertAlmostEqual(&test, F, -1.34361893570375);
    casimir_free(&casimir);

    casimir_init(&casimir, 0.7, 1);
    casimir_set_precision(&casimir, 1e-14);
    casimir_set_lmax(&casimir, 15);
    F = casimir_F(&casimir, NULL);
    AssertAlmostEqual(&test, F, -0.220709222562969);
    casimir_free(&casimir);

    casimir_init(&casimir, 0.85, 2.7);
    casimir_set_precision(&casimir, 1e-14);
    casimir_set_cores(&casimir, 10);
    casimir_set_lmax(&casimir, 30);
    F = casimir_F(&casimir, NULL);
    AssertAlmostEqual(&test, F, -1.34361893570375);
    casimir_free(&casimir);

    return test_results(&test, stderr);
}

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

    /*
    logdet = casimir_logdetD(&casimir, 0, 0, NULL);
    AssertAlmostEqual(&test, logdet, -3.45236396285874);

    logdet = casimir_logdetD(&casimir, 0, 1, NULL);
    AssertAlmostEqual(&test, logdet, -2.63586999367158);

    logdet = casimir_logdetD(&casimir, 0, 10, NULL);
    AssertAlmostEqual(&test, logdet, -0.0276563864490425);
    */

    casimir_mie_cache_init(&cache, 1);
    casimir_mie_cache_alloc(&casimir, &cache);

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
    matrix_edouble_t *M;
    unittest_t test;
    unittest_init(&test, "QR decomposition", "Test QR decomposition using givens rotation");

    {
        M = matrix_edouble_alloc(2);
        matrix_set(M, 0,0, 20e100);
        matrix_set(M, 0,1, 2);
        matrix_set(M, 1,0, 1);
        matrix_set(M, 1,1, 1e-100);

        matrix_edouble_balance(M);

        AssertAlmostEqual(&test, matrix_edouble_logdet(M), log(18));

        matrix_edouble_free(M);
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

    AssertAlmostEqual(&test, casimir_lnLambda(1,1,0,NULL), 0.40546510810816);
    AssertAlmostEqual(&test, casimir_lnLambda(1,1,1,NULL), -0.28768207245178);
    AssertAlmostEqual(&test, casimir_lnLambda(2,1,1,NULL), -1.130881549236895);
    AssertAlmostEqual(&test, casimir_lnLambda(4,5,3,NULL), -10.11921344416614);
    AssertAlmostEqual(&test, casimir_lnLambda(100,1,0,NULL), -1.755760343331055);
    AssertAlmostEqual(&test, casimir_lnLambda(100,1,1,NULL), -6.712479285025704);
    AssertAlmostEqual(&test, casimir_lnLambda(100,33,17,NULL), -140.5615114663231);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,0,NULL),   -3.921875301871158+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,1,NULL),   -11.76572394002363+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,50,NULL),  -367.6612508574347+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(50,20,10,NULL),  -72.40153130583653+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(5,7,3,NULL),     -12.77235621226475+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(16,6,4,NULL),    -20.7139882421892+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(100,6,4,NULL),   -28.88322376001037+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(100,100,0,NULL), -4.61013297533022+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(100,100,50,NULL), -461.1524718729809+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(100,201,10,NULL), -104.0968227550132+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(200,200,0,NULL), -5.300808027860489+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(200,100,70,NULL), -690.4926643211061+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(500,500,0,NULL), -6.215606600751781+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,0,NULL), -6.908254904273569+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,1,NULL), -20.72476496257093+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,2,NULL), -34.54127302286429+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,3,NULL), -48.35777708713769+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,10,NULL), -145.0730258071027+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,20,NULL), -283.2357984094714+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,50,NULL), -697.6921182236243+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,100,NULL), -1388.225291090218+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,499,NULL), -6856.454154873699+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,500,NULL), -6869.983981361874+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,501,NULL), -6883.51247629193+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,999,NULL), -13205.83170295854+0.6931471805599453);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,1000,NULL), -13213.43260541808+0.6931471805599453);

    return test_results(&test, stderr);
}

int test_Xi(void)
{
    int sign;
    unittest_t test;
    unittest_init(&test, "Xi", "Test Xi function for various parameters");

    AssertAlmostEqual(&test, casimir_lnXi(1,1,0,&sign), -3.060270794691562);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, casimir_lnXi(1,1,1,&sign), -3.753417975251507);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, casimir_lnXi(3,2,1,&sign), -1.817138914330164);
    AssertEqual(&test, sign, +1);

    AssertAlmostEqual(&test, casimir_lnXi(4,3,2,&sign), -0.1101206735572);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, casimir_lnXi(4,2,2,&sign), -2.394730234408415);
    AssertEqual(&test, sign, +1);

    AssertAlmostEqual(&test, casimir_lnXi(11,1,0,&sign), 2.60283664295575);
    AssertAlmostEqual(&test, casimir_lnXi(11,7,0,&sign), 19.22557931884024);
    AssertAlmostEqual(&test, casimir_lnXi(11,7,5,&sign), 16.28731202862324);
    AssertAlmostEqual(&test, casimir_lnXi(201,7,5,&sign), 623.3839523251071);

    AssertAlmostEqual(&test, casimir_lnXi(100,10,0,&sign), 269.8159771440838);
    AssertAlmostEqual(&test, casimir_lnXi(100,10,1,&sign), 269.7633468887551);
    AssertAlmostEqual(&test, casimir_lnXi(100,10,10,&sign), 263.2542489687728);

    AssertAlmostEqual(&test, casimir_lnXi(100,100,100,&sign), 587.0039751538028);
    AssertAlmostEqual(&test, casimir_lnXi(100,100,50,&sign),  696.7380895450116);
    AssertAlmostEqual(&test, casimir_lnXi(100,100,0,&sign),  722.7572112350813);
    AssertAlmostEqual(&test, casimir_lnXi(100,100,1,&sign),  722.747260904228);
    AssertAlmostEqual(&test, casimir_lnXi(17,14,10,&sign),    45.8135805997528);

    return test_results(&test, stderr);
}

static double _mie_lna_perf(int l, double arg, int *sign)
{
    casimir_t self;
    casimir_init(&self, 0.5, 2*arg);
    return casimir_lna_perf(&self, l, 1, sign);
}

static double _mie_lnb_perf(int l, double arg, int *sign)
{
    casimir_t self;
    casimir_init(&self, 0.5, 2*arg);
    return casimir_lnb_perf(&self, l, 1, sign);
}

int test_mie(void)
{
    int sign;
    casimir_t self;
    unittest_t test;
    unittest_init(&test, "Mie", "Test Mie functions al,bl for various parameters");

    casimir_init(&self, 0.5,2);

    /* b_l */
    AssertAlmostEqual(&test, _mie_lnb_perf(5,3,&sign), -3.206110089012862);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, _mie_lnb_perf(6,3,&sign), -6.093433624873396);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lnb_perf(50,1,&sign),   -365.8137152817732);
    AssertAlmostEqual(&test, _mie_lnb_perf(50,100,&sign),  174.3104974165916);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,2,&sign),  -726.3166073149845);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,100,&sign), 104.9919945452843);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,200,&sign), 349.7964954441692);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,300,&sign), 565.9447715085943);
    AssertAlmostEqual(&test, _mie_lnb_perf(40,0.01,&sign), -648.6664276814638);
    AssertAlmostEqual(&test, _mie_lnb_perf(4,0.01,&sign),  -52.95166526324419);

    /* a_l */
    AssertAlmostEqual(&test, _mie_lna_perf(3,3,&sign), 1.692450306201961);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, _mie_lna_perf(4,3,&sign), -0.50863950281017);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, _mie_lna_perf(70,1,&sign),     -557.4493819729695);
    AssertAlmostEqual(&test, _mie_lna_perf(70,2,&sign),     -459.6943641467319);
    AssertAlmostEqual(&test, _mie_lna_perf(70,70,&sign),     73.02602649528605);
    AssertAlmostEqual(&test, _mie_lna_perf(70,100,&sign),    151.4135590544529);
    AssertAlmostEqual(&test, _mie_lna_perf(7,0.2,&sign),    -50.34157726932342);
    AssertAlmostEqual(&test, _mie_lna_perf(20,0.1,&sign),   -206.3146872637107);
    AssertAlmostEqual(&test, _mie_lna_perf(20,0.01,&sign),  -300.7209163862779);
    AssertAlmostEqual(&test, _mie_lna_perf(30,0.01,&sign),  -471.3445070668955);
    AssertAlmostEqual(&test, _mie_lna_perf(30,0.001,&sign), -611.8021993589887);

    return test_results(&test, stderr);
}

int test_mie_drude(void)
{
    int sign_a, sign_b;
    double lna, lnb;
    casimir_t casimir;
    unittest_t test;
    unittest_init(&test, "Mie Drude", "Test Mie functions al,bl for various parameters");

    casimir_init(&casimir, 0.85, 2.7);
    casimir_set_omegap_sphere(&casimir, 1);
    casimir_set_gamma_sphere(&casimir, 1);

    // void casimir_lnab(casimir_t *self, const int n, const int l, double *lna, double *lnb, int *sign_a, int *sign_b);

    casimir_lnab(&casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -3.4553021173541371702);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -5.8688410248499200961);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 2, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 1.7793937447527983231);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, 0.47334767009093381418);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 2, 7, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -5.8888749734115481971);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -8.4332679619305466190);
    AssertEqual(&test, sign_b, +1);

    casimir_free(&casimir);


    casimir_init(&casimir, 0.95, 0.1);
    casimir_set_omegap_sphere(&casimir, 0.1);
    casimir_set_gamma_sphere(&casimir, 1.4);

    casimir_lnab(&casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -11.294081746210160588);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -18.282872530968926839);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 1, 150, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -2128.7634085682857047);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -2144.8884593165140032);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 100, 15, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -13.098866648025603965);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -15.636661776735083507);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 200, 20, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 5.8518266664014060296);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, 3.9865380993256289699);
    AssertEqual(&test, sign_b, +1);

    casimir_free(&casimir);


    casimir_init(&casimir, 0.5, 1);
    casimir_set_omegap_sphere(&casimir, 1e-4);
    casimir_set_gamma_sphere(&casimir,  1e-4);

    casimir_lnab(&casimir, 1, 7, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -55.208389583053445904);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -62.068548807408056689);
    AssertEqual(&test, sign_b, +1);

    return test_results(&test, stderr);
}

int test_integration(void)
{
    casimir_integrals_t cint;
    unittest_t test;
    unittest_init(&test, "Integration", "Test integration for various parameters");

    casimir_integrate(&cint, 3, 2, 1, 2);
    AssertAlmostEqual(&test, cint.lnA_TM, -4.094372316589062);
    AssertAlmostEqual(&test, cint.lnB_TM, -1.970116759119433);
    AssertAlmostEqual(&test, cint.lnC_TM, -3.298725852652321);

    casimir_integrate(&cint, 4, 4, 0, 0.005);
    AssertAlmostEqual(&test, cint.lnB_TM, 56.28387814539346+0.6931471805599453);

    casimir_integrate(&cint, 4, 4, 1, 0.005);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +2.4806179125126554e17*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -2.2226323455151368e24*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -6.9457269656680333e20*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +6.9457269656680333e20*-2);

    casimir_integrate(&cint, 40, 40, 1, 0.25);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +1.5754477603435539e159*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -6.3723632215476122e166*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -9.9568222699306801e162*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +9.9568222699306801e162*-2);

    casimir_integrate(&cint, 40, 40, 40, 1);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +6.4140686579381969e91*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -1.0147301906459434e95*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -2.5352219594503741e93*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +2.5352219594503736e93*-2);

    casimir_integrate(&cint, 7, 4, 3, 8.5);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +4.8180365200137397e-9*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -1.3731640166794149e-8*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -6.7659079909128738e-9*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +9.44463292099617e-9*-2);

    casimir_integrate(&cint, 40, 40, 0, 2.5);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -6.0455421304871757e85*-2);

    casimir_integrate(&cint, 100, 41, 0, 2.5);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), 8.8689390374540308e185*-2);

    return test_results(&test, stderr);
}

int test_doublefact()
{
    unittest_t test;

    unittest_init(&test, "doublefact", "Test double factorial");

    AssertEqual(&test, ln_doublefact(0),   0);
    AssertEqual(&test, ln_doublefact(1),   0);

    AssertAlmostEqual(&test, ln_doublefact(2),   0.693147180559945309417232121458176568075500134360255254120680);
    AssertAlmostEqual(&test, ln_doublefact(3),   1.098612288668109691395245236922525704647490557822749451734694);
    AssertAlmostEqual(&test, ln_doublefact(4),   2.079441541679835928251696364374529704226500403080765762362040);
    AssertAlmostEqual(&test, ln_doublefact(5),   2.708050201102210065996004570148713344173091912091267173647342);

    AssertAlmostEqual(&test, ln_doublefact(51),  77.07730784751820516445630791916220149454265081179451116106491);
    AssertAlmostEqual(&test, ln_doublefact(52),  79.28352845556058002961361747099464590876118666235724979722338);
    
    AssertAlmostEqual(&test, ln_doublefact(100), 183.1351259797702975383987999237883518258868824082753118006299);
    AssertAlmostEqual(&test, ln_doublefact(101), 185.2193700926344520565653917127802753588279212418095676769807);
    
    AssertAlmostEqual(&test, ln_doublefact(333), 803.8068699169127948868947011937225159269950980762322648014825);
    AssertAlmostEqual(&test, ln_doublefact(334), 806.9389802679216196222732317043854384535329497966731766081127);

    AssertAlmostEqual(&test, ln_doublefact(499), 1303.998431529316796415491154267014961848434758393410615755047);
    AssertAlmostEqual(&test, ln_doublefact(500), 1307.332026930839287985819861539370764480277316393340117718467);

    return test_results(&test, stderr);
}

int test_plm()
{
    int sign;
    plm_combination_t comb;
    unittest_t test;

    unittest_init(&test, "Plm", "Test associated Legendre polynomials");

    /* Plm */
    AssertAlmostEqual(&test, plm_Plm(2,0,3), 13);
    AssertAlmostEqual(&test, plm_Plm(2,1,3), -25.45584412271572);
    AssertAlmostEqual(&test, plm_Plm(2,2,3), -24);

    AssertAlmostEqual(&test, plm_Plm(3,0,2), 17);
    AssertAlmostEqual(&test, plm_Plm(3,1,2), -49.36344801571300286553222073);
    AssertAlmostEqual(&test, plm_Plm(3,2,2), -90);
    AssertAlmostEqual(&test, plm_Plm(3,3,2), 77.94228634059947);

    AssertAlmostEqual(&test, plm_Plm(3,2,200), -119997000);

    AssertAlmostEqual(&test, plm_Plm(7,0,2), 2199.125);
    AssertAlmostEqual(&test, plm_Plm(7,1,2), -15209.24639443778);
    AssertAlmostEqual(&test, plm_Plm(7,2,2), -88026.75);
    AssertAlmostEqual(&test, plm_Plm(7,3,2), 414721.1628325371);
    AssertAlmostEqual(&test, plm_Plm(7,4,2), 1528065.0);
    AssertAlmostEqual(&test, plm_Plm(7,5,2), -4132071.339203708);
    AssertAlmostEqual(&test, plm_Plm(7,6,2), -7297290.0);
    AssertAlmostEqual(&test, plm_Plm(7,7,2), 6319638.518782143);

    AssertAlmostEqual(&test, plm_Plm(8,0,2), 7691.1484375);
    AssertAlmostEqual(&test, plm_Plm(8,1,2), -60890.46264643482);
    AssertAlmostEqual(&test, plm_Plm(8,2,2), -413142.1875);
    AssertAlmostEqual(&test, plm_Plm(8,3,2), 2354110.359916711);
    AssertAlmostEqual(&test, plm_Plm(8,4,2), 1.0957629375e7);
    AssertAlmostEqual(&test, plm_Plm(8,5,2), -4.0024377285620242e7);
    AssertAlmostEqual(&test, plm_Plm(8,6,2), -1.076350275e8);
    AssertAlmostEqual(&test, plm_Plm(8,7,2), 1.8958915556346425e8);
    AssertAlmostEqual(&test, plm_Plm(8,8,2), 1.64189025e8);

    AssertAlmostEqual(&test, plm_Plm(30,30,10), -2.512712634569675964e70);
    AssertAlmostEqual(&test, plm_Plm(30,0,10), 1.022858298005579512e38);

    AssertAlmostEqual(&test, plm_lnPlm(300,0,100, &sign),    1586.0630493580697);
    AssertAlmostEqual(&test, plm_lnPlm(300,200,100, &sign),  2637.2261846173691);

    AssertAlmostEqual(&test, plm_lnPlm(300,1,100, &sign),    1591.7668317492472);
    AssertAlmostEqual(&test, plm_lnPlm(300,201,100, &sign),  2641.8313213287660);
    AssertAlmostEqual(&test, plm_lnPlm(300,201,1000, &sign), 3332.6176009928417);
    AssertAlmostEqual(&test, plm_lnPlm(300,201,5000, &sign), 3815.4490789776808);

    /* dPlm */
    AssertAlmostEqual(&test, plm_dPlm(3,0,3), 66);
    AssertAlmostEqual(&test, plm_dPlm(3,1,3), -197.2827919510467);
    AssertAlmostEqual(&test, plm_dPlm(3,2,3), -390);
    AssertAlmostEqual(&test, plm_dPlm(3,3,3), 381.8376618407357);


    plm_PlmPlm(4, 3, 2, 2, &comb);
    AssertAlmostEqual(&test, comb.lnPl1mPl2m, log(54675));
    AssertEqual(&test, comb.sign_Pl1mPl2m, +1);

    AssertAlmostEqual(&test, comb.lnPl1mdPl2m, log(100237.5));
    AssertEqual(&test, comb.sign_Pl1mdPl2m, +1);

    AssertAlmostEqual(&test, comb.lndPl1mPl2m, log(129600.0));
    AssertEqual(&test, comb.sign_dPl1mPl2m, +1);

    AssertAlmostEqual(&test, comb.lndPl1mdPl2m, log(237600.0));
    AssertEqual(&test, comb.sign_dPl1mdPl2m, +1);


    plm_PlmPlm(4, 3, 1, 2, &comb);
    AssertAlmostEqual(&test, comb.lnPl1mPl2m, log(10687.5));
    AssertEqual(&test, comb.sign_Pl1mPl2m, -1);

    AssertAlmostEqual(&test, comb.lnPl1mdPl2m, log(18375.0));
    AssertEqual(&test, comb.sign_Pl1mdPl2m, -1);

    AssertAlmostEqual(&test, comb.lndPl1mPl2m, log(24438.75));
    AssertEqual(&test, comb.sign_dPl1mPl2m, -1);

    AssertAlmostEqual(&test, comb.lndPl1mdPl2m, log(42017.5));
    AssertEqual(&test, comb.sign_dPl1mdPl2m, -1);


    plm_PlmPlm(7, 5, 3, 2, &comb);
    AssertAlmostEqual(&test, comb.lnPl1mPl2m, 22.09944134068912);
    AssertEqual(&test, comb.sign_Pl1mPl2m, -1);

    AssertAlmostEqual(&test, comb.lnPl1mdPl2m, 23.20753237331177);
    AssertEqual(&test, comb.sign_Pl1mdPl2m, -1);

    AssertAlmostEqual(&test, comb.lndPl1mPl2m, 23.51706034623926);
    AssertEqual(&test, comb.sign_dPl1mPl2m, -1);

    AssertAlmostEqual(&test, comb.lndPl1mdPl2m, 24.62515137886192);
    AssertEqual(&test, comb.sign_dPl1mdPl2m, -1);

    return test_results(&test, stderr);
}

int main(int argc, char *argv[])
{
    /*
    test_plm();
    test_doublefact();
    test_Lambda();
    test_Xi();
    test_integration();
    test_mie();
    test_mie_drude();
    test_besselI();
    test_besselK();
    test_givens();
    test_logadd();
    */
    test_logdet();
    //test_casimirF();
    
    return 0;
}

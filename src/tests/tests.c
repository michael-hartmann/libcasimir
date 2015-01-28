#include <math.h>
#include <stdio.h>

#include "libcasimir.h"
#include "integration.h"
#include "sfunc.h"
#include "matrix.h"
#include "unittest.h"


/* prototypes */
int test_casimirF(void);
int test_logdet(void);
int test_Lambda(void);
int test_Xi(void);
int test_integration(void);
int test_integration_drude(void);
int test_mie(void);
int test_besselI(void);
int test_besselK(void);
int test_givens(void);
int test_mie_drude(void);
int test_doublefact(void);
int test_plm(void);
int test_epsilon(void);
int test_fresnel(void);

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

    AssertAlmostEqual(&test, casimir_lnLambda(1,1,0,NULL),           0.40546510810816438197801311546434913657);
    AssertAlmostEqual(&test, casimir_lnLambda(1,1,1,NULL),          -0.2876820724517809274392190059938274315);

    AssertAlmostEqual(&test, casimir_lnLambda(2,1,1,NULL),          -1.1308815492368952772317071947645221687);
    AssertAlmostEqual(&test, casimir_lnLambda(4,5,3,NULL),          -10.119213444166137830580776265774926608);
    AssertAlmostEqual(&test, casimir_lnLambda(5,7,3,NULL),          -12.079209031704799296645550329264124926);
    AssertAlmostEqual(&test, casimir_lnLambda(16,6,4,NULL),         -20.020841061629258759138303870262239174);

    AssertAlmostEqual(&test, casimir_lnLambda(50,50,0,NULL),        -3.2287281213112123793793323757932149304);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,1,NULL),        -11.072576759463684209642863018499419017);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,50,NULL),       -366.96810367687470252345932574544885295);
    AssertAlmostEqual(&test, casimir_lnLambda(50,20,10,NULL),       -71.708384125276581706972600369949945466);

    AssertAlmostEqual(&test, casimir_lnLambda(100,1,0,NULL),        -1.7557603433310553429384254137655599901);
    AssertAlmostEqual(&test, casimir_lnLambda(100,1,1,NULL),        -6.7124792850257034071071320626355070602);
    AssertAlmostEqual(&test, casimir_lnLambda(100,6,4,NULL),        -28.190076579450425813590378183217220139);
    AssertAlmostEqual(&test, casimir_lnLambda(100,33,17,NULL),      -140.56151146632312845247188373793277075);
    AssertAlmostEqual(&test, casimir_lnLambda(100,100,0,NULL),      -3.9169857947702750678548639429954691168);
    AssertAlmostEqual(&test, casimir_lnLambda(100,100,50,NULL),     -460.45932469242092685829036419158558034);
    AssertAlmostEqual(&test, casimir_lnLambda(100,201,10,NULL),     -103.40367557445323498326245865777798370);

    AssertAlmostEqual(&test, casimir_lnLambda(200,200,0,NULL),      -4.6076608473005432426678531767638037076);
    AssertAlmostEqual(&test, casimir_lnLambda(200,100,70,NULL),     -689.79951714054617706762753098998776201);
    AssertAlmostEqual(&test, casimir_lnLambda(500,500,0,NULL),      -5.5224594201918359560709544932885934909);

    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,0,NULL),    -6.2151077237136242278894880259320868284);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,1,NULL),    -20.031617782010981865164246152958807085);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,2,NULL),    -33.848125842304345491790305669883313001);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,3,NULL),    -47.664629906577744997358213920918706752);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,10,NULL),   -144.37987862654270719988337448491110499);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,20,NULL),   -282.54265122891145618177303437612939551);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,50,NULL),   -696.99897104306439045438129052648809049);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,100,NULL),  -1387.5321439096580157639308668460135449);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,499,NULL),  -6855.7610076931390806110616109722327157);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,500,NULL),  -6869.2908341813142468437485944654174625);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,501,NULL),  -6882.8193291113699005405673096813816430);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,999,NULL),  -13205.138555777978298175816431569216187);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,1000,NULL), -13212.739458237520380537287638054727456);

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
    double T, RbyScriptL, omegap, gamma;;
    int sign_a, sign_b;
    double lna, lnb;
    casimir_t casimir;
    unittest_t test;
    unittest_init(&test, "Mie Drude", "Test Mie functions al,bl for various parameters");

    RbyScriptL = 0.85;
    T          = 2.7;
    omegap     = 1;
    gamma      = 1;
    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_omegap_sphere(&casimir, omegap);
    casimir_set_gamma_sphere(&casimir, gamma);

    // void casimir_lnab(casimir_t *self, const int n, const int l, double *lna, double *lnb, int *sign_a, int *sign_b);

    casimir_lnab(&casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -3.4553021173541333720476206556765874126);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -5.8688410248499158590684691465517966950);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 2, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 1.77939374475276907218961292198964700248);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, 0.47334767009104316463610995593888380542);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 2, 7, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -5.8888749734115470918000948427901346308);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -8.4332679619305632421065168698778165934);
    AssertEqual(&test, sign_b, +1);

    casimir_free(&casimir);


    RbyScriptL = 0.95;
    T          = 0.1;
    omegap     = 0.1;
    gamma      = 1.4;
    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_omegap_sphere(&casimir, omegap);
    casimir_set_gamma_sphere(&casimir, gamma);

    casimir_lnab(&casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -11.294081746210154580677515252479982979);
    AssertAlmostEqual(&test, lnb, -18.282872530972731935336084568403897476);

    casimir_lnab(&casimir, 1, 150, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -2128.7634085682857996403030143536620910);
    AssertAlmostEqual(&test, lnb, -2144.8884593589829670487681397791770641);

    casimir_lnab(&casimir, 100, 15, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -13.098866648021450586862321076066204163);
    AssertAlmostEqual(&test, lnb, -15.636661776753788796856285766794114116);

    casimir_lnab(&casimir, 200, 20, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 5.85182666636940197335162042117627870470);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, 3.98653809954659615487558488229399701526);
    AssertEqual(&test, sign_b, +1);

    casimir_free(&casimir);


    RbyScriptL = 0.5;
    T          = 1;
    omegap     = 1e-4;
    gamma      = 1e-4;
    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_omegap_sphere(&casimir, omegap);
    casimir_set_gamma_sphere(&casimir,  gamma);

    casimir_lnab(&casimir, 1, 7, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -55.208389631622379721466940745498367474);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -62.068501413438155604254192987355717706);
    AssertEqual(&test, sign_b, +1);

    casimir_free(&casimir);


    RbyScriptL = 0.5;
    T          = 1;
    omegap     = 1;
    gamma      = 1e-4;
    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_omegap_sphere(&casimir, omegap);
    casimir_set_gamma_sphere(&casimir,  gamma);

    casimir_lnab(&casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -3.8540239333215465024615528927141405724);
    AssertAlmostEqual(&test, lnb, -7.2588880177134734386242016983196388100);

    casimir_lnab(&casimir, 1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -7.6570049740166300703828407893184117471);
    AssertAlmostEqual(&test, lnb, -12.197169758793283148024686910238050626);

    casimir_lnab(&casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -12.482146722495979280724692478599547992);
    AssertAlmostEqual(&test, lnb, -17.727243821435497530589658157489862574);

    casimir_lnab(&casimir, 1, 4, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -17.949100704846886740879599953685442860);
    AssertAlmostEqual(&test, lnb, -23.709990656446819196291385264231700122);

    casimir_lnab(&casimir, 1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -23.891599427762872566501942669866619999);
    AssertAlmostEqual(&test, lnb, -30.060470310939164883544475705953387986);

    casimir_lnab(&casimir, 1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -58.285652601186559317984617617473593988);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -65.757440056322685520981556114527343549);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 5, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -0.7747752591093184484501684513821649046);
    AssertAlmostEqual(&test, lnb, -1.5892101743013427987322480836938505974);

    casimir_lnab(&casimir, 5, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -1.8145275883758257249878124562886200521);
    AssertAlmostEqual(&test, lnb, -3.4950263771356228454682283405268413444);

    casimir_lnab(&casimir, 5, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -3.6528777367390642776262971550796544139);
    AssertAlmostEqual(&test, lnb, -5.9236362448960108785813307258565097790);

    casimir_lnab(&casimir, 5, 4, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -6.0449566310035356718674254329839576890);
    AssertAlmostEqual(&test, lnb, -8.7689157304390758571133982892467814487);

    casimir_lnab(&casimir, 5, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -8.8673827346012321060116358013846191791);
    AssertAlmostEqual(&test, lnb, -11.960448909722605823940749038683258544);

    casimir_lnab(&casimir, 5, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -27.402076826754931142375278917970072361);
    AssertAlmostEqual(&test, lnb, -31.720170045619969451321686623861774853);

    casimir_lnab(&casimir, 15, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 7.27860937372782962645683152471467237745);
    AssertAlmostEqual(&test, lnb, 7.19219783208384282132829359294763960883);

    casimir_lnab(&casimir, 15, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 7.02334016320151547642830308466323532804);
    AssertAlmostEqual(&test, lnb, 6.58397773056963416767295911583683756757);

    casimir_lnab(&casimir, 15, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 6.45115961754533443069229264656379407834);
    AssertAlmostEqual(&test, lnb, 5.69616122563189712145054916633333964902);

    casimir_lnab(&casimir, 15, 4, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 5.58804392430302555818524298325189488817);
    AssertAlmostEqual(&test, lnb, 4.54964202124786251853718659079643076704);

    casimir_lnab(&casimir, 15, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 4.45913124438624191784643966629578158975);
    AssertAlmostEqual(&test, lnb, 3.16500182122169513370075408520449339632);

    casimir_lnab(&casimir, 15, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -4.4496934916915767898295519231813122808);
    AssertAlmostEqual(&test, lnb, -6.7270716450698327878555128172661915061);

    casimir_lnab(&casimir, 15, 20, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -33.497220794437611205903569828087297055);
    AssertAlmostEqual(&test, lnb, -36.973404266432968593645620958704189729);

    casimir_lnab(&casimir, 15, 50, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -167.19982979934473385634246969571044075);
    AssertAlmostEqual(&test, lnb, -172.42020385124701110811733746273093103);

    casimir_free(&casimir);

    return test_results(&test, stderr);
}

int test_integration(void)
{
    casimir_integrals_t cint;
    unittest_t test;
    unittest_init(&test, "Integration", "Test integration for various parameters");

    casimir_integrate_perf(&cint, 3, 2, 1, 2);
    AssertAlmostEqual(&test, cint.lnA_TM, -4.094372316589062);
    AssertAlmostEqual(&test, cint.lnB_TM, -1.970116759119433);
    AssertAlmostEqual(&test, cint.lnC_TM, -3.298725852652321);

    casimir_integrate_perf(&cint, 4, 4, 0, 0.005);
    AssertAlmostEqual(&test, cint.lnB_TM, 56.28387814539346+0.6931471805599453);

    casimir_integrate_perf(&cint, 4, 4, 1, 0.005);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +2.4806179125126554e17*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -2.2226323455151368e24*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -6.9457269656680333e20*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +6.9457269656680333e20*-2);

    casimir_integrate_perf(&cint, 40, 40, 1, 0.25);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +1.5754477603435539e159*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -6.3723632215476122e166*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -9.9568222699306801e162*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +9.9568222699306801e162*-2);

    casimir_integrate_perf(&cint, 40, 40, 40, 1);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +6.4140686579381969e91*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -1.0147301906459434e95*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -2.5352219594503741e93*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +2.5352219594503736e93*-2);

    casimir_integrate_perf(&cint, 7, 4, 3, 8.5);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +4.8180365200137397e-9*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -1.3731640166794149e-8*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -6.7659079909128738e-9*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +9.44463292099617e-9*-2);

    casimir_integrate_perf(&cint, 40, 40, 0, 2.5);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -6.0455421304871757e85*-2);

    casimir_integrate_perf(&cint, 100, 41, 0, 2.5);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), 8.8689390374540308e185*-2);

    return test_results(&test, stderr);
}

int test_integration_drude(void)
{
    casimir_integrals_t cint;
    casimir_t casimir;
    double omegap, gamma_;
    unittest_t test;

    unittest_init(&test, "Integration Drude", "Test integration for various parameters for Drude");

    omegap = 1.32e2;
    gamma_ = 6.9e-1;
    casimir_init(&casimir, 0.5, 1);
    casimir_set_omegap_plane(&casimir, omegap);
    casimir_set_gamma_plane(&casimir, gamma_);

    {
        casimir_integrate_drude(&casimir, &cint, 3, 2, 1, 1);

        AssertAlmostEqual(&test, cint.lnA_TE, -0.62981145199252068602408);
        AssertEqual(&test, cint.signA_TE, +1);
        AssertAlmostEqual(&test, cint.lnA_TM, -0.59589434712666196879313);
        AssertEqual(&test, cint.signA_TM, -1);

        AssertAlmostEqual(&test, cint.lnB_TE, 2.7383266248198112347328);
        AssertEqual(&test, cint.signB_TE, -1);
        AssertAlmostEqual(&test, cint.lnB_TM, 2.7945369735963442460566);
        AssertEqual(&test, cint.signB_TM, +1);

        AssertAlmostEqual(&test, cint.lnC_TE, 0.74158885587407677484842);
        AssertEqual(&test, cint.signC_TE, +1);
        AssertAlmostEqual(&test, cint.lnC_TM, 0.78654850992186874297884);
        AssertEqual(&test, cint.signC_TM, -1);

        AssertAlmostEqual(&test, cint.lnD_TE, 1.1256865694785886272321);
        AssertEqual(&test, cint.signD_TE, -1);
        AssertAlmostEqual(&test, cint.lnD_TM, 1.1711803659444285938208);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 3, 2, 2, 1);

        AssertAlmostEqual(&test, cint.lnA_TE, -0.7132786835392315505014);
        AssertEqual(&test, cint.signA_TE, +1);
        AssertAlmostEqual(&test, cint.lnA_TM, -0.67405454481613061169669);
        AssertEqual(&test, cint.signA_TM, -1);

        AssertAlmostEqual(&test, cint.lnB_TE, 1.5622997116874727152691);
        AssertEqual(&test, cint.signB_TE, -1);
        AssertAlmostEqual(&test, cint.lnB_TM, 1.6191154482067796685377);
        AssertEqual(&test, cint.signB_TM, +1);

        AssertAlmostEqual(&test, cint.lnC_TE, 0.18064782184885164701636);
        AssertEqual(&test, cint.signC_TE, +1);
        AssertAlmostEqual(&test, cint.lnC_TM, 0.22781988160375950931957);
        AssertEqual(&test, cint.signC_TM, -1);

        AssertAlmostEqual(&test, cint.lnD_TE, 0.52072721231985447793985);
        AssertEqual(&test, cint.signD_TE, -1);
        AssertAlmostEqual(&test, cint.lnD_TM, 0.56889321108023160164297);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 3, 2, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, -4.1459191747317624709052);
        AssertEqual(&test, cint.signA_TE, +1);
        AssertAlmostEqual(&test, cint.lnA_TM, -4.1197945729671869295841);
        AssertEqual(&test, cint.signA_TM, -1);

        AssertAlmostEqual(&test, cint.lnB_TE, -2.0359815567492122711267);
        AssertEqual(&test, cint.signB_TE, -1);
        AssertAlmostEqual(&test, cint.lnB_TM, -1.9903481981146134279244);
        AssertEqual(&test, cint.signB_TM, +1);

        AssertAlmostEqual(&test, cint.lnC_TE, -3.356496453224571521287);
        AssertEqual(&test, cint.signC_TE, +1);
        AssertAlmostEqual(&test, cint.lnC_TM, -3.3216912696961575288389);
        AssertEqual(&test, cint.signC_TM, -1);

        AssertAlmostEqual(&test, cint.lnD_TE, -3.0252864555803092122331);
        AssertEqual(&test, cint.signD_TE, -1);
        AssertAlmostEqual(&test, cint.lnD_TM, -2.9891216277980328174213);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 4, 2, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, -3.4410543260111500276103);
        AssertEqual(&test, cint.signA_TE, +1);
        AssertAlmostEqual(&test, cint.lnA_TM, -3.4078967115562010093212);
        AssertEqual(&test, cint.signA_TM, -1);

        AssertAlmostEqual(&test, cint.lnB_TE, -0.74990261955450098567219);
        AssertEqual(&test, cint.signB_TE, -1);
        AssertAlmostEqual(&test, cint.lnB_TM, -0.69545756637556400940693);
        AssertEqual(&test, cint.signB_TM, +1);

        AssertAlmostEqual(&test, cint.lnC_TE, -2.5076062999924391950641);
        AssertEqual(&test, cint.signC_TE, +1);
        AssertAlmostEqual(&test, cint.lnC_TM, -2.4646103122548521433085);
        AssertEqual(&test, cint.signC_TM, -1);

        AssertAlmostEqual(&test, cint.lnD_TE, -1.883875993266507818876);
        AssertEqual(&test, cint.signD_TE, -1);
        AssertAlmostEqual(&test, cint.lnD_TM, -1.8392724297362885494622);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 4, 3, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, -2.6626325469377015011493);
        AssertEqual(&test, cint.signA_TE, -1);
        AssertAlmostEqual(&test, cint.lnA_TM, -2.6217256954918375062335);
        AssertEqual(&test, cint.signA_TM, +1);

        AssertAlmostEqual(&test, cint.lnB_TE, 0.69662250934949833938181);
        AssertEqual(&test, cint.signB_TE, +1);
        AssertAlmostEqual(&test, cint.lnB_TM, 0.76016712598231202383933);
        AssertEqual(&test, cint.signB_TM, -1);

        AssertAlmostEqual(&test, cint.lnC_TE, -1.2220766351698472482591);
        AssertEqual(&test, cint.signC_TE, -1);
        AssertAlmostEqual(&test, cint.lnC_TM, -1.1696039758614060731918);
        AssertEqual(&test, cint.signC_TM, +1);

        AssertAlmostEqual(&test, cint.lnD_TE, -0.94454349314724581653695);
        AssertEqual(&test, cint.signD_TE, +1);
        AssertAlmostEqual(&test, cint.lnD_TM, -0.89174332617421184874047);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 4, 3, 1, 20);

        AssertAlmostEqual(&test, cint.lnA_TE, -43.137947441949791493356);
        AssertAlmostEqual(&test, cint.lnA_TM, -43.120372668298594245786);
        AssertAlmostEqual(&test, cint.lnB_TE, -42.712836955286222313897);
        AssertAlmostEqual(&test, cint.lnB_TM, -42.688853217438179634812);
        AssertAlmostEqual(&test, cint.lnC_TE, -42.981679865536465578416);
        AssertAlmostEqual(&test, cint.lnC_TM, -42.961747317757136470438);
        AssertAlmostEqual(&test, cint.lnD_TE, -42.896313862496345284576);
        AssertAlmostEqual(&test, cint.lnD_TM, -42.875331151064598797301);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 4, 3, 1, 0.01);

        AssertAlmostEqual(&test, cint.lnA_TE, 29.361786303876121509058);
        AssertAlmostEqual(&test, cint.lnA_TM, 29.727644399301097042065);
        AssertAlmostEqual(&test, cint.lnB_TE, 43.288690290377248702015);
        AssertAlmostEqual(&test, cint.lnB_TM, 43.774264506257246488789);
        AssertAlmostEqual(&test, cint.lnC_TE, 36.104129446007204258962);
        AssertAlmostEqual(&test, cint.lnC_TM, 36.53003707484319590354);
        AssertAlmostEqual(&test, cint.lnD_TE, 36.39181148264103954446);
        AssertAlmostEqual(&test, cint.lnD_TM, 36.817719115542167950645);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 20, 20, 0, 2);

        AssertAlmostEqual(&test, cint.lnB_TE, 80.261829579383622339087);
        AssertAlmostEqual(&test, cint.lnB_TM, 80.616659994373914035408);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 20, 20, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, 69.659355648443184396324);
        AssertAlmostEqual(&test, cint.lnA_TM, 69.996776884623239716192);
        AssertAlmostEqual(&test, cint.lnB_TE, 80.212799262982493059119);
        AssertAlmostEqual(&test, cint.lnB_TM, 80.567633810049211725342);
        AssertAlmostEqual(&test, cint.lnC_TE, 74.92332975118821662198);
        AssertAlmostEqual(&test, cint.lnC_TM, 75.269464691472079467208);
        AssertAlmostEqual(&test, cint.lnD_TE, 74.92332975118821662198);
        AssertAlmostEqual(&test, cint.lnD_TM, 75.269464691472079467208);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 7, 7, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, 6.6971469912051709882689);
        AssertAlmostEqual(&test, cint.lnA_TM, 6.8044855533996484279007);
        AssertAlmostEqual(&test, cint.lnB_TE, 12.986970305176281775914);
        AssertAlmostEqual(&test, cint.lnB_TM, 13.113563515622480422875);
        AssertAlmostEqual(&test, cint.lnC_TE, 9.8017901880975337522829);
        AssertAlmostEqual(&test, cint.lnC_TM, 9.9188754798338178698508);
        AssertAlmostEqual(&test, cint.lnD_TE, 9.8017901880975337522829);
        AssertAlmostEqual(&test, cint.lnD_TM, 9.9188754798338178698508);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 60, 7, 1, 4);

        AssertAlmostEqual(&test, cint.lnA_TE, 111.142136572991446682);
        AssertAlmostEqual(&test, cint.lnA_TM, 111.67030961586370531194);
        AssertAlmostEqual(&test, cint.lnB_TE, 121.40133118469416742774);
        AssertAlmostEqual(&test, cint.lnB_TM, 121.94548929011007049376);
        AssertAlmostEqual(&test, cint.lnC_TE, 115.18974866771901373092);
        AssertAlmostEqual(&test, cint.lnC_TM, 115.72592401117753046786);
        AssertAlmostEqual(&test, cint.lnD_TE, 117.33853608403136877599);
        AssertAlmostEqual(&test, cint.lnD_TM, 117.87470585363882880168);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 60, 50, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, 317.54496357191706674104);
        AssertAlmostEqual(&test, cint.lnA_TM, 317.54496357191706674104);
        AssertAlmostEqual(&test, cint.lnB_TE, 332.17054025252939692152);
        AssertAlmostEqual(&test, cint.lnB_TM, 332.17054025252939692152);
        AssertAlmostEqual(&test, cint.lnC_TE, 324.76202098855843751046);
        AssertAlmostEqual(&test, cint.lnC_TM, 324.76202098855843751046);
        AssertAlmostEqual(&test, cint.lnD_TE, 324.94434361445935795232);
        AssertAlmostEqual(&test, cint.lnD_TM, 324.94434361445935795232);
    }

    casimir_free(&casimir);

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
    AssertAlmostEqual(&test, plm_Plm(2,1,3), -25.455844122715710878430397035774565414);
    AssertAlmostEqual(&test, plm_Plm(2,2,3), -24);

    AssertAlmostEqual(&test, plm_Plm(3,0,2), 17);
    AssertAlmostEqual(&test, plm_Plm(3,1,2), -49.363448015713002865532220732917362457);
    AssertAlmostEqual(&test, plm_Plm(3,2,2), -90);
    AssertAlmostEqual(&test, plm_Plm(3,3,2), 77.9422863405994782087350853677642565124);


    AssertAlmostEqual(&test, plm_Plm(7,0,2), 2199.125);
    AssertAlmostEqual(&test, plm_Plm(7,1,2), -15209.246394437784569295351790159409371);
    AssertAlmostEqual(&test, plm_Plm(7,2,2), -88026.75);
    AssertAlmostEqual(&test, plm_Plm(7,3,2), 414721.162832537248613903297356202638370);
    AssertAlmostEqual(&test, plm_Plm(7,4,2), 1528065.0);
    AssertAlmostEqual(&test, plm_Plm(7,5,2), -4132071.3392037110374969861832293381768);
    AssertAlmostEqual(&test, plm_Plm(7,6,2), -7297290.0);
    AssertAlmostEqual(&test, plm_Plm(7,7,2), 6319638.51878214629264244945670369368228);

    AssertAlmostEqual(&test, plm_Plm(8,0,2), 7691.1484375);
    AssertAlmostEqual(&test, plm_Plm(8,1,2), -60890.462646434827363619067066431631293);
    AssertAlmostEqual(&test, plm_Plm(8,2,2), -413142.1875);
    AssertAlmostEqual(&test, plm_Plm(8,3,2), 2354110.35991671119020797796713392578050);
    AssertAlmostEqual(&test, plm_Plm(8,4,2), 1.0957629375e7);
    AssertAlmostEqual(&test, plm_Plm(8,5,2), -40024377.285620259853402179892456726654);
    AssertAlmostEqual(&test, plm_Plm(8,6,2), -1.076350275e8);
    AssertAlmostEqual(&test, plm_Plm(8,7,2), 189589155.563464388779273483701110810468);
    AssertAlmostEqual(&test, plm_Plm(8,8,2), 1.64189025e8);

    AssertAlmostEqual(&test, plm_Plm(3,2,200), -119997000);

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

int test_epsilon()
{
    double omegap, gamma_;
    unittest_t test;

    unittest_init(&test, "Epsilon", "Test dielectric function");

    /* parameters of gold */
    omegap = 1.32e16;
    gamma_ = 6.9e13;

    AssertAlmostEqual(&test, casimir_epsilon(1e-10, omegap, gamma_), 2.5252173913043478e28);
    AssertAlmostEqual(&test, casimir_epsilon(1e-9,  omegap, gamma_), 2.5252173913043477e27);
    AssertAlmostEqual(&test, casimir_epsilon(1e-8,  omegap, gamma_), 2.5252173913043479e26);
    AssertAlmostEqual(&test, casimir_epsilon(1e-7,  omegap, gamma_), 2.5252173913043479e25);
    AssertAlmostEqual(&test, casimir_epsilon(1e-6,  omegap, gamma_), 2.5252173913043478e24);
    AssertAlmostEqual(&test, casimir_epsilon(1e-5,  omegap, gamma_), 2.5252173913043478e23);
    AssertAlmostEqual(&test, casimir_epsilon(1e-4,  omegap, gamma_), 2.525217391304348e22);
    AssertAlmostEqual(&test, casimir_epsilon(1e-3,  omegap, gamma_), 2.525217391304348e21);
    AssertAlmostEqual(&test, casimir_epsilon(1e-2,  omegap, gamma_), 2.5252173913043475e20);
    AssertAlmostEqual(&test, casimir_epsilon(1e-1,  omegap, gamma_), 2.525217391304344e19);
    AssertAlmostEqual(&test, casimir_epsilon(1e+0,  omegap, gamma_), 2.5252173913043113e18);
    AssertAlmostEqual(&test, casimir_epsilon(1e+1,  omegap, gamma_), 2.5252173913039818e17);
    AssertAlmostEqual(&test, casimir_epsilon(1e+2,  omegap, gamma_), 2.525217391300688e16);
    AssertAlmostEqual(&test, casimir_epsilon(1e+3,  omegap, gamma_), 2.5252173912677515e15);
    AssertAlmostEqual(&test, casimir_epsilon(1e+4,  omegap, gamma_), 2.5252173909383844e14);
    AssertAlmostEqual(&test, casimir_epsilon(1e+5,  omegap, gamma_), 2.5252173876447125e13);
    AssertAlmostEqual(&test, casimir_epsilon(1e+6,  omegap, gamma_), 2.5252173547079951e12);
    AssertAlmostEqual(&test, casimir_epsilon(1e+7,  omegap, gamma_), 2.5252170253408661e11);
    AssertAlmostEqual(&test, casimir_epsilon(1e+8,  omegap, gamma_), 2.5252137316743023e10);
    AssertAlmostEqual(&test, casimir_epsilon(1e+9,  omegap, gamma_), 2.5251807954812393e9);
    AssertAlmostEqual(&test, casimir_epsilon(1e+10, omegap, gamma_), 2.5248514808013332e8);
    AssertAlmostEqual(&test, casimir_epsilon(1e+11, omegap, gamma_), 2.5215630522431258e7);
    AssertAlmostEqual(&test, casimir_epsilon(1e+12, omegap, gamma_), 2489143.857142857);
    AssertAlmostEqual(&test, casimir_epsilon(1e+13, omegap, gamma_), 220557.9620253165);
    AssertAlmostEqual(&test, casimir_epsilon(1e+14, omegap, gamma_), 10311.05917159763);
    AssertAlmostEqual(&test, casimir_epsilon(1e+15, omegap, gamma_), 163.9934518241347);
    AssertAlmostEqual(&test, casimir_epsilon(1e+16, omegap, gamma_), 2.730459827192373);
    AssertAlmostEqual(&test, casimir_epsilon(1e+17, omegap, gamma_), 1.017411985729847);
    AssertAlmostEqual(&test, casimir_epsilon(1e+18, omegap, gamma_), 1.00017422797827);

    AssertAlmostEqual(&test, casimir_lnepsilon(1e-10, omegap, gamma_), 65.39870975842067);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e-8,  omegap, gamma_), 60.79353957243258);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e-6,  omegap, gamma_), 56.18836938644449);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e-2,  omegap, gamma_), 46.97802901446831);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+0,  omegap, gamma_), 42.3728588284802);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+2,  omegap, gamma_), 37.76768864249068);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+4,  omegap, gamma_), 33.16251845635911);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+6,  omegap, gamma_), 28.55734825602358);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+8,  omegap, gamma_), 23.95217663529314);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+10, omegap, gamma_), 19.34686298546513);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+12, omegap, gamma_), 14.7274493768442);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+14, omegap, gamma_), 9.240972304188087);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+16, omegap, gamma_), 1.00447002988231);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+18, omegap, gamma_), 1.7421280233797525e-4);

    return test_results(&test, stderr);
}

int test_fresnel()
{
    //const double c = 299792458;
    edouble r_TE, r_TM, T;
    double omegap, gamma_;
    unittest_t test;
    casimir_t casimir;

    unittest_init(&test, "Fresnel", "Test Fresnel coefficients");

    T = 1;
    omegap = 1.32e2;
    gamma_ = 6.9e-1;
    casimir_init(&casimir, 0.5, T);

    casimir_set_omegap_plane(&casimir, omegap);
    casimir_set_gamma_plane(&casimir,  gamma_);

    casimir_rp(&casimir, 1*T, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.97252954726278);
    AssertAlmostEqual(&test, r_TM, +0.98616846109802);

    casimir_rp(&casimir, 10*T, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.85446954163643);
    AssertAlmostEqual(&test, r_TM, +0.85579839473205);

    casimir_rp(&casimir, 100*T, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.24595396364878);
    AssertAlmostEqual(&test, r_TM, +0.24598373253191);

    casimir_free(&casimir);

    return test_results(&test, stderr);
}

int main(int argc, char *argv[])
{
    test_Lambda();
    test_mie_drude();

    test_fresnel();
    test_integration_drude();
    test_plm();
    test_doublefact();
    test_epsilon();
    test_Xi();
    test_integration();
    test_mie();
    test_besselI();
    test_besselK();
    test_givens();
    test_logdet();
    test_casimirF();
    
    return 0;
}

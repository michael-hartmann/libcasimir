#include <stdio.h>
#include <math.h>
#include "unittest.h"

void unittest_init(unittest_t *test, const char *func, const char *desc)
{
    test->passed = test->failed = 0;
    test->func = func;
    test->desc = desc;
}

int test_results(unittest_t *test, FILE *stream)
{
    fprintf(stream, "[%d/%d] %s (%s)", test->passed-test->failed, test->passed, test->func, test->desc);
    if(test->failed == 0)
        fprintf(stream, " PASSED\n");
    else
        fprintf(stream, " FAILED\n");

    return test->failed;
}

int _AssertEqual(int line, unittest_t *test, int x, int y)
{
    if(x == y)
    {
        test->passed++;
        return 0;
    }
    else
    {
        test->failed++;
        fprintf(stderr, "FAILED: %d != %d on line %d\n", x, y, line);
        return 1;
    }
}

int _AssertAlmostEqual(int line, unittest_t *test, double x, double y)
{
    if(fabs(1-x/y) < EPS)
    {
        test->passed++;
        return 0;
    }
    else
    {
        test->failed++;
        fprintf(stderr, "FAILED: %g != %g on line %d\n", x, y, line);
        return 1;
    }
}

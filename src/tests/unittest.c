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
    fprintf(stream, "[%3d/%3d]\t%-20s\t%-50s", test->passed, test->passed+test->failed, test->func, test->desc);
    if(test->failed == 0)
        fprintf(stream, " [PASSED]\n");
    else
        fprintf(stream, " [FAILED]\n");

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

int _Assert(int line, unittest_t *test, int boolean)
{
    if(boolean)
    {
        test->passed++;
        return 0;
    }
    else
    {
        test->failed++;
        fprintf(stderr, "FAILED: on line %d\n", line);
        return 1;
    }
}

int _AssertAlmostEqual(int line, unittest_t *test, double x, double y, double eps)
{
    if(fabs(1-x/y) < eps)
    {
        test->passed++;
        return 0;
    }
    else
    {
        test->failed++;
        fprintf(stderr, "FAILED: %.20g != %.20g (%g) on line %d\n", x, y, fabs(1-x/y), line);
        return 1;
    }
}

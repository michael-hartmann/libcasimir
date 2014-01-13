#ifndef UNITTESTS__
#define UNITTESTS__

#define EPS 1e-9
#define AssertAlmostEqual(t,x,y) _AssertAlmostEqual(__LINE__, t, x, y)

typedef struct {
    int passed;
    int failed;
    const char *func;
    const char *desc;
} unittest_t;

void unittest_init(unittest_t *test, const char *func, const char *desc);
int test_results(unittest_t *test, FILE *stream);
int _AssertAlmostEqual(int line, unittest_t *test, double x, double y);

#endif

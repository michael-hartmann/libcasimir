#ifndef UNITTESTS__
#define UNITTESTS__

#define EPS 1e-10
#define Assert(t,boolean) _Assert(__LINE__, t, boolean)
#define AssertEqual(t,x,y) _AssertEqual(__LINE__, t, x, y)
#define AssertAlmostEqual(t,x,y) _AssertAlmostEqual(__LINE__, t, x, y)

typedef struct {
    int passed;
    int failed;
    const char *func;
    const char *desc;
} unittest_t;

void unittest_init(unittest_t *test, const char *func, const char *desc);
int test_results(unittest_t *test, FILE *stream);
int _Assert(int line, unittest_t *test, int boolean);
int _AssertAlmostEqual(int line, unittest_t *test, double x, double y);
int _AssertEqual(int line, unittest_t *test, int x, int y);

#endif

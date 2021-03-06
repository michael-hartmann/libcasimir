#ifndef __UTILS_H
#define __UTILS_H

double now(void);

void set_defaut_error_handler(void (*f)(const char *));
void default_error_handler(const char *str);

void *xmalloc(size_t len);
void *xrealloc(void *p, size_t len);
void xfree(void *p);

int cinstr(const char *str, char c);
const char *indexn(const char *str, char c, int n);

void swap(double *a, double *b);

#endif

.PHONY: casimir tests
all: casimir tests

CC=gcc
CFLAGS=-Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -fshort-enums -fno-common -g
LDFLAGS=-lm -lquadmath
OPTIMIZATION=3
SOURCES=sfunc.c integration.c libcasimir.c givens.c

casimir:
	$(CC) -O$(OPTIMIZATION) $(CFLAGS) $(SOURCES) casimir.c $(LDFLAGS) -o casimir

tests:
	$(CC) -O$(OPTIMIZATION) $(CFLAGS) $(SOURCES) unittest.c tests.c $(LDFLAGS) -o casimir

clean:
	rm -f casimir tests

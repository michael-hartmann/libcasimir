.PHONY: casimir tests
all: casimir tests

CC=gcc
CFLAGS=-Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -fshort-enums -fno-common -g
LDFLAGS=-lm -lquadmath
OPTIMIZATION=3
SOURCES=sfunc.c integration.c libcasimir.c givens.c

casimir:
	$(CC) -O$(OPTIMIZATION) $(SOURCES) $(CFLAGS) casimir.c $(LDFLAGS) -o casimir

tests:
	$(CC) -O$(OPTIMIZATION) $(SOURCES) $(CFLAGS) unittest.c tests.c $(LDFLAGS) -o tests

clean:
	rm -f casimir tests

.PHONY: casimir tests
all: casimir tests

CC=gcc
CFLAGS=-Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -fshort-enums -fno-common
ifeq ($(CC),gcc)
    LDFLAGS=-lquadmath
endif
LDFLAGS+= -lm
OPTIMIZATION=3
SOURCES=sfunc.c integration.c libcasimir.c givens.c

casimir:
	$(CC) -O$(OPTIMIZATION) $(SOURCES) $(CFLAGS) casimir.c $(LDFLAGS) -o casimir

n0:
	$(CC) -O$(OPTIMIZATION) $(SOURCES) $(CFLAGS) n0.c $(LDFLAGS) -o n0

tests:
	$(CC) -O$(OPTIMIZATION) $(SOURCES) $(CFLAGS) unittest.c tests.c $(LDFLAGS) -o tests

clean:
	rm -f casimir tests

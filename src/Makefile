.PHONY: casimir tests
all: casimir tests

CC = gcc
OPTIMIZATION = 3
SOURCES = sfunc.c integration.c libcasimir.c givens.c
CFLAGS += -Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -fshort-enums -fno-common
LDFLAGS += -lm
ifeq ($(CC),gcc)
    LDFLAGS += -lquadmath
endif
ifeq ($(CC),icc)
    LDFLAGS += -Qoption,c,--extended_float_type
endif


casimir:
	$(CC) -O$(OPTIMIZATION) $(SOURCES) $(CFLAGS) casimir.c $(LDFLAGS) -o casimir

n0:
	$(CC) -O$(OPTIMIZATION) $(SOURCES) $(CFLAGS) n0.c $(LDFLAGS) -o n0

tests:
	$(CC) -O$(OPTIMIZATION) $(SOURCES) $(CFLAGS) unittest.c tests.c $(LDFLAGS) -o tests

clean:
	rm -f casimir tests

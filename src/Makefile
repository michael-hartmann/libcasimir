.PHONY: casimir tests n0
all: casimir

CC = gcc
SOURCES = sfunc.c integration.c libcasimir.c matrix.c
CFLAGS += -Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -fshort-enums -fno-common
#CFLAGS += -DNDEBUG
#CFLAGS += -DMATRIX_QUAD
#CFLAGS += -DINTEGRATION_QUAD
LDFLAGS += -lm -llapacke
ifeq ($(CC),gcc)
    OPTIMIZATION = -O3
    LDFLAGS += -lquadmath
endif
ifeq ($(CC),icc)
    OPTIMIZATION = -fast
    LDFLAGS += -Qoption,c,--extended_float_type
endif


casimir:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) casimir.c $(LDFLAGS) -o casimir

n0:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) n0.c $(LDFLAGS) -o casimir_n0

tests:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) unittest.c tests.c $(LDFLAGS) -o casimir_tests

clean:
	rm -f casimir casimir_tests casimir_n0

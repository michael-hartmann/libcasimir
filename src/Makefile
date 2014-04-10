.PHONY: casimir casimir_tests casimir_hiT casimir_logdetD
all: casimir

EXTENDED_DOUBLE = 1

CC = gcc
SOURCES = sfunc.c integration.c libcasimir.c matrix.c
CFLAGS += -Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -fshort-enums -fno-common
#CFLAGS += -DNDEBUG
CFLAGS += -DMATRIX_QUAD
CFLAGS += -DINTEGRATION_QUAD
ifeq ($(CC),gcc)
    OPTIMIZATION = -O3
endif
ifeq ($(CC),icc)
    OPTIMIZATION = -fast
endif

ifdef EXTENDED_DOUBLE
    CFLAGS += -DEXTENDED_DOUBLE
else
    ifeq ($(CC),gcc)
        LDFLAGS += -lquadmath
    endif
    ifeq ($(CC),icc)
        LDFLAGS += -Qoption,c,--extended_float_type
    endif
endif

LDFLAGS += -lm


casimir:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) casimir.c $(LDFLAGS) -o casimir

casimir_logdetD:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) casimir_logdetD.c $(LDFLAGS) -o casimir_logdetD

casimir_hiT:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) casimir_hiT.c $(LDFLAGS) -pthread -o casimir_hiT

casimir_tests:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) unittest.c tests.c $(LDFLAGS) -o casimir_tests

clean:
	rm -f casimir casimir_tests casimir_hiT casimir_logdetD

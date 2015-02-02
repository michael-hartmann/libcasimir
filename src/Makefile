.PHONY: casimir casimir_tests casimir_hiT casimir_logdetD
all: casimir casimir_hiT casimir_logdetD casimir_tests

OBJECTS=$(SOURCES:.c=.o)

# Use extended doubles for calculations. gcc and icc implement long double as
# 80-bit extended precision type by x86 hardware This is less precise than
# quadrupule precision, but should be much faster.
EXTENDED_DOUBLE = 1

# use fftw3
FFTW3 = 1

# compile. You may use gcc or icc.
CC = gcc
#CC = icc

SOURCES = sfunc.c integration.c libcasimir.c matrix.c utils.c tests/unittest.c
CFLAGS += -Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -fshort-enums -fno-common

# -DNDEBUG will disable asserts. This hardly increases performance and in case
# of errors the program won't abort.
#CFLAGS += -DNDEBUG

LDFLAGS += -lm -pthread

# use optimization
ifeq ($(CC),gcc)
    OPTIMIZATION = -O3
endif
ifeq ($(CC),icc)
    OPTIMIZATION = -fast
endif

ifdef FFTW3
    CFLAGS += -DFFT_POLYMULT
    LDFLAGS += -lfftw3
endif

# if we need quadrupole precision
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



casimir: $(OBJECTS)
	$(CC) $(OPTIMIZATION) $(CFLAGS) $(OBJECTS) casimir.c $(LDFLAGS) -o casimir

casimir_logdetD: $(OBJECTS)
	$(CC) $(OPTIMIZATION) $(CFLAGS) $(OBJECTS) casimir_logdetD.c $(LDFLAGS) -o casimir_logdetD

casimir_hiT: $(OBJECTS)
	$(CC) $(OPTIMIZATION) $(CFLAGS) $(OBJECTS) casimir_hiT.c $(LDFLAGS) -o casimir_hiT

casimir_tests: $(OBJECTS)
	$(CC) $(OPTIMIZATION) -I . $(CFLAGS) $(OBJECTS) tests/tests.c $(LDFLAGS) -o casimir_tests
	./casimir_tests

clean:
	rm -f casimir casimir_tests casimir_hiT casimir_logdetD *.o

doc:
	doxygen doxygen.conf

.PHONY: clean doc

.PHONY: casimir casimir_tests casimir_hiT casimir_logdetD
all: casimir casimir_hiT casimir_logdetD casimir_tests

# Use extended doubles for calculations. gcc and icc implement long double as
# 80-bit extended precision type by x86 hardware This is less precise than
# quadrupule precision, but should be much faster.
EXTENDED_DOUBLE = 1

# compile. You may use gcc or icc.
CC = gcc
#CC = icc

SOURCES = sfunc.c integration.c libcasimir.c matrix.c utils.c
CFLAGS += -Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -fshort-enums -fno-common

# -DNDEBUG will disable asserts. This hardly increases performance and in case
# of errors the program won't abort.
#CFLAGS += -DNDEBUG

# use optimization
ifeq ($(CC),gcc)
    OPTIMIZATION = -O3
endif
ifeq ($(CC),icc)
    OPTIMIZATION = -fast
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

LDFLAGS += -lm -pthread


casimir:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) casimir.c $(LDFLAGS) -o casimir

casimir_logdetD:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) casimir_logdetD.c $(LDFLAGS) -o casimir_logdetD

casimir_hiT:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) casimir_hiT.c $(LDFLAGS) -o casimir_hiT

casimir_tests:
	$(CC) $(OPTIMIZATION) $(SOURCES) $(CFLAGS) unittest.c tests.c $(LDFLAGS) -o casimir_tests
	./casimir_tests

clean:
	rm -f casimir casimir_tests casimir_hiT casimir_logdetD

doc:
	doxygen doxygen.conf

all: gcc

gcc:
	gcc -O3 -Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -fshort-enums -fno-common -g la.c logdet1m.c casimir.c integration.c libcasimir.c `pkg-config --cflags --libs gsl` -lm -o casimir

tests:
	gcc -O3 -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -fshort-enums -fno-common -g la.c logdet1m.c integration.c libcasimir.c unittest.c tests.c `pkg-config --cflags --libs gsl` -lm -o tests

clean:
	rm -f casimir

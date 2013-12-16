all: gcc

icc:
	icc -O3 -Wall la.c logdet1m.c gausslaguerre.c casimir.c plm_fast.c libcasimir.c `pkg-config --cflags --libs gsl` -lm -o casimir

gcc:
	gcc -O3 -Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -fshort-enums -fno-common -g la.c logdet1m.c gausslaguerre.c casimir.c plm_fast.c libcasimir.c `pkg-config --cflags --libs gsl` -lm -o casimir

clean:
	rm -f casimir

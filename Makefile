all: plantri stcq

clean:

plantri: plantri.c
	cc -o plantri -O4 plantri.c

stcq: stcq_sa.c
	cc -o stcq -O4 stcq_sa.c liblpsolve55.a -lm -ldl
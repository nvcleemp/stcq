all: plantri stcq_concave

clean:

plantri: plantri.c
	cc -o plantri -O4 plantri.c

stcq_concave: stcq_concave.c
	cc -o stcq_concave -O4 stcq_concave.c liblpsolve55.a -lm -ldl
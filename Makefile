all: plantri stcq

old: plantri_fkt plantri_stcq

clean:

plantri_fkt: plantri.c fkt.c
	cc -o plantri_fkt -O4 '-DPLUGIN="fkt.c"' plantri.c -lm
	
plantri_stcq: plantri.c stcq.c
	cc -o plantri_stcq -O4 '-DPLUGIN="stcq.c"' plantri.c liblpsolve55.a
	
plantri_stcq_debug: plantri.c stcq.c
	cc -o plantri_stcq -O4 '-DPLUGIN="stcq.c"' -D_DEBUG plantri.c liblpsolve55.a

plantri: plantri.c
	cc -o plantri -O4 plantri.c

stcq: stcq_sa.c
	cc -o stcq -O4 stcq_sa.c liblpsolve55.a -lm -ldl
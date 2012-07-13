all: plantri_fkt plantri_stcq

clean:

plantri_fkt: plantri.c fkt.c
	cc -o plantri_fkt -O4 '-DPLUGIN="fkt.c"' plantri.c
	
plantri_stcq: plantri.c stcq.c
	cc -o plantri_stcq -O4 '-DPLUGIN="stcq.c"' plantri.c

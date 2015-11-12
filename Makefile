all: build/plantri build/stcq_concave

clean:
	rm -rf build

build/plantri: plantri.c
	mkdir -p build
	cc -o build/plantri -O4 plantri.c

build/stcq_concave: stcq_concave.c
	mkdir -p build
	cc -o build/stcq_concave -O4 stcq_concave.c liblpsolve55.a -lm -ldl
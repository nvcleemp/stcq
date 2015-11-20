all: build/plantri build/stcq_RTC

clean:
	rm -rf build

build/plantri: plantri.c
	mkdir -p build
	cc -o build/plantri -O4 plantri.c

build/stcq_RTC: stcq_concave.c
	mkdir -p build
	cc -o build/stcq_RTC -O4 stcq_concave.c liblpsolve55.a -lm -ldl
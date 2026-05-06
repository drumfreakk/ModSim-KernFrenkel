level := g

filename := kernfrenkel

run: out
	./build/out

out: $(filename).o
	gcc -Wall -g3 -O$(level) -lm -I../inc/   ./build/$(filename).o -o ./build/out

$(filename).o: $(filename).c
	mkdir -p build
	mkdir -p out
	gcc -Wall -g3 -O$(level) -lm -I../inc/ -c $(filename).c -o ./build/$(filename).o

dev: clean debug_compile

prod: clean build_compile

debug_compile: run

build_compile: level := 3
build_compile: run

clean:
	rm -f build/*



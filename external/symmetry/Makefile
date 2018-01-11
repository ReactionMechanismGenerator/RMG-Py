target: symmetry

symmetry: symmetry.c
	cc -o symmetry -O3 -ansi -Wall symmetry.c -lm

install: symmetry
	cp symmetry ../../bin/symmetry

test: symmetry check
	./check

clean:
	rm symmetry

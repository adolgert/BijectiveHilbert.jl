skilling: skillingtry.c
	gcc skillingtry.c -o skilling

lawder: lawdercorrected.c
	gcc lawdercorrected.c -o lawder

libskilling.so: skilling.o
	gcc -shared skilling.o -o libskilling.so

skilling.o: skilling.c
	gcc -c -fPIC skilling.c -o skilling.o


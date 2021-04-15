kepler.so: kepler.o 
	gcc -shared -o kepler.so kepler.o

kepler.o: kepler.c
	gcc -c -Wall -Werror -fpic kepler.c

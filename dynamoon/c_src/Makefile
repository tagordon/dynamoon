all: kepler.so phot.so

kepler.so: kepler.o 
	gcc -shared -o kepler.so kepler.o

phot.so: phot.o
	gcc -shared -o phot.so phot.o

kepler.o: kepler.c
	gcc -c -Wall -Werror -fpic kepler.c

phot.o: phot.c
	gcc -c -Wall -Werror -fpic phot.c

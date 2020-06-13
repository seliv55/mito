all: mod.a
.cpp.o:
	g++ -c -O3 -I../include -I../integrT $*.cpp

OBJECTS = lab.o solve.o parray.o distr.o coord.o

mod.a: $(OBJECTS)
	rm -f *.a
	ar r mod.a *.o

lab.o: lab.cpp

solve.o: solve.cpp

parray.o: parray.cpp

distr.o: distr.cpp

coord.o: coord.cpp

clean:
	rm -f *.o mod.a

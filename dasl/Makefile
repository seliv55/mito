all: dasl.a
.f.o:
	gfortran -c -O3 $*.f
#	f77 -c -O3 $*.f
	
OBJECTS = dlinpk.o daux.o ddassl.o peter.o

dasl.a: $(OBJECTS)
	rm -f *.a
	ar r dasl.a *.o

dlinpk.o: dlinpk.f

daux.o: daux.f 

ddassl.o: ddassl.f

peter.o: peter.f

clean:
	rm -f *.o *.a *~

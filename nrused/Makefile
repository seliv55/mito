all: nr.a
.cpp.o: 
	g++ -c -O2 -I../include $*.cpp

OBJECTS = odeint.o pythag.o pzextr.o ludcmp.o lubksb.o rkqs.o rkck.o stifbs.o stiff.o simprsparse.o bsstep.o mmid.o sort2a.o svbksb.o svdcmp.o mprove.o

SPARSE = sprsin.o sprsax.o sprstx.o asolve.o atimes.o snrm.o linbcg.o

nr.a: $(OBJECTS) $(SPARSE)
	rm -f nr.a
	ar r nr.a *.o 

bsstep.o: bsstep.cpp 
 
mmid.o: mmid.cpp

odeint.o: odeint.cpp

pythag.o: pythag.cpp

pzextr.o: pzextr.cpp

ludcmp.o: ludcmp.cpp

lubksb.o: lubksb.cpp

rkqs.o: rkqs.cpp

rkck.o: rkck.cpp

stifbs.o: stifbs.cpp
 
stiff.o: stiff.cpp
 
simprsparse.o: simprsparse.cpp
       
sprsin.o: sprsin.cpp

sprsax.o: sprsax.cpp

sprstx.o: sprstx.cpp

asolve.o: asolve.cpp
      
atimes.o: atimes.cpp
       
snrm.o: snrm.cpp
       
linbcg.o: linbcg.cpp

sort2a.o: sort2a.cpp

svdcmp.o: svdcmp.cpp

svbksb.o: svbksb.cpp

mprove.o: mprove.cpp

clean:	
	rm -f *.o nr.a

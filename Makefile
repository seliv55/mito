all: mitodyn.out

mitodyn.out: main.o con512tpl/mod.a nrused/nr.a integrT/integT.a dasl/dasl.a
	g++ -O3 -m64 main.o con512tpl/mod.a nrused/nr.a integrT/integT.a dasl/dasl.a -lgfortran -o mitodyn.out

main.o:	main.cpp
	g++ -Iinclude -IintegrT -O3 -c main.cpp
	
con512tpl/mod.a:
	make -C con512tpl

nrused/nr.a:
	make -C nrused
	
integrT/integT.a:
	make -C integrT
	
dasl/dasl.a:
	make -C dasl
	
clean:
	rm -f *.out *.o
	make clean -C con512tpl


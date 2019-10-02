
all: compile
	
compile: SistLinear.o Gs.o Fronteira.o pdeSolver.o 
	
	gcc pdeSolver.c Fronteira.c Gs.c  SistLinear.c -lm -o  pdeSolver 

clean:
	rm	-f	$(binary)	*.o

doc:
	doxygen solverDoc

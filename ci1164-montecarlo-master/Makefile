.PHONY: clean dist-clean plot3d
.ONESHELL: plot3d

LIKWID = /home/soft/likwid
LIKWID_FLAGS = -I$(LIKWID)/include -DLIKWID_PERFMON
LIKWID_LIBS = -L$(LIKWID)/lib -llikwid
cFLAGS = $(LIKWID_FLAGS) -Wall -O3 -mavx -march=native
lFLAGS = $(LIKWID_LIBS) -lm

all: utils.o
	gcc -o labMonteCarlo utils.o labMonteCarlo.c -lm

utils.o:
	gcc -c utils.c

debug:
	gcc -g -o labMonteCarlo utils.o labMonteCarlo.c

# gnuplot -c "pde3d.plg" "arq.dat" $${NX} $${NY}
plot3d:
	@ARQDAT=$$(cat arqdat)
	@NX=$$(cat nx) && NY=$$(cat ny)
	@gnuplot -c "pde3d.plg" "$${ARQDAT}" $${NX} $${NY}

clean:
	@echo Removendo arquivos-objeto e arquivos temporários...
	@rm -f *.o 
	@rm -f *~ 

dist-clean: clean
	@echo Removendo arquivos auxiliares...
	@rm -f nx 
	@rm -f ny 
	@ARQDAT=$$(cat arqdat) && rm -f $${ARQDAT} 
	@echo Removendo arquivo dos pontos...
	@rm -f arqdat 

purge: clean dist-clean
	@rm -f labMonteCarlo 

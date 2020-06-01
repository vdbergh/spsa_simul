CC=gcc
CFLAGS=-I. -std=gnu99 -O3 -Wall
LIBS= -lm -lpthread -lgsl -lgslcblas
DEPS = spsa_sim.h
OBJ = brentq_gsl.o elo.o gx2.o lf.o prng.o options.o sos.o spsa.o spsa_sim.o util.o
EXE = spsa_sim

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXE): $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

clean:
	rm -f $(EXE)
	rm -f *.o

.PHONY: clean

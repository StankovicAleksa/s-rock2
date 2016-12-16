export SCOREP_EXPERIMENT_DIRECTORY=scorep_profile
CPP= g++ --std=c++11
FAST= -O3
HPCFLAGS= -g -O3
FLAGSGPROF= -pg
SCOREP= scorep g++
EXEC= srock
EXECP= srock_scorep
EXECF= srock_fast
EXECHPC= srock_hpc
EXECG= srock_gprof
OBJ= obj
RUN= run
HEADERS=./srock2/

all:

build: compile

compile: neuron black-scholes test_1

neuron:  $(OBJ)/srock2.o $(OBJ)/simulation.o $(OBJ)/neuron.o 
	$(CPP)  $(OBJ)/srock2.o $(OBJ)/simulation.o $(OBJ)/neuron.o -o $(RUN)/neuron 

black-scholes: $(OBJ)/srock2.o $(OBJ)/simulation.o $(OBJ)/black-scholes.o 
	$(CPP) $(OBJ)/srock2.o $(OBJ)/simulation.o $(OBJ)/black-scholes.o -o $(RUN)/black-scholes

test_1: $(OBJ)/srock2.o $(OBJ)/simulation.o $(OBJ)/test_1.o 
	$(CPP) $(OBJ)/srock2.o $(OBJ)/simulation.o $(OBJ)/test_1.o -o $(RUN)/test_1

$(OBJ)/srock2.o: integrators/srock2.cpp integrators/srock2.h integrators/parameters.h
	$(CPP) -c integrators/srock2.cpp -o $(OBJ)/srock2.o

$(OBJ)/simulation.o: tests/simulation.h tests/simulation.cpp integrators/parameters.h integrators/srock2.h
	$(CPP) -c tests/simulation.cpp -o $(OBJ)/simulation.o

$(OBJ)/neuron.o: tests/neuron.cpp tests/simulation.h integrators/parameters.h 
	$(CPP) -c tests/neuron.cpp -o $(OBJ)/neuron.o

$(OBJ)/black-scholes.o: tests/black-scholes.cpp tests/simulation.h integrators/parameters.h 
	$(CPP) -c tests/black-scholes.cpp -o $(OBJ)/black-scholes.o

$(OBJ)/test_1.o: tests/test_1.cpp tests/simulation.h integrators/parameters.h 
	$(CPP) -c tests/test_1.cpp -o $(OBJ)/test_1.o


run-neuron: build
	./$(RUN)/neuron

run-black-scholes: build
	./$(RUN)/black-scholes

run-test_1: build
	./$(RUN)/test_1

clean:
	-rm -f $(OBJ)/*.o $(EXEC) $(EXECG) $(EXECF) $(EXECHPC) sol.m time.txt mesh.m
	-rm -rf $(SCOREP_EXPERIMENT_DIRECTORY) $(EXECP) *scorep_init*
	-rm -fr hpctoolkit-$(EXECHPC)-database hpctoolkit-$(EXECHPC)-measurements workspace


#BUILD AND RUN PROFILED PROGRAM
buildp: CPP=$(SCOREP)
buildp: EXEC=$(EXECP)
buildp: clean compile

runp:
	./$(EXECP)

showp:
	cube $(SCOREP_EXPERIMENT_DIRECTORY)/profile.cubex &


#BUILD AND RUN OPTIMIZED PROGRAM
buildf: CPP+=$(FAST)
buildf: EXEC=$(EXECF)
buildf: clean compile

runf:
	./$(EXECF)

#BUILD AND RUN PROFILED WITH HPCTOOLS
buildh: CPP+=$(HPCFLAGS)
buildh: EXEC=$(EXECHPC)
buildh: clean compile


# CPUTIME@5000 count effective running time, without idle
# PAPI_TOT_CYC@4000001 count total cycles
# Its recommended to no used times and PAPIs togheter
# IO measures read/write in bytes
# MEMLEAK finds where the programs do not deallocate memory
runh:
	-rm -fr hpctoolkit-$(EXECHPC)* workspace $(EXECHPC).hpc*
	hpcrun -t -e PAPI_TOT_CYC@4000001 -e CPUTIME@5000 -e MEMLEAK ./$(EXECHPC)
	hpcstruct $(EXECHPC)
	hpcprof -S $(EXECHPC).hpcstruct -I ./'*' hpctoolkit-$(EXECHPC)-measurements
	hpcviewer hpctoolkit-$(EXECHPC)-database


#BUILD AND RUN PROFILED PROGRAM WITH GPROF
buildg: CPP+=$(FLAGSGPROF)
buildg: EXEC=$(EXECG)
buildg: compile

rung:
	./$(EXECG)
	gprof $(EXECG) | tee prof.txt | less

# BUILD AND RUN VALGRIND PROFILED CODE
buildv: build

runv:
	valgrind --tool=memcheck --leak-check=full ./srock

	



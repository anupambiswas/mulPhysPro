SOURCE = ../../sources/flowDomain_2.cpp
HEAD = ../../headers/flowDomain_2.h
MAIN = case_t0.cpp
CPPC = g++

run: $(MAIN) $(SOURCE) $(HEAD)
	$(CPPC) $(MAIN) $(SOURCE) -o run

clean:
	rm -f *~

resetf: clean
	rm -f Channel* Extension* res_* fc* run
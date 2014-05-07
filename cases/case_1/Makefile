SOURCE = ../../sources/flowDomain_2.cpp
HEAD = ../../headers/flowDomain_2.h
MAIN = case_t1.cpp
CPPC = g++

run: $(MAIN) $(HEAD) $(SOURCE)
	$(CPPC) $(MAIN) $(SOURCE) -o run

clean:
	rm -f *~

resetf: clean
	rm -f *~ *u *v *vm F*p R*p G*p res_* run
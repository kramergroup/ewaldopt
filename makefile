F90_FILES := $(wildcard src/*.f90)
OBJ_FILES := $(addprefix obj/,$(notdir $(F90_FILES:.f90=.o)))

F90=ifort

FFLAGS=-g -CB -traceback -module mod

main: ewaldopt

ewaldopt: obj/ewaldsum.o obj/structure.o obj/ewaldopt.o
	$(F90) $(FFLAGS) -o ewaldopt obj/ewaldopt.o obj/ewaldsum.o obj/structure.o

test: obj/ewaldsum.o obj/structure.o obj/ewaldtest.o
	$(F90) $(FFLAGS) -o testewald obj/ewaldsum.o obj/ewaldtest.o obj/structure.o

clean:
	rm -rf obj mod *.mod ewaldopt ewaldtest

obj/ewaldsum.o: src/ewaldsum.f90 obj/structure.o obj mod
	$(F90) $(FFLAGS) -c -o $@ $<


obj/ewaldopt.o: src/ewaldopt.f90 obj/ewaldsum.o obj mod
	$(F90) $(FFLAGS) -c -o $@ $<

obj/structure.o: src/structure.f90 obj mod
	$(F90) $(FFLAGS) -c -o $@ $<

obj/%.o: src/%.f90 obj mod
	$(F90) $(FFLAGS) -c -o $@ $<

obj:
	mkdir obj

mod:
	mkdir mod

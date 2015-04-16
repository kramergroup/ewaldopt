F90_FILES := $(wildcard src/*.f90)
OBJ_FILES := $(addprefix obj/,$(notdir $(F90_FILES:.f90=.o)))

F90=ifort

FFLAGS=-g -CB -traceback

main: obj/ewaldsum.o obj/structure.o obj/ewaldopt.o
	$(F90) $(FFLAGS) -o ewaldopt obj/ewaldopt.o obj/ewaldsum.o obj/structure.o

test: obj/ewaldsum.o obj/structure.o obj/ewaldtest.o
	$(F90) $(FFLAGS) -o testewald obj/ewaldsum.o obj/ewaldtest.o obj/structure.o

obj/%.o: src/%.f90 obj
	$(F90) $(FFLAGS) -c -o $@ $<

obj:
	mkdir obj

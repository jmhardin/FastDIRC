CFLAGS_BASE = -O3 -mtune=native -Wno-comment -Wl,--no-as-needed -L./lib/ `root-config --cflags` `root-config --glibs` 
CFLAGS_OPT = -g
CFLAGS_BASE += $(CFLAGS_OPT)
INCLUDE = -I./include/

CFLAGS = $(CFLAGS_BASE) $(goptical_CPPFLAGS)

LIBLOC = ./lib/
OUT = ./dircfit

OBJFILES = dirc_optical_sim.o
OBJFILES += dirc_digitizer.o
OBJFILES += dirc_probability_spread.o
OBJFILES += dirc_spread_relative.o
OBJFILES += dirc_spread_radius.o
OBJFILES += dirc_spread_linear_soft.o
OBJFILES += dirc_spread_gaussian.o
OBJFILES += dirc_probability_separation.o
OBJFILES += dirc_progressive_separation.o

OBJLOC = $(patsubst %,$(LIBLOC)/%,$(OBJFILES))
LIBFILES = $(LIBLOC)
vpath %.o ./lib/
vpath %.cpp ./source/

%.o : %.cpp
	g++ -Wall $(CFLAGS) $(INCLUDE) -g -o $@ -c $<
	mv $@ $(LIBLOC)

.PHONY : all
all: dircfit.cpp $(OBJFILES)
	g++ -Wall dircfit.cpp $(OBJLOC) $(CFLAGS) $(INCLUDE) -o $(OUT)

.PHONY : clean
clean:
	rm lib/*.o
	rm $(OUT)
	
.PHONY : cleanall
cleanall:
	rm lib/*
	rm *.gcda
	rm $(OUT)
	

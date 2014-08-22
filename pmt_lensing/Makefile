CFLAGS_BASE = -O5 -L./lib/ -lgsl -lgslcblas `root-config --cflags` `root-config --glibs` 
INCLUDE = -I./include/ -I${HOME}/srcinstall/dlib-18.7
goptical_BASEDIR = ${HOME}/srcinstall/goptical-1.0
goptical_LIBDIR = $(goptical_BASEDIR)/goptical_core/src/.libs
goptical_LIBS = $(goptical_BASEDIR)/goptical_core/src/.libs/libgoptical-1.0.so
goptical_HEADERS = -I$(goptical_BASEDIR)/goptical_core/src -I$(goptical_BASEDIR)/goptical_design/src

CFLAGS = $(CFLAGS_BASE) $(goptical_CPPFLAGS)

LIBLOC = ./lib/
OUT = ./check_lensing

OBJLOC = $(patsubst %,$(LIBLOC)/%,$(OBJFILES))
LIBFILES = $(LIBLOC)
vpath %.o ./lib/
vpath %.cpp ./source/

%.o : %.cpp
	g++ -Wall $(CFLAGS) $(goptical_HEADERS) $(INCLUDE) -g -o $@ -c $<
	mv $@ $(LIBLOC)

.PHONY : all
all: check_lensing.cpp $(OBJFILES) 
	g++ -Wall check_lensing.cpp $(OBJLOC) $(goptical_LIBS) -L$(goptical_LIBDIR) $(CFLAGS) $(INCLUDE) -o $(OUT)

.PHONY : clean
clean:
	rm $(OUT)

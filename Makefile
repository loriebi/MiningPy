CC = g++
INC_DIR = -I/usr/include/python2.7 -I/home/nkasradze/ccfits/include -I/usr/include/cfitsio/
LIB_DIR = -L/home/nkasradze/ccfits/lib
LIB = -lCCfits -pthread
MODULE = averager


.PHONY: all
all: swig compile build

swig:
	swig -c++ -python $(MODULE).i

compile:
	g++ -O4 -fPIC -c Averager.h $(MODULE)_wrap.cxx $(INC_DIR) -std=gnu++11
	
build:
	g++ -O4 -shared $(MODULE)_wrap.o -o _$(MODULE).so $(LIB_DIR) $(LIB) -std=gnu++11

.PHONY: clean
clean:
	rm -rf *.o
	rm -rf *.so
	rm -rf $(MODULE)_wrap.cxx
	rm -rf $(MODULE).py
	rm -rf $(MODULE).pyc
	rm -rf Averager.h.gch

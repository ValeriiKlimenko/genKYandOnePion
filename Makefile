ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTINCLUDE  := -I$(shell root-config --incdir)


all: eg_ky

eg_ky:
	$(CXX) -O3 $(ROOTINCLUDE) $(ROOTCFLAGS) -o eg_ky eg_ky.cpp $(ROOTLIBS)

clean:
	rm -rf eg_ky
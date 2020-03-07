ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTINCLUDE  := -I$(shell root-config --incdir)


all: genKYandOnePion

eg_ky:
	$(CXX) -O3 $(ROOTINCLUDE) $(ROOTCFLAGS) -o genKYandOnePion genKYandOnePion.cpp $(ROOTLIBS)

clean:
	rm -rf genKYandOnePion

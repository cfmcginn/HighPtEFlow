CXX = g++
#O3 for max optimization (go to 0 for debug)
#c++17 required for use of <filesystem>
CXXFLAGS = -Wall -Werror -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++17 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

INCLUDE=-I$(PWD)

MKDIR_BIN=mkdir -p bin
MKDIR_OUTPUT=mkdir -p output

ROOT=`root-config --cflags --glibs`
FASTJET=`fastjet-config --cxxflags --libs --plugins`

all: mkdirBin mkdirOutput bin/quickAna.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

bin/quickAna.exe: src/quickAna.C
	$(CXX) $(CXXFLAGS) src/quickAna.C -o bin/quickAna.exe $(INCLUDE) $(ROOT)

clean:
	rm -f ./*~
	rm -f ./#*#
	rm -f bin/*.exe
	rm -rf bin
	rm -f src/*~
	rm -f src/#*#

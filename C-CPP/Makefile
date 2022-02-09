# Project: HH_BBT2020_allP

CPP      = g++
OBJ      = HH_BBT2020_allP.o
LINKOBJ  = HH_BBT2020_allP.o
LIBS     = -L"E:/Pesquisa/precisionPaper/boost_1_72_0/libs"
CXXINCS  = -I"E:/Pesquisa/precisionPaper/boost_1_72_0"
BIN      = HH_BBT2020_allP
# -O0 -> no optimization, -O1, -O2, -O3 -> maximum optimization
CXXOPTS  = -std=c++17 -O3
CXXFLAGS = $(CXXINCS) $(CXXOPTS)
 
RM       = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o $(BIN) $(LIBS)

HH_BBT2020_allP.o: HH_BBT2020_allP.cpp
	$(CPP) -c HH_BBT2020_allP.cpp -o HH_BBT2020_allP.o $(CXXFLAGS)

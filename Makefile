SHAPE_GRADIENT_DOMAIN_TARGET=ShapeGradientDomain
SHAPE_GRADIENT_DOMAIN_SOURCE=ShapeGradientDomain/ShapeGradientDomain.cpp

COMPILER ?= gcc
#COMPILER ?= clang

CFLAGS += -std=c++20 -Wno-deprecated -Wno-invalid-offsetof
LFLAGS += -lstdc++
ifeq ($(COMPILER),gcc)
	CFLAGS += -fopenmp
	LFLAGS += -lgomp -lpthread
	CC=gcc
	CXX=g++
else
	CFLAGS += -Wno-dangling-else -Wno-null-dereference
	CC=clang
	CXX=clang++
	ifeq ($(SANITIZATION),none)
		LFLAGS += -static
	else
		LFLAGS += -g -fsanitize=$(SANITIZATION)
		CFLAGS += -O1 -g -fsanitize=$(SANITIZATION)
	endif
endif

CFLAGS += -O3 -DRELEASE -funroll-loops -g
LFLAGS += -O3 -g

BIN = Bin/Linux/
BIN_O = Obj/Linux/
INCLUDE = /usr/include/ -I.

MD=mkdir

SHAPE_GRADIENT_DOMAIN_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(SHAPE_GRADIENT_DOMAIN_SOURCE))))
SHAPE_GRADIENT_DOMAIN_OBJECT_DIR=$(dir $(SHAPE_GRADIENT_DOMAIN_OBJECTS))

all: make_dirs
all: $(BIN)$(SHAPE_GRADIENT_DOMAIN_TARGET)

shapegradientdomain: make_dirs
shapegradientdomain: $(BIN)$(SHAPE_GRADIENT_DOMAIN_TARGET)

clean:
	rm -rf $(BIN)$(SHAPE_GRADIENT_DOMAIN_TARGET)
	rm -rf $(BIN_O)

make_dirs: FORCE
	$(MD) -p $(BIN)
	$(MD) -p $(BIN_O)
	$(MD) -p $(SHAPE_GRADIENT_DOMAIN_OBJECT_DIR)


$(BIN)$(SHAPE_GRADIENT_DOMAIN_TARGET): $(SHAPE_GRADIENT_DOMAIN_OBJECTS)
	$(CXX) -o $@ $(SHAPE_GRADIENT_DOMAIN_OBJECTS) -L$(BIN) $(LFLAGS)


$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

FORCE:
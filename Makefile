## var

CXX	= g++
CFLAGS	= -std=c++11
LIB	= -lpthread -ldivsufsort64 -ldivsufsort -lsdsl -L$(HOME)/lib
INC	= -I$(HOME)/include
OFLAGS	= -O2 -funroll-loops

## files

OBJ	= obj/main.o  obj/tools.o
	  
TGT	= bin/mphf

## building rule

all:	$(TGT)
	@echo "building done"

## linking

$(TGT):	$(OBJ)
	$(CXX) $(CFLAGS) $(OFLAGS) -o $@ $^ $(LIB) 

## compiling

obj/%.o: src/%.cpp
	$(CXX) $(INC) -c $(CFLAGS) $(OFLAGS) -o $@ $< $(LIB)

## cleaning

clean: $(OBJ) $(TGT)
	rm $^



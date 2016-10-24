#CXX=clang++
CXX=g++
CFLAGS=-Wall -std=c++14

OBJ_DIR=obj
SRC_DIR=src
LIB_DIR=lib
INC_DIR=include
TEST_DIR=tests
BIN_DIR=bin

LIB=libParcours.a

SOURCES=*.cpp
LIB_SOURCES=$(SRC_DIR)/$(SOURCES)
testINC=./tests/
testHelpers=$(testINC)test_helpers.cpp
#buildDIR=./build/
#poa_graph_lib=./build/poa_graph_lib.a

.PHONY: all clean lib

pre: 
	mkdir -v -p $(OBJ_DIR)
	mkdir -v -p $(LIB_DIR)
	mkdir -v -p $(INC_DIR)
	mkdir -v -p $(BIN_DIR)

$(LIB_DIR)/$(LIB): pre
	$(CXX) $(CFLAGS) -I $(SRC_DIR) -c $(LIB_SOURCES) -g
	mv *.o $(OBJ_DIR)/
	ar rc $(LIB) $(OBJ_DIR)/*.o
	ranlib $(LIB)
	mv $(LIB) $(LIB_DIR)/
	cp $(SRC_DIR)/*.h $(INC_DIR)

test: lib
	$(CXX) $(CFLAGS) -I $(INC_DIR) -I $(testINC) -o $(BIN_DIR)/ParcoursLibTests $(testHelpers) $(TEST_DIR)/poa_graph_tests.cpp $(LIB_DIR)/$(LIB)
	$(BIN_DIR)/ParcoursLibTests -d yes

#$(poa_graph_lib): $(LIB_SOURCES)
#	mkdir -v -p $(buildDIR)
#	$(CXX) $(CFLAGS) -I $(INC) -c $(LIB_SOURCES) -g 
#	ar rc poa_graph_lib.a *.o
#	ranlib poa_graph_lib.a
#	rm *.o
#	cp $(INC)*.h $(buildDIR)
#	mv poa_graph_lib.a $(buildDIR)

clean:
	rm -rf $(OBJ_DIR)
	rm -rf $(LIB_DIR)
	rm -rf $(INC_DIR)
	rm -rf $(BIN_DIR)

lib: $(LIB_DIR)/$(LIB)

ifeq ($(shell uname -s),Darwin)
	CXX=g++
else
	CXX=g++
endif
#CXX=clang++
#CXX=g++-5
CFLAGS=-Wall -std=c++11

OBJ_DIR=obj
UNITTEST_OBJ_DIR=$(OBJ_DIR)/unittest
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

.PHONY: all clean lib test

pre: 
	mkdir -v -p $(OBJ_DIR)
	mkdir -v -p $(UNITTEST_OBJ_DIR)
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

$(UNITTEST_OBJ_DIR)/hmm_graph_tests.o: $(TEST_DIR)/HmmGraph_tests.cpp 
	+$(CXX) $(CFLAGS) -c $< -I $(testINC) -I $(INC_DIR) $< 
	mv HmmGraph_tests.o $@

$(UNITTEST_OBJ_DIR)/pairwise_aligner_tests.o: $(TEST_DIR)/PairwiseAligner_tests.cpp 
	+$(CXX) $(CFLAGS) -c $< -I $(testINC) -I $(INC_DIR) $< 
	mv PairwiseAligner_tests.o $@

test: $(LIB_DIR)/$(LIB) $(UNITTEST_OBJ_DIR)/hmm_graph_tests.o $(UNITTEST_OBJ_DIR)/pairwise_aligner_tests.o $(TEST_DIR)/allTests.cpp
	$(CXX) $(CFLAGS) -I $(INC_DIR) -I $(testINC) -o $(BIN_DIR)/ParcoursLibTests $(testHelpers) $^ -L$(LIB_DIR) -lParcours
	#$(BIN_DIR)/ParcoursLibTests -d yes
	$(BIN_DIR)/ParcoursLibTests

clean:
	rm -rf $(OBJ_DIR)
	rm -rf $(UNITTEST_OBJ_DIR)
	rm -rf $(LIB_DIR)
	rm -rf $(INC_DIR)
	rm -rf $(BIN_DIR)

lib: $(LIB_DIR)/$(LIB)

dryrun: 
	make clean
	make test
	make clean

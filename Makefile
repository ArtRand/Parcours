CXX=g++
CFLAGS=-Wall -std=c++11
INC=./src/
LIB_SOURCES=./src/*.cpp
testINC=./tests/
testHelpers=$(testINC)test_helpers.cpp
buildDIR=./build/
poa_graph_lib=./build/poa_graph_lib.a


test: lib
	mkdir -v -p ./bin
	$(CXX) $(CFLAGS) -I $(buildDIR) -I $(testINC) -o ./bin/graphsTests $(testHelpers) ./tests/poa_graph_tests.cpp $(poa_graph_lib)
	./bin/graphsTests -d yes

$(poa_graph_lib): $(LIB_SOURCES)
	mkdir -v -p $(buildDIR)
	$(CXX) $(CFLAGS) -I $(INC) -c $(LIB_SOURCES) -g 
	ar rc poa_graph_lib.a *.o
	ranlib poa_graph_lib.a
	rm *.o
	cp $(INC)*.h $(buildDIR)
	mv poa_graph_lib.a $(buildDIR)

clean:
	rm -rf ./bin 
	rm -rf $(buildDIR)

lib: $(poa_graph_lib)
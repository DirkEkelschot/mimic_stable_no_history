TESTBIN = ../bin

SRC_OBJ = ../../src/*.cpp
TES_OBJ = main.cpp
TEST    = ConvertHex2HybridMesh

include ../../module.mk

test:makebin
	$(CC) $(CXXFLAGS) $(SRC_OBJ) $(TES_OBJ) -o $(TESTBIN)/$(TEST) $(LDFLAGS) $(LDLIBS)

makebin:
	mkdir -p $(TESTBIN)

clean:	
	rm -rf testing

CXX		= c++
CXXFLAGS 	= --std=c++11 -Wall


markovsim: markovsim.cpp state.h
	$(CXX) $(CXXFLAGS) markovsim.cpp -o markovsim

clean:
	/bin/rm -f markovsim

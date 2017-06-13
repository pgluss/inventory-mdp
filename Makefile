CXX		= c++
CXXFLAGS 	= --std=c++11 -Wall


markovsim: ConsoleApplication2.cpp state.h
	$(CXX) $(CXXFLAGS) ConsoleApplication2.cpp -o markovsim

clean:
	/bin/rm -f markovsim

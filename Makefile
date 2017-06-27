CXX		= c++
CXXFLAGS 	= --std=c++11 -Wall


mainsim: main.cpp state.hpp space.hpp space.cpp
	$(CXX) $(CXXFLAGS) main.cpp space.cpp -o mainsim

clean:
	/bin/rm -f mainsim

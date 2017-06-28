#ifndef STATE_H
#define STATE_H

struct State {
  int x;
  int y;
  double f; 

  State(int x, int y, double f) : x(x), y(y), f(f) {}
  State() { }
 
};

#endif

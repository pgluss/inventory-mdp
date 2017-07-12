#ifndef STATE_H
#define STATE_H

struct State {
  unsigned int N;
  int x1;
  int x2;
  int y;
  double f; 

  State(unsigned int N, int x1, int x2, int y, double f) : N(N), x1(x1), x2(x2), y(y), f(f) {}
  State() { }
 
};

#endif

#ifndef STATE_H
#define STATE_H

struct State {
  int x1;
  int x2;
  double f; 

  State(int x1, int x2, double f) {this->x1 = x1; this->x2 = x2; this->f = f;}
  State() { }
 
};

#endif

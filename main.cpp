#include <iostream>

#include "space.hpp"
#include "state.hpp"

int main() {

  // Intialize simulation for n to m
  Space* s = new Space(20);

  // Set parameters
  s->mu = 1;
  s->delta = 1;
  s->lambdax = 0.4;
  s->lambday = 0.5;
  s->hx = 3;
  s->hy = 2;
  s->cx = 80;
  s->cy = 100;
  s->alpha = 0.05;

  // Perform value iteration
  double error = s->vi(1e-5, 1e6);

  std::cout << "Error = " << error << std::endl;
  std::cout << s->state(0,0).f << std::endl;
  return 0;
}

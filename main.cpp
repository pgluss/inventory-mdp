#include <iostream>

#include "space.hpp"
#include "state.hpp"

int main() {

  // Intialize simulation for n to m
  Space* s = new Space(100);

  // Set parameters
  s->mu = 1;
  s->delta = 0.8;
  s->lambdax = 0.85;
  s->lambday = 0.85;
  s->hx = 2;
  s->hy = 1;
  s->cx = 11.4;
  s->cy = 5.7;
  s->alpha = 0.05;

  // Perform value iteration
  double error = s->vi(1e-10, 1e10);

  // Get decisions
  s->decide("test1.txt");

  std::cout << "Error = " << error << std::endl;
  std::cout << s->state(0,0).f << std::endl;
  return 0;
}

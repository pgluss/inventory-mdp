#include <iostream>

#include "space.hpp"
#include "state.hpp"

int main() {

  // Intialize simulation for n to m
  Space* s = new Space(-20, 20);

  // Set parameters
  s->mu1 = 1;
  s->mu2 = 1;
  s->lambda1 = 0.4;
  s->lambda2 = 0.5;
  s->h1 = 3;
  s->h2 = 2;
  s->b1 = 80;
  s->b2 = 100;
  s->alpha = 0.05;

  // Perform value iteration
  double error = s->vi(1e-5, 1e6);

  std::cout << "Error = " << error << std::endl;
  std::cout << s->state(0,0).f << std::endl;
  return 0;
}

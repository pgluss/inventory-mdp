#include <iostream>

#include "space.hpp"

int main() {

  // Intialize simulation for n to m
  Space* s = new Space(-20, 20);

  // Set parameters
  s->mu1 = 1;
  s->mu2 = 1;
  s->lambda1 = 0.4;
  s->lambda2 = 0.5;
  s->h1 = 1;
  s->h2 = 1;
  s->b1 = 30;
  s->b2 = 40;
  s->alpha = 0.01;

  // Perform value iteration
  double error = s->vi(1e-5, 1e6);

  s->decide("haTest");

  std::cout << "Error = " << error << std::endl;
  std::cout << s->f(0,0) << std::endl;
  return 0;
}

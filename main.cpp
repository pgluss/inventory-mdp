#include <iostream>

#include "space.hpp"
#include "state.hpp"

int main() {

  // Intialize simulation for n to m
  Space* s = new Space(-20, 20, 500);

  // Set parameters
  s->mu = 1;
  s->gamma = 0.8;
  s->alpha = 0.05;

  s->lambdax1 = 0.85;
  s->lambdax2 = 0.85;
  s->lambday = 0.85;

  s->hx1 = 2;
  s->hx2 = 2;
  s->hy = 1;

  s->bx1 = 11.4;
  s->bx2 = 11.4;
  s->by = 5.7;

  s->Dx1 = 3;
  s->Dx2 = 3;
  s->Dy  = 3;

  s->pn = 0.3;
  s->po = 0.3;

  // Perform value iteration
  double error = s->vi(1e-10, 1e10);

  // Get decisions
  //s->decide("test1.txt");

  std::cout << "Error = " << error << std::endl;
  std::cout << s->state(0,0,0,5).f << std::endl;
  return 0;
}

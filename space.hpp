#ifndef SPACE_H
#define SPACE_H

#include "state.hpp"
#include <vector>

class Space {
private:
  unsigned int m_;
  std::vector<State> states_;

public:
  double lambdax;
  double lambday;
  double mu;
  double delta;
  double cx;
  double cy;
  double hx;
  double hy;
  double alpha;

  /* Creates a new Markov Space object with the given bounds
   *
   * Args: n = minimum inventory to simulate
   *       m = maximum inventory to simulate
   */
  Space(unsigned int m);

  // Gets the State at the given indices (by reference)
  State& state(unsigned int x, unsigned int y);

  // Perform value iteration
  double vi(double thresh, unsigned int maxIter);
};

#endif

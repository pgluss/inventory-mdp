#ifndef SPACE_H
#define SPACE_H

#include <vector>

class Space {
private:
  int n_;
  int m_;
  std::vector<float> fvals_;

public:
  double lambda1;
  double lambda2;
  double mu1;
  double mu2;
  double b1;
  double b2;
  double h1;
  double h2;
  double alpha;

  /* Creates a new Markov Space object with the given bounds
   *
   * Args: n = minimum inventory to simulate
   *       m = maximum inventory to simulate
   */
  Space(int n, int m);

  // Gets the State at the given indices (by reference)
  float& f(int x1, int x2);

  // Calculates inventory cost function
  double h(int x1, int x2);

  // Perform value iteration
  double vi(double thresh, unsigned int maxIter);

  // Get decisions
  void decide(std::string filename);
};

#endif

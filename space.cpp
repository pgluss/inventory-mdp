#include <iostream>
#include <algorithm>
#include <cmath>

#include "space.hpp"

Space::Space(int n, int m) {
  this->m_ = m;
  this->n_ = n;

  // Initialize all states
  for (int i = n; i <= m; i++) {
    for (int j = n; j <= m; j++) {
      this->states_.push_back(State(i,j,0));
    }
  }
}

State& Space::state(int x1, int x2) {
  x1 = x1 - this->n_;
  x2 = x2 - this->n_;
  int index = x1*(m_ - n_ + 1) + x2;
  return this->states_[index];
}

double Space::h(int x1, int x2) {
  double hVal = 0.0;
  hVal += (x1 < 0) ? -x1*b1 : x1*h1;
  hVal += (x2 < 0) ? -x2*b2 : x2*h2;
  return hVal;
}

double Space::vi(double thresh, unsigned int maxIter) {
  double error = thresh + 1;
  int i = 0;
  while (error > thresh && i++ < maxIter) {
    double localError = 0;
    for (State& s: this->states_) {
      int x1 = s.x1;
      int x2 = s.x2;
      double old_f = s.f;

      if (x1 == m_) {
        s.f = this->state(x1-1, x2).f + h1/alpha;
      } else if (x2 == m_) {
        s.f = this->state(x1, x2-1).f + h2/alpha;
      } else if (x1 == n_) {
        s.f = this->state(x1+1, x2).f + b1/alpha;
      } else if (x2 == n_) {
        s.f = this->state(x1, x2+1).f + b2/alpha;
      } else {
        double f1 = this->state(x1-1, x2).f;
        double f2 = this->state(x1, x2-1).f;
        double f3 = this->state(x1+1, x2).f;
        double f4 = this->state(x1, x2+1).f;

        double minimum = std::min(s.f, f3);
        minimum = std::min(minimum, f4);

        s.f = (h(x1,x2) + lambda1*f1 + lambda2*f2 + mu1*minimum) / (mu1 + lambda1 + lambda2 + alpha);
      }

      localError = std::max(localError, std::fabs(s.f - old_f));
    }    
    error = localError;
  }
  return error;
}

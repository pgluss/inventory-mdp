#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>

#include "space.hpp"

Space::Space(int n, int m) {
  this->m_ = m;
  this->n_ = n;

  // Initialize all states
  for (int i = n; i <= m; i++) {
    for (int j = n; j <= m; j++) {
      this->fvals_.push_back(0);
    }
  }
}

float& Space::f(int x1, int x2) {
  x1 = x1 - this->n_;
  x2 = x2 - this->n_;
  int index = x1*(m_ - n_ + 1) + x2;
  return this->fvals_[index];
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

    for (int x1 = this->n_; x1 <= m_; x1++) {
      for (int x2 = this->n_; x2 <= m_; x2++) {
        double old_f = f(x1, x2);

        if (x1 == m_) {
          f(x1,x2) = f(x1-1, x2) + h1/alpha;
        } else if (x2 == m_) {
          f(x1,x2) = f(x1, x2-1) + h2/alpha;
        } else if (x1 == n_) {
          f(x1,x2) = f(x1+1, x2) + b1/alpha;
        } else if (x2 == n_) {
          f(x1,x2) = f(x1, x2+1) + b2/alpha;
        } else {
          float f1 = f(x1-1, x2);
          float f2 = f(x1, x2-1);
          float f3 = f(x1+1, x2);
          float f4 = f(x1, x2+1);

          float minimum = std::min(f(x1,x2), f3);
          minimum = std::min(minimum, f4);

          f(x1,x2) = (h(x1,x2) + lambda1*f1 + lambda2*f2 + mu1*minimum) / (mu1 + lambda1 + lambda2 + alpha);
        }

        localError = std::max(localError, std::fabs(f(x1,x2) - old_f));
      }
    }

    error = localError;
  }
  return error;
}

enum {
  NoProduce = 0,
  ProduceX1 = 1,
  ProduceX2 = 2
};

void Space::decide(std::string filename) {
  std::ofstream file(filename);

  file << "x1\t" << "x2\t" << "f-value\t" << "decision" << std::endl;

  for (int x1 = this->n_; x1 <= m_; x1++) {
    for (int x2 = this->n_; x2 <= m_; x2++) {
      
      if (x1 == m_ || x2 == m_) {
        continue;
      }

      file << x1 << "\t" << x2 << "\t" << f(x1,x2) << "\t";

      double f1 = f(x1,x2);      //NoProduce f-val
      double f2 = f(x1 + 1, x2); //ProduceX1 f-val
      double f3 = f(x1, x2 + 1); //ProduceX2 f-val

      if (f1 < f2 && f1 < f3) {
        file << NoProduce << std::endl;
      } else if (f2 < f1 && f2 < f3) {
        file << ProduceX1 << std::endl;
      } else {
        file << ProduceX2 << std::endl;
      }
    }
  }
  file.close();
}

/*  Markov decision process simulation for UW Foster School of Business Research
    Copyright (C) 2017 Jieling Han, Paula Gluss

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// ConsoleApplication2.cpp : Defines the entry point for the console application.
//


// value iteration.cpp : Defines the entry point for the console application.
//
// COMMENT
// 1. How to simulate (alpha + v = 1)??.
// 2. HOw to initialize values?
// 3. paper p48. Generating sequence of larger state spaces.

// Memo: Jieling Han 11/19/12
// Modified the code for reverse supply chain model with (x, y) where x = new product, y = returned product. 
// In the code, x = x1, y = x2

// Memo: Jieling Han		05/06/2011
// The optimal value now matches Ha's paper!!
// When updating the boundary values, no need to /NORMALIZECOST
// When updating the middle values, ALL parts need to /NORMALIZECOST

// Memo: JL		3/30/2013
// Changed the global parameters to local parameters

// Memo: JL 4/8/2013
// Changed to average cost
// The boundary of the average cost is updated in this way, not similar to Ha's.
// m_ is the upper bound of the space. When the current state is on the boundary,
// on the bounday, then it will bounce back to itself (instead of any of its neighbours).
// See the update of T4v

//#include "stdafx.h"	//for using with visual studio
#include "state.h"

#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <array>
#include <cmath>

using namespace std;

double startingValue;
double avgCost;
int BOUND = 40;

double space::h(int x1, int x2) {
	double currentH = 0;
	if (x1 < 0) {
		x1 = -x1;
		currentH += x1 * B1;
	} else {
		currentH += x1 * H1;
  }

	if (x2 < 0)	{
		x2 = -x2;
		currentH += x2 * B2;
	} else {
		currentH += x2 * H2;
  }

	return currentH;
}

/* Constructs a new space object
 * 
 * Takes two ints n and m
 */
space::space(int n, int m) : n_(n), m_(m) {
	int k = 0;

	for (int i = n; i <= m; ++i ) {
		for (int j = n; j <= m; ++j ) {
			// add boundary conditions here
			state *s = new state( i, j);
			// state *s = new state( i, j, 100);
			states_.push_back(s);
			if ( i == 0 && j == 0 ) {
				initial_ = s;
      }
			++k;
		}
	}
	//std::cout << k << " states were generated\n";
}

long space::index(int x1, int x2) {
  if (x1 < n_) x1 = n_ ;
	if (x1 > m_) x1 = m_ ;
	if (x2 < n_) x2 = n_ ;
	if (x2 > m_) x2 = m_ ;

  // returns index
	return (x2 - n_) + ((x1 - n_) * (m_ - n_ + 1));
}

double space::iteration(int iter, int method, int Sx, int Sy) {
	double error = 0.0;
	double max_diff = 0.0;
	double min_diff = 1000.0;

	std::vector<state*>::iterator si;
	int x1, x2;

	double T1v, T2v, T3v, T4v;
	double best_value_T2 = 0; 
	double best_value_T4 = 0;
	
	//std::cout << "Iteration " << iter << " starts..." << std::endl;

	for (si = states_.begin(); si != states_.end(); ++si) {
		state *s = (*si);
		x1 = s->x1();
		x2 = s->x2();

		if (x1 == 0 && x2 == 0) {
			startingValue = s->f();
    }

		double old_value = s->old_f();
	
		double value;
		long index1 = 0, index2 = 0, index3 = 0, index4 = 0, index5 = 0;
		state *s1, *s2, *s3, *s4, *s5;

		// handle boundary conditions

		//if ( x1 == m_ && x2 == m_ )
		//	continue;

		//if ( x1 == m_ ) {
		//	index1 = index( m_-1, x2 );
		//	s1 = states_[index1];
		//	value = s1->old_f() ; // avg cost
		//	//value = s1->old_f() + H1 / ALPHA; // total cost
		//	s->change_value( value );
		//	continue;
		//}

		//if ( x2 == m_ ) {
		//	index1 = index( x1, m_-1 );
		//	s1 = states_[index1];
		//	value = s1->old_f() ;
		//	//value = s1->old_f() + H2 / ALPHA; // total cost
		//	s->change_value( value );
		//	continue;
		//}


		// find indices of all successors
    index1 = (x1 > n_ ) ? index(x1-1, x2) : index(x1, x2);
    index2 = (x2 > n_ ) ? index(x1, x2-1) : index(x1, x2);
    index3 = (x1 < m_ ) ? index(x1+1, x2) : index(x1, x2);
    index4 = (x2 < m_ ) ? index(x1, x2+1) : index(x1, x2);

		if (x2 > n_ && x1 < m_) {
		  index5 = index(x1+1, x2-1);
		} else if (x2 == n_ && x1 < m_) {
			index5 = index(x1+1, x2);
		} else if (x2 > n_ && x1 == m_) {
			index5 = index(x1, x2);
		} else if (x2 == n_ && x1 == m_) {
			index5 = index(x1, x2);
		}

    // Get successors
    s1 = states_[index1];
    s2 = states_[index2];
		s3 = states_[index3];
		s4 = states_[index4];
		s5 = states_[index5];

		double best_value = 1e50;

		////////
    T1v = (x1 == 0) ? s->old_f() + B1 : s1->old_f();
    T2v = (x2 == 0) ? s->old_f() + B2 : s2->old_f();
		T3v = s4->old_f();

		// Update T4v

		// method = 0: Optimal
		// method = 1: ...
		// method = 2: Sx, Sy two thresholds
		if (method == 0) {
			// First, test when x2 == 0
			if (s->old_f() < s3->old_f())
			{
				T4v = s->old_f();
				s->modify_best_action(0);
			}
			else
			{
				T4v = s3->old_f();
				s->modify_best_action(1); // produce x from raw material
			}

			if (x2 > 0)
			{
				if (s5->old_f() < T4v)
				{
					s->modify_best_action(2); // produce x from y
					T4v = s5->old_f();
				}
			}
		} else if (method == 2) {
			// First, test when x2 == 0
			if (x2 == 0 ) {
				if (x1 < Sx) {
					T4v = s3->old_f();
					s->modify_best_action(1);
				} else {
					T4v = s->old_f();
					s->modify_best_action(0);
				}
			} else if (x2 > 0) {
				if (x1 < Sx) {
					if (x2 < Sy) {
						T4v = s3->old_f();
						s->modify_best_action(1);
					} else if (x2 >= Sy) {
						T4v = s5->old_f();
						s->modify_best_action(2);
					}
				} else {
					T4v = s->old_f();
					s->modify_best_action(0);
				}
			}
		}
		////////

		value = (h( x1, x2 ) + LAMBDA1 * T1v + LAMBDA2 * T2v + DELTA * T3v + MU * T4v)/NORMALIZECOST;

		s->change_value( value );

		double diff = fabs(value - old_value);

		if (diff > max_diff) {
			max_diff = diff;
		}
		if (diff < min_diff) {
			min_diff = diff;
		}

	}
	//std::cout << "Iteration " << iter << " " << initial_->f() << " ends. max diff: " << max_diff << " min diff: " << min_diff << std::endl;

	//replace all the old_f_ with f_
	for ( si = states_.begin(); si != states_.end(); ++si ) {
		state *s = (*si);
		s->change_old_value(s->f());
	}

	avgCost = max_diff;
	return max_diff - min_diff;
}

int 
space::printAll(int fileIdx)
{
	
	ofstream myfile, myfile1;
	string filename1 = "Avg_allResult" + to_string(fileIdx) + ".txt";
	myfile.open (filename1);
	
	myfile << "x1\tx2\tbestAction\t\n";
	
	// output all values
	std::vector<state*>::iterator si;
	int x1, x2;
	int k = 1;
	
	for ( si = states_.begin(); si != states_.end(); ++si, ++k ) {
		state *s = (*si);
		x1 = s->x1();
		x2 = s->x2();
		
		int action = s->best_action();
		int admission_action = s->best_admission_action();
	    //double value = s->f();
	
	
		myfile << x1 << "	" << x2 << "	"<< action /*<< "	| admission = "<< admission_action << " value: " << s->f()*/ <<"\n";
			}

	for ( si = states_.begin(); si != states_.end(); ++si, ++k ) {
		state *s = (*si);
		x1 = s->x1();
		x2 = s->x2();
		
		int action = s->best_action();
		int admission_action = s->best_admission_action();
	}


	myfile.close();

	return 1;
}

int
space::printSwitching(int fileIdx)
{
	
	ofstream myfile;
	string filename;
	filename = "Avg_switchingResult" + to_string(fileIdx) + ".txt";
	myfile.open (filename);
	myfile << "Switching Result... (lambda1 = " << LAMBDA1 << ", lambda2 = " << LAMBDA2 << ") \n";
	
	int ind1;
	int ind2;
	int flag = 0;
	
	state *st;
	int preBA = 1; // previous best action
	int curBA = 1; // current best action
	
	int *swtiching = new int[m_ - n_ + 1];
	int *actionSwi = new int[m_ - n_ + 1];
	int *actionPre = new int[m_ - n_ + 1];
	
	for (ind1 = n_; ind1 <= m_; ++ind1) {
		for (ind2 = n_; ind2 <= m_; ++ind2) {
			st = states_[index(ind1, ind2)];
			curBA = st->best_action();
						
			if (curBA == -1) {
				continue;
			}
			
			if (flag == 0) {
				preBA = curBA;
				flag = 1;
			}
			else {
				if (curBA != preBA) {
					swtiching[ind1-n_] = ind2;
					actionSwi[ind1-n_] = curBA;
					actionPre[ind1-n_] = preBA;
					flag = 0;
					break;
				}
			}			
		
		}
	}
	
	//std::cout <<	"Printing the switching points.." << std::endl;
	
	myfile << "(x1	x2	preAction	curAction) \n";
	
	for (int i = n_; i <= m_; ++i) {
		//std::cout << "x1 = " << i << "; x2 = " << swtiching[i-n_] << std::endl;
		myfile << i << "	" << swtiching[i-n_] << "	"<< actionPre[i-n_] << "	" << actionSwi[i-n_] << "\n";
	}
	
	myfile.close();
	return 1;
}

int
space::vi( int fileIndex, int method, int sx, int sy )
{
	int iter = 0;
	double error = 1.0;

	if (method == 0) {
		// Value iteration for optimal cost function
		while ( error > EPSILON && iter < MAXITER) {
			error = iteration( ++iter, method, 0, 0 );
		}
		printAll(fileIndex);
		printSwitching(fileIndex);
		cout << "Optimal value function. Avg cost is " << avgCost << endl ;
	} else if (method == 2) {
		while ( error > EPSILON && iter < MAXITER) {
			error = iteration( ++iter, method, sx, sy );
		}
		//cout << "(sx,sy) = " << sx << ", " << sy << "; Heuristics avg cost is ";
	}
	//cout << avgCost << endl;

	return iter;
}

//For testing the properties
double space::v(int x1, int x2)
{
	return states_[index(x1, x2)]->f();
}

double space::D1_v(int x1, int x2)
{
	return v(x1+1,x2) - v(x1, x2);
}
double space::D2_v(int x1, int x2)
{
	return v(x1, x2+1) - v(x1, x2);
}

double space::D11_v(int x1, int x2)
{
	return D1_v(x1+1, x2) - D1_v(x1, x2);
}

double space::D22_v(int x1, int x2)
{
	return D2_v(x1, x2+1) - D2_v(x1, x2);
}

double space::D12_v(int x1, int x2)
{
	return D2_v(x1+1, x2) - D2_v(x1, x2);
}

double space::D1_2_v(int x1, int x2)
{
	if (x2 <= 0)
	{
		cout << "D1_2_v has non-positive x2" << endl;
		return -1;
	}

	return v(x1+1, x2-1) - v(x1, x2);
	
}

double space::D11_2_v(int x1, int x2)
{
	return D1_2_v(x1+1, x2) - D1_2_v(x1, x2);
}

double space::D21_2_v(int x1, int x2)
{
	return D1_2_v(x1, x2+1) - D1_2_v(x1, x2);
}

double space::D1_21_2_v(int x1, int x2)
{
	if (x2 <= 1)
	{
		cout << "D1_21_2_v invalid x2" << endl;
		return -1;
	}
	return D1_2_v(x1+1, x2-1) - D1_2_v(x1, x2);
}

////////////////////////////////
double space::T1V(int x1, int x2)
{
	if (x1 == 0)
		return states_[index(x1,x2)]->f();
	else if (x1 > 0)
		return states_[index(x1-1, x2)]->f();
	else
	{
		cout << "T1V error inputs" << endl;
		return -1;
	}
}

double space::T2V(int x1, int x2)
{
	if (x2 == 0)
		return states_[index(x1, x2)]->f() + B2;
	else if (x2 > 0)
		return min(states_[index(x1, x2)]->f() + B2, states_[index(x1, x2-1)]->f());
	else
	{
		cout << "T2V error inputs" << endl;
		return -1;
	}
}

double space::T4V(int x1, int x2)
{
	double curMin = min(states_[index(x1, x2)]->f(), states_[index(x1+1, x2)]->f());
	if (x2 == 0) {
		return curMin;
	} else if (x2 > 0) {
		return  min(curMin, states_[index(x1+1, x2-1)]->f());
	} else {
    throw invalid_argument("x2 is less than zero");
  }

}
///////////////////Test optimal Operator T2

double space::D1_T2v(int x1, int x2)
{
	return T2V(x1+1,x2) - T2V(x1, x2);
}
double space::D2_T2v(int x1, int x2)
{
	return T2V(x1, x2+1) - T2V(x1, x2);
}

double space::D11_T2v(int x1, int x2)
{
	return D1_T2v(x1+1, x2) - D1_T2v(x1, x2);
}

double space::D22_T2v(int x1, int x2)
{
	return D2_T2v(x1, x2+1) - D2_T2v(x1, x2);
}

double space::D12_T2v(int x1, int x2)
{
	return D2_T2v(x1+1, x2) - D2_T2v(x1, x2);
}

double space::D1_2_T2v(int x1, int x2)
{
	if (x2 <= 0)
	{
		cout << "D1_2_v has non-positive x2" << endl;
		return -1;
	}

	return T2V(x1+1, x2-1) - T2V(x1, x2);
	
}

double space::D11_2_T2v(int x1, int x2)
{
	return D1_2_T2v(x1+1, x2) - D1_2_T2v(x1, x2);
}

double space::D21_2_T2v(int x1, int x2)
{
	return D1_2_T2v(x1, x2+1) - D1_2_T2v(x1, x2);
}

double space::D1_21_2_T2v(int x1, int x2)
{
	if (x2 <= 1)
	{
		cout << "D1_21_2_v invalid x2" << endl;
		return -1;
	}
	return D1_2_T2v(x1+1, x2-1) - D1_2_T2v(x1, x2);
}

/////////////////Test Operator T4
double space::D1_T4v(int x1, int x2)
{
	return T4V(x1+1,x2) - T4V(x1, x2);
}
double space::D2_T4v(int x1, int x2)
{
	return T4V(x1, x2+1) - T4V(x1, x2);
}

double space::D11_T4v(int x1, int x2)
{
	return D1_T4v(x1+1, x2) - D1_T4v(x1, x2);
}

double space::D22_T4v(int x1, int x2)
{
	return D2_T4v(x1, x2+1) - D2_T4v(x1, x2);
}

double space::D12_T4v(int x1, int x2)
{
	return D2_T4v(x1+1, x2) - D2_T4v(x1, x2);
}

double space::D1_2_T4v(int x1, int x2)
{
	if (x2 <= 0)
	{
		cout << "D1_2_v has non-positive x2" << endl;
		return -1;
	}

	return T4V(x1+1, x2-1) - T4V(x1, x2);
	
}

double space::D11_2_T4v(int x1, int x2)
{
	return D1_2_T4v(x1+1, x2) - D1_2_T4v(x1, x2);
}

double space::D21_2_T4v(int x1, int x2)
{
	return D1_2_T4v(x1, x2+1) - D1_2_T4v(x1, x2);
}

double space::D1_21_2_T4v(int x1, int x2)
{
	if (x2 <= 1)
	{
		cout << "D1_21_2_v invalid x2" << endl;
		return -1;
	}
	return D1_2_T4v(x1+1, x2-1) - D1_2_T4v(x1, x2);
}
//////////////////////


int space::testProperties()
{
	int ind1;
	int ind2;

	for (ind1 = 0; ind1 <= m_-3; ++ind1)
	{
		for (ind2 = 0; ind2 <= m_-3; ++ind2)
		{
			//// **********Test for the basic v function **************
			if (D1_v(ind1, ind2) < (-1)*B1)
				cout << "A1 failed" << endl; // passed

			if (D11_v(ind1, ind2) < 0)
				cout << "A2 failed" << endl; // passed

			if (D11_2_v(ind1, ind2+1) < 0)
				cout << "A3 failed" << endl; // passed

			if (D11_v(ind1, ind2+1) < D11_2_v(ind1, ind2+1))
				cout << "A3 long form failed at (" << ind1 << ", " << ind2+1 << endl; // failed

			if (D12_v(ind1, ind2) < 0)
				cout << "A4 failed" << endl; // passed

			if (D22_v(ind1, ind2) <= 0)
				cout << "A5 failed" << endl; // passed

			if (D21_2_v(ind1, ind2+1) > 0)
				cout << "A6 failed" << endl; // passed

			if (D1_21_2_v(ind1, ind2+2) < 0)
				cout << "A7 failed" << endl; // passed

			//*************Testing Optimal for T2********
			if (D1_T2v(ind1, ind2) < (-1)*B1)
				cout << "A1 failed" << endl; // passed

			if (D11_T2v(ind1, ind2) < 0)
				cout << "A2 failed" << endl; // passed

			if (D11_2_T2v(ind1, ind2+1) < 0)
				cout << "A3 failed" << endl; // passed

			if (D11_T2v(ind1, ind2+1) < D11_2_T2v(ind1, ind2+1))
				cout << "T2v: A3 long form failed at (" << ind1 << ", " << ind2+1 << endl; // FAILED AT (21, 1)

			if (D12_T2v(ind1, ind2) < 0)
				cout << "T2v: A4 failed (" << ind1 << ", " << ind2 << endl; // FAILED at (X1, 0)

			if (D22_T2v(ind1, ind2) <= 0)
				cout << "A5 failed" << endl; // passed

			if (D21_2_T2v(ind1, ind2+1) > 0)
				cout << "A6 failed" << endl; // passed

			if (D1_21_2_T2v(ind1, ind2+2) < 0)
				cout << "A7 failed" << endl; // passed

			//**********Testing Optimal for T4**********
			if (D1_T4v(ind1, ind2) < (-1)*B1)
				cout << "A1 failed" << endl; // passed

			if (D11_T4v(ind1, ind2) < 0)
				cout << "T4v: A2 failed at (" << ind1 << ", " << ind2 << endl; // passed

			if (D11_2_T4v(ind1, ind2+1) < 0)
				cout << "T4v: A3 failed at (" << ind1 << ", " << ind2 << endl; // passed

			if (D11_T4v(ind1, ind2+1) < D11_2_T4v(ind1, ind2+1))
				cout << "T4v: long A3 form failed at (" << ind1 << ", " << ind2+1 << endl; // passed

			if (D12_T4v(ind1, ind2) < 0)
				cout << "T4v: A4 failed at (" << ind1 << ", " << ind2 << endl; // passed

			if (D22_T4v(ind1, ind2) <= 0)
				cout << "A5 failed at (" << ind1 << ", " << ind2 << endl; // passed

			if (D21_2_T4v(ind1, ind2+1) > 0)
				cout << "A6 failed at (" << ind1 << ", " << ind2+1 << endl; // passed

			if (D1_21_2_T4v(ind1, ind2+2) < 0)
				cout << "A7 failed at (" << ind1 << ", " << ind2+2 << endl; // passed
		}
	}

	return 0;
}


int main (int argc, char* const argv[]) {
	clock_t starting;

	char c;
	
	int method = 2;

	int sx = 0;
	int sy = 0;
	int fileIdx = 0;

	double H1 = 2;
	double H2 = 1;
	double B1 = 500;
	double B2 = 250;
	double P1 = 0;
	double P2 = 0;
	double MU = 1;
	double LAMBDA1 = 0.8;
	double LAMBDA2 = 0.4;
	double DELTA = 0.45; 
	double ALPHA = 0;

	double optimalAvgCost, avgCostHeuristics, costIncrease;
	double avgCostPre, avgCostCur;

  // (Sx, Sy), method = 2 since there are two parameters
  // (Sx, Sy=1) is a special case. Always produce from used, when there is any.
  // (Sx, Sy=+inf) is a special case. Always produce from raw material.

  int sxOptimal;
  int syOptimal;

  fileIdx = 0;

  int sxMin = 1;
  int sxMax = BOUND;
  int syMin, syMax;

  starting = clock();
  avgCostHeuristics = 10000;

  /* Assume the avg cost rate is convex in both x and y. Find sxOptimal and syOptimal by
   * searching for the smallest x, such that f(x+1,y) >= f(x,y)
   * Reset avgCostPre and avgCostCur before finding the optimal sy and sx
   */
  avgCostPre = avgCostCur = 10000; 

  //First, find the optimal sx. With sy = syMin fixed.
  for (sx = sxMin; sx <= sxMax; sx++ ) {
    /* due to the orginal program, the starting point should be -1 such that the real states
     * beginning from 0.
     */
    space *sp = new space(0, BOUND); 

    //H1, H2, MU, B1, B2, P1, P2, LAMBDA1, LAMBDA2, DELTA, 
    sp->setParameters(H1, H2, MU, B1, B2, P1, P2, LAMBDA1, LAMBDA2, DELTA, ALPHA);
    sp->vi(method*10+fileIdx, method, sx, syMin);

    avgCostCur = avgCost;
    delete sp;

    if ((avgCostCur < avgCostPre)) {
      avgCostPre = avgCostCur;
      if (sx == sxMin && sy == syMin) {
        continue;
      } else if (sx == sxMax) {
        sxOptimal = sxMax;
        //avgCostHeuristics = avgCostCur;
        break;
      } else {
        continue;
      }
    } else {
      sxOptimal = sx-1;
      //avgCostHeuristics = avgCostPre;
      break;
    } 

  }

  // Next, find the optimal sy
  avgCostPre = avgCostCur = 10000; 
  for (sy = syMin; sy <= syMax; sy++) {

    /* due to the orginal program, the starting point should be -1 such that the real states 
     * beginning from 0.
     */
    space *sp = new space(0, BOUND);

    //H1, H2, MU, B1, B2, P1, P2, LAMBDA1, LAMBDA2, DELTA, 
    sp->setParameters(H1, H2, MU, B1, B2, P1, P2, LAMBDA1, LAMBDA2, DELTA, ALPHA);

    sp->vi(method*10+fileIdx, method, sxOptimal, sy);						
    avgCostCur = avgCost;
    delete sp;

    if ((avgCostCur < avgCostPre)) {
      avgCostPre = avgCostCur;
      if (sy == syMin) {
        continue;
      } else if(sy == syMax) {
        syOptimal = syMax;
        avgCostHeuristics = avgCostCur;
        break;
      } else {
        continue;
      }
    } else {
      syOptimal = sy-1;
      avgCostHeuristics = avgCostPre;
      break;
    } 
  }			

  costIncrease = (avgCostHeuristics-optimalAvgCost)/optimalAvgCost;
  if (costIncrease < 0)
    //costIncrease = 0;
    costIncrease = (optimalAvgCost-avgCostHeuristics)/avgCostHeuristics;

  //	cout << " (Sx,Sy, AvgCost) = "<< sxOptimal <<", " << syOptimal << "	" << avgCostHeuristics << ".	Time: " << clock()-starting << endl;
  cout << H1 <<", " << H2 << ", " << B1 << ", " << B2;
  cout << ", " << DELTA <<", " << avgCostHeuristics << ", ("  << sxOptimal << ", ";
  cout << syOptimal <<") Time:" << double((clock()-starting)/1000) <<"s" << endl;

  //cout << "Heuristics, optimal (sx, sy) is (" << sxOptimal << ", " << syOptimal << "); Avg cost is " << avgCostOptimal << endl;
  //std::cin >> c;

	cout << "Complete!!" << endl;
	return 0;
}

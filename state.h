#ifndef STATE_H        
#define STATE_H

#include <vector>

class state
{
private:
	int x1_;
	int x2_;
	double f_;
	double old_f_;
	int best_action_;
	int best_admission_action_;

public:
	state( int x1, int x2 ) : x1_(x1), x2_(x2), f_(0.0), old_f_(0.0), best_action_(-1), best_admission_action_(-1) {}
	state( int x1, int x2, double f ) : x1_(x1), x2_(x2), f_(f), old_f_(f), best_action_(-1), best_admission_action_(-1) {}
	int x1( void ) { return x1_; }
	int x2( void ) { return x2_; }
	double f( void ) { return f_; }
	double old_f( void ) { return old_f_;}
	void change_value( double f ) { f_ = f; }
	void change_old_value(double f) { old_f_ = f;}
	int best_action( void ) { return best_action_; }
	int best_admission_action( void ) {return best_admission_action_;}
	void modify_best_action( int a ) { best_action_ = a; }
	void modify_best_admission_action( int b) { best_admission_action_ = b; }
};

class space
{
private:
	int n_;
	int m_;
	state *initial_;
	std::vector<state*> states_;

	//parameters
	double H1; // NEW
	double H2; // OLD
	double MU; // PRODUCTION
	double B1; // LOST SALES FOR NEW
	double B2; // LOST SALES FOR OLD
	double P1; // PRODUCTION COST FOR NEW
	double P2; // PRODUCTION COST FOR OLD
	double LAMBDA1; // NEW
	double LAMBDA2; // OLD
	double DELTA;  // RETURNED
	double ALPHA; // DISCOUNT RATE

	double NORMALIZECOST; 
	double EPSILON; 
	double	MAXITER; 

public:
	space( int, int );
    long index( int, int );
	double iteration( int, int, int, int );
	int vi( int, int, int, int);
	int printAll(int);
	int printSwitching(int);
	int setParameters(double h1, double h2, double mu, double b1, double b2, double p1, double p2, double lambda1, double lambda2,
		double delta, double alpha){H1=h1; H2=h2; MU=mu; B1=b1; B2=b2; P1=p1; P2=p2; LAMBDA1=lambda1; LAMBDA2=lambda2; DELTA=delta; ALPHA=alpha; 
		NORMALIZECOST=LAMBDA1+LAMBDA2+DELTA+MU+ALPHA; EPSILON=1e-6; MAXITER=1e6; return 1;}

	//for cost rate
	double h( int x1, int x2 );

	//for test the properties
	double v(int, int);
	double D1_v(int, int);
	double D2_v(int, int);
	double D11_v(int, int);
	double D22_v(int, int);
	double D12_v(int, int);
	double D1_2_v(int, int);
	double D11_2_v(int, int);
	double D21_2_v(int, int);
	double D1_21_2_v(int, int);

	double T1V(int, int);
	double T2V(int, int);
	double T4V(int, int);

	double D1_T2v(int, int);
	double D2_T2v(int, int);
	double D11_T2v(int, int);
	double D22_T2v(int, int);
	double D12_T2v(int, int);
	double D1_2_T2v(int, int);
	double D11_2_T2v(int, int);
	double D21_2_T2v(int, int);
	double D1_21_2_T2v(int, int);

	double D1_T4v(int, int);
	double D2_T4v(int, int);
	double D11_T4v(int, int);
	double D22_T4v(int, int);
	double D12_T4v(int, int);
	double D1_2_T4v(int, int);
	double D11_2_T4v(int, int);
	double D21_2_T4v(int, int);
	double D1_21_2_T4v(int, int);


	int testProperties(void);
};


#endif

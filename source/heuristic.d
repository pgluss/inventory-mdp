import std.stdio;
import std.math : abs, isNaN;
import std.algorithm.comparison : min,max;
import std.parallelism;
import std.datetime;
import std.range : iota;

enum double lambdax1 = 0.5;
enum double lambdax2 = 0.5;
enum double lambday  = 0.5; 

enum double mu = 1;
enum double gamma = 0.5 / 1000;
enum double alpha = 0.05;

enum double bx1 = 80;
enum double bx2 = 100;
enum double bx  = bx1*(lambdax1_hat / lambdax_hat) + bx2*(lambdax2_hat / lambdax_hat);
enum double by  = 50;

enum double hx1 = 2;
enum double hx2 = 2;
enum double hx  = hx1*(lambdax1_hat / lambdax_hat) + hx2*(lambdax2_hat / lambdax_hat);
enum double hy  = 1;

enum int Dx1 = 3;
enum int Dx2 = 3;
enum int Dy  = 3;

enum int EDx1 = (Dx1 + 1) / 2;
enum int EDx2 = (Dx2 + 1) / 2;
enum int EDy  = (Dy  + 1) / 2;

enum double lambdax1_hat = lambdax1*EDx1;
enum double lambdax2_hat = lambdax2*EDx2;
enum double lambdax_hat  = lambdax1_hat + lambdax2_hat;
enum double lambday_hat  = lambday*EDy;

enum double pn = 0.3;
enum double po = 0.3;

enum int n = -100;
enum int m = 100;

enum int EN = cast(int) ((lambdax_hat + lambday_hat) / (gamma*(1 - pn)));
enum double delta = EN*gamma;

__gshared float[] fvals;

// The main function. Runs the value iteration function.
void main() {
    fvals = new float[(m-n+1)^^2];

    for (int i = 0; i < fvals.length; i++) {
      fvals[i] = 0;
    }
    
    auto startTime = Clock.currTime(UTC());
    int numIter = vi(1e-2);
    auto elapsedTime = Clock.currTime(UTC())- startTime;

    decide("test-h.csv");

    stdout.writefln("Num States: %d, \tEN: %d, \tTotal Time: %s", fvals.length, EN, elapsedTime);
}

// Returns the f-value for a given (x,y) state
ref float f(int x, int y) {
    if (x < n) {
        x = n;
    } else if (x > m) {
        x = m;
    }

    if (y < n) {
        y = n;
    } else if (y > m) {
        y = m;
    }

    x = x - n;
    y = y - n;
    int i = x*(m-n+1) + y;
    return fvals[i];
}

// Returns the x value for a given index in the array of f-values
int x(int index) {
    return index / (m-n+1) + n;
}

// Returns the y value for a given index in the array of f-values
int y(int index) {
    return index % (m-n+1) + n;
}

// Computes the holding or backorder costs for the given inventory (x,y)
float h(int x, int y) {
    float hVal = 0.0;
    hVal += (x < 0) ? -x*bx : x*hx;
    hVal += (y < 0) ? -y*by : y*hy;
    return hVal;
}

// Performs value iterations until the error is below the specified threshold
// or the max number of iterations have ocurred
int vi(float thresh, uint maxIter = cast(int)1e6) {
    float maxError = thresh + 1;
    float avgError = 0;
    int i = 0;

    while (maxError > thresh && i++ < maxIter) {
      
        auto startTime = Clock.currTime(UTC());

        // Sequential
        //float localError = 0;
        //for (int index = 0; index < fvals.length; index++) {
        //  localError = max(viTask(index), localError);
        //}
        //maxError = localError;


        // Parallel
        auto localErrors = taskPool.amap!viTask(iota(cast(int)fvals.length));
        maxError = taskPool.reduce!max(localErrors);
        avgError = taskPool.reduce!"a+b"(localErrors)/localErrors.length;


        auto elapsedTime = Clock.currTime(UTC()) - startTime;

        stdout.writefln("%d \t\t %15f %15f \t\t %s", i, avgError, maxError, elapsedTime);
    }
    return i;
}

// Does value iteration for a single state
float viTask(int index) {
    int x = x(index);
    int y = y(index);

    float old_f = f(x,y);

    float T1v = f(x-1,y);
    float T2v = f(x,y-1);
    float T3v = f(x,y+1);
    float T4v = 0;

    if (y <= 0) {
        T4v = min(f(x,y), f(x+1,y));
    } else {
        T4v = min(f(x,y), f(x+1,y), f(x+1,y-1));
    }

    f(x,y) = (h(x, y) + lambdax_hat*T1v + lambday_hat*T2v + delta*T3v + mu*T4v) /
                                      (alpha + lambdax_hat + lambday_hat + mu + EN*gamma);

    return abs(f(x, y) - old_f);
}

enum int NoProduce = 0;   // Idle
enum int ProduceRaw = 1;  // Produce new product from raw materials
enum int ProduceUsed = 2; // Produce new product by refurbishing used product

// Get decisions for each state and print to a file
void decide(string filename) {
    File file = File(filename, "w");
    file.writefln("x, y, f-value, decision");

    for (int i = 0; i < fvals.length; i++) {
        int x = x(i);
        int y = y(i);

        double fval = f(x,y);
        int d;

        if (x == m) {
            d = NoProduce;
        } else if (y <= 0) {
            if (f(x,y) < f(x+1,y)) {
                d = NoProduce;
            } else {
                d = ProduceRaw;
            }
        } else if (y > 0) {
            if (f(x,y) < f(x+1,y) && f(x,y) < f(x+1,y-1)) {
                d = NoProduce;
            } else if (f(x+1,y) < f(x,y) && f(x+1,y) < f(x+1,y-1)) {
                d = ProduceRaw;
            } else {
                d = ProduceUsed;
            }
        } else {
            assert(false, "Error: Check x and y!");
        }

        file.writefln("%d, %d, %f, %d", x, y, fval, d);
    }
}


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

enum int n = -50;
enum int m = 50;

enum int EN = cast(int) ((lambdax_hat + lambday_hat) / (gamma*(1 - pn)));
enum double delta = EN*gamma;

__gshared float[] fvals;

void main() {
    fvals = new float[(m-n+1)^^2];

    for (int i = 0; i < fvals.length; i++) {
      fvals[i] = 0;
    }
    
    auto startTime = Clock.currTime(UTC());
    int numIter = vi(1e-2);
    auto elapsedTime = Clock.currTime(UTC())- startTime;

    stdout.writefln("Num States: %d, \tEN: %d, \tTotal Time: %s", fvals.length, EN, elapsedTime);
}


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


int x(int index) {
    return index / (m-n+1) + n;
}

int y(int index) {
    return index % (m-n+1) + n;
}

float h(int x, int y) {
    float hVal = 0.0;
    hVal += (x < 0) ? -x*bx : x*hx;
    hVal += (y < 0) ? -y*by : y*hy;
    return hVal;
}

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



/+
double calcEN(double thresh, int maxIter = 100) {
    double ENFrac = (lambdax_hat + lambday_hat) / (gamma*(1 - pn));
    double maxError = thresh + 1;

    double inside = 1;
    double en = 0;
    double fact = 1;

    for (int j = 1; j <= maxIter; j++) {
        fact *= j;
        inside += (1.0/fact)*(ENFrac^^j);
    }

    maxError = thresh + 1;
    fact = 1;

    for (int i = 1; i <= maxIter; i++) {
        fact *= i;
        double pi = ((1.0/fact) * (ENFrac^^i)) / inside;

        double old_val = en;
        en += i*pi;

        maxError = abs(en - old_val);
    }

    /+while (maxError > thresh && i++ < maxIter) {
        fact *= i;
        double pi = ((1.0/fact) * (ENFrac^^i)) / inside;

        double old_val = en;
        en += i*pi;

        maxError = abs(en - old_val);
    } +/

    //writeln("Average N is: ", en);

    return en;
}

+/

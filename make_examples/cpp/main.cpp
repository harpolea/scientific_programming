#include <math.h>

#include <iostream>
#include <fstream>

#define PI 3.14159265

using namespace std;

int main() {

    ofstream outfile ("data/data.csv");

    int n = 100;
    float xs[n];
    float ys[n];

    outfile << "x, y\n";

    for (int i = 0; i < n; i++) {
        float x = 2.0 * PI * i / n;
        float y = sin(x);

        outfile << x << ',' << y <<  '\n';
    }

    outfile.close();

}

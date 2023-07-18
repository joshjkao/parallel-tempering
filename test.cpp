#include <random>
#include <iostream>
#include <chrono> 

using namespace std;

int main(int argc, char** argv) {
    int x = 0;

    auto start = chrono::high_resolution_clock::now();

    mt19937_64 mt(0);
    uniform_int_distribution<int> distr(0, 1000);

    for (int i = 0; i < 1000000; ++i) {
        for (int j = 0; j < 1000; ++j) {
            x += distr(mt);
            int* y = new int;
            *y = x;
            x *= *y;
            delete y;
        }
    }
    cout << x << endl;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

    cout << duration.count()/1e6 << endl;

    return 0;
}
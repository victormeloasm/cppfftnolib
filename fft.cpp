#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <fstream>
#include <thread>
#include <windows.h>

using namespace std;
using namespace chrono;

using cd = complex<double>;
const double PI = acos(-1);

// Optimized iterative FFT with parallel blocks
void fft(vector<cd>& a, bool invert, int num_threads) {
    int n = a.size();
    int lg_n = 0;
    while ((1 << lg_n) < n) ++lg_n;

    vector<int> rev(n);
    for (int i = 0; i < n; ++i) {
        rev[i] = 0;
        for (int j = 0; j < lg_n; ++j)
            if (i & (1 << j))
                rev[i] |= 1 << (lg_n - 1 - j);
        if (i < rev[i])
            swap(a[i], a[rev[i]]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cos(ang), sin(ang));

        int chunk_size = (n / len + num_threads - 1) / num_threads;
        vector<thread> threads;

        for (int t = 0; t < num_threads; ++t) {
            int chunk_start = t * chunk_size;
            int chunk_end = min(chunk_start + chunk_size, n / len);
            if (chunk_start >= chunk_end) break;

            threads.emplace_back([&, chunk_start, chunk_end]() {
                for (int k = chunk_start; k < chunk_end; ++k) {
                    int i = k * len;
                    cd w(1);
                    for (int j = 0; j < len / 2; ++j) {
                        cd u = a[i + j], v = a[i + j + len / 2] * w;
                        a[i + j] = u + v;
                        a[i + j + len / 2] = u - v;
                        w *= wlen;
                    }
                }
            });
        }

        for (auto& th : threads) th.join();
    }

    if (invert)
        for (cd& x : a)
            x /= n;
}

// Multiply two big integer digit arrays using FFT
vector<int> multiply(const vector<int>& a, const vector<int>& b, int num_threads,
                     double& fft_time, double& carry_time) {
    vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < int(a.size() + b.size())) n <<= 1;
    fa.resize(n);
    fb.resize(n);

    auto fft_start = high_resolution_clock::now();
    fft(fa, false, num_threads);
    fft(fb, false, num_threads);
    for (int i = 0; i < n; ++i)
        fa[i] *= fb[i];
    fft(fa, true, num_threads);
    auto fft_end = high_resolution_clock::now();
    fft_time = duration<double>(fft_end - fft_start).count();

    // Carry correction
    auto carry_start = high_resolution_clock::now();
    vector<int> result(n);
    long long carry = 0;
    for (int i = 0; i < n; ++i) {
        long long t = (long long)(round(fa[i].real())) + carry;
        result[i] = t % 10;
        carry = t / 10;
    }
    while (carry) {
        result.push_back(carry % 10);
        carry /= 10;
    }
    while (result.size() > 1 && result.back() == 0)
        result.pop_back();
    auto carry_end = high_resolution_clock::now();
    carry_time = duration<double>(carry_end - carry_start).count();

    return result;
}

// Generate random big integer with 'size' digits
vector<int> generateRandom(int size) {
    vector<int> num(size);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 9);
    for (int& d : num)
        d = dis(gen);
    return num;
}

// Save big integer to file
void saveNumber(const string& filename, const vector<int>& num) {
    ofstream fout(filename);
    if (!fout) {
        cerr << "Error opening " << filename << " for writing!" << endl;
        return;
    }
    for (int i = num.size() - 1; i >= 0; --i)
        fout << num[i];
    fout.close();
}

int main() {
    const int SIZE = 1000000;

    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    int num_threads = sysinfo.dwNumberOfProcessors;

    cout << "Detected CPU cores: " << num_threads << endl;

    cout << "Generating random numbers..." << endl;
    auto num1 = generateRandom(SIZE);
    auto num2 = generateRandom(SIZE);

    saveNumber("num1.txt", num1);
    saveNumber("num2.txt", num2);

    double fft_time = 0, carry_time = 0;

    auto start = high_resolution_clock::now();
    auto result = multiply(num1, num2, num_threads, fft_time, carry_time);
    auto end = high_resolution_clock::now();
    double total_time = duration<double>(end - start).count();

    saveNumber("result.txt", result);

    cout << "Numbers saved to num1.txt, num2.txt" << endl;
    cout << "Result saved to result.txt" << endl;
    cout << "Total time: " << total_time << " seconds" << endl;
    cout << "  FFT time: " << fft_time << " seconds" << endl;
    cout << "  Carry time: " << carry_time << " seconds" << endl;

    return 0;
}

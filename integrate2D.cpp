#include <iostream>
#include <fstream>
#include "integrate2D.h"
#include "abscissae.h"
#include "weights.h"

using namespace std;

#define PRECISION 64

long double integral(long double (Node::*ref)(long double, long double, int, int), int alpha, int beta, long double a1, long double b1, long double a2, long double b2, Node& node) {
    long double diff1 = (b1 - a1) / 2;
    long double avg1 = (b1 + a1) / 2;
    long double diff2 = (b2 - a2) / 2;
    long double avg2 = (b2 + a2) / 2;

    long double ans = 0;

    for (int i = 0; i < PRECISION; i++) {
        for (int j = 0; j < PRECISION; j++) {
            ans += weights[i] * weights[j] * (node.*ref)(diff1 * abscissae[i] + avg1, diff2 * abscissae[j] + avg2, alpha, beta);
        }
    }

    return diff1 * diff2 * ans;
}
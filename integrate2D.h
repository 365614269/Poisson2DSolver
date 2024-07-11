#pragma once

#include "Node.h"

using namespace std;

long double integral(long double (Node::*)(long double, long double, int, int), int, int, long double, long double, long double, long double, Node&);
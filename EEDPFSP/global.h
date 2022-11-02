#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <windows.h>
#include<assert.h>
#include<algorithm>
#include<numeric>
#define numObjectives 2
#define jobsize 20
#define machinesize 8
#define factorysize 2
#define popsize 30
#define URAND (rand()/(RAND_MAX+1.0))

using namespace std;

vector<vector<int>>Time(jobsize,vector<int>(machinesize,0));
static double Speed[] = { 1,1.3,1.55,1.75,2.10 };
static double SP = 1.0;

void ReadFile(int k)
{
	ifstream in;
	if (k == 1) {
		in.open("./instances/I_20_8_2_1.txt");
	}
	else if (k == 2) {
		in.open("./instances/I_20_16_2_2.txt");
	}
	else if (k == 3) {
		in.open("./instances/I_20_16_2_3.txt");
	}
	else if (k == 4) {
		in.open("./instances/I_20_16_2_4.txt");
	}
	else if (k == 5) {
		in.open("./instances/I_20_16_2_5.txt");
	}
	else if (k == 6) {
		in.open("./instances/I_20_16_2_6.txt");
	}
	else if (k == 7) {
		in.open("./instances/I_20_16_2_7.txt");
	}
	else if (k == 8) {
		in.open("./instances/I_20_16_2_8.txt");
	}
	else if (k == 9) {
		in.open("./instances/I_20_16_2_9.txt");
	}
	else if (k == 10) {
		in.open("./instances/I_20_16_2_10.txt");
	}
	assert(in.is_open());
	int number;
	vector<int> Arr;
	vector<int>::iterator itA;
	while (in >> number, !in.eof()) {
		if (in.fail()) {
			in.clear();
			in.ignore();
			continue;
		}
		Arr.push_back(number);
	}
	itA = Arr.begin() + 2;
	for (int i = 0; i < jobsize; i++) {
		for (int j = 0; j < machinesize; j++) {
			itA = itA + 2;
			Time[i][j] = *itA;
		}
	}
}

double cmp(pair<int, double>a, pair<int, double>b)
{
	return a.second < b.second;
}

int rand_int(int low, int high)
{
	return int((high - low + 1) * URAND) + low;
}

double rand_real(double low, double high)
{
	double h;
	h = (high - low) * URAND + low + 0.001;
	if (h >= high)
		h = high - 0.001;
	return h;
}

bool is_dominated(const vector<double>& a, const vector<double>& b) {
	if (a[0] <= b[0] && a[1] <= b[1]) {
		if (a[0] == b[0] && a[1] == b[1])
			return false;
		else
			return true;
	}
	else
		return false;
}
#endif

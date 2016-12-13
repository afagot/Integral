#ifndef FONCTIONS_H
#define FONCTIONS_H

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <sstream>
#include <assert.h>

#include "TFile.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"

using namespace std;

//Function to cast an integer into a string
string intTostring(int value){
	string word;
	stringstream ss;
	ss << setw(2) << setfill('0') << value;
	ss >> word;

	return word;
};

//Function to cast a float into a string
string floatTostring(float value){
	string word;
	stringstream ss;
	ss.setf(ios::fixed, ios::floatfield);
	ss.precision(3);
	ss << value;
	ss >> word;

	return word;
}

#endif // FONCTIONS_H

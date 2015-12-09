#include <iostream>
#include <fstream>
#include <cstdlib>
#include <list>
#include <sstream>
#include <string>
#include "implicant.h"

using namespace std;

int main (int argc, char ** argv) {

    ifstream input;
    input.open("in_tmp.txt");

    int nVars, nMinterms;
    string word;

    input >> nVars >> nMinterms;
    list<Implicant> minterms;
    while(input >> word){
        minterms.push_back(Implicant(word));
    }

    Implicant i1(string("0101"));
    Implicant i2(string("1101"));
    cout << Implicant::consensus(i1,i2) << endl;

    return 0;
}


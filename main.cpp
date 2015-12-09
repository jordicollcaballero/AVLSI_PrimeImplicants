#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <cmath>
#include <bitset>

using namespace std;

typedef dynamic_bitset minterm;

void print(int * x, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            cout << (x[i*n+j]>0 ? "X" : ".");
        }
        cout << endl;
    }
}

minterm parseMinterm(const string &word){
    char c;
    minterm m;
    for(int i = 0; i < word.length(); i++){
        m[i] = c == '0' ? 0 : 1;
    }
}

int main (int argc, char ** argv) {

    ifstream input;
    input.open("in_tmp.txt");

    int nVars, nMinterms;
    string word;

    input >> nVars >> nMinerms;
    list<minterm> minterms;
    while(input >> word){
        minterms << parseMinterm(word);
    }


    return 0;
}


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <list>
#include <sstream>
#include <string>
#include "implicant.h"
#include "matrix.h"

using namespace std;
using namespace boost;

void printTable(list<Implicant> ** table, int nVars){
    for(int col = 0; col <= nVars; col++){
        cout << "n var in literal = " << nVars-col << endl;
        cout << "=====================" << endl;
        for(int row = 0; row <= nVars-col; row++){
            cout << "n negated = " << row << endl;
            cout << "--------------" << endl;
            for(Implicant &im : table[col][row])
                cout << im << endl;
            cout << endl;
        }
        cout << "=====================" << endl;
    }
}

list<Implicant> tabularMethod(const list<Implicant>& minterms){
    list<Implicant> implicants;
    int nVars = implicants.front().getNVars();
    list<Implicant> ** table = new list<Implicant> * [nVars+1];
    for(int i = 0; i <= nVars; i++)
        table[i] = new list<Implicant> [nVars-i+1];
    for(const Implicant& im : minterms)
        table[0][im.countNegatedVars()].push_back(im);

    for(int col = 0; col < nVars; col++){
        for(int row = 0; row < nVars-col; row++){
            for(Implicant &im1 : table[col][row+1]){
                for(Implicant &im2 : table[col][row]){
                    cout << im1 << " " << im2 << endl;
                    Implicant consensus = Implicant::distance1Merging(im1,im2);
                    if(consensus.isValid()){
                        //Vigilar al afegir, pot ser que ja hi sigui!
                        table[col+1][row].push_back(consensus);
                    }
                }
            }
        }
    }

    printTable(table, nVars);

    for(int i = 0; i < nVars; i++)
        delete [] table[i];
    delete [] table;

    return implicants;
}


list<Implicant> iteratedConsensus(const list<Implicant> implicants){
    list<Implicant> prime = implicants;
    list<Implicant>::iterator it1 = prime.begin();
    it1++;
    while(it1!=prime.end()){
        list<Implicant>::iterator it2 = prime.begin();
        while(it2!=it1){
            if((*it1).covers(*it2)){
                prime.erase(it2++);
            }
            else if((*it2).covers(*it1)){
                prime.erase(it1++);
                break;
            }
        }
        if(it1==it2) it1++;
    }

    it1 = prime.begin();
    it1++;
    while(it1!=prime.end()){
        list<Implicant>::iterator it2 = prime.begin();
        while(it2!=it1){
            Implicant consensus = Implicant::consensus(*it1,*it2);
            if(consensus.isValid()){
                prime.push_back(consensus);
                if(consensus.covers(*it2))
                    prime.erase(it2++);
                if(consensus.covers(*it1)){
                    prime.erase(it1++);
                    break;
                }
            }
            it2++;
        }
        if(it1==it2) it1++;
    }
}

list<Implicant> iteratedConsensus2(const list<Implicant> implicants){
    list<Implicant> prime = implicants;
    list<Implicant>::iterator it1 = prime.begin();
    it1++;
    while(it1!=prime.end()){
        list<Implicant>::iterator it2 = prime.begin();
        while(it2!=it1){
            if((*it1).covers(*it2)){
                prime.erase(it2++);
            }
            else if((*it2).covers(*it1)){
                prime.erase(it1++);
                break;
            }
            else{
                Implicant consensus = Implicant::consensus(*it1,*it2);
                if(consensus.isValid())
                    prime.push_back(consensus);
                it2++;
            }
        }
        if(it1==it2) it1++;
    }
}

dynamic_bitset<> exactCover(Matrix &A, dynamic_bitset<> x, dynamic_bitset<> b){
    int c;

    A.reduce(x);
    if(x.count() >= b.count()) return b;
    if(A.empty()) return x;

    c = A.selectBranchingColumn();

    Matrix A2 = A;
    x.set(A2.removeColumnAndRows(c));
    dynamic_bitset<> x2 = exactCover(A2,x,b);
    if(x2.count() < b.count()) b = x2;

    A2 = A;
    x.reset(A2.removeColumn(c));
    x2 = exactCover(A2,x,b);
    if(x2.count() < b.count()) b = x2;

    return b;
}


int main (int argc, char ** argv) {

    ifstream input;
    input.open("instances/test2.txt");

    string word;

    list<Implicant> minterms;
    while(input >> word)
        minterms.push_back(Implicant(word));

    input.close();

    /*
    Implicant im1(string("0101"));
    Implicant im2(string("1101"));
    Implicant im3(string("0001"));
    Implicant im4(string("1001"));
    Implicant consensus1 = Implicant::distance1Merging(im1,im2);
    Implicant consensus2 = Implicant::distance1Merging(im3,im4);
    Implicant consensus = Implicant::distance1Merging(consensus1,consensus2);
    cout << consensus1 << endl << consensus2 << endl;
    cout << (consensus.isValid() ? "valid" : "invalid") << endl;
    cout << consensus << endl;
    */

    //tabularMethod(minterms);


    return 0;
}


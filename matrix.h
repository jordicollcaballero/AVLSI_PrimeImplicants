#ifndef MATRIX_H
#define MATRIX_H

#include <boost/dynamic_bitset.hpp> //sudo apt install libboost1.55-all-dev
#include <list>
#include <set>
#include <string>

#include "implicant.h"

using namespace boost;

class Matrix
{
private:
    int nRows;
    int nColumns;
    int * idx;
    vector<bool> ** matrix;

    list<int> removedRowsStack;
    list<int> nRemovedRowsStack;
    int nRemovedRows;
    list<int> removedColumnsStack;
    list<int> nRemovedColumnsStack;
    int nRemovedColumns;


    int selectFirst() const;
    int selectRandom(std::default_random_engine &rnd_eng) const;
    int selectMaxOnes() const;
    int selectMinOnes() const;

    inline vector<bool>::reference mat(int row, int col);
    inline vector<bool>::const_reference const_mat(int row, int col) const;

    int essentialColumn(int row);
    bool columnDominance(int i, int j);
    bool rowDominance(int i, int j);


public:

    static int colSelectCriteria;

    static const int FIRST = 0;
    static const int RANDOM = 1;
    static const int MAX_ONES = 2;
    static const int MIN_ONES = 3;

    Matrix(bool ** vals, int nRow, int nCol);
    Matrix(const list<Implicant> &f, const list<Implicant> &primes);
    void reduce(dynamic_bitset<> &x);
    int selectBranchingColumn(std::default_random_engine &rnd_eng) const;
    void removeRow(int row);
    int removeColumn(int col);
    int removeColumnAndRows(int col);
    bool empty() const;
    void saveState();
    void restoreState();
    void print() const;
    ~Matrix();
};

#endif // MATRIX_H

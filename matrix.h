#ifndef MATRIX_H
#define MATRIX_H

#include <boost/dynamic_bitset.hpp> //sudo apt install libboost1.55-all-dev
#include <list>
#include <set>
#include "implicant.h"

using namespace boost;

class Matrix
{
private:
    int nRows;
    int nColumns;
    int * idx;
    vector<bool> ** matrix;
    list<list<int> > removedRowsStack;
    list<list<int> > removedColumnsStack;

    int essentialColumn(int row);
    bool columnDominance(int i, int j);
    bool rowDominance(int i, int j);

public:
    Matrix(bool ** vals, int nRow, int nCol);
    Matrix(const list<Implicant> &minterms, const set<Implicant> &idx);
    void reduce(dynamic_bitset<> &x);
    int selectBranchingColumn() const;
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

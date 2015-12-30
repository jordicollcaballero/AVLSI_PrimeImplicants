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

    list<int> removedRowsStack;
    list<int> nRemovedRowsStack;
    int nRemovedRows;
    list<int> removedColumnsStack;
    list<int> nRemovedColumnsStack;
    int nRemovedColumns;

    int colSelectCriteria;
    int selectFirst() const;
    int selectRandom() const;
    int selectMaxOnes() const;
    int selectMinOnes() const;

    inline vector<bool>::reference mat(int row, int col);

    int essentialColumn(int row);
    bool columnDominance(int i, int j);
    bool rowDominance(int i, int j);

    const static int FIRST = 0;
    const static int RANDOM = 1;
    const static int MAX_ONES = 2;
    const static int MIN_ONES = 3;

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
    void setColumnSelectionCriteria(int c);
    void print() const;
    ~Matrix();
};

#endif // MATRIX_H

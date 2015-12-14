#ifndef MATRIX_H
#define MATRIX_H

#include<list>
#include <boost/dynamic_bitset.hpp> //sudo apt install libboost1.55-all-dev
#include "implicant.h"

using namespace boost;

class Matrix
{
private:
    int nRows;
    int nColumns;
    vector<int> implicants;
    vector<vector<bool> > matrix;

    int essentialColumn(int row);
    bool columnDominance(int i, int j);
    bool rowDominance(int i, int j);

public:
    Matrix(const list<Implicant> &minterms, const list<Implicant> &implicants);
    void reduce(dynamic_bitset<> &x);
    int selectBranchingColumn();
    void removeRow(int row);
    int removeColumn(int col);
    int removeColumnAndRows(int col);
    bool empty();
};

#endif // MATRIX_H

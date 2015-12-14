#ifndef MATRIX_H
#define MATRIX_H

#include "implicant.h"

using namespace boost;

class Matrix
{
private:
    int nRows;
    int nColumns;

    vector<int> implicants;
    vector<vector<bool> > matrix;

public:
    Matrix(const list<Implicant> &minterms, const list<Implicant> &implicants);
    void reduce(dynamic_bitset x);
    int selectBranchingColumn();
    int removeColumn(int col);
    int removeColumnAndRows(int col);
    bool empty();
};

#endif // MATRIX_H

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
    int * idx; //Indexes of the columns
    vector<bool> ** matrix;

    //Keep track of the changes made to the matrix and be able to revert it
    list<int> removedRowsStack;
    list<int> nRemovedRowsStack;
    int nRemovedRows;
    list<int> removedColumnsStack;
    list<int> nRemovedColumnsStack;
    int nRemovedColumns;


    //Methods of selecting branching column
    int selectFirst() const;
    int selectRandom(std::default_random_engine &rnd_eng) const;
    int selectMaxOnes() const;
    int selectMinOnes() const;

    //Aux methods to index the 'matrix' data structure
    inline vector<bool>::reference mat(int row, int col);
    inline vector<bool>::const_reference const_mat(int row, int col) const;

    //Return the first essential column for this 'row' (if any)
    int essentialColumn(int row);
    //Check if column i dominates column j
    bool columnDominance(int i, int j);
    //Check if row i dominates row j
    bool rowDominance(int i, int j);


public:

    //Criteria to select the branhing column
    static int colSelectCriteria;

    static const int FIRST = 0;
    static const int RANDOM = 1;
    static const int MAX_ONES = 2;
    static const int MIN_ONES = 3;

    Matrix(bool ** vals, int nRow, int nCol);
    Matrix(const list<Implicant> &f, const list<Implicant> &primes);

    //Reduce the matrix applying column and row dominance, and essential columns
    void reduce(dynamic_bitset<> &x);
    //Select the branching column according to 'colSelectCriteria'
    int selectBranchingColumn(std::default_random_engine &rnd_eng) const;
    //Remove row 'row'
    void removeRow(int row);
    //Remove column 'col' and return its corresponding implicant
    int removeColumn(int col);
    //Remove column 'col', the rows that it satisfies, and return its corresponding implicant
    int removeColumnAndRows(int col);
    //Return true iff the matrix is empty
    bool empty() const;
    //Save the current state of the matrix.
    void saveState();
    //Return the matrix to the last saved state
    void restoreState();
    //print the matrix
    void print() const;

    ~Matrix();
};

#endif // MATRIX_H

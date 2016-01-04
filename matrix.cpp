#include "matrix.h"
#include <boost/range/adaptor/reversed.hpp>
#include <set>

int Matrix::colSelectCriteria = FIRST;

Matrix::Matrix(bool **vals, int nRow, int nCol)
{
    nRows = nRow;
    nColumns = nCol;
    idx = new int[nColumns];
    matrix = new vector<bool> *[nRows];
    for(int j = 0; j < nColumns; j++)
        this->idx[j] = j;
    for(int i = 0; i < nRows; i++){
        matrix[i] = new vector<bool>(nColumns);
        for(int j = 0; j < nColumns; j++)
            mat(i,j) = vals[i][j];
    }

    nRemovedRows = 0;
    nRemovedColumns = 0;
}

Matrix::Matrix(const list<Implicant> &minterms, const set<Implicant> &implicants)
{
    nRows = minterms.size();
    nColumns = implicants.size();
    idx = new int[nColumns];
    matrix = new vector<bool> *[nRows];
    for(int j = 0; j < nColumns; j++)
        this->idx[j] = j;

    int i = 0;
    for(const Implicant &minterm : minterms){
        int j = 0;
        matrix[i] = new vector<bool>(nColumns);
        for(const Implicant &implicant : implicants){
            mat(i,j) = implicant.covers(minterm);
            j++;
        }
        i++;
    }

    nRemovedRows = 0;
    nRemovedColumns = 0;
}

void Matrix::reduce(dynamic_bitset<> &x)
{

    set<int> to_remove;
    bool changes;

    do{
        changes = false;
        //Essential columns
        for(int i = 0; i < nRows; i++){
            int essential = essentialColumn(i);
            if(essential != -1)
                to_remove.insert(essential);
        }
        if(!to_remove.empty()){
            changes = true;
            for(const int &col : adaptors::reverse(to_remove))
                x.set(removeColumnAndRows(col));
            to_remove.clear();
        }

        //Column dominance
        for(int i = 1; i < nColumns; i++){
            for(int j = 0; j < i; j++){
                if(columnDominance(i,j))
                    to_remove.insert(j);
                else if(columnDominance(j,i))
                    to_remove.insert(i);
            }
        }
        if(!to_remove.empty()){
            changes = true;
            for(const int &col : adaptors::reverse(to_remove))
                removeColumn(col);
            to_remove.clear();
        }

        //Row dominance
        for(int i = 1; i< nRows; i++){
            for(int j = 0; j < i; j++){
                if(rowDominance(i,j))
                    to_remove.insert(i);
                else if(rowDominance(j,i))
                    to_remove.insert(j);
            }
        }
        if(!to_remove.empty()){
            changes = true;
            for(const int &row : adaptors::reverse(to_remove))
                removeRow(row);
            to_remove.clear();
        }

    }while(changes);

}

int Matrix::selectBranchingColumn(std::default_random_engine &rnd_eng) const
{
    switch(colSelectCriteria){
        case MAX_ONES:
            return selectMaxOnes();
        case MIN_ONES:
            return selectMinOnes();
        case RANDOM:
            return selectRandom(rnd_eng);
        case FIRST:
        default:
            return selectFirst();
    }
}

void Matrix::removeRow(int row)
{
    nRows--;
    vector<bool> *aux = matrix[row];
    matrix[row] = matrix[nRows];
    matrix[nRows] = aux;

    removedRowsStack.push_front(row);
    nRemovedRows++;
}

int Matrix::removeColumn(int col)
{
    nColumns--;
    int implicant = idx[col];
    idx[col] = idx[nColumns];
    idx[nColumns] = implicant;

    removedColumnsStack.push_front(col);
    nRemovedColumns++;
    return implicant;
}

int Matrix::removeColumnAndRows(int col)
{
    for(int i = 0; i < nRows; i++){
        if(mat(i,col)){
            removeRow(i);
            i--;
        }
    }
    return removeColumn(col);
}

bool Matrix::empty() const
{
    return nColumns == 0 || nRows == 0;
}

void Matrix::saveState()
{
    nRemovedColumnsStack.push_front(nRemovedColumns);
    nRemovedColumns = 0;
    nRemovedRowsStack.push_front(nRemovedRows);
    nRemovedRows = 0;
}

void Matrix::restoreState()
{
    for(int i = 0; i < nRemovedRows; i++){
        int row = removedRowsStack.front();
        removedRowsStack.pop_front();
        vector<bool> *aux = matrix[row];
        matrix[row] = matrix[nRows];
        matrix[nRows] = aux;
        nRows++;
    }
    nRemovedRows = nRemovedRowsStack.front();
    nRemovedRowsStack.pop_front();

    for(int i = 0; i < nRemovedColumns; i++){
        int col = removedColumnsStack.front();
        removedColumnsStack.pop_front();
        int aux = idx[col];
        idx[col] = idx[nColumns];
        idx[nColumns] = aux;
        nColumns++;
    }
    nRemovedColumns = nRemovedColumnsStack.front();
    nRemovedColumnsStack.pop_front();
}

void Matrix::print() const
{
    for(int j = 0; j < nColumns; j++)
        cout << idx[j] << " ";
    cout << endl << endl;
    for(int i = 0; i < nRows; i++){
        for(int j = 0; j < nColumns; j++)
            cout << (const_mat(i,j) ? 1 : 0) << " ";
        cout << endl;
    }
}

Matrix::~Matrix()
{
    delete [] idx;
    for(int i = 0; i < nRows; i++)
        delete matrix[i];
    delete [] matrix;
}

int Matrix::selectFirst() const
{
    return 0;
}

int Matrix::selectRandom(std::default_random_engine &rnd_eng) const
{
    std::uniform_int_distribution<int> uir(0, nColumns-1);
    return uir(rnd_eng);
}

int Matrix::selectMaxOnes() const
{
    int maxCol = 0;
    int maxOnes = 0;
    for(int col = 0; col < nColumns; col++){
        int nOnes = 0;
        for(int row = 0; row < nRows; row++){
            if(const_mat(row,col))
                nOnes++;
        }
        if(nOnes > maxOnes){
            maxOnes = nOnes;
            maxCol = col;
        }
    }
    return maxCol;
}

int Matrix::selectMinOnes() const
{
    int minCol = 0;
    int minOnes = nRows;
    for(int col = 0; col < nColumns; col++){
        int nOnes = 0;
        for(int row = 0; row < nRows; row++){
            if(const_mat(row,col))
                nOnes++;
        }
        if(nOnes < minOnes){
            minOnes = nOnes;
            minCol = col;
        }
    }
    return minCol;
}

vector<bool>::reference Matrix::mat(int row, int col)
{
    return (*(matrix[row]))[idx[col]];
}

vector<bool>::const_reference Matrix::const_mat(int row, int col) const
{
    return (*(matrix[row]))[idx[col]];
}

int Matrix::essentialColumn(int row){
    int nImp = 0;
    int selectedCol = -1;
    for(int col = 0; nImp <= 1 && col < nColumns ; col++){
        if(mat(row,col)){
            nImp++;
            selectedCol = col;
        }
    }
    return nImp == 1 ? selectedCol : -1;
}

bool Matrix::columnDominance(int i, int j){
    for(int row = 0; row < nRows; row++)
        if(mat(row,j) && !mat(row,i))
            return false;
    return true;
}

bool Matrix::rowDominance(int i, int j){
    for(int col = 0; col < nColumns; col++)
        if(mat(j,col) && !mat(i,col))
            return false;
    return true;
}


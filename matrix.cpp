#include "matrix.h"

#include <boost/range/adaptor/reversed.hpp>
#include <set>

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
        for(int j = 0; j < nColumns; j++){
            matrix[i]->resize(nColumns);
            (*(matrix[i]))[j] = vals[i][j];
        }
    }
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
        for(const Implicant &implicant : implicants){
            matrix[i] = new vector<bool>(nColumns);
            (*(matrix[i]))[j] = implicant.covers(minterm);
            j++;
        }
        i++;
    }
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
            for(const int &col : adaptors::reverse(to_remove)){
                x.set(removeColumnAndRows(col));
               // cout << "column and rows" << endl;
               // print();
            }
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
            for(const int &col : adaptors::reverse(to_remove)){
                removeColumn(col);
             //   cout << "column" << endl;
             //   print();
            }
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
            for(const int &row : adaptors::reverse(to_remove)){
                removeRow(row);
           //     cout << "row" << endl;
           //     print();
            }
            to_remove.clear();
        }

    }while(changes);

}

int Matrix::selectBranchingColumn() const
{
    return 0;
}

void Matrix::removeRow(int row)
{
    nRows--;
    vector<bool> *aux = matrix[row];
    matrix[row] = matrix[nRows];
    matrix[nRows] = aux;
    removedRowsStack.front().push_front(row);
}

int Matrix::removeColumn(int col)
{
    nColumns--;
    int implicant = idx[col];
    idx[col] = idx[nColumns];
    idx[nColumns] = implicant;
    removedColumnsStack.front().push_front(col);
    return implicant;
}

int Matrix::removeColumnAndRows(int col)
{
    for(int i = 0; i < nRows; i++){
        if((*(matrix[i]))[idx[col]]){
            removeRow(i);
            i--;
        }
    }

    nColumns--;
    int implicant = idx[col];
    idx[col] = idx[nColumns];
    idx[nColumns] = implicant;

    removedColumnsStack.front().push_front(col);
    return implicant;
}

bool Matrix::empty() const
{
    return nRows == 0;
}

void Matrix::saveState()
{
    removedColumnsStack.push_front(list<int>());
    removedRowsStack.push_front(list<int>());
}

void Matrix::restoreState()
{
    for(const int &row : removedRowsStack.front()){
        vector<bool> *aux = matrix[row];
        matrix[row] = matrix[nRows];
        matrix[nRows] = aux;
        nRows++;
    }
    removedRowsStack.pop_front();

    for(const int &col : removedColumnsStack.front()){
        int aux = idx[col];
        idx[col] = idx[nColumns];
        idx[nColumns] = aux;
        nColumns++;
    }
    removedColumnsStack.pop_front();
}

void Matrix::print() const
{
    for(int j = 0; j < nColumns; j++)
        cout << idx[j] << " ";
    cout << endl << endl;
    for(int i = 0; i < nRows; i++){
        for(int j = 0; j < nColumns; j++)
            cout << ((*(matrix[i]))[idx[j]] ? 1 : 0) << " ";
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

int Matrix::essentialColumn(int row){
    int nImp = 0;
    int selectedCol = -1;
    for(int col = 0; nImp <= 1 && col < nColumns ; col++){
        if((*(matrix[row]))[idx[col]]){
            nImp++;
            selectedCol = col;
        }
    }
    return nImp == 1 ? selectedCol : -1;
}

bool Matrix::columnDominance(int i, int j){
    for(int row = 0; row < nRows; row++)
        if((*(matrix[row]))[idx[j]] && !(*(matrix[row]))[idx[i]])
            return false;
    return true;
}

bool Matrix::rowDominance(int i, int j){
    for(int col = 0; col < nColumns; col++)
        if((*(matrix[j]))[idx[col]] && !(*(matrix[i]))[idx[col]])
            return false;
    return true;
}


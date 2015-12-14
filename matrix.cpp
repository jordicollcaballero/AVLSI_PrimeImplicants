#include "matrix.h"

Matrix::Matrix(const list<Implicant> &minterms, const list<Implicant> &implicants)
{
    nRows = minterms.size();
    nColumns = implicants.size();
    this->implicants.resize(nColumns);
    matrix.resize(nRows);
    for(int j = 0; j < nColumns; j++)
        this->implicants[j] = j;
    int i = 0;
    int j = 0;
    for(const Implicant &minterm : minterms){
        for(const Implicant &implicant : implicants){
            matrix[i].resize(nColumns);
            matrix[i][j] = implicant.covers(minterm);
            j++;
        }
        i++;
    }
}

void Matrix::reduce(dynamic_bitset<> &x)
{
    list<int> to_remove;
    bool changes;

    do{
        changes = false;

        //Essential columns
        for(int i = 0; i < nRows; i++){
            int essential = essentialColumn(i);
            if(essential != -1)
                to_remove.push_back(essential);
        }

        if(!to_remove.empty()){
            //sort to_remove decreasing
            changes = true;
            int lastChanged = nColumns;
            for(int &col : to_remove){
                if(col != lastChanged){ //may contain repetitions
                    x.set(removeColumnAndRows(col));
                    lastChanged = col;
                }
            }
        }

        //Column dominance
        for(int i = 1; i < nColumns; i++){
            for(int j = 0; j < i; j++){
                if(columnDominance(i,j))
                    to_remove.push_back(j);
                else if(columnDominance(j,i))
                    to_remove.push_back(i);
            }
        }

        if(!to_remove.empty()){
            //sort to_remove decreasing
            changes = true;
            int lastChanged = nColumns;
            for(int &col : to_remove){
                if(col != lastChanged){ //may contain repetitions
                    removeColumn(col);
                    lastChanged = col;
                }
            }
        }

        //Row dominance
        for(int i = 1; i< nRows; i++){
            for(int j = 0; j < i; j++){
                if(rowDominance(i,j))
                    to_remove.push_back(i);
                else if(rowDominance(j,i))
                    to_remove.push_back(j);
            }
        }

        if(!to_remove.empty()){
            //sort to_remove decreasing
            changes = true;
            int lastChanged = nRows;
            for(int &row : to_remove){
                if(row != lastChanged){
                    removeRow(row);
                    lastChanged = row;
                }
            }
        }

    }while(changes);

}

int Matrix::selectBranchingColumn()
{
    return 0;
}

void Matrix::removeRow(int row)
{
    nRows--;
    if(row < nRows)
        matrix[row] = matrix[nRows];
}

int Matrix::removeColumn(int col)
{
    nColumns--;
    int implicant = implicants[col];
    if(col < nColumns)
    {
        implicants[col]=implicants[nColumns];
        for(int i = 0; i < nRows; i++)
            matrix[i][col] = matrix[i][nColumns];
    }
    return implicant;
}

int Matrix::removeColumnAndRows(int col)
{
    nColumns--;
    int implicant = implicants[col];
    if(col < nColumns)
    {
        implicants[col]=implicants[nColumns];
        for(int i = 0; i < nRows; i++){
            if(matrix[i][col]){
                nRows--;
                if(i < nRows) matrix[i] = matrix[nRows];
                i--;
            }
            else matrix[i][col] = matrix[i][nColumns];
        }
    }
    return implicant;
}

bool Matrix::empty()
{
    return nRows == 0;
}

int Matrix::essentialColumn(int row){
    int nImp = 0;
    int selectedCol = -1;
    for(int col = 0; nImp <= 1 && col < nColumns ; col++){
        if(matrix[row][col]){
            nImp++;
            selectedCol = col;
        }
    }
    return nImp == 1 ? selectedCol : -1;
}

bool Matrix::columnDominance(int i, int j){
    for(int row = 0; row < nRows; row++)
        if(matrix[row][j] && !matrix[row][i])
            return false;
    return true;
}

bool Matrix::rowDominance(int i, int j){
    for(int col = 0; col < nColumns; col++)
        if(matrix[col][j] && !matrix[col][i])
            return false;
    return true;
}


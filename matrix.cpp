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

void Matrix::reduce(dynamic_bitset &x)
{
    do{
        bool changes = false;
        //Essential columns
        for(int j = nColumns-1; j >= 0; j--)
            if(isEssential(j))
                x.set(removeColumn(j));



        //Column cover


        //Row cover

    }while(changes);

}

int Matrix::selectBranchingColumn()
{
    return 0;
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



#include "implicant.h"

Implicant::Implicant()
{
}

Implicant::Implicant(const Implicant &i)
{
    mask = i.mask;
    variables = i.variables;
}

Implicant::Implicant(const string &s)
{
    mask = dynamic_bitset<>(s.size());
    variables = dynamic_bitset<>(s.size());
    int i = 0;
    for(const char& c : s) {
        if(c=='x')
            mask[i]=0;
        else{
            mask[i] = 1;
            variables[i] = c=='1' ? 1 : 0;
        }
        i++;
    }
}

dynamic_bitset<> Implicant::maskedVars() const{
    return mask & variables;
}

ostream& operator<<(ostream& os, const Implicant& imp)
{
    for(int i = 0; i < imp.mask.size(); i++){
        if(imp.mask[i])
            os << imp.variables[i];
        else
            os << 'x';
    }
    return os;
}

Implicant Implicant::consensus(const Implicant &i1, const Implicant &i2){
    Implicant i(i1);
    dynamic_bitset<> uncommon_variables = i1.mask ^ i2.mask;
    if(uncommon_variables.count()==0){
        i.mask = i1.mask;
        int disjoint_var = (i1.maskedVars() ^ i2.maskedVars()).find_first();
        i.mask[disjoint_var]=false;
    }
    else
        i.mask.reset();
    return i;
}


#include "implicant.h"
#include <cstdlib>
#include <iostream>


Implicant::Implicant()
{
    valid = false;
}

Implicant::Implicant(const Implicant &i)
{
    mask = i.mask;
    variables = i.variables;
    valid = i.valid;
}

Implicant::Implicant(const dynamic_bitset<> &mask, const dynamic_bitset<> &variables)
{
    this->mask = mask;
    this->variables = variables;
    valid = true;
}

Implicant::Implicant(const string &s)
{
    mask = dynamic_bitset<>(s.size());
    variables = dynamic_bitset<>(s.size());
    int i = 0;
    for(const char& c : s) {
        if(c=='x'){
            mask[i]=0;
            variables[i]=0;
        }
        else{
            mask[i] = 1;
            variables[i] = c=='1' ? 1 : 0;
        }
        i++;
    }
    valid = true;
}

int Implicant::getNVars() const
{
    return mask.size();
}

int Implicant::countNegatedVars() const
{
    dynamic_bitset<> aux(mask.size());
    aux.set();
    return (mask & (variables^aux)).count();
}

dynamic_bitset<> Implicant::maskedVars(const dynamic_bitset<>  & mask) const{
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

bool Implicant::covers(const Implicant &i) const{
    return this->mask.is_subset_of(i.mask) &&
            this->variables==i.maskedVars(this->mask);
}

bool Implicant::isValid() const
{
    return valid;
}

Implicant Implicant::distance1Merging(const Implicant &i1, const Implicant &i2){
    Implicant i;
    if(i1.mask == i2.mask){
        dynamic_bitset<> uncommonVariables = i1.variables ^ i2.variables;
        int distance = uncommonVariables.count();
        if(distance==1){
            i.mask = i1.mask;
            i.variables = i1.variables;
            i.valid = true;
            int uncommonVar = uncommonVariables.find_first();
            i.mask[uncommonVar]=false;
            i.variables[uncommonVar]=false;
        }
        else
            i.valid = false;
    }
    else
        i.valid = false;
    return i;
}


Implicant Implicant::consensus(const Implicant &i1, const Implicant &i2){
    Implicant i;
    dynamic_bitset<> commonVariables = i1.mask & i2.mask;
    dynamic_bitset<> uncommonVariables = i1.maskedVars(commonVariables) ^ i2.maskedVars(commonVariables);
    int distance = uncommonVariables.count();
    if(distance==1){
        int uncommonVar = uncommonVariables.find_first();
        i.mask = i1.mask|i2.mask;
        i.variables = i1.variables|i2.variables;
        i.mask[uncommonVar]=false;
        i.variables[uncommonVar]=false;
        i.valid = true;
    }
    else
        i.valid = false;
    return i;
}


Implicant Implicant::trueImplicant(int nVars)
{
    dynamic_bitset<> mask(nVars);
    dynamic_bitset<> variables(nVars);
    mask.reset();
    variables.reset();
    Implicant i(mask,variables);
    return i;
}


bool Implicant::operator<(const Implicant &rhs) const
{
    unsigned long maskulong = mask.to_ulong();
    unsigned long rhsmaskulong = rhs.mask.to_ulong();
    if(maskulong < rhsmaskulong) return true;
    else if(maskulong > rhsmaskulong) return false;
    else return variables.to_ulong() < rhs.variables.to_ulong();
}

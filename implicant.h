#ifndef IMPLICANT_H
#define IMPLICANT_H

#include <boost/dynamic_bitset.hpp> //sudo apt install libboost1.55-all-dev

using namespace std;
using namespace boost;

class Implicant
{
private:
    dynamic_bitset<> mask;
    dynamic_bitset<> variables;
    bool valid;

public:
    Implicant();
    Implicant(const Implicant &i);
    Implicant(const string &s);
    int getNVars() const;
    int countNegatedVars() const;
    bool covers(const Implicant &i) const;
    bool isValid() const;
    dynamic_bitset<> maskedVars() const;
    dynamic_bitset<> maskedVars(const dynamic_bitset<> &mask) const;
    friend ostream& operator<< (ostream &out, const Implicant& imp);
    static Implicant distance1Merging(const Implicant &i1, const Implicant &i2);
    static Implicant consensus(const Implicant &i1, const Implicant &i2);
    static Implicant trueImplicant(int nVars);

};


#endif // IMPLICANT_H

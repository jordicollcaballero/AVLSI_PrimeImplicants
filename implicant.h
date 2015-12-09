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

public:
    Implicant();
    Implicant(const Implicant &i);
    Implicant(const string &s);
    dynamic_bitset<> maskedVars() const;
    friend ostream& operator<< (ostream &out, const Implicant& imp);
    static Implicant consensus(const Implicant &i1, const Implicant &i2);
};


#endif // IMPLICANT_H

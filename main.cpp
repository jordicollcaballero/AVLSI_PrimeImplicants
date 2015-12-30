#include <iostream>
#include <fstream>
#include <cstdlib>
#include <list>
#include <set>
#include <sstream>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>
#include "time.h"
#include "implicant.h"
#include "matrix.h"

using namespace std;
using namespace boost;

void printTable(set<Implicant> ** table, int nVars){
    for(int col = 0; col <= nVars; col++){
        cout << "n var in literal = " << nVars-col << endl;
        cout << "=====================" << endl;
        for(int row = 0; row <= nVars-col; row++){
            cout << "n negated = " << row << endl;
            cout << "--------------" << endl;
            for(const Implicant &im : table[col][row])
                cout << im << endl;
            cout << endl;
        }
        cout << "=====================" << endl;
    }
}

void tabularMethod(const list<Implicant>& minterms, set<Implicant> &implicants){
    if(minterms.size() == 0)
        return;
    int nVars = minterms.front().getNVars();
    set<Implicant> ** table = new set<Implicant> * [nVars+1];
    for(int i = 0; i <= nVars; i++)
        table[i] = new set<Implicant> [nVars-i+1];
    for(const Implicant& im : minterms){
        table[0][im.countNegatedVars()].insert(im);
        implicants.insert(im);
    }

    for(int col = 0; col < nVars; col++){
        for(int row = 0; row < nVars-col; row++){
            for(const Implicant &im1 : table[col][row+1]){
                for(const Implicant &im2 : table[col][row]){
                    Implicant consensus = Implicant::distance1Merging(im1,im2);
                    if(consensus.isValid()){
                        table[col+1][row].insert(consensus);
                        implicants.insert(consensus);
                        implicants.erase(im1);
                        implicants.erase(im2);
                    }
                }
            }
        }
    }

    //printTable(table, nVars);

    for(int i = 0; i <= nVars; i++)
        delete [] table[i];
    delete [] table;
}


list<Implicant> iteratedConsensus(const list<Implicant> implicants){
    list<Implicant> prime = implicants;
    list<Implicant>::iterator it1 = prime.begin();
    it1++;
    while(it1!=prime.end()){
        list<Implicant>::iterator it2 = prime.begin();
        while(it2!=it1){
            if((*it1).covers(*it2)){
                prime.erase(it2++);
            }
            else if((*it2).covers(*it1)){
                prime.erase(it1++);
                break;
            }
        }
        if(it1==it2) it1++;
    }

    it1 = prime.begin();
    it1++;
    while(it1!=prime.end()){
        list<Implicant>::iterator it2 = prime.begin();
        while(it2!=it1){
            Implicant consensus = Implicant::consensus(*it1,*it2);
            if(consensus.isValid()){
                prime.push_back(consensus);
                if(consensus.covers(*it2))
                    prime.erase(it2++);
                if(consensus.covers(*it1)){
                    prime.erase(it1++);
                    break;
                }
            }
            it2++;
        }
        if(it1==it2) it1++;
    }
}

list<Implicant> iteratedConsensus2(const list<Implicant> implicants){
    list<Implicant> prime = implicants;
    list<Implicant>::iterator it1 = prime.begin();
    it1++;
    while(it1!=prime.end()){
        list<Implicant>::iterator it2 = prime.begin();
        while(it2!=it1){
            if((*it1).covers(*it2)){
                prime.erase(it2++);
            }
            else if((*it2).covers(*it1)){
                prime.erase(it1++);
                break;
            }
            else{
                Implicant consensus = Implicant::consensus(*it1,*it2);
                if(consensus.isValid())
                    prime.push_back(consensus);
                it2++;
            }
        }
        if(it1==it2) it1++;
    }
}

/*
 * A: matrix of covering of implicants on minterms
 * x: partial solution
 * b: best solution found until the moment
 */
dynamic_bitset<> exactCover(Matrix *A, dynamic_bitset<> x, dynamic_bitset<> b){
    int c;

    A->reduce(x); //Reduce the matrix using essentials and dominance
    if(x.count() >= b.count()) return b; //Bound
    if(A->empty()) return x; //Solution completed

    c = A->selectBranchingColumn();
    A->saveState();
    x.set(A->removeColumnAndRows(c)); //Include implicant of column c to the solution
    dynamic_bitset<> x2 = exactCover(A,x,b); //Branch
    A->restoreState();
    if(x2.count() < b.count()) b = x2;

    A->saveState();
    x.reset(A->removeColumn(c)); //Exclude implicant of column c to the solution
    x2 = exactCover(A,x,b); //Branch
    A->restoreState();
    if(x2.count() < b.count()) b = x2;

    return b;
}

void QuineMcCluskey(const list<Implicant> &minterms, list<Implicant> &result){
    set<Implicant> prime;
    tabularMethod(minterms,prime);

    dynamic_bitset <>x(prime.size()); x.reset();
    dynamic_bitset <>b(prime.size()); b.set();
    Matrix * m = new Matrix(minterms,prime);
    dynamic_bitset <>minset = exactCover(m,x,b);
    delete m;

    int i = 0;
    for(const Implicant &im : prime){
        if(minset[i])
            result.push_back(im);
        i++;
    }
}

void generateUniformRandom(int nvars, double q, default_random_engine &rnd_eng, list<Implicant> &minterms){
    dynamic_bitset<> mask(nvars,(1<<nvars)-1);
    uniform_real_distribution<double> uir(0.0, 1.0);
    for(int i = 0; i < 1<<nvars; i++){
        double p = uir(rnd_eng);
        if(p <= q) minterms.push_back(Implicant(mask,dynamic_bitset<> (nvars,i)));
    }
}

/*
bool verifyResult(const list<Implicant> minterms, const list<Implicant> implicants){
    set<Implicant> mintermsset;
    for(const Implicant & min : minterms)
        mintermsset.insert(min);
    for(const Implicant & im : implicants){
        for(const Implicant &min : im.getImplicants()){
            if(mintermsset.find(min) == mintermsset.end() && im.covers(min))
                return false;
            else if(mintermsset.find(min) != mintermsset.end() && !im.covers(min))
                return false;
        }
    }
    return true;
}*/

void experiments1(int argc, char ** argv){
    std::default_random_engine rnd_eng(1234);

    double fraction = 0.3;
    int min_nvars = 4;
    int max_nvars = 10;
    int nreps = 5;

    for(int nvars = min_nvars; nvars <= max_nvars; nvars++){
        double runtime = 0;

        for(int rep = 0; rep < nreps; rep++){
            list<Implicant> minterms;
            list<Implicant> result;
            generateUniformRandom(nvars,fraction,rnd_eng,minterms);

            clock_t initemps=clock();
            QuineMcCluskey(minterms,result);
            runtime += (clock()-initemps)/(double)CLOCKS_PER_SEC;

        }
        cout << nvars << ";" << (1<<nvars)*fraction << ";" << runtime/nreps << endl;
    }
}

void experiments11(int argc, char ** argv){
    std::default_random_engine rnd_eng(123);

    double fraction = 0.3;
    int min_nvars = 4;
    int max_nvars = 10;
    int nreps = 5;

    for(int nvars = min_nvars; nvars <= max_nvars; nvars++){
        double runtime = 0;

        for(int rep = 0; rep < nreps; rep++){
            list<Implicant> minterms;
            generateUniformRandom(nvars,fraction,rnd_eng,minterms);
            clock_t initemps=clock();
            set<Implicant> prime;

            tabularMethod(minterms,prime);
            runtime += (clock()-initemps)/(double)CLOCKS_PER_SEC;

        }
        cout << nvars << ";" << (1<<nvars)*fraction << ";" << runtime/nreps << endl;
    }
}

void experiments2(int argc, char ** argv){
    std::default_random_engine rnd_eng(1234);

    double min_fraction = 0.05;
    double max_fraction = 0.85;
    int nvars = 8;
    int nreps = 5;

    for(double fraction = min_fraction; fraction <= max_fraction; fraction+=0.05){
        double runtime = 0;
        for(int rep = 0; rep < nreps; rep++){
            list<Implicant> minterms;
            list<Implicant> result;
            generateUniformRandom(nvars,fraction,rnd_eng,minterms);
            clock_t initemps=clock();
            QuineMcCluskey(minterms,result);
            runtime += (clock()-initemps)/(double)CLOCKS_PER_SEC;
        }
        cout << fraction << ";" << (1<<nvars)*fraction << ";" << runtime/nreps << endl;
    }
}
void experiments22(int argc, char ** argv){
    std::default_random_engine rnd_eng(1234);

    double min_fraction = 0.05;
    double max_fraction = 0.85;
    int nvars = 8;
    int nreps = 5;

    for(double fraction = min_fraction; fraction <= max_fraction; fraction+=0.05){
        double runtime = 0;
        for(int rep = 0; rep < nreps; rep++){
            list<Implicant> minterms;
            generateUniformRandom(nvars,fraction,rnd_eng,minterms);
            clock_t initemps=clock();
            set<Implicant> prime;

            tabularMethod(minterms,prime);
            runtime += (clock()-initemps)/(double)CLOCKS_PER_SEC;
        }
        cout << fraction << ";" << (1<<nvars)*fraction << ";" << runtime/nreps << endl;
    }
}

void experiments3(int argc, char ** argv){
    std::default_random_engine rnd_eng(1234);

    double min_fraction = 0.05;
    double max_fraction = 1.01;
    int nvars = 6;
    int nreps = 10;

    for(double fraction = min_fraction; fraction <= max_fraction; fraction+=0.05){
        double nimp = 0;
        for(int rep = 0; rep < nreps; rep++){
            list<Implicant> minterms;
            list<Implicant> result;
            generateUniformRandom(nvars,fraction,rnd_eng,minterms);
            QuineMcCluskey(minterms,result);
            cout << minterms.size() << endl;
            cout << result.size() << endl;
            //cout << result.front() << endl;
            nimp+=result.size();
        }
        cout << fraction << ";" << (1<<nvars)*fraction << ";" << nimp/nreps << endl;
    }
}

int main (int argc, char ** argv) {

    // Parse input arguments for options
    /*// http://www.boost.org/doc/libs/1_41_0/doc/html/program_options/tutorial.html
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("n", po::value<double>(), "number of variables")
        ("seed,s", po::value<int>(), "random seed")
        ("width,w", po::value<int>(), "width of the lattice (only for lattice graphs)")
        ("height,h", po::value<int>(), "height of the lattice (only for lattice graphs)")
        ("p0", po::value<float>(), "fraction of initially infected nodes")
        ("p", po::value<float>(), "connection probability in Erdos-Renyi graph")
        ("n,n", po::value<int>(), "number of nodes (for non-lattice graphs)")
        ("beta", po::value<double>(), "probability of infecting each neighbour at a time")
        ("gamma", po::value<double>(), "probability of recovering at a time"
         "vertex")
        ("tmax,t", po::value<int>(), "amount of time of simulation")
        ("model,m", po::value<string>(), "tye of graph: "
         "complete, lattice, star, er, tree")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }


    if (vm.count("seed")) seed = vm["seed"].as<int>();
    if (vm.count("p0")) p0 = vm["p0"].as<float>();
    if (vm.count("p")) p_er = vm["p"].as<float>();
    if (vm.count("tmax")) tmax = vm["tmax"].as<int>();
    if (vm.count("n")) n = vm["n"].as<int>();
    if (vm.count("width")) w = vm["width"].as<int>();
    if (vm.count("height")) h = vm["height"].as<int>();
    if (vm.count("model")) model = vm["model"].as<string>();
    if (vm.count("beta")) beta = vm["beta"].as<double>();
    if (vm.count("gamma")) gamma = vm["gamma"].as<double>();


    //Check errors
    if(model == "lattice" && (!vm.count("width") || !vm.count("height")) ){
        cout << "w and h must be specifiec for lattice graph" << endl;
        return -1;
    }
    if(model != "lattice" && !vm.count("n")){
        cout << "n must be specifiec for " << model << " graph" << endl;
        return -1;
    }
*/
    //experiments3(argc,argv);

    ifstream input;
    string instance("instances/test.txt");

    input.open(instance);

    string word;

    list<Implicant> minterms;
    while(input >> word)
        minterms.push_back(Implicant(word));

    input.close();

    list<Implicant> result;
    QuineMcCluskey(minterms,result);

    ofstream output;
    output.open(instance + ".res");
    for(const Implicant & im : result)
        output << im << endl;
    output.close();


    /*input.open("instances/testmatriu.txt");

    int nRow, nCol, val;
    input >> nRow >> nCol;

    bool ** mat = new bool *[nRow];
    for(int i = 0; i < nRow; i++){
        mat[i] = new bool[nCol];
        for(int j = 0; j < nCol; j++){
            input >> val;
            mat[i][j] = val ? true : false;
        }
    }


    Matrix m(mat,nRow,nCol);
    m.saveState();
    m.print();
    dynamic_bitset <>x(nCol); x.reset();
    //m.reduce(x);
    m.removeColumnAndRows(9);
    m.print();
    //cout << endl << endl << x << endl;
    m.restoreState();
    m.print();*/
/*
    std::default_random_engine rnd_eng(1234);
    double fraction = 0.3;
    int nvars = 10;
    int nminterms = (int)((1<<nvars)*fraction);

    generateBinomialRandom(nvars,fraction,rnd_eng);
    */
    return 0;
}


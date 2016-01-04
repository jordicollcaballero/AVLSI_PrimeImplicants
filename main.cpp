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
namespace po = boost::program_options;

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
                if(consensus.covers(*it2)){
                    prime.erase(it2++);
                    list<Implicant>::iterator it3 = it1;
                    it3++;
                }
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
dynamic_bitset<> exactCover(Matrix *A, dynamic_bitset<> x, dynamic_bitset<> b, std::default_random_engine & rnd_eng){
    int c;

    A->reduce(x); //Reduce the matrix using essentials and dominance
    if(x.count() >= b.count()) return b; //Bound
    if(A->empty()) return x; //Solution completed

    c = A->selectBranchingColumn(rnd_eng);
    A->saveState();
    x.set(A->removeColumnAndRows(c)); //Include implicant of column c to the solution
    dynamic_bitset<> x2 = exactCover(A,x,b,rnd_eng); //Branch
    A->restoreState();
    if(x2.count() < b.count()) b = x2;

    A->saveState();
    x.reset(A->removeColumn(c)); //Exclude implicant of column c to the solution
    x2 = exactCover(A,x,b,rnd_eng); //Branch
    A->restoreState();
    if(x2.count() < b.count()) b = x2;

    return b;
}

void QuineMcCluskey(const list<Implicant> &minterms, list<Implicant> &result, std::default_random_engine & rnd_eng,
                    double *prime_time, double *cover_time){

    //Prime implicants computation
    clock_t initemps=clock();

    set<Implicant> prime;
    tabularMethod(minterms,prime);

    *prime_time = (clock()-initemps)/(double)CLOCKS_PER_SEC;

    //Minimum cover computation
    initemps=clock();

    dynamic_bitset <>x(prime.size()); x.reset();
    dynamic_bitset <>b(prime.size()); b.set();
    Matrix * m = new Matrix(minterms,prime);
    dynamic_bitset <>minset = exactCover(m,x,b,rnd_eng);
    delete m;

    int i = 0;
    for(const Implicant &im : prime){
        if(minset[i])
            result.push_back(im);
        i++;
    }

    *cover_time = (clock()-initemps)/(double)CLOCKS_PER_SEC;
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

void runIncProb(int nvars, int nreps, int seed){
    std::default_random_engine rnd_eng(seed);
    double min_fraction = 0.05;
    double max_fraction = 1.01;

    double *prime_time = new double(0);
    double *cover_time = new double(0);
    cout << "prob ; n minterms ; prime comp. time ; min cover comp. time ; total time ; sol. size" << endl << endl;

    for(double fraction = min_fraction; fraction <= max_fraction; fraction+=0.05){
        double sum_prime_time = 0;
        double sum_cover_time = 0;
        double sum_sizes = 0;

        for(int rep = 0; rep < nreps; rep++){
            list<Implicant> minterms;
            list<Implicant> result;
            generateUniformRandom(nvars,fraction,rnd_eng,minterms);

            QuineMcCluskey(minterms,result,rnd_eng, prime_time, cover_time);
            sum_prime_time += *prime_time;
            sum_cover_time += *cover_time;
            sum_sizes += result.size();
        }
        cout << fraction << " ; ";
        cout << (1<<nvars)*fraction << " ; ";
        cout << sum_prime_time/nreps << " ; ";
        cout << sum_cover_time/nreps << " ; ";
        cout << sum_prime_time/nreps + sum_cover_time/nreps<< " ; ";
        cout << sum_sizes / nreps << endl;
    }

    delete prime_time;
    delete cover_time;
}

void runIncVars(double fraction, int nreps, int seed){
    std::default_random_engine rnd_eng(seed);

    int min_nvars = 4;
    int max_nvars = 10;

    double *prime_time = new double(0);
    double *cover_time = new double(0);
    cout << "nvars ; n minterms ; prime comp. time ; min cover comp. time ; total time ; sol. size" << endl << endl;

    for(int nvars = min_nvars; nvars <= max_nvars; nvars++){
        double sum_prime_time = 0;
        double sum_cover_time = 0;
        double sum_sizes = 0;

        for(int rep = 0; rep < nreps; rep++){
            list<Implicant> minterms;
            list<Implicant> result;
            generateUniformRandom(nvars,fraction,rnd_eng,minterms);

            QuineMcCluskey(minterms,result,rnd_eng, prime_time, cover_time);
            sum_prime_time += *prime_time;
            sum_cover_time += *cover_time;
            sum_sizes += result.size();
        }
        cout << nvars << " ; ";
        cout << (1<<nvars)*fraction << " ; ";
        cout << sum_prime_time/nreps << " ; ";
        cout << sum_cover_time/nreps << " ; ";
        cout << sum_prime_time/nreps + sum_cover_time/nreps<< " ; ";
        cout << sum_sizes / nreps << endl;
    }

    delete prime_time;
    delete cover_time;
}

void runSingle(int nvars, double fraction, int nreps, int seed){
    std::default_random_engine rnd_eng(seed);

    double *prime_time = new double(0);
    double *cover_time = new double(0);

    double sum_prime_time = 0;
    double sum_cover_time = 0;
    double sum_sizes = 0;
    for(int rep = 0; rep < nreps; rep++){
        list<Implicant> minterms;
        list<Implicant> result;
        generateUniformRandom(nvars,fraction,rnd_eng,minterms);
        QuineMcCluskey(minterms,result,rnd_eng, prime_time, cover_time);
        sum_prime_time += *prime_time;
        sum_cover_time += *cover_time;
        sum_sizes += result.size();
    }
    cout << "N minterms = " << (1<<nvars)*fraction << endl;
    cout << "Primes computation time = " << sum_prime_time/nreps << endl;
    cout << "Min cover computation time = " << sum_cover_time/nreps << endl;
    cout << "Total computation time = " << sum_prime_time/nreps + sum_cover_time/nreps << endl;
    cout << "Min cover size = " << sum_sizes/nreps << endl;

    delete prime_time;
    delete cover_time;
}

int main (int argc, char ** argv) {

    int n;
    double p;
    int seed = 0;
    int r = 1;
    string b = "first";
    string e = "single";

    // Parse input arguments for options
    // http://www.boost.org/doc/libs/1_41_0/doc/html/program_options/tutorial.html
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("seed,s",po::value<int>(), "seed of the random generator. DEFAULT 0")
        ("numvar,n", po::value<int>(), "number of variables of the minterms (only for 'single' and 'incprob')")
        ("prob,p", po::value<double>(), "probability for a minterm of being in the formula (only for 'single' and 'incvars')")
        ("reps,r", po::value<int>(), "number of repetitions. DEFAULT 1")
        ("branch,b", po::value<string>(), "selection strategy of the branching column"
         "\nfirst: first column. DEFAULT"
         "\nrandom: random column"
         "\nmax: column with maximum number of 1s"
         "\nmin: column with minimum number of 1s")
        ("execution,e", po::value<string>(), "type of execution"
         "\nsingle: average of 'r' executions for defined 'n' and 'p'. DEFAULT"
         "\nincvars: average of 'r' repetitions for defined 'p', and 'n'' from 5 to 12"
         "\nincprob: average of 5 repetitions for defined 'n', and 'p'' from 0.05 to 1 step 0.05")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    if(vm.count("seed")) seed = vm["seed"].as<int>();
    if(vm.count("numvar")) n = vm["numvar"].as<int>();
    if(vm.count("prob")) p = vm["prob"].as<double>();
    if(vm.count("reps")) r = vm["reps"].as<int>();
    if(vm.count("branch")) b = vm["branch"].as<string>();
    if(vm.count("execution")) e = vm["execution"].as<string>();

    if(b == "random") Matrix::colSelectCriteria = Matrix::RANDOM;
    else if(b=="max") Matrix::colSelectCriteria = Matrix::MAX_ONES;
    else if(b=="min") Matrix::colSelectCriteria = Matrix::MIN_ONES;
    else Matrix::colSelectCriteria = Matrix::FIRST;

    //Check errors
    if(e == "incprob")
        runIncProb(n,r,seed);
    else if(e == "incvars")
        runIncVars(p,r,seed);
    else
        runSingle(n,p,r,seed);

/*
    ifstream input;
    string instance("instances/test.txt");

    input.open(instance);

    string word;

    list<Implicant> minterms;
    while(input >> word)
        minterms.push_back(Implicant(word));

    input.close();

    list<Implicant> result;
    //QuineMcCluskey(minterms,result);

    ofstream output;
    output.open(instance + ".res");
    for(const Implicant & im : result)
        output << im << endl;
    output.close();
*/

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


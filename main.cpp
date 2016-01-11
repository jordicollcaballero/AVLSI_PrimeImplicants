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

void testResult(const list<Implicant> & f, const list<Implicant>& d, const list<Implicant> &cover);

//Auxiliar method used in developing and testing. Prints the table in tabular method.
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

//Auxiliar method used in developing and testing. Prints the temporal list of prime implicants in the iterated consensus
void printList(const list<Implicant> & primes){
    for(const Implicant & im : primes)
        cout << im << endl;
    cout << endl;
}

/*
 * Tabular method to compute the prime implicants of a function
 * f: function (list of minterms)
 * d: don't care set (list of minterms)
 * primes: out param. List of prime implicants
 */
void tabularMethod(const list<Implicant>& f, const list<Implicant>& d,list<Implicant> &primes){
    if(f.size() == 0)
        return;
    set<Implicant> implicantset;
    int nVars = f.front().getNVars();

    //Construct the distance-1 table
    set<Implicant> ** table = new set<Implicant> * [nVars+1];
    for(int i = 0; i <= nVars; i++)
        table[i] = new set<Implicant> [nVars-i+1];
    for(const Implicant& im : f){
        table[0][im.countNegatedVars()].insert(im);
        implicantset.insert(im);
    }
    for(const Implicant& im : d)
        table[0][im.countNegatedVars()].insert(im);

    //Make distance 1 mergings
    for(int col = 0; col < nVars; col++){
        for(int row = 0; row < nVars-col; row++){
            for(const Implicant &im1 : table[col][row+1]){
                for(const Implicant &im2 : table[col][row]){
                    Implicant consensus = Implicant::distance1Merging(im1,im2);
                    if(consensus.isValid()){
                        table[col+1][row].insert(consensus);
                        if(!consensus.isDontCare()) implicantset.insert(consensus);
                        if(!im1.isDontCare()) implicantset.erase(im1);
                        if(!im2.isDontCare()) implicantset.erase(im2);
                    }
                }
            }
        }
    }

    //Put the solution in the return
    primes.insert(primes.end(),implicantset.begin(),implicantset.end());

    for(int i = 0; i <= nVars; i++)
        delete [] table[i];
    delete [] table;
}


/*
 * Iterated consensus to compute the prime implicants of a function
 * f: function (list of implicants)
 * d: don't care set (list of implicants)
 * primes: out param. List of prime implicants.
 */
void iteratedConsensus(const list<Implicant> &f, const list<Implicant> &d, list<Implicant> & primes){
    if(f.empty())
        return;

    //Insert f and d in the list
    primes = f;
    primes.insert(primes.end(),d.begin(),d.end());

    //Ensure that there aren't implicants in the list covered by other implicants in the list
    list<Implicant>::iterator it1 = primes.begin();
    it1++;
    while(it1!=primes.end()){
        list<Implicant>::iterator it2 = primes.begin();
        while(it2!=it1){
            if((*it1).covers(*it2)){
                it2 = primes.erase(it2);
            }
            else if((*it2).covers(*it1)){
                it1 = primes.erase(it1);
                break;
            }
            else it2++;
        }
        if(it1==it2) it1++;
    }

    //Iterate over all pairs and make the consensus
    it1 = primes.begin();
    it1++;
    while(it1!=primes.end()){
        Implicant imp1 = *it1;
        bool it1displaced = false;
        list<Implicant>::iterator it2 = primes.begin();
        while(it2!=it1){
            Implicant imp2 = *it2;
            Implicant consensus = Implicant::consensus(imp1,imp2);
            bool toAdd = true;
            bool it1deleted = false;
            if(consensus.isValid()){ //If exists consensus
                list<Implicant>::iterator it3 = primes.begin();

                //Ensure that the consensus is not covered, and erase the implicants covered by the consensus
                while(it3 != primes.end()){
                    if((*it3).covers(consensus)){
                        toAdd = false;
                        break;
                    }
                    else if(consensus.covers(*it3)){
                        if(it3 == it1){
                            it1--;
                            it1deleted = true;
                            it1displaced = true;
                        }
                        if(it3 == it2) it2--;
                        it3 = primes.erase(it3);
                    }
                    else it3++;
                }
                //If the consensus is not covered, add it.
                if(toAdd) primes.push_back(consensus);
            }
            if(it1deleted) it1++;
            it2++;
        }
        if(!it1displaced)it1++;
    }

    //Detect don't care implicants. Contrarilly to the tabular method, there are many consensus that give a same implicant.
        //Thus, it is not ensured that all the don't care implicants will be detected and deleted
    it1 = primes.begin();
    while(it1!=primes.end()){
        if((*it1).isDontCare())
            it1 = primes.erase(it1);
        else it1++;
    }
}

/*
 * A: matrix of covering of implicants on minterms
 * x: partial solution
 * b: best solution found until the moment
 * rnd_eng: random engine used for random operations
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

/*
 * Quine McCluskey method.
 * f: function
 * d: don't care set
 * cover: out param with the solution
 * rnd_eng: random engine used for random operations
 * prime_time, cover_time, nprime: out params to get statistics
 * tabular: Method to find the prime implicants. true -> tabular method. false -> iterated consensus.
 */
void QuineMcCluskey(const list<Implicant> &f, const list<Implicant> &d, list<Implicant> &cover,
                    std::default_random_engine & rnd_eng,
                    double *prime_time, double *cover_time, int *nprime, bool tabular){

    //Prime implicants computation
    clock_t initemps=clock();

    list<Implicant> prime;

    if(tabular) tabularMethod(f,d,prime);
    else iteratedConsensus(f,d,prime);

    *prime_time = (clock()-initemps)/(double)CLOCKS_PER_SEC;
    *nprime = prime.size();

    //Minimum cover computation
    initemps=clock();

    dynamic_bitset <>x(prime.size()); x.reset();
    dynamic_bitset <>b(prime.size()); b.set();
    Matrix * m = new Matrix(f,prime);
    dynamic_bitset <>minset = exactCover(m,x,b,rnd_eng);
    delete m;

    int i = 0;
    for(const Implicant &im : prime){
        if(minset[i])
            cover.push_back(im);
        i++;
    }

    *cover_time = (clock()-initemps)/(double)CLOCKS_PER_SEC;
}

/*
 * Generate a random function
 * nvars: number of variables of the function
 * fp: ratio of minterms that are in the function
 * dp: ratio of minterms that are in the don't care set
 * rnd_eng: random engine used for random operations
 * f: out param. function
 * d: out param. don't care set
 */
void generateUniformRandom(int nvars, double fp, double dp, default_random_engine &rnd_eng, list<Implicant> &f, list<Implicant> &d){
    dynamic_bitset<> mask(nvars,(1<<nvars)-1);
    uniform_real_distribution<double> uir(0.0, 1.0);
    for(int i = 0; i < 1<<nvars; i++){
        double q = uir(rnd_eng);
        if(q < fp)
            f.push_back(Implicant(mask,dynamic_bitset<> (nvars,i),false));
        else if(q < fp + dp)
            d.push_back(Implicant(mask,dynamic_bitset<> (nvars,i),true));
    }
}

/*
 * Experiment with increasing ratio of minterms in the function for a fixed number of vars
 * nvars: number of variables of the function
 * dontcare: ratio of minterms that are in the don't care set
 * nreps: number of repetitions of each experiment
 * seed: seed for the random generator
 * tabular: Method to find the prime implicants. true -> tabular method. false -> iterated consensus.
 */
void runIncProb(int nvars, double dontcare, int nreps, int seed, bool tabular){
    std::default_random_engine rnd_eng(seed);
    double min_fraction = 0.05;
    double max_fraction = 1.01;

    double *prime_time = new double(0);
    double *cover_time = new double(0);
    int *nprimes = new int(0);
    cout << "prob ; n minterms f ; n minterms d ; prime comp. time ; min cover comp. time ; total time ; n primes; sol size" << endl;

    for(double fraction = min_fraction; fraction <= max_fraction; fraction+=0.05){
        double sum_prime_time = 0;
        double sum_cover_time = 0;
        int sum_primes = 0;
        double sum_sizes = 0;

        for(int rep = 0; rep < nreps; rep++){
            list<Implicant> f;
            list<Implicant> d;
            list<Implicant> cover;
            generateUniformRandom(nvars,fraction,dontcare,rnd_eng,f,d);

            QuineMcCluskey(f,d,cover,rnd_eng, prime_time, cover_time, nprimes, tabular);
            sum_prime_time += *prime_time;
            sum_cover_time += *cover_time;
            sum_primes += *nprimes;
            sum_sizes += cover.size();
        }
        cout << fraction << " ; ";
        cout << (1<<nvars)*fraction << " ; ";
        cout << (1<<nvars)*dontcare << " ; ";
        cout << sum_prime_time/nreps << " ; ";
        cout << sum_cover_time/nreps << " ; ";
        cout << sum_prime_time/nreps + sum_cover_time/nreps<< " ; ";
        cout << sum_primes/nreps << " ; ";
        cout << sum_sizes / nreps << endl;
    }

    delete prime_time;
    delete cover_time;
}

/*
 * Experiment with increasing number of variables for a fixed ratio of minterms in the function
 * fraction: ratio of minterms that are in the function
 * dontcare: ratio of minterms that are in the don't care set
 * nreps: number of repetitions of each experiment
 * seed: seed for the random generator
 * tabular: Method to find the prime implicants. true -> tabular method. false -> iterated consensus.
 */
void runIncVars(double fraction, double dontcare, int nreps, int seed, bool tabular){
    std::default_random_engine rnd_eng(seed);

    int min_nvars = 4;
    int max_nvars = 10;

    double *prime_time = new double(0);
    double *cover_time = new double(0);
    int *nprimes = new int(0);
    cout << "nvars ; n minterms f ; n minterms d ; prime comp time ; min cover comp time ; total time ; n primes; sol size" << endl;

    for(int nvars = min_nvars; nvars <= max_nvars; nvars++){
        double sum_prime_time = 0;
        double sum_cover_time = 0;
        int sum_primes = 0;
        double sum_sizes = 0;

        for(int rep = 0; rep < nreps; rep++){
            list<Implicant> f;
            list<Implicant> d;
            list<Implicant> cover;
            generateUniformRandom(nvars,fraction,dontcare,rnd_eng,f,d);

            QuineMcCluskey(f,d,cover,rnd_eng, prime_time, cover_time, nprimes, tabular);
            sum_prime_time += *prime_time;
            sum_cover_time += *cover_time;
            sum_primes += *nprimes;
            sum_sizes += cover.size();
        }
        cout << nvars << " ; ";
        cout << (1<<nvars)*fraction << " ; ";
        cout << (1<<nvars)*dontcare << " ; ";
        cout << sum_prime_time/nreps << " ; ";
        cout << sum_cover_time/nreps << " ; ";
        cout << sum_prime_time/nreps + sum_cover_time/nreps<< " ; ";
        cout << sum_primes/nreps << " ; ";
        cout << sum_sizes / nreps << endl;
    }

    delete prime_time;
    delete cover_time;
}

/*
 * Experiment with fixed number of variables and fixed ratio of minterms in the function
 * nvars: number of variables of the function
 * fraction: ratio of minterms that are in the function
 * dontcare: ratio of minterms that are in the don't care set
 * nreps: number of repetitions of each experiment
 * seed: seed for the random generator
 * tabular: Method to find the prime implicants. true -> tabular method. false -> iterated consensus.
 */
void runSingle(int nvars, double fraction, double dontcare, int nreps, int seed, bool tabular){
    std::default_random_engine rnd_eng(seed);

    double *prime_time = new double(0);
    double *cover_time = new double(0);
    int *nprimes = new int(0);

    double sum_prime_time = 0;
    double sum_cover_time = 0;
    int sum_primes = 0;
    double sum_sizes = 0;
    for(int rep = 0; rep < nreps; rep++){
        list<Implicant> f;
        list<Implicant> d;
        list<Implicant> cover;
        generateUniformRandom(nvars,fraction,dontcare,rnd_eng,f,d);
        QuineMcCluskey(f,d,cover,rnd_eng, prime_time, cover_time, nprimes, tabular);
        sum_prime_time += *prime_time;
        sum_cover_time += *cover_time;
        sum_primes = *nprimes;
        sum_sizes += cover.size();
    }
    cout << "N minterms in f = " << (1<<nvars)*fraction << endl;
    cout << "N minterms in d = " << (1<<nvars)*dontcare << endl;
    cout << "Primes computation time = " << sum_prime_time/nreps << endl;
    cout << "Min cover computation time = " << sum_cover_time/nreps << endl;
    cout << "Total computation time = " << sum_prime_time/nreps + sum_cover_time/nreps << endl;
    cout << "Number of prime implicants = " << sum_primes/nreps << endl;
    cout << "Min cover size = " << sum_sizes/nreps << endl;

    delete prime_time;
    delete cover_time;
}

//Auxiliar method used in developing and testing. Checks the validity of the solutions.
void testResult(const list<Implicant> & f, const list<Implicant>& d, const list<Implicant> &cover){
    set<Implicant> fset; fset.insert(f.begin(),f.end());
    set<Implicant> dset; dset.insert(d.begin(),d.end());

    bool allCovered = true;
    bool onlyNeededCovered = true;
    int nvars = f.front().getNVars();
    dynamic_bitset<> mask(nvars,(1<<nvars)-1);

    for(int i = 0; i < 1<<nvars; i++){
        Implicant im1(Implicant(mask,dynamic_bitset<> (nvars,i)));
        if(fset.find(im1)!=fset.end()){
            bool covered = false;
            for(const Implicant & im2 : cover){
                if(im2.covers(im1)){
                    covered = true;
                    break;
                }
            }
            if(!covered) allCovered = false;
        }
        else if(dset.find(im1)== dset.end()){
            bool covered = false;
            for(const Implicant & im2 : cover){
                if(im2.covers(im1)){
                    covered = true;
                    break;
                }
            }
            if(covered) onlyNeededCovered = false;
        }
    }
    cout << "All the function is covered: " << (allCovered ? "yes" : "no") << endl;
    cout << "Only the function and don't care set are covered: " << (onlyNeededCovered ? "yes" : "no") << endl;
}

//Auxiliar method used in developing and testing. Checks the results of tabular method and iterated consensus.
//The result should be the same if there isn't don't care set
//If there is don't care set, the iterated consensus should give all the prime implicants that gives the tabular method,
    //and could also give prime implicants that only cover members of the don't care set
void testTabularVsIterated(int nvars, double prob, double dontcare, int seed){

    list<Implicant> f,d;
    std::default_random_engine rnd_eng(seed);
    generateUniformRandom(nvars,prob,dontcare,rnd_eng,f,d);

    list<Implicant> primesTabular;
    tabularMethod(f,d,primesTabular);
    cout << "With tabular method: " << primesTabular.size() << " prime implicants" << endl;
    for(const Implicant & im : primesTabular)
        cout << im << endl;
    cout << endl;

    list<Implicant> primesConsensus;
    iteratedConsensus(f,d,primesConsensus);
    primesConsensus.sort();
    cout << "With iterated consensus: " << primesConsensus.size() << " prime implicants" << endl;
    for(const Implicant & im : primesConsensus)
        cout << im << endl;
    cout << endl;
}

//Main methods. Arguments can be checked with '-h' / '--help'
int main (int argc, char ** argv) {

    int n;
    double p;
    double d = 0;
    int seed = 0;
    int r = 1;
    string a = "tabular";
    string b = "first";
    string e = "single";

    // Parse input arguments for options
    // http://www.boost.org/doc/libs/1_41_0/doc/html/program_options/tutorial.html
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("seed,s",po::value<int>(), "seed of the random generator. DEFAULT 0")
        ("numvar,n", po::value<int>(), "number of variables of the minterms. Required in: single, incprob")
        ("prob,p", po::value<double>(), "proportion of the minterms in the formula. Required in: single, incvars")
        ("dontcare,d", po::value<double>(), "proportion of the minterms in the don't care set. DEFAULT 0. Used in: single, incvars, incprob.")
        ("reps,r", po::value<int>(), "number of repetitions. DEFAULT 1. Used in: single, incvars, incprob")
        ("primealg,a", po::value<string>(), "algorithm for computing prime implicants"
         "\ntabular: tabular method. DEFAULT"
         "\niterated: iterated consensus")
        ("branch,b", po::value<string>(), "selection strategy of the branching column. Used in: single, incvars, incprob."
         "\nfirst: first column. DEFAULT"
         "\nrandom: random column"
         "\nmax: column with maximum number of 1s"
         "\nmin: column with minimum number of 1s")
        ("execution,e", po::value<string>(), "type of execution"
         "\nsingle: average of 'r' executions for defined 'n','p' and 'd'. DEFAULT"
         "\nincvars: average of 'r' repetitions for defined 'p' and 'd', and 'n' from 4 to 10"
         "\nincprob: average of 5 repetitions for defined 'n' and 'd', and 'p' from 0.05 to 1 step 0.05")
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
    if(vm.count("dontcare")) d = vm["dontcare"].as<double>();
    if(vm.count("reps")) r = vm["reps"].as<int>();
    if(vm.count("primealg")) a = vm["primealg"].as<string>();
    if(vm.count("branch")) b = vm["branch"].as<string>();
    if(vm.count("execution")) e = vm["execution"].as<string>();

    if(b == "random") Matrix::colSelectCriteria = Matrix::RANDOM;
    else if(b=="max") Matrix::colSelectCriteria = Matrix::MAX_ONES;
    else if(b=="min") Matrix::colSelectCriteria = Matrix::MIN_ONES;
    else Matrix::colSelectCriteria = Matrix::FIRST;

    if(e == "incprob")
        runIncProb(n,d,r,seed,a == "iterated" ? false : true);
    else if(e == "incvars")
        runIncVars(p,d,r,seed,a == "iterated" ? false : true);
    else
        runSingle(n,p,d,r,seed,a == "iterated" ? false : true);

    return 0;
}

/*ifstream input;
string instance("instances/test.txt");

input.open(instance);

string word;

list<Implicant> minterms;
while(input >> word)
    minterms.push_back(Implicant(word));

input.close();*/

/*   ofstream output;
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

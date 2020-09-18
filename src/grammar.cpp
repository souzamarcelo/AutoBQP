#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace boost;
using namespace std;

void search();
string determineTenure();
string determineMaxSteps();
string determineMaxStagnate();

string preamble = "";
string body = "";
unsigned level = 0;
unsigned level_perturbation = 1000;

mt19937 randomNumberGenerator;
unsigned setupRandom(unsigned seed) {
    if(seed == 0) {
        seed = time(0);
        ifstream f("/dev/urandom");
        if(f.good()) {
            f.read((char *) (&seed), sizeof(unsigned int));
        }
    }
    randomNumberGenerator.seed(seed);
    srand48(seed);
    return seed;
}

int getRandomUnsigned(unsigned range) {
    return rand() % range;
}

enum {
    ALG = 0,
    SEARCH,
    MODIFICATION,
    TS,
    PERT,
    STEP,
    CONSTRUCTION,
    CONSTRUCTOR
};

struct CodonDecider {
    vector<unsigned> codons;
    unsigned actual;
    
    vector<unsigned> cAlg;
    vector<unsigned> cSearch;
    vector<unsigned> cModification;
    vector<unsigned> cTs;
    vector<unsigned> cPert;
    vector<unsigned> cStep;
    vector<unsigned> cConstruction;
    vector<unsigned> cConstructor;
    
    unsigned nextAlg = 0;
    unsigned nextSearch = 0;
    unsigned nextModification = 0;
    unsigned nextTs = 0;
    unsigned nextPert = 0;
    unsigned nextStep = 0;
    unsigned nextConstruction = 0;
    unsigned nextConstructor = 0;
    
    unsigned decide(unsigned options, unsigned from) {
        
        if(options == 1) return 0;
        
        switch (from) {
            case ALG:
                if(nextAlg >= cAlg.size())
                    return invalidGrammar(options);
                nextAlg++;
                return cAlg[nextAlg - 1] % options;
            
            case SEARCH:
                if(nextSearch >= cSearch.size())
                    return invalidGrammar(options);
                nextSearch++;
                return cSearch[nextSearch - 1] % options;
            
            case MODIFICATION:
                if(nextModification >= cModification.size())
                    return invalidGrammar(options);
                nextModification++;
                return cModification[nextModification - 1] % options;
            
            case TS:
                if(nextTs >= cTs.size())
                    return invalidGrammar(options);
                nextTs++;
                return cTs[nextTs - 1] % options;
            
            case PERT:
                if(nextPert >= cPert.size())
                    return invalidGrammar(options);
                nextPert++;
                return cPert[nextPert - 1] % options;
            
            case STEP:
                if(nextStep >= cStep.size())
                    return invalidGrammar(options);
                nextStep++;
                return cStep[nextStep - 1] % options;
            
            case CONSTRUCTION:
                if(nextConstruction >= cConstruction.size())
                    return invalidGrammar(options);
                nextConstruction++;
                return cConstruction[nextConstruction - 1] % options;
            
            case CONSTRUCTOR:
                if(nextConstructor >= cConstructor.size())
                    return invalidGrammar(options);
                nextConstructor++;
                return cConstructor[nextConstructor - 1] % options;
            
            default:
                cout << "Error on CODON DECIDER" << endl;
                return invalidGrammar(options);
        }
    }
    
    unsigned invalidGrammar(unsigned options) {
        srand(setupRandom(0));
        return getRandomUnsigned(options);
    }
    
    unsigned invalidGrammar() {
        return 0;
    }
    
    void setCodons(vector<unsigned> input_codons) {
        codons = input_codons;
        actual = -1;
    }
    
    void setSpecificCodons(vector<unsigned> inputAlg,
                           vector<unsigned> inputSearch,
                           vector<unsigned> inputModification,
                           vector<unsigned> inputTs,
                           vector<unsigned> inputPert,
                           vector<unsigned> inputStep,
                           vector<unsigned> inputConstruction,
                           vector<unsigned> inputConstructor) {
        
        cAlg = inputAlg;
        cSearch = inputSearch;
        cModification = inputModification;
        cTs = inputTs;
        cPert = inputPert;
        cStep = inputStep;
        cConstruction = inputConstruction;
        cConstructor = inputConstructor;
    }
    
} codon_decider;

string pTenureType = "<not_defined>";
int pTenureValue = -1;
int pTenurePercent = -1;
int pTenureDivision = -1;
int pTenureVariation = -1;
int pMum = -1;
int pCm = -1;

string pMaxStagnateType = "<not_defined>";
int pMaxStagnateValue = -1;
int pIc = -1;

string pMaxStepsType = "<not_defined>";
int pMaxStepsValue = -1;

double pP = -1;
int pD1 = -1;
int pD2 = -1;
int pB = -1;
int pR = -1;
double pGamma = -1;
int pGammam = -1;
double pBeta = -1;
double pAlpha = -1;
int pN = -1;
int pBsize = -1;
double pLambda = -1;
int pSiPartialSize = -1;

void modification_ls() {
    unsigned option = codon_decider.decide(7, MODIFICATION);
    
    switch (option) {
        
        case 0:
            body += "fi";
            break;
        
        case 1:
            preamble += "FiRR firr;\n";
            body += "firr";
            break;
        
        case 2:
            body += "bi";
            break;
        
        case 3:
            body += "ri";
            break;
        
        case 4:
            body += "si";
            break;
        
        case 5:
            preamble += "SiPartial sipartial(false, " + to_string(pSiPartialSize) + ");\n";
            body += "sipartial";
            break;
        
        case 6:
            preamble += "SiPartial sipartial(true, " + to_string(pSiPartialSize) + ");\n";
            body += "sipartial";
            break;
        
        default:
            cout << "Error on LOCAL SEARCH MODIFICATION" << endl;
            break;
    }
}

void modification_recombination() {
    
    unsigned option = codon_decider.decide(7, MODIFICATION);
    switch (option) {
        
        case 0:
            preamble += "recombine::Fn fn(false);\n";
            body += "fn";
            break;
        
        case 1:
            preamble += "recombine::Fn fn(true);\n";
            body += "fn";
            break;
        
        case 2:
            body += "recombine::bn";
            break;
        
        case 3:
            body += "recombine::rn";
            break;
        
        case 4:
            body += "recombine::si";
            break;
        
        case 5:
            preamble += "recombine::SiPartial sipartial_recombiner(false, " + to_string(pSiPartialSize) + ");\n";
            body += "sipartial_recombiner";
            break;
        
        case 6:
            preamble += "recombine::SiPartial sipartial_recombiner(true, " + to_string(pSiPartialSize) + ");\n";
            body += "sipartial_recombiner";
            break;
        
        default:
            cout << "Error on RECOMBINATION MODIFICATION" << endl;
            break;
    }
}

enum {
    ls = 0, recomb
};

void modification(unsigned from) {
    
    if(from == ls)
        modification_ls();
    if(from == recomb)
        modification_recombination();
}

void local_search() {
    body += "localsearch(S, ";
    modification(ls);
    body += ", startTime, target);\n";
}

void non_monotone_local_search() {
    preamble += "double p = " + to_string(pP) + ";\n";
    preamble += "RII rii(p);\n";
    
    body += "localsearch(S, [&](const Solution& S) { return rii.blmr(S, ";
    modification(ls);
    body += "); }, startTime, target);\n";
}

void tabu_search() {
    
    preamble += "unsigned maxstagnate = " + determineMaxStagnate() + ";\n";
    preamble += "unsigned maxsteps = " + determineMaxSteps() + ";\n";
    
    unsigned option = codon_decider.decide(2, TS);
    switch (option) {
        
        case 0:
            preamble += "BTS bt;\n";
            body += "tabusearch(S, bt, [&]() { return " + determineTenure() +
                    "; }, startTime, target, maxstagnate, maxsteps);\n";
            break;
        
        case 1:
            preamble += "double p = " + to_string(pP) + ";\n";
            preamble += "BTR bt(p);\n";
            body += "tabusearch(S, bt, [&]() { return " + determineTenure() +
                    "; }, startTime, target, maxstagnate, maxsteps);\n";
            break;
        
        default:
            cout << "Error on TABU SEARCH" << endl;
            break;
    }
}

void step_generator() {
    unsigned option = codon_decider.decide(4, STEP);
    
    switch (option) {
        case 0: {
            preamble +=
                    "StepGenerator sg" + to_string(level) + "(" + to_string(pD1) + ", I.n>=" + to_string(pD1) + "*" +
                    to_string(pD2) + "?I.n/" + to_string(pD2) + ":" + to_string(pD1) + ");\n";
            body += "sg" + to_string(level);
            break;
        }
        
        case 1: {
            preamble += "StepGaussian sg" + to_string(level) + "(" + to_string(pD1) + ", I.n>=" + to_string(pD1) + "*" +
                        to_string(pD2) + "?I.n/" + to_string(pD2) + ":" + to_string(pD1) + ", I.n);\n";
            body += "sg" + to_string(level);
            break;
        }
        
        case 2: {
            preamble +=
                    "StepExponential sg" + to_string(level) + "(" + to_string(pD1) + ", I.n>=" + to_string(pD1) + "*" +
                    to_string(pD2) + "?I.n/" + to_string(pD2) + ":" + to_string(pD1) + ", I.n);\n";
            body += "sg" + to_string(level);
            break;
        }
        
        case 3:
            body += "[&]() { return I.n / " + to_string(pGammam) + "; }";
            break;
        
        default:
            cout << "Error on STEP GENERATOR" << endl;
            break;
    }
}

void perturbation(bool ilse) {
    body += "[&](Solution& S) { return ";
    
    unsigned option = codon_decider.decide(3, PERT);
    switch (option) {
        case 0:
            body += "perturbRandom(S, ";
            step_generator();
            body += ", " + to_string(pB) + "); },\n";
            break;
        
        case 1:
            body += "perturbLeastLoss(S, ";
            step_generator();
            body += ", " + to_string(pB) + "); },\n";
            break;
        
        case 2:
            unsigned level_final;
            
            if(!ilse) {
                level_perturbation++;
                level_final = level_perturbation;
                preamble += "Elite e" + to_string(level_final) + "(I," + to_string(pR) + ");\n";
            } else {
                level_final = level;
            }
            
            if(preamble.find("rank_proportional") == std::string::npos) {
                preamble += "vector<double> w(I.n,0);";
                preamble += "for(unsigned i=0; i<I.n; i++)";
                preamble += "  w[i]=pow(i+1," + to_string(pLambda) + ");";
                preamble += "rank_proportional = std::discrete_distribution<>(w.begin(),w.end());";
            }
            
            body += "perturbDiversity(S, ";
            step_generator();
            body += ", e" + to_string(level_final) + ", " + to_string(pBeta) + "); },\n";
            break;
        
        default:
            cout << "Error on PERTURBATION" << endl;
            break;
    }
}

void iterated_local_search() {
    level++;
    
    body += "iterated_localsearch(S, [&](Solution& S) { return \n";
    search();
    body += " },\n";
    perturbation(false);
    body += "startTime, target);\n";
    
    level--;
}

void iterated_local_search_elite() {
    level++;
    
    preamble += "Elite e" + to_string(level) + "(I, " + to_string(pR) + ");\n";
    
    body += "iterated_localsearch_elite(S, [&](Solution& S) { return \n";
    search();
    body += " },\n";
    perturbation(true);
    body += "e" + to_string(level) + ", startTime, target);\n";
    
    level--;
}

void recombination() {
    preamble += "double gamma = " + to_string(pGamma) + ";\n";
    
    body += "[&](const Instance& I, Solution& S, const Solution& T) { return recombine::recombine(I, S, T,";
    modification(recomb);
    body += ", gamma); },";
}

void repeated_elite_recombination() {
    preamble += "unsigned bsize = " + to_string(pBsize) + ";\n";
    
    body += "recombine::recombiner(I, S, bsize, ";
    recombination();
    body += "[&](Solution& S) { return ";
    search();
    body += "},";
    body += "startTime, target);";
}

void constructor() {
    unsigned option = codon_decider.decide(2, CONSTRUCTOR);
    
    switch (option) {
        case 0:
            body += "construct0";
            break;
        case 1:
            body += "constructHalf";
            break;
        default:
            cout << "Error on CONSTRUCTOR" << endl;
            break;
    }
}

void construction() {
    preamble += "double alpha = " + to_string(pAlpha) + ";\n";
    preamble += "unsigned N = " + to_string(pN) + ";\n";
    
    unsigned option = codon_decider.decide(2, CONSTRUCTION);
    switch (option) {
        
        case 0:
            preamble += "ofstream o(\"\");";
            body += "gra(o, I, S, ";
            constructor();
            body += ", N, alpha, startTime);\n";
            break;
        
        case 1:
            body += "grasp(I, S, ";
            constructor();
            body += ", N, alpha, ";
            body += "[&](Solution& S) { return ";
            local_search();
            body += " }, startTime);\n";
            break;
        
        default:
            cout << "Error on CONSTRUCT" << endl;
            break;
    }
}

void search() {
    unsigned option = codon_decider.decide(5, SEARCH);
    
    switch (option) {
        
        case 0:
            local_search();
            break;
        
        case 1:
            non_monotone_local_search();
            break;
        
        case 2:
            tabu_search();
            break;
        
        case 3:
            iterated_local_search();
            break;
        
        case 4:
            iterated_local_search_elite();
            break;
        
        default:
            cout << "Error on SEARCH" << endl;
            break;
    }
}

void algorithm() {
    unsigned option = codon_decider.decide(3, ALG);
    
    switch (option) {
        case 0:
            search();
            break;
        
        case 1:
            construction();
            break;
        
        case 2:
            repeated_elite_recombination();
            break;
        
        default:
            cout << "Error on ALGORITHM" << endl;
            break;
    }
}

void grammar() {
    preamble += "void run(Solution& S, Instance& I, double target, chrono::system_clock::time_point startTime) {\n";
    algorithm();
    body += "}\n";
    
    cout << preamble << endl << body << endl;
}

int main(int argc, char *argv[]) {
    
    po::options_description desc("Options");
    desc.add_options()
            ("help", "show help")
            ("seed", po::value<unsigned>(), "random seed")
            ("codon", po::value<vector<unsigned> >(), "general codon values")
            ("cAlg", po::value<vector<unsigned> >(), "ALG codon values")
            ("cSearch", po::value<vector<unsigned> >(), "SEARCH codon values")
            ("cModification", po::value<vector<unsigned> >(), "MODIFICATION codon values")
            ("cTs", po::value<vector<unsigned> >(), "TS codon values")
            ("cPert", po::value<vector<unsigned> >(), "PERTURBATION codon values")
            ("cStep", po::value<vector<unsigned> >(), "STEP GENERATOR codon values")
            ("cConstruction", po::value<vector<unsigned> >(), "CONSTRUCTION codon values")
            ("cConstructor", po::value<vector<unsigned> >(), "CONSTRUCTOR codon values")
            ("pTenureType", po::value<string>(), "Parameter TENURE type")
            ("pTenureValue", po::value<int>(), "Parameter TENURE value")
            ("pTenurePercent", po::value<int>(), "Parameter TENURE percent")
            ("pTenureDivision", po::value<int>(), "Parameter TENURE division")
            ("pTenureVariation", po::value<int>(), "Parameter TENURE variation")
            ("pCm", po::value<int>(), "Parameter CM")
            ("pMaxStagnateType", po::value<string>(), "Parameter MAX STAGNATE type")
            ("pMaxStagnateValue", po::value<int>(), "Parameter MAX STAGNATE value")
            ("pMum", po::value<int>(), "Parameter MUM")
            ("pIc", po::value<int>(), "Parameter IC")
            ("pMaxStepsType", po::value<string>(), "Parameter MAX STEPS type")
            ("pMaxStepsValue", po::value<int>(), "Parameter MAX STEPS value")
            ("pP", po::value<double>(), "Parameter P")
            ("pD1", po::value<int>(), "Parameter D1")
            ("pD2", po::value<int>(), "Parameter D2")
            ("pB", po::value<int>(), "Parameter B")
            ("pR", po::value<int>(), "Parameter R")
            ("pGamma", po::value<double>(), "Parameter GAMMA")
            ("pGammam", po::value<int>(), "Parameter GAMMAM")
            ("pBeta", po::value<double>(), "Parameter BETA")
            ("pAlpha", po::value<double>(), "Parameter ALPHA")
            ("pN", po::value<int>(), "Parameter N")
            ("pBsize", po::value<int>(), "Parameter BSIZE")
            ("pLambda", po::value<double>(), "Parameter LAMBDA")
            ("pSiPartialSize", po::value<int>(), "Parameter SI PARTIAL SIZE");
    
    po::positional_options_description pod;
    pod.add("codon", -1);
    pod.add("cAlg", -1);
    pod.add("cSearch", -1);
    pod.add("cModification", -1);
    pod.add("cTs", -1);
    pod.add("cPert", -1);
    pod.add("cStep", -1);
    pod.add("cConstruction", -1);
    pod.add("cConstructor", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
    po::notify(vm);
    
    if(vm.count("help")) {
        cout << desc << endl;
        return 0;
    }
    
    if(vm.count("seed")) {
        setupRandom(vm["seed"].as<unsigned>());
    } else {
        setupRandom(0);
    }
    
    vector<unsigned> inputAlg;
    vector<unsigned> inputSearch;
    vector<unsigned> inputModification;
    vector<unsigned> inputTs;
    vector<unsigned> inputPert;
    vector<unsigned> inputStep;
    vector<unsigned> inputConstruction;
    vector<unsigned> inputConstructor;
    
    if(vm.count("cAlg")) inputAlg = vm["cAlg"].as<vector<unsigned> >();
    if(vm.count("cSearch")) inputSearch = vm["cSearch"].as<vector<unsigned> >();
    if(vm.count("cModification")) inputModification = vm["cModification"].as<vector<unsigned> >();
    if(vm.count("cTs")) inputTs = vm["cTs"].as<vector<unsigned> >();
    if(vm.count("cPert")) inputPert = vm["cPert"].as<vector<unsigned> >();
    if(vm.count("cStep")) inputStep = vm["cStep"].as<vector<unsigned> >();
    if(vm.count("cConstruction")) inputConstruction = vm["cConstruction"].as<vector<unsigned> >();
    if(vm.count("cConstructor")) inputConstructor = vm["cConstructor"].as<vector<unsigned> >();
    
    codon_decider.setSpecificCodons(inputAlg,
                                    inputSearch,
                                    inputModification,
                                    inputTs,
                                    inputPert,
                                    inputStep,
                                    inputConstruction,
                                    inputConstructor);
    
    if(vm.count("pTenureType")) pTenureType = vm["pTenureType"].as<string>();
    if(vm.count("pTenureValue")) pTenureValue = vm["pTenureValue"].as<int>();
    if(vm.count("pTenurePercent")) pTenurePercent = vm["pTenurePercent"].as<int>();
    if(vm.count("pTenureDivision")) pTenureDivision = vm["pTenureDivision"].as<int>();
    if(vm.count("pTenureVariation")) pTenureVariation = vm["pTenureVariation"].as<int>();
    if(vm.count("pCm")) pCm = vm["pCm"].as<int>();
    if(vm.count("pMaxStagnateType")) pMaxStagnateType = vm["pMaxStagnateType"].as<string>();
    if(vm.count("pMaxStagnateValue")) pMaxStagnateValue = vm["pMaxStagnateValue"].as<int>();
    if(vm.count("pMum")) pMum = vm["pMum"].as<int>();
    if(vm.count("pIc")) pIc = vm["pIc"].as<int>();
    if(vm.count("pMaxStepsType")) pMaxStepsType = vm["pMaxStepsType"].as<string>();
    if(vm.count("pMaxStepsValue")) pMaxStepsValue = vm["pMaxStepsValue"].as<int>();
    if(vm.count("pP")) pP = vm["pP"].as<double>();
    if(vm.count("pD1")) pD1 = vm["pD1"].as<int>();
    if(vm.count("pD2")) pD2 = vm["pD2"].as<int>();
    if(vm.count("pB")) pB = vm["pB"].as<int>();
    if(vm.count("pR")) pR = vm["pR"].as<int>();
    if(vm.count("pGamma")) pGamma = vm["pGamma"].as<double>();
    if(vm.count("pGammam")) pGammam = vm["pGammam"].as<int>();
    if(vm.count("pBeta")) pBeta = vm["pBeta"].as<double>();
    if(vm.count("pTenureVariation")) pTenureVariation = vm["pTenureVariation"].as<int>();
    if(vm.count("pAlpha")) pAlpha = vm["pAlpha"].as<double>();
    if(vm.count("pN")) pN = vm["pN"].as<int>();
    if(vm.count("pBsize")) pBsize = vm["pBsize"].as<int>();
    if(vm.count("pLambda")) pLambda = vm["pLambda"].as<double>();
    if(vm.count("pSiPartialSize")) pSiPartialSize = vm["pSiPartialSize"].as<int>();
    
    grammar();
    
    return 0;
}

string determineTenure() {
    string value;

    switch (pTenureType.back()) {

        case 'a':
            preamble += "double davg = double(accumulate(I.deg.begin(),I.deg.end(),0))/I.n;";
            value = "max(1u, unsigned(davg+1))";
            break;

        case 'v':
            value = "I.n";
            break;

        case 'p': {
            string t = to_string(pTenurePercent);
            t.pop_back();
            value = "max(1.0," + t + " * I.n / 100.0)";
            break;
        }

        case 'c': {
            string t = to_string(pTenureDivision);
            t.pop_back();
            value = "I.n / " + t;
            break;
        }

        case 'r': {
            value = "(I.n / " + to_string(pCm) + ") + uniform_int_distribution<>(1, " + to_string(pTenureVariation) +
                    ")(rng)";
            break;
        }

        case 'n':
            value = "max(1, " + to_string(pTenureValue) + ")";
            break;

        default:
            value = "<undefined>";
            break;
    }

    return value;
}

string determineMaxStagnate() {
    string value;

    switch (pMaxStagnateType.back()) {

        case 'a':
            value = to_string(pIc) + " * I.n";
            break;

        case 'i':
            value = "numeric_limits<unsigned>::max()";
            break;

        case 'm':
            value = to_string(pMum) + " * I.n";
            break;

        case 'n':
            value = "max(1," + to_string(pMaxStagnateValue) + ")";
            break;

        default:
            value = "<undefined>";
            break;
    }

    return value;
}

string determineMaxSteps() {
    string value;

    switch (pMaxStepsType.back()) {

        case 'a':
            value = "I.n>5000?15000:(I.n>3000?12000:10000)";
            break;

        case 'i':
            value = "numeric_limits<unsigned>::max()";
            break;

        case 'n':
            value = "max(1," + to_string(pMaxStepsValue) + ")";
            break;

        default:
            value = "<undefined>";
            break;
    }

    return value;
}
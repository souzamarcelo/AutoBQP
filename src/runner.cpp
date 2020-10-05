#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <chrono>
#include <iomanip>
#include <random>
#include <boost/multi_array.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;
using namespace boost;
using namespace std;

unsigned timeLimit;
string format;
chrono::system_clock::time_point startTime;

#include "rng.hpp"
#include "instance.hpp"
#include "solution.hpp"
#include "report.hpp"
Report report;

#include "library/elite.hpp"
#include "library/localsearch.hpp"
#include "library/tabusearch.hpp"
#include "library/construct.hpp"
#include "library/recombine.hpp"
#include "algorithm.cpp"

int main(int argc, char *argv[]) {
    po::options_description desc("General options");
    desc.add_options()
            ("help", "show help")
            ("ins", po::value<string>(), "instance")
            ("seed", po::value<unsigned>()->default_value(0), "random seed")
            ("runtime", po::value<unsigned>()->default_value(300), "running time limit (s)")
            ("maximize", po::value<bool>()->default_value(false), "multiply weight by -1")
            ("format", po::value<string>()->default_value("sparse"), "instance format (sparse, pgen, maxcut, tap...)")
            ("target", po::value<double>()->default_value(numeric_limits<int>::min()), "target value")
            ;

    po::positional_options_description pod;
    pod.add("ins", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
    po::notify(vm);

    cout << fixed;

    if (vm.count("help") || !vm.count("ins")) {
        cout << desc << endl;
        return 0;
    }

    unsigned seed = vm["seed"].as<unsigned>();
    setupRandom(seed);

    Instance I;
    string instanceName = vm["ins"].as<string>();
    ifstream ins(instanceName);
    if (!ins.good()) {
        cout << "Can't open instance file " << instanceName << endl;
        return 1;
    }
    double target = vm["target"].as<double>();
    timeLimit = vm["runtime"].as<unsigned>();
    bool maximize = vm["maximize"].as<bool>();
    string format = vm["format"].as<string>();

    I.readInstance(ins, maximize, format);
    Solution S(I);
    startTime = chrono::system_clock::now();

    run(S, I, target, startTime);
    if(format == "tap") S.value += (I.btb * I.P);
    cout << S.value << endl;
}

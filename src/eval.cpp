#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include <cassert>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <chrono>
#include <iterator>
#include <iomanip>
#include <boost/multi_array.hpp>
#include <boost/program_options.hpp>
#include <sys/ioctl.h>

namespace po = boost::program_options;
using namespace boost;
using namespace std;

#include "instance.hpp"
#include "solution.hpp"
#include "instancedata.hpp"

int main(int argc, char *argv[]) {
    
    po::options_description desc("General options");
    desc.add_options()
            ("help", "show help")
            ("sol", po::value<string>(), "solution, format: 10110101011...")
            ("ins", po::value<string>(), "instance");
    
    po::positional_options_description pod;
    pod.add("ins", 1);
    
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
    po::notify(vm);
    
    if(vm.count("help") || !vm.count("sol") || !vm.count("ins")) {
        cout << desc << endl;
        return 0;
    }

    string instanceName = vm["ins"].as<string>();
    ifstream ins(instanceName);
    if (!ins.good()) {
        cout << "Can't open instance file " << instanceName << endl;
        return 1;
    }

    InstanceData instanceData;
    bool maximize = instanceData.maximize(instanceName);
    string instanceFormat = instanceData.format(instanceName);

    I.readInstance(ins, maximize, instanceFormat);

    string solution = vm["sol"].as<string>();
    Solution S(I, solution);
    cout << S.value << endl;
    
}

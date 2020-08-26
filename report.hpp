#pragma once
#include <iostream>
#include <chrono>

struct Report {
    double bestSoFar;

    Report() {
        bestSoFar = std::numeric_limits<int>::max();
    }

    void newBestKnownValue(double value, std::chrono::system_clock::time_point timePoint = chrono::system_clock::now()) {
        if(value < bestSoFar) {
            bestSoFar = value;
            std::cout << bestSoFar << std::endl;
        }
    }
};
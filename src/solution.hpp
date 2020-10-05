#pragma once

// solution, variable indices 0-based
struct Solution {
    const Instance &I; // corresponding instance
    vector<char> x;    // x-variables
    double value;         // value
    vector<double> d;     // deltas
    std::chrono::system_clock::time_point time; // time found
    
    // data for diversity-driven search
    unsigned ones;
    vector<unsigned> freq; // number of flips;
    
    // data for recombiner
    bool novel;
    
    Solution(const Instance &_I) : I(_I), x(I.n), value(0), d(I.n), ones(0), freq(I.n, 0), novel(true) {
        // generate random solution
        for(char &c : x)
            c = random01(rng) < 0.5 ? 0 : 1;
        
        //code for the test assignment problem
        //(generate a feasible solution)
        if(format == "tap") {
            for(char& c : x)
                c = 0;
            //fill F desks with null color
            for(unsigned i = 0; i < I.uncolored; i++) {
                unsigned desk = rand() % I.desks;
                unsigned index = I.colors * desk;
                if(!x[index]) {
                    x[index] = 1;
                } else {
                    i--;
                }
            }
            //fill each table with exactly one color
            for(unsigned i = 0; i < I.desks; i++) {
                unsigned start = I.colors * i;
                if(!x[start]) {
                    unsigned index = 0;
                    while(index == 0) {
                        index = rand() % I.colors;
                    }
                    index += start;
                    x[index] = 1;
                }
            }
        }
        
        evaluate();
    }

    Solution(const Instance &_I, string s) : I(_I), x(I.n), value(0), d(I.n), ones(0), freq(I.n, 0), novel(true) {
        // generate random solution
        assert(x.size() == I.n);
        for(unsigned i = 0; i < x.size(); i++)
            x[i] = (s[i] == '1' ? 1 : 0);
        evaluate();
    }
    
    Solution &operator=(const Solution &other) {
        x = other.x;
        value = other.value;
        d = other.d;
        time = other.time;
        ones = other.ones;
        freq = other.freq;
        novel = other.novel;
        return *this;
    }
    
    bool operator==(const Solution &other) {
        return x == other.x;
    }
    
    // full evalation in O(n^2)
    void evaluate() {
        value = 0;
        ones = 0;
        for(unsigned i = 0; i < x.size(); i++) {
            d[i] = 0;
            if(x[i])
                ones++;
            for(unsigned j = 0; j < x.size(); j++) {
                if(j == i)
                    d[i] += I[i][i];
                else
                    d[i] += 2 * I[i][j] * x[j];
                value += I[i][j] * x[i] * x[j];
            }
            d[i] *= (1 - 2 * x[i]);
        }
    }
    
    // get delta & update delta's in O(n)
    double adelta(unsigned i) {
        double result = d[i];
        const int sign = x[i] ? -1 : 1;
        const unsigned n = x.size();
        
        for(unsigned j = 0; j < n; j++)
            if(j != i)
                d[j] += 2 * I[j][i] * sign * (1 - 2 * x[j]);
            else
                d[j] = -d[j];
        return result;
    }
    
    double delta(unsigned i) const {
        return d[i];
    }
    
    // update: set a variable in O(n)
    void set(unsigned i) {
        if(x[i])
            return;
        value += adelta(i);
        x[i] = 1;
        ones++;
        freq[i]++;
    }
    
    // update: unset a variable in O(n)
    void reset(unsigned i) {
        if(!x[i])
            return;
        value += adelta(i);
        x[i] = 0;
        ones--;
        freq[i]++;
    }
    
    // flip variable
    void flip(unsigned i) {
        if(x[i])
            reset(i);
        else
            set(i);
    }
    
    // get value of flip
    double flipvalue(unsigned i) const {
        return delta(i);
    }
    
    // distance from given solution
    unsigned distance(const Solution &other) const {
        unsigned result = 0;
        for(unsigned i = 0; i < x.size(); i++)
            if(x[i] != other.x[i])
                result++;
        return result;
    }
};

// a move is simply an integer `i`; by convention i==n is an invalid move
bool isValid(unsigned i, const Solution &S) {
    return i < S.x.size();
}

#pragma once

// BT: Tabu search
//   Random step with probability p, some non-tabu improving move, if it exists, random otherwise.
//   Next variable "i" select in order 0,1,2,...
struct BT {
    double p;
    
    BT(double _p) : p(_p) {}


    // some improvement w/ reservoir sampling
    template<typename Excluder>
    unsigned sie(const Solution &S, Excluder exclude, int bkv) {
        unsigned imp = S.x.size();
        unsigned nimp = 0;
        for(unsigned i = 0; i < S.x.size(); i++) {
            // aspiration?
            if(S.flipvalue(i) + S.value < bkv)
                return i;
            if(!exclude(i) && S.flipvalue(i) < 0) {
                nimp++;
                if(nimp * random01(rng) <= 1)
                    imp = i;
            }
        }
        return imp;
    }
    
    template<typename Excluder>
    unsigned blmre(const Solution &S, Excluder exclude, int bkv) {
        if(random01(rng) < p)
            return randomInt(rng) % S.x.size();
        else {
            unsigned r = sie(S, exclude, bkv);
            if(isValid(r, S))
                return r;
            else
                return randomInt(rng) % S.x.size();
        }
    }
};

// BTR: randomized tabu search
//   Random step with probability p, some non-tabu improving move, if
//   it exists, some non-tabu move of equal value, if that exists,
//   random otherwise.
//   Next variable "i" in round-robin manner.
struct BTR {
    double p;
    unsigned i;
    
    BTR(double _p) : p(_p), i(0) {}
    
    // some improvement w/ reservoir sampling
    template<typename Excluder>
    unsigned sie(const Solution &S, Excluder exclude, int bkv) {
        unsigned imp = S.x.size(), emp = S.x.size();
        unsigned nimp = 0, nemp = 0;
        int nt = (S.x.size() + 9) / 10;
        unsigned tt = S.x.size();
        while (nt-- > 0 || (tt-- > 0 && nimp == 0 && nemp == 0)) {
            // aspiration?
            if(S.flipvalue(i) + S.value < bkv)
                return i;
            if(!exclude(i) && S.flipvalue(i) < 0) {
                nimp++;
                if(nimp * random01(rng) <= 1)
                    imp = i;
            }
            if(!exclude(i) && S.flipvalue(i) <= 0) {
                nemp++;
                if(nemp * random01(rng) <= 1)
                    emp = i;
            }
            i = (i + 1) % S.x.size();
        }
        if(nimp > 0)
            return imp;
        else return emp;
    }
    
    template<typename Excluder>
    unsigned step(const Solution &S, Excluder exclude, int bkv) {
        if(random01(rng) < p)
            return randomInt(rng) % S.x.size();
        else {
            unsigned r = sie(S, exclude, bkv);
            if(isValid(r, S))
                return r;
            else
                return randomInt(rng) % S.x.size();
        }
    }
};

// BTS: Simple tabu search
//   as specified in list 2, exercise 2, 2014/1
//
struct BTS {
    // some improvement w/ reservoir sampling
    BTS() {}
    
    // bi with reservoir sampling
    template<typename Excluder>
    unsigned bi(const Solution &S, Excluder exclude, int bkv) {
        pair<int, unsigned> best = make_pair(numeric_limits<int>::max(), S.x.size()); // delta, index
        unsigned nbest = 0;
        for(unsigned i = 0; i < S.x.size(); i++) {
            // aspiration
            if(S.flipvalue(i) + S.value < bkv)
                return i;
            if(exclude(i))
                continue;
            if(S.flipvalue(i) < get<0>(best)) {
                best = make_pair(S.flipvalue(i), i);
                nbest = 1;
            } else if(get<1>(best) != S.x.size() && S.flipvalue(i) == get<0>(best)) {
                nbest++;
                if(nbest * random01(rng) <= 1)
                    best = make_pair(S.flipvalue(i), i);
            }
        }
        return get<1>(best);
    }
    
    template<typename Excluder>
    unsigned step(const Solution &S, Excluder exclude, int bkv) {
        return bi(S, exclude, bkv);
    }
};

template<typename Improvement, typename Duration>
unsigned tabusearch(Solution &S, Improvement improve, Duration dgen, chrono::system_clock::time_point start, int target,
                    const unsigned maxstagnate = numeric_limits<unsigned>::max(),
                    const unsigned maxsteps = numeric_limits<unsigned>::max()) {
    // local searches
    //chrono::system_clock::time_point last_report = chrono::system_clock::now();
    unsigned i, steps = 0, notimproved = 0;
    Solution B(S.I);
    B = S;
    report.newBestKnownValue(B.value);

    vector<unsigned> tabu(S.x.size() + 1, 0);
    
    do {
        // stop on target value or time limit or non-improvement
        if(S.value <= target)
            break;
        if(chrono::system_clock::now() - start > chrono::duration<int>(timeLimit))
            break;
        if(notimproved > maxstagnate)
            break;
        if(steps > maxsteps)
            break;
        // one step
        steps++;
        i = improve.step(S, [&steps, &tabu](unsigned j) { return steps < tabu[j]; }, B.value);
        // no move: stop
        if(!isValid(i, S))
            continue;
        S.flip(i);
        
        // set tabu, S.x.size() indicates that we want the degree
        unsigned d = dgen();
        if(d == S.x.size())
            tabu[i] = steps + S.I.deg[i] + 1;
        else if(d < S.x.size())
            tabu[i] = steps + d + 1; //+randomInt(rng)%10+1; // hack for list 2, ex 2
        //tabu[i]=steps+d+1+10000*2*double(numdiff[i])/double(steps);  // exp
        
        if(S.value < B.value) {
            //for(auto &t : tabu)
            //t = 0;
            B = S;
            B.time = chrono::system_clock::now();
            notimproved = 0;
            report.newBestKnownValue(B.value);
        } else
            notimproved++;
    } while (true);
    steps--;
    S = B;
    return steps;
}
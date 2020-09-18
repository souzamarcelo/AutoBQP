#pragma once

// semi-greedy construction
//
//  Starts from 0, repeatedly sets variables to 1
//  `alpha`: greediness, 0=none, 1=full
void construct0(const Instance &I, Solution &S, double alpha) {
    unordered_set<unsigned> N; // list of unset variables
    
    // (1) initial solution 0, all variables unfixed
    for(unsigned i = 0; i < I.n; i++) {
        N.insert(i);
        S.reset(i);
    }
    
    // (2) go over free variavles
    while (N.size() > 0) {
        // (2.1) determine cutoff value
        double mv = numeric_limits<double>::max(), Mv = numeric_limits<double>::lowest();
        for(unsigned v : N)
            if(S.flipvalue(v) < 0) {
                mv = min(mv, S.flipvalue(v));
                Mv = max(Mv, S.flipvalue(v));
            }
        if(mv > Mv)
            break;
        double maxvalue = ceil(alpha * mv + (1 - alpha) * Mv);
        
        // (2.2) select
        unsigned cv;
        unsigned ce = 1;
        for(unsigned v : N)
            if(S.flipvalue(v) <= maxvalue && lrand48() % ce == 0) {
                cv = v;
                ce++;
            }
        
        // (2.3) set
        //cout << "Flipping " << cv << endl;
        S.flip(cv);
        N.erase(cv);
    }
}

// semi-greedy construction
//
//  Starts from 1/2, repeatedly fixes variables at 0 or 1 (Merz &
//  Freisleben 2000). `alpha`: greediness, 0=none, 1=full
void constructHalf(const Instance &I, Solution &S, double alpha) {
    // (1) compute initial gains in O(n^2); gains are multiplied by 4 to void fractions
    unordered_set<unsigned> N; // list of unset variables
    vector<int> g0(I.n), g1(I.n);
    for(unsigned i = 0; i < I.n; i++) {
        g1[i] = 3 * I[i][i];
        g0[i] = -1 * I[i][i];
        for(unsigned j = 0; j < I.n; j++) {
            if(j == i)
                continue;
            g1[i] += 2 * I[i][j];
            g0[i] -= 2 * I[i][j];
        }
        N.insert(i);
    }
    
    // (2) go over free variables
    typedef tuple<int, int> Cel;
    
    while (N.size() > 0) {
        // (2.1) determine cutoff value
        int mv = numeric_limits<int>::max(), Mv = numeric_limits<int>::lowest();
        for(unsigned v : N) {
            mv = min(mv, g0[v]);
            mv = min(mv, g1[v]);
            Mv = max(Mv, g0[v]);
            Mv = max(Mv, g1[v]);
        }
        int maxvalue = ceil(alpha * mv + (1 - alpha) * Mv);
        
        // (2.2) select
        Cel sc;
        unsigned ce = 1;
        for(unsigned v : N) {
            if(g0[v] <= maxvalue && lrand48() % ce == 0) {
                sc = make_tuple(v, 0);
                ce++;
            }
            if(g1[v] <= maxvalue && lrand48() % ce == 0) {
                sc = make_tuple(v, 1);
                ce++;
            }
        }
        
        // (2.3) set
        // TBD: could go over non-zeros only
        const unsigned cv = get<0>(sc);
        if(get<1>(sc) == 1) {
            S.set(cv);
            for(unsigned i = 0; i < I.n; i++) {
                g1[i] += 2 * I[i][cv];
                g0[i] -= 2 * I[i][cv];
            }
        } else {
            S.reset(cv);
            for(unsigned i = 0; i < I.n; i++) {
                g1[i] -= 2 * I[i][cv];
                g0[i] += 2 * I[i][cv];
            }
        }
        N.erase(cv);
    }
}

// repeated greedy randomized construction
template<typename Constructor>
void gra(ostream &o, const Instance &I, Solution &S, Constructor construct, unsigned rep, double alpha,
         chrono::system_clock::time_point start) {
    Solution B(I);
    B = S;
    report.newBestKnownValue(B.value);

    o << "# repeated greedy randomized construction" << endl;
    o << "value" << endl;
    
    while (rep > 0) {
        
        if(chrono::system_clock::now() - start > chrono::duration<int>(timeLimit))
            break;
        
        construct(I, S, alpha);
        
        o << S.value << endl;
        
        if(S.value < B.value) {
            B = S;
            B.time = chrono::system_clock::now();
            report.newBestKnownValue(B.value);
        }
        rep--;
    }
    S = B;
}

template<typename Constructor, typename Improver>
void grasp(const Instance &I, Solution &S, Constructor construct, unsigned rep, double alpha, Improver improve,
           chrono::system_clock::time_point start) {
    Solution B(I);
    B = S;
    report.newBestKnownValue(B.value);
    
    while (rep > 0) {
        
        if(chrono::system_clock::now() - start > chrono::duration<int>(timeLimit))
            break;
        
        construct(I, S, alpha);
        improve(S);
        
        if(S.value < B.value) {
            B = S;
            B.time = chrono::system_clock::now();
            report.newBestKnownValue(B.value);
        }
        rep--;
    }
    S = B;
}

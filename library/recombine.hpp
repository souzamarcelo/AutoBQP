#pragma once

namespace recombine {
    // rn: random neighbor
    unsigned rn(const Solution &S, const vector<unsigned> &NC) {
        return uniform_int_distribution<>(0, NC.size() - 1)(rng);
    }
    
    // bn: best neighbor with reservoir sampling
    unsigned bn(const Solution &S, const vector<unsigned> &NC) {
        pair<int, unsigned> best = make_pair(numeric_limits<int>::max(), S.x.size()); // delta, index
        unsigned nbest = 0;
        for(unsigned i = 0; i < NC.size(); i++)
            if(S.flipvalue(NC[i]) < get<0>(best)) {
                best = make_pair(S.flipvalue(NC[i]), i);
                nbest = 1;
            } else if(S.flipvalue(NC[i]) == get<0>(best)) {
                assert(get<1>(best) != S.x.size());
                nbest++;
                if((nbest) * random01(rng) <= 1.0)
                    best = make_pair(S.flipvalue(NC[i]), i);
            }
        return get<1>(best);
    }
    
    unsigned si(const Solution &S, const vector<unsigned> &NC) {
        unsigned imp = 0;
        unsigned nimp = 0;
        for(unsigned i = 0; i < NC.size(); i++) {
            if(S.flipvalue(NC[i]) < 0) {
                nimp++;
                if(nimp * random01(rng) <= 1) {
                    imp = i;
                }
            }
        }
        
        if(nimp == 0)
            return rn(S, NC);
        return imp;
    }
    
    struct SiPartial {
        unsigned i;
        bool rr;
        double partial_size;
        
        SiPartial(bool _rr, double _partial_size) : i(0), rr(_rr), partial_size(_partial_size) {}
        
        unsigned operator()(const Solution &S, const vector<unsigned> &NC) {
            unsigned imp = 0, emp = 0;
            unsigned nimp = 0, nemp = 0;
            int nt = (NC.size() + 9) * (partial_size / 100);
            unsigned tt = NC.size();
            while (nt-- > 0 || (tt-- > 0 && nimp == 0)) { //  && nemp==0
                if(S.flipvalue(NC[i]) < 0) {
                    nimp++;
                    if(nimp * random01(rng) <= 1)
                        imp = i;
                }
                if(S.flipvalue(NC[i]) <= 0) {
                    nemp++;
                    if(nemp * random01(rng) <= 1)
                        emp = i;
                }
                i = (i + 1) % NC.size();
            }
            
            if(!rr) i = 0;
            
            if(nimp > 0)
                return imp;
            else if(nemp > 0)
                return emp;
            else
                return rn(S, NC);
        }
    };
    
    // Fn: first neighbor with round-robin
    struct Fn {
        unsigned i;
        bool rr;
        
        Fn(bool _rr) : i(0), rr(_rr) {}
        
        unsigned operator()(const Solution &S, const vector<unsigned> &NC) {
            pair<int, unsigned> best = make_pair(numeric_limits<int>::max(), S.x.size()); // delta, index
            unsigned nbest = 0;
            unsigned ncand = NC.size();
            i = i % NC.size();
            while (ncand-- > 0 && S.value <= get<0>(best)) {
                if(S.flipvalue(NC[i]) < get<0>(best)) {
                    best = make_pair(S.flipvalue(NC[i]), i);
                    nbest = 1;
                } else if(S.flipvalue(NC[i]) == get<0>(best)) {
                    assert(get<1>(best) != S.x.size());
                    nbest++;
                    if((nbest) * random01(rng) <= 1.0)
                        best = make_pair(S.flipvalue(NC[i]), i);
                }
                i = (i + 1) % NC.size();
            }
            
            if(!rr) i = 0;
            
            return get<1>(best);
        }
    };
    
    // path-relinking from `S` to `T`, return the best solution in `S` with minimum
    // distance gamma*H(S,T) for Hamming distance H between `S` and `T`
    template<typename Improver>
    void recombine(const Instance &I, Solution &S, const Solution &T, Improver improve, double gamma) {
        const unsigned n = S.x.size();
        assert(n == T.x.size());
        assert(0 <= gamma && gamma < 0.5);
        
        Solution B(S.I);
        B.value = numeric_limits<int>::max();
        
        // (1) generate the difference set
        vector<unsigned> NC;
        for(unsigned i = 0; i < n; i++)
            if(S.x[i] != T.x[i])
                NC.push_back(i);
        
        unsigned dmin = gamma * NC.size();
        unsigned dmax = NC.size() - dmin;
        //unsigned dmax = NC.size()-1;
        
        // (2) repeatedly let "improve" select one of the different bits
        unsigned d = 0;
        while (NC.size() > 0) {
            // (2.1) select a variable, flip it, remove it
            unsigned i = improve(S, NC);
            //cout << "Selected " << i << " from "; copy(NC.begin(),NC.end(),ostream_iterator<unsigned>(cout," ")); cout << endl;
            S.flip(NC[i]);
            d++;
            
            swap(NC[i], NC.back());
            NC.pop_back();
            
            if(S.value < B.value && dmin <= d && d <= dmax) {
                B = S;
                B.time = chrono::system_clock::now();
            }
        }
        
        S = B;
    }
    
    // simple recombiner, here to reproduce Wang et al. (2012)
    //
    //   maintains `b` solutions and systematicially recombines them
    //   until
    template<typename Improver, typename Recombiner>
    void recombiner(const Instance &I, Solution &S, unsigned b, Recombiner recombine, Improver improve,
                    chrono::system_clock::time_point start, int target) {
        chrono::system_clock::time_point last_report = chrono::system_clock::now();
        unsigned steps = 0;
        vector<Solution> current;
        Solution B(S.I);
        B = S;
        
        Elite e(I, b);
        
        while (true) {
            // (1) create b random solutions, improve them
            e.clear();
            for(unsigned i = 0; i < b; i++) {
                Solution S(I);
                improve(S);
                
                if(i == 0 || S.value < e[0].value)
                    report.newBestKnownValue(S.value);

                e.add(S);
                if(termination(e[0].value, e.back().value, steps, target, start, last_report))
                    goto done;
            }
            
            // (2) extract current solutions
            current = e;
            for(Solution &s : e)
                s.novel = false;
            
            unsigned num_novel = b;
            
            while (num_novel > 0) {
                num_novel = 0;
                
                // (3) recombine until stagnation
                for(unsigned i = 0; i < current.size(); i++)
                    for(unsigned j = i + 1; j < current.size(); j++) {
                        if(!current[i].novel && !current[j].novel)
                            continue;
                        steps++;
                        
                        // forward
                        Solution S = current[i];
                        recombine(I, S, current[j]);
                        improve(S);
                        S.novel = true;

                        if(S.value < e[0].value)
                            report.newBestKnownValue(S.value);

                        if(e.add(S))
                            num_novel++;
                        
                        if(termination(e[0].value, e.back().value, steps, target, start, last_report))
                            goto done;
                        
                        // backward
                        S = current[j];
                        recombine(I, S, current[i]);
                        improve(S);
                        S.novel = true;

                        if(S.value < e[0].value)
                            report.newBestKnownValue(S.value);

                        if(e.add(S))
                            num_novel++;
                        
                        if(termination(e[0].value, e.back().value, steps, target, start, last_report))
                            goto done;
                    }
                
                // (4) recreate pairset
                current = e;
                for(Solution &s : e)
                    s.novel = false;
            } // not converged (i.e. num_novel>0)
            if(e[0].value < B.value)
                B = e[0];
        } // while "restart"
        done:
        if(e[0].value < B.value)
            B = e[0];
        S = B;
    }
    
} // namespace recombine

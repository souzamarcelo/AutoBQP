void run(Solution& S, Instance& I, double target, chrono::system_clock::time_point startTime) {
    Elite e1(I, 17);
    double p = 0.200000;
    RII rii(p);
    vector<double> w(I.n,0);for(unsigned i=0; i<I.n; i++)  w[i]=pow(i+1,3.000000);rank_proportional = std::discrete_distribution<>(w.begin(),w.end());StepGaussian sg1(51, I.n>=51*60?I.n/60:51, I.n);

    iterated_localsearch_elite(S, [&](Solution& S) { return
                                       localsearch(S, [&](const Solution& S) { return rii.blmr(S, si); }, startTime, target);
                               },
                               [&](Solution& S) { return perturbDiversity(S, sg1, e1, 0.200000); },
                               e1, startTime, target);
}

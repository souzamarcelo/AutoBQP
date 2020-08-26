#pragma once

struct Instance : public boost::multi_array<double, 2> {
    unsigned n, m;
    std::vector<unsigned> deg;
    unsigned nnz;
    
    unsigned desks, colors, uncolored;
    double P;
    double btb;
    
    unsigned deskDistances;
    boost::multi_array<double,2> penaltiesColors;
    boost::multi_array<double,2> distancesDesks;
    boost::multi_array<double,2> q;
    boost::multi_array<double,2> A;
    boost::multi_array<double,2> At;
    boost::multi_array<double,2> AtA;
    vector<double> b;
    vector<double> twoDiagAtb;
    
    void setSize(unsigned _n) {
        n = _n;
        checkPlausible();
        resize(extents[n][n]);
        deg.resize(n);
        nnz = 0;
    }
    
    void checkPlausible() {
        if(n > 15000) {
            cerr << "Instance has " << n << " variables. Refused." << endl;
            exit(1);
        }
    }

    void readInstance(istream &ins, bool maximize, string format) {
        if(format == "sparse")
            readSparse(ins, maximize);
        if(format == "dense")
            readDense(ins, maximize);
        if(format == "pgen")
            readPalubeckis(ins, maximize, false);
        if(format == "lgen")
            readPalubeckis(ins, maximize, true);
        if(format == "maxcut")
            readMaxcut(ins);
        if(format == "clique")
            readClique(ins);
        if(format == "tap")
            readTap(ins);
    }
    
    void readSparse(istream &ins, bool maximize = false) {
        ins >> n;
        setSize(n);
        
        ins >> m;
        this->nnz = 2 * m; // TBD: wrong? diagonals!
        unsigned nnz = m;
        while (nnz > 0) {
            unsigned i, j;
            ins >> i >> j;
            ins >> (*this)[i - 1][j - 1];
            if(maximize)
                (*this)[i - 1][j - 1] = -(*this)[i - 1][j - 1];
            (*this)[j - 1][i - 1] = (*this)[i - 1][j - 1];
            if((*this)[j - 1][i - 1] != 0) {
                deg[i - 1]++;
                deg[j - 1]++;
            }
            nnz--;
        }
    }
    
    void readDense(istream &ins, bool maximize = false) {
        ins >> n; // dummy: 1
        ins >> n;
        setSize(n);
        
        for(unsigned i = 0; i < n; i++)
            for(unsigned j = 0; j < n; j++) {
                ins >> (*this)[i][j];
                if(maximize)
                    (*this)[i][j] = -(*this)[i][j];
                (*this)[j][i] = (*this)[i][j];
                if((*this)[j][i] != 0) {
                    deg[i]++;
                    deg[j]++;
                    nnz++;
                }
            }
    }
    
    void readMaxcut(istream &ins) {
        ins >> n;
        setSize(n);
        
        ins >> m;
        this->nnz = 2 * m + n;
        unsigned nedge = m;
        while (nedge > 0) {
            unsigned i, j;
            int wij;
            ins >> i >> j;
            ins >> wij;
            (*this)[i - 1][i - 1] -= wij;
            (*this)[j - 1][j - 1] -= wij;
            (*this)[i - 1][j - 1] += wij;
            (*this)[j - 1][i - 1] += wij;
            nedge--;
            if(wij != 0) {
                deg[i - 1]++;
                deg[j - 1]++;
            }
        }
    }

    void readClique(istream &ins) {
        double p = 1000;
        string line;
        string value;
        while(getline(ins, line)) {
            ins >> value;
            if(value != "c")
                break;
        }
        ins >> value;
        ins >> n;
        setSize(n);

        for(unsigned i = 0; i < n; i++) {
            for(unsigned j = 0; j < n; j++) {
                (*this)[i][j] = p;
                if(i == j) {
                    (*this)[i][j] = (double) (-(((i + 1) % 200) + 1));
                }
            }
        }

        ins >> m;
        for(unsigned a = 0; a < m; a++) {
            ins >> value;
            unsigned i, j;
            ins >> i;
            ins >> j;

            (*this)[i - 1][j - 1] = 0;
            (*this)[j - 1][i - 1] = 0;
        }
    }
    
    void genPalubeckis(int density, int lb_linear, int ub_linear, int lb_quadr, int ub_quadr, int clumpiness, double seed,
                  bool negate = true, bool clumpy = false) {
        double coef = 2048;
        coef *= 1024;
        coef *= 1024;
        coef -= 1;
        
        for(unsigned i = 0; i < n; i++) {
            int r = lehmer::random(&seed, coef) * (ub_linear - lb_linear + 1);
            // check bqpdata/data/Readme.org for an explanation of setting `j` here to n+1
            unsigned j = n + 1;
            if(clumpy && (i + 1 + j + 1) % clumpiness == 0)
                r -= ub_linear;        // here's the perturbation which could create a small diagonal Q(i,i)
            (*this)[i][i] = r + lb_linear;
            if(negate)
                (*this)[i][i] = -(*this)[i][i];
            if((*this)[i][i] != 0)
                nnz++;
            for(unsigned j = i + 1; j < n; j++) {
                float fl = lehmer::random(&seed, coef) * 100;
                if(fl <= density) {
                    r = lehmer::random(&seed, coef) * (ub_quadr - lb_quadr + 1);
                    if(clumpy && (i + 1 + j + 1) % clumpiness == 0)
                        r += ub_quadr;    // here's the perturbation which could create a large quadratic Q(i,j)
                    (*this)[i][j] = r + lb_quadr;
                    if(negate)
                        (*this)[i][j] = -(*this)[i][j];
                    if((*this)[i][j] != 0) {
                        nnz += 2;
                        deg[i]++;
                        deg[j]++;
                    }
                } else
                    (*this)[i][j] = 0;
                (*this)[j][i] = (*this)[i][j];
            }
        }
    }
    
    // read generator format of Palubeckis
    //   for clumpy==true: generator of Lewis with extra `clumpiness` parameter
    //   NOTE: Lewis' variant needs more tests.
    void readPalubeckis(istream &ins, bool negate = true, bool clumpy = false) {
        ins >> n;
        setSize(n);
        
        int density, lb_linear, ub_linear, lb_quadr, ub_quadr;
        int clumpiness = 0;
        double seed;
        ins >> density >> lb_linear >> ub_linear >> lb_quadr >> ub_quadr >> seed;
        if(clumpy) {
            ins >> clumpiness;
            if(clumpiness == 0)
                clumpiness = 2 * n;
        }
        genPalubeckis(density, lb_linear, ub_linear, lb_quadr, ub_quadr, clumpiness, seed, negate, clumpy);
    }
    
    void readTap(istream &ins) {
        ins >> desks;
        ins >> deskDistances;
        ins >> colors;
        ins >> uncolored;
        
        nnz = 0;
        setSize(desks * colors);
        q.resize(extents[desks * colors][desks * colors]);
        penaltiesColors.resize(extents[colors][colors]);
        distancesDesks.resize(extents[desks][desks]);
        
        vector<string> deskNames;
        for(unsigned i = 0; i < desks; i++) {
            string deskName;
            ins >> deskName;
            deskNames.push_back(deskName);
        }
        
        for(unsigned i = 0; i < desks; i++) {
            for(unsigned j = 0; j < desks; j++) {
                distancesDesks[i][j] = 0;
            }
        }
        
        for(unsigned i = 0; i < deskDistances; i++) {
            string name1, name2;
            double value;
            ins >> name1;
            ins >> name2;
            ins >> value;
        
            unsigned index1 = -1, index2 = -1;
            for(unsigned j = 0; j < deskNames.size(); j++) {
                if(deskNames[j] == name1) index1 = j;
                if(deskNames[j] == name2) index2 = j;
            }
        
            distancesDesks[index1][index2] = value;
            distancesDesks[index2][index1] = value;
        }
        
        for(unsigned i = 0; i < colors; i++) {
            for(unsigned j = 0; j < colors; j++) {
                penaltiesColors[i][j] = 0;
            }
        }
        
        unsigned penalties = (colors) * (colors - 1);
        for(unsigned i = 0; i < penalties; i++) {
            unsigned color1, color2;
            double penalty;
            ins >> color1;
            ins >> color2;
            ins >> penalty;
        
            penaltiesColors[color1][color2] = penalty;
            penaltiesColors[color2][color1] = penalty;
        }
        
        for(unsigned i = 0; i < n; i++) {
            for(unsigned j = 0; j < n; j++) {
                q[i][j] = 0;
            }
        }
        
        for(unsigned c1 = 0; c1 < colors; c1++) {
            for(unsigned d1 = 0; d1 < desks; d1++) {
                for(unsigned c2 = 0; c2 < colors; c2++) {
                    for(unsigned d2 = 0; d2 < desks; d2++) {
                        unsigned index1 = (colors * d1) + c1;
                        unsigned index2 = (colors * d2) + c2;
                        q[index1][index2] = distancesDesks[d1][d2] * penaltiesColors[c1][c2];
                        deg[index1]++;
                        deg[index2]++;
                        if(q[index1][index2] != 0) nnz++;
                    }
                }
            }
        }
        
        for(unsigned i = 0; i < n; i++) {
            for(unsigned j = 0; j < n; j++) {
                q[i][j] = q[i][j] / 2;
            }
        }
    
        for(unsigned i = 0; i < n; i++) {
            for(unsigned j = 0; j < n; j++) {
                (*this)[i][j] = q[i][j];
            }
        }
        
        computeTransformation();
    }
    
    void writeDense(ostream &out) {
        out << 1 << " " << n << endl;
        
        for(unsigned i = 0; i < n; i++) {
            for(unsigned j = 0; j < n; j++) {
                out << setw(5) << (*this)[i][j];
                if((j + 1) % 15 == 0)
                    out << endl;
            }
            if(n < 14)
                out << endl;
        }
        out << endl;
    }
    
    void writeSparse(ostream &out) {
        unsigned ndiag = 0;
        for(unsigned i = 0; i < n; i++)
            if((*this)[i][i] != 0)
                ndiag++;
        out << n << " " << ndiag + (nnz - ndiag) / 2 << endl;
        for(unsigned i = 0; i < n; i++)
            for(unsigned j = i; j < n; j++)
                if((*this)[i][j] != 0)
                    out << i + 1 << " " << j + 1 << " " << (*this)[i][j] << " "
                        << endl; // space is to match Lewis' generator
    }
    
    // show an instance, positive values are red, negative ones green
    void writePPM(ostream &out, unsigned bf = 1) {
        vector<unsigned> order(n, 0);
        iota(order.begin(), order.end(), 0);
        writePPM(out, order, bf);
    }
    
    // show an instance, positive values are red, negative ones green
    void writePPM(ostream &out, const vector<unsigned> &order, unsigned bf = 1) {
        assert(order.size() == n);
        out << "P3 " << " " << n / bf << " " << n / bf << endl;
        out << "255" << endl;
        
        // (1) go over all blocks, find maximum absolute sum
        int bmax = 0;
        for(unsigned i = 0; i < n / bf; i++) {
            for(unsigned j = 0; j < n / bf; j++) {
                int v = 0;
                for(unsigned bi = 0; bi < bf; bi++)
                    for(unsigned bj = 0; bj < bf; bj++)
                        v += (*this)[order[i * bf + bi]][order[j * bf + bj]];
                if(v > bmax)
                    bmax = v;
                if(-v > bmax)
                    bmax = -v;
            }
        }
        
        // (2) go over them again and create an image
        for(unsigned i = 0; i < n / bf; i++) {
            for(unsigned j = 0; j < n / bf; j++) {
                //int v = (*this)[order[i]][order[j]];
                int v = 0;
                for(unsigned bi = 0; bi < bf; bi++)
                    for(unsigned bj = 0; bj < bf; bj++)
                        v += (*this)[order[i * bf + bi]][order[j * bf + bj]];
                if(v > 0)
                    out << int(v * 255 / bmax) << " 0 0 ";
                else if(v == 0)
                    out << "0 0 0 ";
                else
                    out << "0 " << int(-v * 255 / bmax) << " 0 ";
            }
            out << endl;
        }
    }
    
    void computeTransformation() {

        //Create matrix A
        A.resize(extents[desks + 1][n]);
        for(unsigned d = 0; d < desks; d++) {
            unsigned start = d * colors;
            unsigned finish = start + colors;
            for(unsigned i = 0; i < n; i++) {
                if(i >= start && i < finish) {
                    A[d][i] = 1;
                } else {
                    A[d][i] = 0;
                }
            }
        }
        unsigned last = desks;
        unsigned color = 0;
        for(unsigned i = 0; i < n; i++) {
            color++;
            if(color == 1) {
                A[last][i] = 1;
            } else {
                A[last][i] = 0;
            }
            if(color == colors) {
                color = 0;
            }
        }
        
        //Create matrix At
        At.resize(extents[n][desks + 1]);
        for(unsigned i = 0; i < desks + 1; i++) {
            for(unsigned j = 0; j < n; j++) {
                At[j][i] = A[i][j];
            }
        }
        
        //Create AtA
        AtA.resize(extents[n][n]);
        for(unsigned i = 0; i < n; i ++) {
            for(unsigned j = 0; j < n; j++) {
                unsigned result = 0;
                for(unsigned k = 0; k < desks + 1; k++) {
                    result += At[i][k] * A[k][j];
                }
                AtA[i][j] = result;
            }
        }
        
        //Create vector b
        b.clear();
        for(unsigned i = 0; i < desks; i++) {
            b.push_back(1);
        }
        b.push_back(uncolored);
        
        //Create vector twoDiagAtb
        vector<unsigned> twoDiagAtb;
        twoDiagAtb.clear();
        for(unsigned i = 0; i < n; i++) {
            unsigned value = 0;
            for(unsigned j = 0; j < desks + 1; j++) {
                value += At[i][j] * b[j];
            }
            twoDiagAtb.push_back(2 * value);
        }

        //Create penalty P
        P = 0;
        for(unsigned i = 0; i < n; i++) {
            for(unsigned j = 0; j < n; j++) {
                P += q[i][j];
            }
        }
        //Compute AtA - twoDiagAtb
        for(unsigned i = 0; i < n; i++) {
            for(unsigned j = 0; j < n; j++) {
                if(i == j) {
                    AtA[i][j] = AtA[i][j] - twoDiagAtb[i];
                }
            }
        }
        
        //Compute P * (AtA - twoDiagAtb)
        for(unsigned i = 0; i < n; i++) {
            for(unsigned j = 0; j < n; j++) {
                AtA[i][j] = P * AtA[i][j];
            }
        }
        
        //Compute Q_hat
        for(unsigned i = 0; i < n; i++) {
            for(unsigned j = 0; j < n; j++) {
                (*this)[i][j] = q[i][j] + AtA[i][j];
            }
        }
        
        //Compute btb
        btb = desks + (uncolored * uncolored);
    }
};

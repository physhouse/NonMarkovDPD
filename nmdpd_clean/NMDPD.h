#ifndef _FITGLE_H_
#define _FITGLE_H_

#include <cstdlib>
#include <cstdio>
#include <vector>
#include "Frame.h"
#include <memory>
#include <gsl/gsl_bspline.h>

namespace NMDPD_NS {

struct InputParameters
{
    double start;  // start distance r0
    double end;    // end distance r1
    double interval;
    unsigned tableLength;
    double boxLength; 
    int    steps;
    int    decaySteps;
    int    numParticles;
    double temperature;
    double dt;
    double mass;
};  //Structure to store input parameters

class NMDPD
{
public:
    NMDPD(int argc, char** argv);
    ~NMDPD();
    void exec();
    void init_one();

    //helper functions
    void updateVelocities();
    void updatePositions();
    void pairCompute();
    void dumpAndThermo(int);
    void endOfStep();
    
    //helper functions
    std::vector<double> parallelVector(std::vector<double>, std::vector<double>);
    void pbc(std::vector<double> &);
    double rand_gauss();

private:
    // trajFrame stores the previous frames of simulation
    std::shared_ptr<class Frame> trajFrame;
    //class Frame* trajFrame;
    std::shared_ptr<struct InputParameters> info; 
    std::vector<std::vector<double> > table;
    std::vector<double> tapCoeffParallel;
    std::vector<double> tapCoeffVertical;

    std::vector<std::vector<double> > x;
    std::vector<std::vector<double> > v;
    std::vector<std::vector<double> > f;
    std::vector<std::vector<double> > randoms;
    std::vector<std::vector<double> > randomsV1;
    std::vector<std::vector<double> > randomsV2;

    FILE* output;
    double dt;
    std::vector<double> dtfm;
    double dtv;
    double kT;
    double Ek;
};

}

#endif

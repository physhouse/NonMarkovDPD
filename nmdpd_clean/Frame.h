#ifndef _FRAME_H_
#define _FRAME_H_

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <deque>
#include "NMDPD.h"

namespace NMDPD_NS {

class Frame
{
public:
    friend class NMDPD;
    Frame(int n, char* filename, int decaySteps);
    ~Frame();
    void readFrame(std::vector<std::vector<double> >, std::vector<std::vector<double> >, std::vector<std::vector<double> >, std::vector<std::vector<double> >, std::vector<std::vector<double> >, std::vector<std::vector<double> >);
    void readOne(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);
    int get();
    double rand_gauss();

private:
    FILE* trajectory;
    int numParticles;
    int lagStep;
    std::deque<std::vector<std::vector<double> > > positions;
    std::deque<std::vector<std::vector<double> > > residualForces;
    std::deque<std::vector<std::vector<double> > > velocities;
    std::deque<std::vector<std::vector<double> > > pairRandoms;
    std::deque<std::vector<std::vector<double> > > pairRandomsV1;
    std::deque<std::vector<std::vector<double> > > pairRandomsV2;
};
}

#endif

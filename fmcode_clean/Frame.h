#ifndef _FRAME_H_
#define _FRAME_H_

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <deque>
#include "FitGLE.h"

namespace FITGLE_NS {

class Frame
{
public:
    friend class FitGLE;
    Frame(int n, char* filename, int decaySteps);
    ~Frame();
    void readFrame();
    int get();
private:
    FILE* trajectory;
    int numParticles;
    int lagStep;
    std::deque<std::vector<std::vector<double> > > positions;
    std::deque<std::vector<std::vector<double> > > residualForces;
    std::deque<std::vector<std::vector<double> > > velocities;
};
}

#endif

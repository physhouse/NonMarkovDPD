#include "Frame.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>
using namespace FITGLE_NS;

Frame::Frame(int n, char* fileName, int decayStep)
{
    trajectory = fopen(fileName, "r");
    numParticles = n;
    lagStep = 5 * decayStep;
    printf("lagStep = %d\n", lagStep);

    positions.resize(lagStep);
    velocities.resize(lagStep);
    residualForces.resize(lagStep);

    for (int step=0; step<lagStep; step++)
    {
       positions[step].resize(numParticles);
       residualForces[step].resize(numParticles);
       velocities[step].resize(numParticles);
       for (int i=0; i<numParticles; i++)
       {
          positions[step][i].resize(3);
          residualForces[step][i].resize(3);
          velocities[step][i].resize(3);
       }

       ssize_t read;
       size_t len;
       char*  line = NULL;
       int    lineID;
       for (int iline = 0; iline < numParticles+9; iline++)
       {
           read = getline(&line, &len, trajectory);
           if (iline >= 9)
           {
              char* pch = strtok(line, " \t");
              int atomID = atoi(pch) - 1;
              pch = strtok(NULL, " \t");
              int type = atoi(pch);
              pch = strtok(NULL, " \t");
              double atomMass = atof(pch);
              pch = strtok(NULL, " \t");
              positions[step][atomID][0] = atof(pch);
              pch = strtok(NULL, " \t");
              positions[step][atomID][1] = atof(pch);
              pch = strtok(NULL, " \t");
              positions[step][atomID][2] = atof(pch);
              pch = strtok(NULL, " \t");
              velocities[step][atomID][0] = atof(pch);
              pch = strtok(NULL, " \t");
              velocities[step][atomID][1] = atof(pch);
              pch = strtok(NULL, " \t");
              velocities[step][atomID][2] = atof(pch);
              pch = strtok(NULL, " \t");
              residualForces[step][atomID][0] = atof(pch);
              pch = strtok(NULL, " \t");
              residualForces[step][atomID][1] = atof(pch);
              pch = strtok(NULL, " \t");
              residualForces[step][atomID][2] = atof(pch);
           }
       }
    }
    printf("finishing initializing trajectory frames, %d Particles\n", numParticles);
}

Frame::~Frame()
{
    fclose(trajectory);
    printf("cleaning up Frame Information\n");
}

int Frame::get()
{
  return numParticles;
}

void Frame::readFrame()
{
    ssize_t read;
    size_t len;
    char*  line = NULL;
    int    lineID;

    positions.pop_front();
    velocities.pop_front();
    residualForces.pop_front();

    auto thisPosition = std::vector<std::vector<double> > (numParticles, std::vector<double>(3));
    auto thisVelocity = std::vector<std::vector<double> > (numParticles, std::vector<double>(3));
    auto thisResidualForce = std::vector<std::vector<double> > (numParticles, std::vector<double>(3));

    for (int iline = 0; iline < numParticles+9; iline++)
    {
        read = getline(&line, &len, trajectory);
        if (iline >= 9)
        {
           char* pch = strtok(line, " \t");
           int atomID = atoi(pch) - 1;
           pch = strtok(NULL, " \t");
           int type = atoi(pch);
           pch = strtok(NULL, " \t");
           double atomMass = atof(pch);
           pch = strtok(NULL, " \t");
           thisPosition[atomID][0] = atof(pch);
           pch = strtok(NULL, " \t");
           thisPosition[atomID][1] = atof(pch);
           pch = strtok(NULL, " \t");
           thisPosition[atomID][2] = atof(pch);
           pch = strtok(NULL, " \t");
           thisVelocity[atomID][0] = atof(pch);
           pch = strtok(NULL, " \t");
           thisVelocity[atomID][1] = atof(pch);
           pch = strtok(NULL, " \t");
           thisVelocity[atomID][2] = atof(pch);
           pch = strtok(NULL, " \t");
           thisResidualForce[atomID][0] = atof(pch);
           pch = strtok(NULL, " \t");
           thisResidualForce[atomID][1] = atof(pch);
           pch = strtok(NULL, " \t");
           thisResidualForce[atomID][2] = atof(pch);
        }
    }
    
    positions.push_back(thisPosition);
    velocities.push_back(thisVelocity);
    residualForces.push_back(thisResidualForce);
}

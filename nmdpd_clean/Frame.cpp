#include "Frame.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>
using namespace NMDPD_NS;

Frame::Frame(int n, char* fileName, int decayStep)
{
    trajectory = fopen(fileName, "r");
    numParticles = n;
    lagStep = decayStep;
    printf("lagStep = %d\n", lagStep);

    positions.resize(lagStep);
    velocities.resize(lagStep);
    residualForces.resize(lagStep);
    pairRandoms.resize(lagStep);
    pairRandomsV1.resize(lagStep);
    pairRandomsV2.resize(lagStep);

    for (int step=0; step<lagStep; step++)
    {
       positions[step].resize(numParticles);
       residualForces[step].resize(numParticles);
       velocities[step].resize(numParticles);
       pairRandoms[step].resize(numParticles);
       pairRandomsV1[step].resize(numParticles);
       pairRandomsV2[step].resize(numParticles);
       for (int i=0; i<numParticles; i++)
       {
          positions[step][i].resize(3);
          residualForces[step][i].resize(3);
          velocities[step][i].resize(3);
          pairRandoms[step][i].resize(numParticles);
          pairRandomsV1[step][i].resize(numParticles);
          pairRandomsV2[step][i].resize(numParticles);
       }

       ssize_t read;
       size_t len;
       char*  line = NULL;
       int    lineID;
       for (int iline = 0; iline < numParticles+9; iline++)
       {
           read = getline(&line, &len, trajectory);
           //printf("%s", line);
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
              for (int j=atomID+1; j<numParticles; j++)
              {
                  double value = rand_gauss();
                  pairRandoms[step][atomID][j] = value;
                  pairRandoms[step][j][atomID] = -value;
                  value = rand_gauss();
                  pairRandomsV1[step][atomID][j] = value;
                  pairRandomsV1[step][j][atomID] = -value;
                  value = rand_gauss();
                  pairRandomsV2[step][atomID][j] = value;
                  pairRandomsV2[step][j][atomID] = -value;
              }
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

void Frame::readOne(std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &v, std::vector<std::vector<double> > &f)
{
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
           x[atomID][0] = atof(pch);
           pch = strtok(NULL, " \t");
           x[atomID][1] = atof(pch);
           pch = strtok(NULL, " \t");
           x[atomID][2] = atof(pch);
           pch = strtok(NULL, " \t");
           v[atomID][0] = atof(pch);
           pch = strtok(NULL, " \t");
           v[atomID][1] = atof(pch);
           pch = strtok(NULL, " \t");
           v[atomID][2] = atof(pch);
           pch = strtok(NULL, " \t");
           f[atomID][0] = atof(pch);
           pch = strtok(NULL, " \t");
           f[atomID][1] = atof(pch);
           pch = strtok(NULL, " \t");
           f[atomID][2] = atof(pch);
        }
    }
    printf("finishing initialize one\n");
}

int Frame::get()
{
  return numParticles;
}

void Frame::readFrame(std::vector<std::vector<double> > x, 
                      std::vector<std::vector<double> > v, 
                      std::vector<std::vector<double> > f, 
                      std::vector<std::vector<double> > r, 
                      std::vector<std::vector<double> > v1, 
                      std::vector<std::vector<double> > v2)
{
    ssize_t read;
    size_t len;
    char*  line = NULL;
    int    lineID;

    positions.pop_front();
    velocities.pop_front();
    residualForces.pop_front();
    pairRandoms.pop_front();
    pairRandomsV1.pop_front();
    pairRandomsV2.pop_front();
    
    positions.push_back(x);
    velocities.push_back(v);
    residualForces.push_back(f);
    pairRandoms.push_back(r);
    pairRandomsV1.push_back(v1);
    pairRandomsV2.push_back(v2);
}

double Frame::rand_gauss()
{
    static double n2 = 0.0;
    static int n2_cached = 0;
    double mean = 0.0;
    double stddev = 1.0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

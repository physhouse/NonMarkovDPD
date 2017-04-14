#include "NMDPD.h"
#include "Frame.h"
#include "comm.h"
#include <numeric>
#include <functional>
#include <algorithm>
#include <cassert>

using namespace NMDPD_NS;

NMDPD::NMDPD(int argc, char** argv)
{
    
    if (argc != 3)
    {
        printf("./FitGLE.x [Bootrap Filename] [number of Particles]\n");
    }

    assert(argc == 3);
    printf("Initializing FitGLE parameters...\n");

    // parsing the configuration parameters
    info = std::make_shared<InputParameters>();
    VAR_BEGIN
      GET_REAL(info->start)
      GET_REAL(info->end)
      GET_REAL(info->interval)
      GET_INT(info->tableLength)
      GET_REAL(info->boxLength)
      GET_INT(info->steps)
      GET_INT(info->decaySteps)
      GET_INT(info->numParticles)
      GET_REAL(info->temperature)
      GET_REAL(info->dt)
      GET_REAL(info->mass)
    VAR_END

    // Dealing With the Table
    table.resize(info->tableLength);
    printf("start reading tables...\n");
    FILE* fp = fopen("table.dat", "r");
    for (int i=0; i<info->tableLength; i++)
    {
       table[i].resize(4);
       fscanf(fp, "%lf %lf %lf %lf\n", &table[i][0], &table[i][1], &table[i][2], &table[i][3]);
       printf("%lf %lf %lf %lf\n", table[i][0], table[i][1], table[i][2], table[i][3]);
    }
    fclose(fp);
    printf("finishing reading tables...\n");

    // Read in the tap coefficients for calculating the random forces
    FILE* ff = fopen("random.dat", "r");
    tapCoeffParallel.resize(info->decaySteps + 1);
    tapCoeffVertical.resize(info->decaySteps + 1);
    for (int i=0; i<info->decaySteps + 1; i++)
    {
        fscanf(ff, "%lf %lf\n", &tapCoeffParallel[i], &tapCoeffVertical[i]);
        printf("%lf %lf\n", tapCoeffParallel[i], tapCoeffVertical[i]);
    }
    fclose(ff);
    printf("finishing reading random coeffs...\n");

    // Initialize the arrays
    x.resize(info->numParticles);
    v.resize(info->numParticles);
    f.resize(info->numParticles);
    randoms.resize(info->numParticles);
    randomsV1.resize(info->numParticles);
    randomsV2.resize(info->numParticles);
    for (int i=0; i<info->numParticles; i++)
    {
        x[i].resize(3);
        v[i].resize(3);
        f[i].resize(3);
        randoms[i].resize(info->numParticles);
        randomsV1[i].resize(info->numParticles);
        randomsV2[i].resize(info->numParticles);
    }

    // Initialize the first several frames by reading trajectories, meaning to bootstrap the simulation
    printf("set up frames files\n");
    trajFrame = std::make_shared<Frame>(atoi(argv[2]), argv[1], info->decaySteps);    

    // Some Local Parameters
    dt = info->dt;
    dtv = dt;
    dtfm.resize(info->numParticles);
    for (int i=0; i<info->numParticles; i++) dtfm[i] = 0.5 * dt / info->mass;

    kT = 1.0 * info->temperature;

    // OUput pointer
    output = fopen("out.trj", "w");
    init_one();
}

NMDPD::~NMDPD()
{
    fclose(output);
    printf("Finishing successfully!\n");
}

void NMDPD::exec()
{
    for (int i=0; i<info->steps; i++)
    {
       updateVelocities();
       updatePositions();
       pairCompute();
       updateVelocities();
       dumpAndThermo(i);
       endOfStep();
    }
}

// Velocity-Verlet Algorithm
void NMDPD::updateVelocities()
{
    Ek = 0.0;
    for (int i=0; i<info->numParticles; i++)
    {
       v[i][0] += dtfm[i] * f[i][0];
       v[i][1] += dtfm[i] * f[i][1];
       v[i][2] += dtfm[i] * f[i][2];
       Ek += 0.5 * info->mass * (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    }
    Ek = 2.0 * Ek / (3.0 * info->numParticles);
}

void NMDPD::updatePositions()
{
    for (int i=0; i<info->numParticles; i++)
    {
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
        pbc(x[i]);
    }
}

void NMDPD::pairCompute()
{
    for (int i=0; i<info->numParticles; i++)
    {
        f[i][0] = 0.0;
        f[i][1] = 0.0;
        f[i][2] = 0.0;
    }
    for (int i=0; i<info->numParticles-1; i++)
    {
        for (int j=i+1; j<info->numParticles; j++)
        {
             //printf("new pair\n");
             auto pvector = parallelVector(x[i],x[j]);
             int index = (pvector[0] - info->start) / info->interval;
             double rand = rand_gauss();
             randoms[i][j] = rand;
             randoms[j][i] = -rand;
             double randV1 = rand_gauss();
             randomsV1[i][j] = randV1;
             randomsV1[j][i] = -randV1;
             double randV2 = rand_gauss();
             randomsV2[i][j] = randV2;
             randomsV2[j][i] = -randV2;
             double fij, vfij, vfx, vfy, vfz, vvx, vvy, vvz, e1x, e1y, e1z, e2x, e2y, e2z;
             
             if (index < info->tableLength && index>=0)
             {
                 fij = vfij = vvx = vvy = vvz = e1x = e1y = e1z = e2x = e2y = e2z = vfx = vfy = vfz = 0.0;
                 fij += table[index][1];
                 double paraFriction = table[index][2];
                 double projection = (v[i][0] - v[j][0]) * pvector[1] + (v[i][1]-v[j][1]) * pvector[2] + (v[i][2]-v[j][2]) * pvector[3];
                 fij += -dt * paraFriction * projection + sqrt(kT * paraFriction) * tapCoeffParallel[0] * rand;
                 // vertical friction forces
                                  
                 double verticalFriction = table[index][3];
                 vfij = -dt * verticalFriction;
                 vvx = (v[i][0] - v[j][0]) - projection * pvector[1]; 
                 vvy = (v[i][1] - v[j][1]) - projection * pvector[2]; 
                 vvz = (v[i][2] - v[j][2]) - projection * pvector[3]; 
         
                 // vertical random forces (hard part)
                 e1x = pvector[2] / sqrt(pvector[1] * pvector[1] + pvector[2] * pvector[2]);
                 e1y = - (pvector[1] / pvector[2]) * e1x;
                 e1z = 0.0;
                 //double p1 = (v[i][0] - v[j][0]) * e1x + (v[i][1] - v[j][1]) * e1y + (v[i][2] - v[j][2]) * e1z;
                 e2x = -pvector[3] * e1y;
                 e2y = pvector[3] * e1x;
                 e2z = pvector[1] * e1y - pvector[2] * e1x;
                 //double p2 = (v[i][0] - v[j][0]) * e2x + (v[i][1] - v[j][1]) * e2y + (v[i][2] - v[j][2]) * e2z;
                 double magnitude = tapCoeffVertical[0] * sqrt(kT * verticalFriction);
         
                 vfx = magnitude * (randV1 * e1x + randV2 * e2x);
                 vfy = magnitude * (randV1 * e1y + randV2 * e2y);
                 vfz = magnitude * (randV1 * e1z + randV2 * e2z);

                 //f[i][0] += fij * pvector[1] + vfij * p1 * e1x + vfij * p2 * e2x + vfx;
                 //f[i][1] += fij * pvector[2] + vfij * p1 * e1y + vfij * p2 * e2y + vfy;
                 //f[i][2] += fij * pvector[3] + vfij * p1 * e1z + vfij * p2 * e2z + vfz;
                 //f[j][0] -= (fij * pvector[1] + vfij * p1 * e1x + vfij * p2 * e2x + vfx);
                 //f[j][1] -= (fij * pvector[2] + vfij * p1 * e1y + vfij * p2 * e2y + vfy);
                 //f[j][2] -= (fij * pvector[3] + vfij * p1 * e1z + vfij * p2 * e2z + vfz);
                 f[i][0] += fij * pvector[1] + vfij * vvx + vfx;
                 f[i][1] += fij * pvector[2] + vfij * vvy + vfy;
                 f[i][2] += fij * pvector[3] + vfij * vvz + vfz;
                 f[j][0] -= (fij * pvector[1] + vfij * vvx + vfx);
                 f[j][1] -= (fij * pvector[2] + vfij * vvy + vfy);
                 f[j][2] -= (fij * pvector[3] + vfij * vvz + vfz);
             }

       

             for (int istep = 0; istep<info->decaySteps; istep++)
             {
                 fij = vfij = vvx = vvy = vvz = e1x = e1y = e1z = e2x = e2y = e2z = vfx = vfy = vfz = 0.0;
                 double deltaT = dt * (info->decaySteps - istep);
                 double tau = dt * (info->decaySteps / 5.0);
                 auto eij = parallelVector(trajFrame->positions[istep][i], trajFrame->positions[istep][j]);
                 int id = (eij[0] - info->start) / info->interval;
                 if (id < info->tableLength && id>=0)
                 {
                     double paraFriction = table[id][2];
                     double projection = (trajFrame->velocities[istep][i][0] - trajFrame->velocities[istep][j][0]) * eij[1]
                                       + (trajFrame->velocities[istep][i][1] - trajFrame->velocities[istep][j][1]) * eij[2]
                                       + (trajFrame->velocities[istep][i][2] - trajFrame->velocities[istep][j][2]) * eij[3];
                     fij = -dt * paraFriction * exp(-deltaT / tau) * projection + sqrt(kT * paraFriction) * tapCoeffParallel[info->decaySteps-istep] * trajFrame->pairRandoms[istep][i][j];

                     // vertical friction forces
                     double verticalFriction = table[id][3]; 
                     vfij = -dt * verticalFriction * exp(- deltaT / tau);
                     vvx = (trajFrame->velocities[istep][i][0] - trajFrame->velocities[istep][j][0]) - projection * eij[1]; 
                     vvy = (trajFrame->velocities[istep][i][1] - trajFrame->velocities[istep][j][1]) - projection * eij[2]; 
                     vvz = (trajFrame->velocities[istep][i][2] - trajFrame->velocities[istep][j][2]) - projection * eij[3]; 
         
                     // vertical random forces (hard part)
                     e1x = eij[2] / sqrt(eij[1] * eij[1] + eij[2] * eij[2]);
                     e1y = - (eij[1] / eij[2]) * e1x;
                     e1z = 0.0;
                     //double p1 = (trajFrame->velocities[istep][i][0] - trajFrame->velocities[istep][j][0]) * e1x
                     //                  + (trajFrame->velocities[istep][i][1] - trajFrame->velocities[istep][j][1]) * e1y
                     //                  + (trajFrame->velocities[istep][i][2] - trajFrame->velocities[istep][j][2]) * e1z;
                     e2x = -eij[3] * e1y;
                     e2y = eij[3] * e1x;
                     e2z = eij[1] * e1y - eij[2] * e1x;
                     //double p2 = (trajFrame->velocities[istep][i][0] - trajFrame->velocities[istep][j][0]) * e2x
                     //                  + (trajFrame->velocities[istep][i][1] - trajFrame->velocities[istep][j][1]) * e2y
                     //                  + (trajFrame->velocities[istep][i][2] - trajFrame->velocities[istep][j][2]) * e2z;
                     double magnitude = tapCoeffVertical[info->decaySteps-istep] * sqrt(kT * verticalFriction);
                     //printf("e1 %lf e2 %lf eij %lf\n", e1x*e1x+e1y*e1y+e1z*e1z,e2x*e2x+e2y*e2y+e2z*e2z,eij[3]*eij[3]+eij[1]*eij[1]+eij[2]*eij[2]);
                     //printf("e1 %lf e2 %lf eij %lf\n", e1x*e2x+e1y*e2y+e1z*e2z,e2x*eij[1]+e2y*eij[2]+e2z*eij[3],e1x*eij[1]+e1y*eij[2]+e1z*eij[3]);
                     //printf("vvx %lf alternative %lf\n",vvx,p1*e1x+p2*e2x);
                     //printf("step %d magnitude = %lf verticalFriction = %lf\n", istep, magnitude, verticalFriction);
         
                     vfx = magnitude * (trajFrame->pairRandomsV1[istep][i][j] * e1x + trajFrame->pairRandomsV2[istep][i][j] * e2x);
                     vfy = magnitude * (trajFrame->pairRandomsV1[istep][i][j] * e1y + trajFrame->pairRandomsV2[istep][i][j] * e2y);
                     vfz = magnitude * (trajFrame->pairRandomsV1[istep][i][j] * e1z + trajFrame->pairRandomsV2[istep][i][j] * e2z);

                     //f[i][0] += fij * eij[1] + vfij * p1 * e1x + vfij * p2 * e2x + vfx;
                     //f[i][1] += fij * eij[2] + vfij * p1 * e1y + vfij * p2 * e2y + vfy;
                     //f[i][2] += fij * eij[3] + vfij * p1 * e1z + vfij * p2 * e2z + vfz;
                     //f[j][0] -= (fij * eij[1] + vfij * p1 * e1x + vfij * p2 * e2x + vfx);
                     //f[j][1] -= (fij * eij[2] + vfij * p1 * e1y + vfij * p2 * e2y + vfy);
                     //f[j][2] -= (fij * eij[3] + vfij * p1 * e1z + vfij * p2 * e2z + vfz);
                     f[i][0] += fij * eij[1] + vfij * vvx + vfx;
                     f[i][1] += fij * eij[2] + vfij * vvy + vfy;
                     f[i][2] += fij * eij[3] + vfij * vvz + vfz;
                     f[j][0] -= (fij * eij[1] + vfij * vvx + vfx);
                     f[j][1] -= (fij * eij[2] + vfij * vvy + vfy);
                     f[j][2] -= (fij * eij[3] + vfij * vvz + vfz);
                 }
                 //if (i==0 && j==31)
                   //printf("rij = %lf id %d fij = %lf random = %lf pp = %lf\n", eij[0], id, fij, trajFrame->pairRandoms[istep][i][j], trajFrame->positions[istep][i][0]);
             }
        }
    }
}

void NMDPD::endOfStep()
{
    //pop the first frame of trajFrame, add the last one in
    trajFrame->readFrame(x,v,f,randoms,randomsV1,randomsV2);
}

void NMDPD::dumpAndThermo(int step)
{
    fprintf(output, "ITEM: TIMESTEP\n");
    fprintf(output, "%d\n", step);
    fprintf(output, "ITEM: NUMBER OF ATOMS\n");
    fprintf(output, "%d\n", info->numParticles);
    fprintf(output, "ITEM: BOX BOUNDS pp pp pp\n");
    for (int i=0; i<3; i++) fprintf(output, "0.00 %lf\n", info->boxLength); 
    fprintf(output, "ITEM: ATOMS id type x y z vx vy vz fx fy fz\n");
    for (int i=0; i<info->numParticles; i++)
    {
       fprintf(output, "%d 1 %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", i+1, x[i][0],x[i][1],x[i][2],
								   v[i][0], v[i][1], v[i][2], f[i][0], f[i][1], f[i][2]);
    }
    fflush(output);
 
    printf("Step: %d Ek = %lf\n", step, Ek);
}

void NMDPD::init_one()
{
    trajFrame->readOne(x,v,f);   
}

inline std::vector<double>  NMDPD::parallelVector(std::vector<double> x, std::vector<double> y)
{
    double L = info->boxLength;
    double dx = x[0] - y[0];
    if (dx > 0.5*L) dx -= L;
    else if (dx < -0.5*L) dx += L;
    double dy = x[1] - y[1];
    if (dy > 0.5*L) dy -= L;
    else if (dy < -0.5*L) dy += L;
    double dz = x[2] - y[2];
    if (dz > 0.5*L) dz -= L;
    else if (dz < -0.5*L) dz += L;

    double r = sqrt(dx*dx + dy*dy + dz*dz);
    std::vector<double> result(4);
    result[0] = r;
    result[1] = dx / r;
    result[2] = dy / r;
    result[3] = dz / r;
    //printf("%lf %lf %lf <--> %lf %lf %lf %lf\n", dx, dy, dz, result[0], result[1], result[2], result[3]);
    return result;
}

inline void NMDPD::pbc(std::vector<double> &x)
{
   double L = info->boxLength;
   for (int i=0; i<3; i++)
   {
       if (x[i] > L) x[i] -= L;
       else if (x[i] < 0) x[i] += L;
   }
}

inline double NMDPD::rand_gauss()
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
        while (r > 1.0);
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

#include "FitGLE.h"
#include "comm.h"
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_bspline.h>
#include <lapacke.h>

using namespace FITGLE_NS;

FitGLE::FitGLE(int argc, char** argv)
{
    if (argc != 3)
    {
        printf("./FitGLE.x [trajectory Filename] [number of Particles]\n");
    }

    assert(argc == 3);
    printf("Initializing FitGLE parameters...\n");

    // parsing the configuration parameters
    info = std::make_shared<InputParameters>();
    VAR_BEGIN
      GET_REAL(info->start)
      GET_REAL(info->end)
      GET_REAL(info->boxLength)
      GET_REAL(info->outputPrecision)
      GET_REAL(info->dt)
      GET_INT(info->splineOrder)
      GET_INT(info->numSplines)
      GET_INT(info->steps)
      GET_INT(info->decayStep)
    VAR_END

    printf("set up trajectory files\n");
    
    trajFrame = std::make_shared<Frame>(atoi(argv[2]), argv[1], info->decayStep);
    //trajFrame = new Frame(atoi(argv[2]), argv[1]);
    //
    // Initialize the Normal Equation matrix and vectors
    printf("set up containers\n");
    info->totalBasisSize = 2 * info->numSplines + 1; 
    normalVector.resize(info->totalBasisSize);
    splineCoefficients.resize(info->totalBasisSize);
    normalMatrix.resize(info->totalBasisSize);
    for (auto&& i : normalMatrix)
    {
        i.resize(info->totalBasisSize);
    }

    // Set up the size of splines according to order and numbers
    // Initialize the spline set up
    printf("set up b-spline data structures\n");
    int numBreaks = info->numSplines + 2 - info->splineOrder;
    bw = gsl_bspline_alloc(info->splineOrder, numBreaks);
    splineValue = gsl_vector_alloc(info->numSplines);
    gsl_bspline_knots_uniform(info->start, info->end, bw);
    printf("finishing configuration, entering normal equation accumulation\n");
}

FitGLE::~FitGLE()
{
    gsl_bspline_free(bw);
    gsl_vector_free(splineValue);
    printf("Exiting the Fitting GLE process...\n");
}

inline double FitGLE::distance(std::vector<double> & A, std::vector<double> & B)
{
    double dx = A[0] - B[0];
    if (dx > 0.5 * info->boxLength) dx = dx - info->boxLength;
    if (dx < -0.5 * info->boxLength) dx = dx + info->boxLength;
    double dy = A[1] - B[1];
    if (dy > 0.5 * info->boxLength) dy = dy - info->boxLength;
    if (dy < -0.5 * info->boxLength) dy = dy + info->boxLength;
    double dz = A[2] - B[2];
    if (dz > 0.5 * info->boxLength) dz = dz - info->boxLength;
    if (dz < -0.5 * info->boxLength) dz = dz + info->boxLength;
  
    return sqrt(dx*dx + dy*dy + dz*dz);
}

inline std::vector<double> FitGLE::parallelVelocity(int step, int i, int j)
{
    double dx = trajFrame->positions[step][i][0] - trajFrame->positions[step][j][0];
    if (dx > 0.5 * info->boxLength) dx = dx - info->boxLength;
    if (dx < -0.5 * info->boxLength) dx = dx + info->boxLength;
    double dy = trajFrame->positions[step][i][1] - trajFrame->positions[step][j][1];
    if (dy > 0.5 * info->boxLength) dy = dy - info->boxLength;
    if (dy < -0.5 * info->boxLength) dy = dy + info->boxLength;
    double dz = trajFrame->positions[step][i][2] - trajFrame->positions[step][j][2];
    if (dz > 0.5 * info->boxLength) dz = dz - info->boxLength;
    if (dz < -0.5 * info->boxLength) dz = dz + info->boxLength;

    double rij = sqrt(dx*dx + dy*dy + dz*dz);
    double eij[] = {dx/rij, dy/rij, dz/rij};
    //printf("%lf %lf %lf %lf %lf d %lf %lf %lf\n", eij[0], rij, dx, dy, dz, trajFrame->positions[i][0], trajFrame->positions[j][0], info->boxLength);
    std::vector<double> vij;
    std::transform(trajFrame->velocities[step][i].begin(), trajFrame->velocities[step][i].end(), trajFrame->velocities[step][j].begin(), std::back_inserter(vij), std::minus<double>());
    //printf("%lf %lf delta = %lf\n", trajFrame->velocities[i][0], trajFrame->velocities[j][0], vij[0]); 
    double projection = vij[0] * eij[0] + vij[1] * eij[1] + vij[2] * eij[2];
    vij[0] = projection * eij[0];
    vij[1] = projection * eij[1];
    vij[2] = projection * eij[2];
    return vij;
}

inline std::vector<double> FitGLE::paraVerticalVelocity(int step, int i, int j)
{
    double dx = trajFrame->positions[step][i][0] - trajFrame->positions[step][j][0];
    if (dx > 0.5 * info->boxLength) dx = dx - info->boxLength;
    if (dx < -0.5 * info->boxLength) dx = dx + info->boxLength;
    double dy = trajFrame->positions[step][i][1] - trajFrame->positions[step][j][1];
    if (dy > 0.5 * info->boxLength) dy = dy - info->boxLength;
    if (dy < -0.5 * info->boxLength) dy = dy + info->boxLength;
    double dz = trajFrame->positions[step][i][2] - trajFrame->positions[step][j][2];
    if (dz > 0.5 * info->boxLength) dz = dz - info->boxLength;
    if (dz < -0.5 * info->boxLength) dz = dz + info->boxLength;

    double rij = sqrt(dx*dx + dy*dy + dz*dz);
    double eij[] = {dx/rij, dy/rij, dz/rij};
    //printf("%lf %lf %lf %lf %lf d %lf %lf %lf\n", eij[0], rij, dx, dy, dz, trajFrame->positions[i][0], trajFrame->positions[j][0], info->boxLength);
    std::vector<double> vij;
    std::transform(trajFrame->velocities[step][i].begin(), trajFrame->velocities[step][i].end(), trajFrame->velocities[step][j].begin(), std::back_inserter(vij), std::minus<double>());
    //printf("%lf %lf delta = %lf\n", trajFrame->velocities[i][0], trajFrame->velocities[j][0], vij[0]); 
    double projection = vij[0] * eij[0] + vij[1] * eij[1] + vij[2] * eij[2];
    std::vector<double> result (6, 0.0);
    result[0] = projection * eij[0];
    result[1] = projection * eij[1];
    result[2] = projection * eij[2];
    result[3] = vij[0] - result[0];
    result[4] = vij[1] - result[1];
    result[5] = vij[2] - result[2];
    return result;
}

inline std::vector<double> FitGLE::parallelUnitVector(int step, int i, int j)
{
    double dx = trajFrame->positions[step][i][0] - trajFrame->positions[step][j][0];
    if (dx > 0.5 * info->boxLength) dx = dx - info->boxLength;
    if (dx < -0.5 * info->boxLength) dx = dx + info->boxLength;
    double dy = trajFrame->positions[step][i][1] - trajFrame->positions[step][j][1];
    if (dy > 0.5 * info->boxLength) dy = dy - info->boxLength;
    if (dy < -0.5 * info->boxLength) dy = dy + info->boxLength;
    double dz = trajFrame->positions[step][i][2] - trajFrame->positions[step][j][2];
    if (dz > 0.5 * info->boxLength) dz = dz - info->boxLength;
    if (dz < -0.5 * info->boxLength) dz = dz + info->boxLength;

    double rij = sqrt(dx*dx + dy*dy + dz*dz);
    std::vector<double> eij;
    eij.push_back(dx / rij);
    eij.push_back(dy / rij);
    eij.push_back(dz / rij);
    //printf("%lf %lf %lf d %lf %lf %lf\n", dx, dy, dz, trajFrame->positions[i][0], trajFrame->positions[j][0], info->boxLength);
    return eij;
}

// Accumulate the normal equation for this particular frame
void FitGLE::accumulateNormalEquation()
{
    int nall = trajFrame->numParticles;
    int nSplines = info->numSplines;
    std::vector<std::vector<double> > frameMatrix(info->totalBasisSize, std::vector<double>(3 * nall));
    double normalFactor = 1.0 / info->steps;
    double decays = (double)info->decayStep;
   
    // Computing Matrix F_km 
    for (int step = 0; step<5*info->decayStep; step++)
    {
        double dtexp_lagtime = info->dt * exp(-(5*decays - step -1) * info->dt / decays);
        for (int i = 0; i < nall; i++)
        {
            for (int j = i + 1; j < nall; j++)
            {
                double rij = distance(trajFrame->positions[step][i], trajFrame->positions[step][j]);
                if (rij < info->end && rij > info->start)
                {
                    gsl_bspline_eval(rij, splineValue, bw);
                    //size_t istart, iend;
                    //gsl_bspline_eval_nonezero(rij, Bk, &istart, &iend, bw);
                    //printf("rij = %lf, %d %d\n", rij, i, j);    
                    std::vector<double> dv = paraVerticalVelocity(step, i, j);
                    //Check if force matching works fine
                    //std::vector<double> dv = parallelUnitVector(i, j);
            
                    for (int m=0; m<nSplines; m++)
                    {
                         double phim = gsl_vector_get(splineValue, m);
                         if (phim < 1e-20)
		             continue;
                         // For all three dimensions
                         // Do the parallel part
                         frameMatrix[m][3*i] += -phim * dv[0] * dtexp_lagtime;
                         frameMatrix[m][3*i + 1] += -phim * dv[1] * dtexp_lagtime;
                         frameMatrix[m][3*i + 2] += -phim * dv[2] * dtexp_lagtime;
                         frameMatrix[m][3*j] -= -phim * dv[0] * dtexp_lagtime;
                         frameMatrix[m][3*j + 1] -= -phim * dv[1] * dtexp_lagtime;
                         frameMatrix[m][3*j + 2] -= -phim * dv[2] * dtexp_lagtime;

                         // Do the perpendicular part
                         frameMatrix[m + nSplines][3*i] += -phim * dv[3] * dtexp_lagtime;
                         frameMatrix[m + nSplines][3*i + 1] += -phim * dv[4] * dtexp_lagtime;
                         frameMatrix[m + nSplines][3*i + 2] += -phim * dv[5] * dtexp_lagtime;
                         frameMatrix[m + nSplines][3*j] -= -phim * dv[3] * dtexp_lagtime;
                         frameMatrix[m + nSplines][3*j + 1] -= -phim * dv[4] * dtexp_lagtime;
                         frameMatrix[m + nSplines][3*j + 2] -= -phim * dv[5] * dtexp_lagtime;
                    }
                }  
            }
        }
    }
        
    // Constructing the normal Matrix and normal Vector
    for (int m = 0; m < info->totalBasisSize; m++)
    {
        for (int n = 0; n < info->totalBasisSize; n++)
        {
            double sum = 0.0;
            for (int k = 0; k < 3 * nall; k++)
                sum += frameMatrix[m][k] * frameMatrix[n][k];
            normalMatrix[m][n] += sum * normalFactor;
        }
   
        double sum_b = 0.0; 
        for (int k = 0; k < 3 * nall; k++)
            sum_b += frameMatrix[m][k] * trajFrame->residualForces.back()[k/3][k%3];
        normalVector[m] += sum_b * normalFactor;
    }
}

void FitGLE::leastSquareSolver()
{
    // Solving the least square normal equation G*phi = b
    double* G = new double[info->totalBasisSize * info->totalBasisSize];
    double* b = new double[info->totalBasisSize];

    
    // Preconditioning the Normal Matrix
    std::vector<double> h(info->totalBasisSize, 0.0);
    
    for (int i = 0; i < info->totalBasisSize; i++) {
        for (int j = 0; j < info->totalBasisSize; j++) {
            h[j] = h[j] + normalMatrix[j][i] * normalMatrix[j][i];  //mat->dense_fm_matrix[j * mat->accumulation_matrix_rows + i] * mat->dense_fm_matrix[j * mat->accumulation_matrix_rows + i];
        }
    }

    for (int i = 0; i < info->totalBasisSize; i++) {
        if (h[i] < 1.0E-20) {
            h[i] = 1.0;
            printf("Row %d has negligible nonzero entries; normal equations are ill-formed.", i);
        }
        else {
            h[i] = 1.0 / sqrt(h[i]);
        }
    }
    for (int i = 0; i < info->totalBasisSize; i++)
    {
        for (int j = 0; j < info->totalBasisSize; j++)
           normalMatrix[i][j] *= h[j];
    }


    // Store the normalMatrix in container 
    for (int i = 0; i < info->totalBasisSize; i++)
    {
        for (int j = 0; j < info->totalBasisSize; j++)
        {
            G[ i * info->totalBasisSize + j] = normalMatrix[i][j];
        }
        b[i] = normalVector[i];
        printf("m %d %lf\n", i, b[i]);
    }
    

    // Solving the linear system using SVD decomposition

    int m = info->totalBasisSize;
    int n = info->totalBasisSize;
    int nrhs = 1;
    int lda = info->totalBasisSize;
    int ldb = 1;
    double rcond = -1.0;
    int irank;
    double* singularValue = new double[info->totalBasisSize];
    int solverInfo = LAPACKE_dgelss(LAPACK_ROW_MAJOR, m, n, nrhs, G, lda, b, ldb, singularValue, rcond, &irank); 
   
    printf("LSQ Solver Info: %d\n", solverInfo);

    for (int m = 0; m < info->totalBasisSize; m++)
        splineCoefficients[m] = b[m] * h[m];

    delete[] G;
    delete[] b;
}

// Output function for gamma(R)
void FitGLE::output()
{
    printf("output\n");
    double start = info->start;
    double end = info->end;
    double precision = info->outputPrecision;

    FILE* fb = fopen("spline_coeff.dat", "w");
    for (int m=0; m < info->totalBasisSize; m++)
    {
        fprintf(fb, "%lf\n", splineCoefficients[m]);
    }
    fclose(fb);

    FILE* fp = fopen("gamma_out.dat", "w");
    FILE* fv = fopen("gamma_vet.dat", "w");
    while (start < end)
    {
        double gamma_r = 0.0;
        double gamma_v = 0.0;
        gsl_bspline_eval(start, splineValue, bw);
        for (int m = 0; m < info->numSplines; m++)
        {
           gamma_r += splineCoefficients[m] * gsl_vector_get(splineValue, m);
           gamma_v += splineCoefficients[m + info->numSplines] * gsl_vector_get(splineValue, m);
        }
        fprintf(fp, "%lf\t%lf\n", start, gamma_r);
        fprintf(fv, "%lf\t%lf\n", start, gamma_v);
        start = start + precision;
    }

    // Print the Langevin friction value from fit
    // ATTN YINING
    printf("%lf", splineCoefficients[info->totalBasisSize - 1]);
    
    fclose(fp);
    fclose(fv);
}

// Execution Process
void FitGLE::exec()
{
    printf("Accumulating the LSQ normal Matrix\n");
    for (int i=0; i < info->steps-1; i++)
    {
        accumulateNormalEquation();
        trajFrame->readFrame();
        printf("finishing step %d (total %d)\r", i+1, info->steps);
    }
    printf("\n");
    leastSquareSolver();
    output();
}

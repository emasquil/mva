// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure
// Normalization factor for points
static const float NORMALIZATION_FACTOR = 0.0001f;

struct Match
{
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color, 2> I1, Image<Color, 2> I2,
              vector<Match> &matches)
{
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0, 0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(), 0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0 * 100.0;
    for (size_t i = 0; i < feats1.size(); i++)
    {
        SIFTDetector::Feature f1 = feats1[i];
        for (size_t j = 0; j < feats2.size(); j++)
        {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if (d < MAX_DISTANCE)
            {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

// Compute F for with exactly 8 points
FMatrix<float, 3, 3> solveFsquare(vector<Match> &matches)
{
    FMatrix<float, 3, 3> NORMALIZATION_MATRIX;
    NORMALIZATION_MATRIX.fill(0);
    NORMALIZATION_MATRIX(0, 0) = NORMALIZATION_FACTOR;
    NORMALIZATION_MATRIX(1, 1) = NORMALIZATION_FACTOR;
    NORMALIZATION_MATRIX(2, 2) = 1;

    FMatrix<float, 3, 3> F;
    // Define A matrix
    FMatrix<float, 9, 9> A;
    // Fill A matrix
    A.fill(0);
    for (int i = 0; i < 9; ++i)
    {
        Match match = matches[i];
        float x1 = match.x1 * NORMALIZATION_FACTOR;
        float x2 = match.x2 * NORMALIZATION_FACTOR;
        float y1 = match.y1 * NORMALIZATION_FACTOR;
        float y2 = match.y2 * NORMALIZATION_FACTOR;
        A(i, 0) = x1 * x2;
        A(i, 1) = x1 * y2;
        A(i, 2) = x1;
        A(i, 3) = y1 * x2;
        A(i, 4) = y1 * y2;
        A(i, 5) = y1;
        A(i, 6) = x2;
        A(i, 7) = y2;
        A(i, 8) = 1;
    }
    // Linear solution
    FVector<float, 9> S;
    FMatrix<float, 9, 9> U, Vt;
    svd(A, U, S, Vt);
    F(0, 0) = Vt(8, 0);
    F(0, 1) = Vt(8, 1);
    F(0, 2) = Vt(8, 2);
    F(1, 0) = Vt(8, 3);
    F(1, 1) = Vt(8, 4);
    F(1, 2) = Vt(8, 5);
    F(2, 0) = Vt(8, 6);
    F(2, 1) = Vt(8, 7);
    F(2, 2) = Vt(8, 8);
    // Constraint enforcement (SVD and setting last sv to 0)
    FVector<float, 3> S2;
    FMatrix<float, 3, 3> U2, Vt2;
    svd(F, U2, S2, Vt2);
    FMatrix<float, 3, 3> S3;
    S3.fill(0);
    S3(0, 0) = S2[0];
    S3(1, 1) = S2[1];
    F = NORMALIZATION_MATRIX * U2 * S3 * Vt2 * NORMALIZATION_MATRIX;
    return F;
}

// Compute F for with more points
FMatrix<float, 3, 3> solveFrectangular(vector<Match> &matches)
{
    FMatrix<float, 3, 3> NORMALIZATION_MATRIX;
    NORMALIZATION_MATRIX.fill(0);
    NORMALIZATION_MATRIX(0, 0) = NORMALIZATION_FACTOR;
    NORMALIZATION_MATRIX(1, 1) = NORMALIZATION_FACTOR;
    NORMALIZATION_MATRIX(2, 2) = 1;

    FMatrix<float, 3, 3> F;
    // Define A matrix
    Matrix<float> A(matches.size(), 9);
    // Fill A matrix
    for (int i = 0; i < matches.size(); ++i)
    {
        Match match = matches[i];
        float x1 = match.x1 * NORMALIZATION_FACTOR;
        float x2 = match.x2 * NORMALIZATION_FACTOR;
        float y1 = match.y1 * NORMALIZATION_FACTOR;
        float y2 = match.y2 * NORMALIZATION_FACTOR;
        A(i, 0) = x1 * x2;
        A(i, 1) = x1 * y2;
        A(i, 2) = x1;
        A(i, 3) = y1 * x2;
        A(i, 4) = y1 * y2;
        A(i, 5) = y1;
        A(i, 6) = x2;
        A(i, 7) = y2;
        A(i, 8) = 1;
    }
    // Linear solution
    Vector<float> S(9);
    Matrix<float> U(matches.size(), matches.size());
    Matrix<float> Vt(9, 9);
    svd(A, U, S, Vt, false);
    F(0, 0) = Vt(8, 0);
    F(0, 1) = Vt(8, 1);
    F(0, 2) = Vt(8, 2);
    F(1, 0) = Vt(8, 3);
    F(1, 1) = Vt(8, 4);
    F(1, 2) = Vt(8, 5);
    F(2, 0) = Vt(8, 6);
    F(2, 1) = Vt(8, 7);
    F(2, 2) = Vt(8, 8);
    // Constraint enforcement (SVD and setting last sv to 0)
    FVector<float, 3> S2;
    FMatrix<float, 3, 3> U2, Vt2;
    svd(F, U2, S2, Vt2);
    FMatrix<float, 3, 3> S3;
    S3.fill(0);
    S3(0, 0) = S2[0];
    S3(1, 1) = S2[1];
    F = U2 * S3 * Vt2;
    return NORMALIZATION_MATRIX * F * NORMALIZATION_MATRIX;
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float, 3, 3> computeF(vector<Match> &matches)
{
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter = 100000;         // Adjusted dynamically
    FMatrix<float, 3, 3> bestF;
    vector<int> bestInliers;
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF POINTS
    // RANSAC loop
    int ransacIter = 0;
    while (ransacIter < Niter)
    {
        ransacIter += 1;
        cout << "Iteration: " << ransacIter << endl;
        // Sampling 8 points for RANSAC
        vector<Match> selectedMatches;
        vector<int> selectedIndexes;
        while (selectedMatches.size() < 8)
        {
            int index = intRandom(0, matches.size());
            Match candidateMatch = matches[index];
            if (find(selectedIndexes.begin(), selectedIndexes.end(), index) == selectedIndexes.end())
            {
                selectedMatches.push_back(candidateMatch);
                selectedIndexes.push_back(index);
            }
        }
        // Compute F
        FMatrix<float, 3, 3> currentF = solveFsquare(selectedMatches);
        vector<int> currentInLiers;
        // Compute the distance of all points and select inliers
        for (int i = 0; i < matches.size(); ++i)
        {
            Match match = matches[i];
            FVector<float, 3> X2{match.x2, match.y2, 1};
            FVector<float, 3> X1{match.x1, match.y1, 1};
            FVector<float, 3> ftX1 = transpose(currentF) * X1;
            auto numerator = abs(ftX1 * X2);
            auto denominator = sqrt(ftX1[0] * ftX1[0] + ftX1[1] * ftX1[1]);
            auto distance = numerator / denominator;
            if (distance <= distMax)
            {
                currentInLiers.push_back(i);
            }
        }
        if (currentInLiers.size() > bestInliers.size())
        {
            bestInliers = currentInLiers;
            bestF = currentF;
            // Update Niter
            float m = bestInliers.size();
            float n = matches.size();
            long int newNiter = (log(BETA) / log(1 - pow(m / n, 8)));
            // Only updating Niter if having enough inliers (to avoid numerical problems with Niter)
            if (m > 80)
            {
                Niter = newNiter;
            }
            cout << m << endl;
            cout << n << endl;
            cout << "Updated best F" << endl;
            cout << "# inliers: " << bestInliers.size() << endl;
            cout << "Niter: " << Niter << endl;
        }
        cout << "----------------" << endl;
    }
    cout << "Best F obtained with RANSAC: " << endl;
    cout << bestF << endl;

    // Updating matches with inliers only
    vector<Match> all = matches;
    matches.clear();
    for (size_t i = 0; i < bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    // One last refinement using all inliers
    bestF = solveFrectangular(matches);
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float, 3, 3> &F)
{
    while (true)
    {
        int x, y;
        if (getMouse(x, y) == 3)
            break;
        // --------------- TODO ------------
        IntPoint2 point{x, y};
        drawCircle(point, 2, RED, 2);
        FVector<float, 3> X;
        FVector<float, 3> epipolarLine;
        float xLine0, xLine1, yLine0, yLine1;
        // Decide if it's left or right image
        if (x <= I1.width())
        {
            // User clicked on left image
            cout << "User clicked on left image" << endl;
            X[0] = x;
            X[1] = y;
            X[2] = 1;
            epipolarLine = transpose(F) * X;
            // Points passing through the line
            xLine0 = 0;
            xLine1 = I2.width();
            yLine0 = (-epipolarLine[2] - epipolarLine[0] * xLine0) / epipolarLine[1];
            yLine1 = (-epipolarLine[2] - epipolarLine[0] * xLine1) / epipolarLine[1];
            // Traslating x coordinates to the second image
            drawLine(xLine0 + I1.width(), yLine0, xLine1 + I1.width(), yLine1, BLUE, 1);
        }
        else
        {
            cout << "User clicked on right image" << endl;
            // Traslating X to the first image
            X[0] = x - I1.width();
            X[1] = y;
            X[2] = 1;
            epipolarLine = F * X;
            // Points passing through the line
            xLine0 = 0;
            xLine1 = I1.width();
            yLine0 = (-epipolarLine[2] - epipolarLine[0] * xLine0) / epipolarLine[1];
            yLine1 = (-epipolarLine[2] - epipolarLine[0] * xLine1) / epipolarLine[1];
            drawLine(xLine0, yLine0, xLine1, yLine1, BLUE, 1);
        }
    }
}

int main(int argc, char *argv[])
{
    srand((unsigned int)time(0));

    const char *s1 = argc > 1 ? argv[1] : srcPath("im1.jpg");
    const char *s2 = argc > 2 ? argv[2] : srcPath("im2.jpg");

    // Load and display images
    Image<Color, 2> I1, I2;
    if (!load(I1, s1) ||
        !load(I2, s2))
    {
        cerr << "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2 * w, I1.height());
    display(I1, 0, 0);
    display(I2, w, 0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    const int n = (int)matches.size();
    cout << " matches: " << n << endl;
    drawString(100, 20, std::to_string(n) + " matches", RED);
    click();

    FMatrix<float, 3, 3> F = computeF(matches);
    cout << "F=" << endl
         << F;

    // Redisplay with matches
    display(I1, 0, 0);
    display(I2, w, 0);
    for (size_t i = 0; i < matches.size(); i++)
    {
        Color c(rand() % 256, rand() % 256, rand() % 256);
        fillCircle(matches[i].x1 + 0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2 + w, matches[i].y2, 2, c);
    }
    drawString(100, 20, to_string(matches.size()) + "/" + to_string(n) + " inliers", RED);
    click();

    // Redisplay without SIFT points
    display(I1, 0, 0);
    display(I2, w, 0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}

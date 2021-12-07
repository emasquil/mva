// Imagine++ project
// Project:  Panorama
// Author:   Pascal Monasse
// Date:     2013/10/08

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <sstream>
using namespace Imagine;
using namespace std;

// Record clicks in two images, until right button click
void getClicks(Window w1, Window w2,
               vector<IntPoint2> &pts1, vector<IntPoint2> &pts2)
{
    // ------------- TODO/A completer ----------
    int button = 1;
    int sw;
    // Initialize variable to know which window was clicked
    Window clicked_window;
    // While don't right click, record points
    while (button != 3)
    {
        IntPoint2 clicked_point;
        button = anyGetMouse(clicked_point, clicked_window, sw);

        if (clicked_window == w1)
        {
            pts1.push_back(clicked_point);
        }
        else if (clicked_window == w2)
        {
            pts2.push_back(clicked_point);
        }
    }
}

// Return homography compatible with point matches
Matrix<float> getHomography(const vector<IntPoint2> &pts1,
                            const vector<IntPoint2> &pts2)
{
    size_t n = min(pts1.size(), pts2.size());
    if (n < 4)
    {
        cout << "Not enough correspondences: " << n << endl;
        return Matrix<float>::Identity(3);
    }
    Matrix<double> A(2 * n, 8);
    Vector<double> B(2 * n);
    // ------------- TODO/A completer ----------
    // For each pair of points, i'll add a pair of equations for H
    for (int i = 0; i < n; ++i)
    {
        double x = pts1[i][0];
        double y = pts1[i][1];
        double x_tilde = pts2[i][0];
        double y_tilde = pts2[i][1];
        // A even rows
        A(2 * i, 0) = x;
        A(2 * i, 1) = y;
        A(2 * i, 2) = 1;
        A(2 * i, 3) = 0;
        A(2 * i, 4) = 0;
        A(2 * i, 5) = 0;
        A(2 * i, 6) = -x * x_tilde;
        A(2 * i, 7) = -x_tilde * y;
        // A odd rows
        A(2 * i + 1, 0) = 0;
        A(2 * i + 1, 1) = 0;
        A(2 * i + 1, 2) = 0;
        A(2 * i + 1, 3) = x;
        A(2 * i + 1, 4) = y;
        A(2 * i + 1, 5) = 1;
        A(2 * i + 1, 6) = -y_tilde * x;
        A(2 * i + 1, 7) = -y_tilde * y;
        // B even entries
        B[2 * i] = x_tilde;
        // B odd entries
        B[2 * i + 1] = y_tilde;
    }

    B = linSolve(A, B);
    Matrix<float> H(3, 3);
    H(0, 0) = B[0];
    H(0, 1) = B[1];
    H(0, 2) = B[2];
    H(1, 0) = B[3];
    H(1, 1) = B[4];
    H(1, 2) = B[5];
    H(2, 0) = B[6];
    H(2, 1) = B[7];
    H(2, 2) = 1;

    // Sanity check
    for (size_t i = 0; i < n; i++)
    {
        float v1[] = {(float)pts1[i].x(), (float)pts1[i].y(), 1.0f};
        float v2[] = {(float)pts2[i].x(), (float)pts2[i].y(), 1.0f};
        Vector<float> x1(v1, 3);
        Vector<float> x2(v2, 3);
        x1 = H * x1;
        cout << x1[1] * x2[2] - x1[2] * x2[1] << ' '
             << x1[2] * x2[0] - x1[0] * x2[2] << ' '
             << x1[0] * x2[1] - x1[1] * x2[0] << endl;
    }
    return H;
}

// Grow rectangle of corners (x0,y0) and (x1,y1) to include (x,y)
void growTo(float &x0, float &y0, float &x1, float &y1, float x, float y)
{
    if (x < x0)
        x0 = x;
    if (x > x1)
        x1 = x;
    if (y < y0)
        y0 = y;
    if (y > y1)
        y1 = y;
}

// Panorama construction
void panorama(const Image<Color, 2> &I1, const Image<Color, 2> &I2,
              Matrix<float> H)
{
    Vector<float> v(3);
    float x0 = 0, y0 = 0, x1 = I2.width(), y1 = I2.height();

    v[0] = 0;
    v[1] = 0;
    v[2] = 1;
    v = H * v;
    v /= v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0] = I1.width();
    v[1] = 0;
    v[2] = 1;
    v = H * v;
    v /= v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0] = I1.width();
    v[1] = I1.height();
    v[2] = 1;
    v = H * v;
    v /= v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0] = 0;
    v[1] = I1.height();
    v[2] = 1;
    v = H * v;
    v /= v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    cout << "x0 x1 y0 y1=" << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1 << endl;

    Image<Color> I(int(x1 - x0), int(y1 - y0));
    setActiveWindow(openWindow(I.width(), I.height()));
    I.fill(WHITE);
    // ------------- TODO/A completer ----------

    /*Compute inverse matrix because we're going to apply the inverse transformation
    on the panorama image and then pull the pixels from the original */
    Matrix<float> H_inv = inverse(H);

    // Traverse all the destination image and pull the pixels from the originals
    for (int i = 0; i < I.width(); ++i)
    {
        for (int j = 0; j < I.height(); ++j)
        {
            // Homogeneus coordinates
            Vector<float> t(3);
            // Traslation of the coordinate center due to growTo
            t[0] = i + x0;
            t[1] = j + y0;
            t[2] = 1;
            // Map v to I1
            Vector<float> v1(3);
            v1 = H_inv * t;
            // Back to inhomogeneus coordinates
            v1 /= v1[2];
            // Check if the point lays in I1
            bool in_t1 = false;
            if ((v1[0] > 0) && (v1[1] > 0) && v1[0] < I1.width() && v1[1] < I1.height())
            {
                in_t1 = true;
                I(i, j) = I1.interpolate(v1[0], v1[1]);
            }
            // Check if the point lays in I2
            if ((t[0] > 0) && (t[1] > 0) && t[0] < I2.width() && t[1] < I2.height())
            {
                // If also lays in T1, take the average of the pixels from both images
                if (in_t1)
                {
                    I(i, j) /= 2.;
                    auto aux = I2(t[0], t[1]);
                    aux /= 2.;
                    I(i, j) += aux;
                }
                // If it's only in T2, just assign that value
                else
                {
                    I(i, j) = I2(t[0], t[1]);
                }
            }
        }
    }

    display(I, 0, 0);
}

// Main function
int main(int argc, char *argv[])
{
    const char *s1 = argc > 1 ? argv[1] : srcPath("image0006.jpg");
    const char *s2 = argc > 2 ? argv[2] : srcPath("image0007.jpg");

    // Load and display images
    Image<Color> I1, I2;
    if (!load(I1, s1) ||
        !load(I2, s2))
    {
        cerr << "Unable to load the images" << endl;
        return 1;
    }
    Window w1 = openWindow(I1.width(), I1.height(), s1);
    display(I1, 0, 0);
    Window w2 = openWindow(I2.width(), I2.height(), s2);
    setActiveWindow(w2);
    display(I2, 0, 0);

    // Get user's clicks in images
    vector<IntPoint2> pts1, pts2;
    getClicks(w1, w2, pts1, pts2);

    vector<IntPoint2>::const_iterator it;
    cout << "pts1=" << endl;
    for (it = pts1.begin(); it != pts1.end(); it++)
        cout << *it << endl;
    cout << "pts2=" << endl;
    for (it = pts2.begin(); it != pts2.end(); it++)
        cout << *it << endl;

    // Compute homography
    Matrix<float> H = getHomography(pts1, pts2);
    cout << "H=" << H / H(2, 2);

    // Apply homography
    panorama(I1, I2, H);

    endGraphics();
    return 0;
}

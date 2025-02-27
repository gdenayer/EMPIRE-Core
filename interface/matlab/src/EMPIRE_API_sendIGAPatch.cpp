#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "HelperFunctions.h"

#include <iostream>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==10);
    assert(nlhs==0);

#define P_DEGREE_IN                 prhs[0]
#define U_NUM_KNOTS_IN              prhs[1]
#define U_KNOT_VECTOR_IN            prhs[2]
#define Q_DEGREE_IN                 prhs[3]
#define V_NUM_KNOTS_IN              prhs[4]
#define V_KNOT_VECTOR_IN            prhs[5]
#define U_NUM_CONTROL_PTS_IN        prhs[6]
#define V_NUM_CONTROL_PTS_IN        prhs[7]
#define CP_NET_IN                   prhs[8]
#define NODE_NET_IN                  prhs[9]


    // p degree
    assert(mxIsDouble(P_DEGREE_IN));
    assert(mxGetNumberOfElements(P_DEGREE_IN) == 1);
    int pDegree = (int) mxGetPr(P_DEGREE_IN)[0]; // cast from double to int

    // u number of knots
    assert(mxIsDouble(U_NUM_KNOTS_IN));
    assert(mxGetNumberOfElements(U_NUM_KNOTS_IN) == 1);
    int uNumKnots = (int) mxGetPr(U_NUM_KNOTS_IN)[0]; // cast from double to int

    // u knot vector
    assert(mxIsDouble(U_KNOT_VECTOR_IN));
    assert(mxGetNumberOfElements(U_KNOT_VECTOR_IN) == uNumKnots);
    double *uKnotVector = mxGetPr(U_KNOT_VECTOR_IN);

    // q degree
    assert(mxIsDouble(Q_DEGREE_IN));
    assert(mxGetNumberOfElements(Q_DEGREE_IN) == 1);
    int qDegree = (int) mxGetPr(Q_DEGREE_IN)[0]; // cast from double to int

    // v number of knots
    assert(mxIsDouble(V_NUM_KNOTS_IN));
    assert(mxGetNumberOfElements(V_NUM_KNOTS_IN) == 1);
    int vNumKnots = (int) mxGetPr(V_NUM_KNOTS_IN)[0]; // cast from double to int

    // v knot vector
    assert(mxIsDouble(V_KNOT_VECTOR_IN));
    assert(mxGetNumberOfElements(V_KNOT_VECTOR_IN) == vNumKnots);
    double *vKnotVector = mxGetPr(V_KNOT_VECTOR_IN);

    // u number of control points
    assert(mxIsDouble(U_NUM_CONTROL_PTS_IN));
    assert(mxGetNumberOfElements(U_NUM_CONTROL_PTS_IN) == 1);
    int uNumControlPoints = (int) mxGetPr(U_NUM_CONTROL_PTS_IN)[0]; // cast from double to int

    // v number of control points
    assert(mxIsDouble(V_NUM_CONTROL_PTS_IN));
    assert(mxGetNumberOfElements(V_NUM_CONTROL_PTS_IN) == 1);
    int vNumControlPoints = (int) mxGetPr(V_NUM_CONTROL_PTS_IN)[0]; // cast from double to int

    // coordinates of control points
    assert(mxIsDouble(CP_NET_IN));
    assert(mxGetNumberOfElements(CP_NET_IN) == uNumControlPoints * vNumControlPoints * 4);
    double *cpNet = mxGetPr(CP_NET_IN);

    // DOF ids
    assert(mxIsDouble(NODE_NET_IN));
    assert(mxGetNumberOfElements(NODE_NET_IN) == uNumControlPoints * vNumControlPoints);
    int *nodeNet = doubleArrayToIntArray(mxGetPr(NODE_NET_IN), uNumControlPoints * vNumControlPoints);

	// Send the patch through API to Empire
    EMPIRE_API_sendIGAPatch(pDegree,  uNumKnots, uKnotVector, qDegree, vNumKnots,
            vKnotVector, uNumControlPoints, vNumControlPoints, cpNet, nodeNet);

    delete[] nodeNet;

#undef P_DEGREE_IN
#undef U_NUM_KNOTS_IN
#undef U_KNOT_VECTOR_IN
#undef Q_DEGREE_IN
#undef V_NUM_KNOTS_IN
#undef V_KNOT_VECTOR_IN
#undef U_NUM_CONTROL_PTS_IN
#undef V_NUM_CONTROL_PTS_IN
#undef CP_NET_IN
#undef NODE_NET_IN
}

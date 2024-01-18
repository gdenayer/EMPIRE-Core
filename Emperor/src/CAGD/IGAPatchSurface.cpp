/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EMPIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
 */
// Inclusion of standard libraries
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits>
#include <algorithm>

// Inclusion of user defined libraries
#include "IGAPatchSurface.h"
#include "MathLibrary.h"
#include "MatrixVectorMath.h"
#include "DataField.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

int IGAPatchSurface::MAX_NUM_ITERATIONS = 20;
int IGAPatchSurface::REL_MAX_NUM_ITERATIONS = 50;

double IGAPatchSurface::TOL_ORTHOGONALITY = 1e-9;
double IGAPatchSurface::REL_TOL_ORTHOGONALITY = 1e-1;

double IGAPatchSurface::TOL_DISTANCE = 1e-5;
double IGAPatchSurface::REL_TOL_DISTANCE = 1e-1;

double IGAPatchSurface::TOL_PARAMETER_ON_STRAIGHT_EDGE = 1e-3;
double IGAPatchSurface::TOL_PARAMETER_ON_PATCH_PARAMETER_SPACE = 1e-6;

const char IGAPatchSurface::EDGE_U0 = 1 << 0;
const char IGAPatchSurface::EDGE_UN = 1 << 1;
const char IGAPatchSurface::EDGE_V0 = 1 << 2;
const char IGAPatchSurface::EDGE_VN = 1 << 3;
const char IGAPatchSurface::EDGES[4] = {EDGE_U0, EDGE_UN, EDGE_V0, EDGE_VN};

IGAPatchSurface::IGAPatchSurface(int _IDBasis, int _pDegree, int _uNoKnots, double* _uKnotVector,
        int _qDegree, int _vNoKnots, double* _vKnotVector, int _uNoControlPoints,
        int _vNoControlPoints, IGAControlPoint** _controlPointNet) :
        uNoControlPoints(_uNoControlPoints), vNoControlPoints(_vNoControlPoints) {

    // Read input
    bool ucondition = _uNoControlPoints != _uNoKnots - _pDegree - 1;
    bool vcondition = _vNoControlPoints != _vNoKnots - _qDegree - 1;

    if (ucondition || vcondition) {
        ERROR_OUT() << " in IGAPatchSurface::IGAPatchSurface" << endl;
        ERROR_OUT()
                << "Number of Control Points, number of knots and polynomial degree do not match!"
                << endl;
        exit(-1);
    }

    // Figure out whether the patch has a B-Spline or a NURBS underlying basis
    int isNurbs = 0;
    int counter = 0;
    for (int j = 0; j < vNoControlPoints; j++) {
        for (int i = 0; i < uNoControlPoints; i++) {
            if (_controlPointNet[counter]->getW() != 1.0) {
                isNurbs = 1;
                break;
            }
            // Update the counter
            counter++;
        }
    }

    // Create the NURBS or the B-Spline underlying basis
    if (!isNurbs) {
        IGABasis = new BSplineBasis2D(_IDBasis, _pDegree, _uNoKnots, _uKnotVector, _qDegree,
                _vNoKnots, _vKnotVector);
    } else {
        double* controlPointWeights = new double[uNoControlPoints * vNoControlPoints];
        for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
            controlPointWeights[i] = _controlPointNet[i]->getW();
        IGABasis = new NurbsBasis2D(_IDBasis, _pDegree, _uNoKnots, _uKnotVector, _qDegree,
                _vNoKnots, _vKnotVector, _uNoControlPoints, _vNoControlPoints, controlPointWeights);
    }

    // On the Control Point net
    assert(_controlPointNet != NULL);
    ControlPointNet = _controlPointNet;
}

IGAPatchSurface::~IGAPatchSurface() {
    delete IGABasis;
    delete[] ControlPointNet;
}

void IGAPatchSurface::computeBoundingBox() {
    if (boundingBox.isComputed())
        return;
    boundingBox[0] = getControlPointNet()[0]->getX();
    boundingBox[1] = getControlPointNet()[0]->getX();
    boundingBox[2] = getControlPointNet()[0]->getY();
    boundingBox[3] = getControlPointNet()[0]->getY();
    boundingBox[4] = getControlPointNet()[0]->getZ();
    boundingBox[5] = getControlPointNet()[0]->getZ();
    for (int cpCount = 0; cpCount < getNoControlPoints(); cpCount++) {
        double x = getControlPointNet()[cpCount]->getX();
        double y = getControlPointNet()[cpCount]->getY();
        double z = getControlPointNet()[cpCount]->getZ();
        if (x < boundingBox[0])
            boundingBox[0] = x;
        else if (x > boundingBox[1])
            boundingBox[1] = x;
        if (y < boundingBox[2])
            boundingBox[2] = y;
        else if (y > boundingBox[3])
            boundingBox[3] = y;
        if (z < boundingBox[4])
            boundingBox[4] = z;
        else if (z > boundingBox[5])
            boundingBox[5] = z;
    }
    boundingBox.isComputed(true);
}

double IGAPatchSurface::computePostprocessingScalarValue(double _u, double _v,
        double* _valuesOnCP) {
    /*
     *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
     *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
     */
    // Read input
    assert(_valuesOnCP != NULL);

    int spanU = IGABasis->getUBSplineBasis1D()->findKnotSpan(_u);
    int spanV = IGABasis->getVBSplineBasis1D()->findKnotSpan(_v);

    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);
    double* localBasisFunctions = new double[noLocalBasisFunctions];

    IGABasis->computeLocalBasisFunctions(localBasisFunctions, _u, spanU, _v, spanV);

    // Initialize the Control Point index
    int CPindex = 0;

    // Initialize a basis functions counter
    int counter_basis = 0;

    double result = 0.0;

    // Loop over all the non-zero contributions
    for (int j = 0; j <= qDegree; j++) {
        for (int i = 0; i <= pDegree; i++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            CPindex = (spanV - qDegree + j) * uNoControlPoints + (spanU - pDegree + i);

            //
            result += localBasisFunctions[counter_basis] * _valuesOnCP[CPindex];

            // Update basis function's counter
            counter_basis++;
        }
    }

    // Free the memory from the heap
    delete[] localBasisFunctions;

    return result;
}

void IGAPatchSurface::addTrimLoop(int inner, int numCurves) {
    Trimming.addTrimLoop(inner, numCurves);
}

void IGAPatchSurface::addTrimCurve(int _direction, int _pDegree, int _uNoKnots, double* _uKnotVector,
                  int _uNoControlPoints, double* _controlPointNet) {
    int IDBasis = 0; ///???

    Trimming.addTrimCurve(_direction, IDBasis, _pDegree, _uNoKnots, _uKnotVector,
                                               _uNoControlPoints, _controlPointNet);
}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates, double _uPrm,
        int _uKnotSpanIndex, double _vPrm, int _vKnotSpanIndex) const{

    /*
     *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
     *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
     */

    // Read input
    assert(_cartesianCoordinates != NULL);

    // Initialize the coordinates of the point
    for (int i = 0; i < 3; i++)
        _cartesianCoordinates[i] = 0;

    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);
    double* localBasisFunctions = new double[noLocalBasisFunctions];

    IGABasis->computeLocalBasisFunctions(localBasisFunctions, _uPrm, _uKnotSpanIndex, _vPrm,
            _vKnotSpanIndex);

    // Initialize the Control Point index
    int CPindex = 0;

    // Initialize a basis functions counter
    int counter_basis = 0;

    // Loop over all the non-zero contributions
    for (int j = 0; j <= qDegree; j++) {
        for (int i = 0; i <= pDegree; i++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            CPindex = (_vKnotSpanIndex - qDegree + j) * uNoControlPoints
                    + (_uKnotSpanIndex - pDegree + i);

            // Compute iteratively the x-coordinate of the point
            _cartesianCoordinates[0] += localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getX();
            // Compute iteratively the y-coordinate of the point
            _cartesianCoordinates[1] += localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getY();
            // Compute iteratively the z-coordinate of the point
            _cartesianCoordinates[2] += localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getZ();

            // Update basis function's counter
            counter_basis++;
        }
    }

    // Free the memory from the heap
    delete[] localBasisFunctions;
}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates,
        double* _localBasisFunctions, int _uKnotSpanIndex, int _vKnotSpanIndex) const{

    /*
     *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
     *  It is also expected that the local basis functions have been precomputed outside the scope of this function and are given as arguments.
     *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
     */

    // Read input
    assert(_cartesianCoordinates != NULL);
    assert(_localBasisFunctions != NULL);

    // Initialize the coordinates of the point
    for (int i = 0; i < 3; i++)
        _cartesianCoordinates[i] = 0;

    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);

    // Initialize the Control Point index
    int CPindex = 0;

    // Initialize a basis functions counter
    int counter_basis = 0;

    // Loop over all the non-zero contributions
    for (int j = 0; j <= qDegree; j++) {
        for (int i = 0; i <= pDegree; i++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            CPindex = (_vKnotSpanIndex - qDegree + j) * uNoControlPoints
                    + (_uKnotSpanIndex - pDegree + i);

            // Compute iteratively the x-coordinate of the point
            _cartesianCoordinates[0] += _localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getX();
            // Compute iteratively the y-coordinate of the point
            _cartesianCoordinates[1] += _localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getY();
            // Compute iteratively the z-coordinate of the point
            _cartesianCoordinates[2] += _localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getZ();

            // Update basis function's counter
            counter_basis++;
        }
    }
}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates,
        double* _localCoordinates) const {

    int _uKnotSpanIndex = IGABasis->getUBSplineBasis1D()->findKnotSpan(_localCoordinates[0]);
    int _vKnotSpanIndex = IGABasis->getVBSplineBasis1D()->findKnotSpan(_localCoordinates[1]);
    IGAPatchSurface::computeCartesianCoordinates(_cartesianCoordinates, _localCoordinates[0],
            _uKnotSpanIndex, _localCoordinates[1], _vKnotSpanIndex);
}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates,
        double* _localBasisFctsAndDerivs, int _derivDegree, int _uKnotSpanIndex,
        int _vKnotSpanIndex) const {

    /*
     *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
     *  It is also expected that the local basis functions have been precomputed outside the scope of this function and are given as arguments.
     *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
     */

    // Read input
    assert(_cartesianCoordinates != NULL);
    assert(_localBasisFctsAndDerivs != NULL);

    // Initialize the coordinates of the point
    for (int i = 0; i < 3; i++)
        _cartesianCoordinates[i] = 0.0;

    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);

    // Initialize the Control Point index
    int CPindex = 0;

    // The basis function index
    int indexBasis = 0;

    // Initialize a basis functions counter
    int counter_basis = 0;

    // The derivative index so that it is exploited only the basis functions themselves
    int derivIndex = 0;

    // Loop over all the non-zero contributions
    for (int j = 0; j <= qDegree; j++) {
        for (int i = 0; i <= pDegree; i++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            CPindex = (_vKnotSpanIndex - qDegree + j) * uNoControlPoints
                    + (_uKnotSpanIndex - pDegree + i);

            // Update the basis function index
            indexBasis = IGABasis->indexDerivativeBasisFunction(_derivDegree, derivIndex,
                    derivIndex, counter_basis);
            // Compute iteratively the x-coordinate of the point
            _cartesianCoordinates[0] += _localBasisFctsAndDerivs[indexBasis]
                    * ControlPointNet[CPindex]->getX();
            // Compute iteratively the y-coordinate of the point
            _cartesianCoordinates[1] += _localBasisFctsAndDerivs[indexBasis]
                    * ControlPointNet[CPindex]->getY();
            // Compute iteratively the z-coordinate of the point
            _cartesianCoordinates[2] += _localBasisFctsAndDerivs[indexBasis]
                    * ControlPointNet[CPindex]->getZ();
            // Update basis function's counter
            counter_basis++;
        }
    }
}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates, double _uPrm, int _patchBLIndex, int _patchBLTrCurveIndex) const {

    /*
     * This function gets the curve parameter on a trimming curve and
     * computes the physical coordinates of the point that lies on the designated
     * trimming curve
     */

    const int noCoordParam = 2;

    double uv[noCoordParam];

    // Compute the patch parameters (u, v)
    Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).computeCartesianCoordinates(uv, _uPrm);

    // Compute the physical coordinates (x, y, z)
    computeCartesianCoordinates(_cartesianCoordinates, uv);

}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates, double _uPrm, IGAPatchCurve* _curve) const {

    /*
     * This function gets the curve parameter on a trimming curve and
     * computes the physical coordinates of the point that lies on the designated
     * trimming curve
     */

    const int noCoordParam = 2;

    double uv[noCoordParam];

    // Compute the patch parameters (u, v)
    _curve->computeCartesianCoordinates(uv, _uPrm);

    // Compute the physical coordinates (x, y, z)
    computeCartesianCoordinates(_cartesianCoordinates, uv);

}

void IGAPatchSurface::computeBaseVectors(double* _baseVectors,
        double* _localBasisFunctionsAndDerivatives, int _uKnotSpanIndex, int _vKnotSpanIndex) const {
    /*
     *  Returns the base vectors at a given pair of surface parameters on the NURBS 2D patch. The function expects the computation of the
     *  basis functions and their derivatives to have been carried out outside the function and given as arguments.
     *
     *  The sorting idea is:
     *
     *                 | g1x g1y g1z |
     *  _baseVectors = | g2x g2y g2z |
     *
     *  The output is sorted in an 1D pointer array as follows:
     *
     *  _baseVectors = [g1x g1y g1z g2x g2y g2z]
     *
     *  E.g. _baseVectors[i * 3 + j] returns the j-th coordinate, j=0,…,2 (x,y,z respectively) of the i-th base vector, i=0,…,1 (with respect to
     *  u-parametric line and v-parametric line respectively)
     */

    // Read input
    assert(_baseVectors != NULL);
    assert(_localBasisFunctionsAndDerivatives != NULL);

    // Number of tangent base vectors
    int noBaseVec = 2;

    // Number of coordinates
    int noCoordinates = 3;

    // Initialize the coordinates of the base vectors
    for (int i = 0; i < noBaseVec * noCoordinates; i++)
        _baseVectors[i] = 0;

    // Initialize the index of the derivative functional value
    int indexBasis = 0;

    // Order of derivatives for the basis functions needed for the computation of the base vectors
    int derivDegree = 1;

    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // Initialize the Control Point index
    int CPindex = 0;

    // Initialize a basis functions counter
    int counterBasis = 0;

    // Initialize the derivative counter
    int counterBaseVector = 0;

    // Loop over all the non-zero contributions
    for (int j = 0; j <= qDegree; j++) {
        for (int i = 0; i <= pDegree; i++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            CPindex = (_vKnotSpanIndex - qDegree + j) * uNoControlPoints
                    + (_uKnotSpanIndex - pDegree + i);

            // Reset the base vector's counter
            counterBaseVector = 0;

            // Loop over all the first derivatives of the basis functions and sum up the contributions
            for (int l = 0; l <= derivDegree; l++) {
                for (int k = 0; k <= derivDegree - l; k++) {

                    // Skip the case when k=l=0 i.e. the case where the basis functions themselves will be issued
                    if (k == l)
                        continue;

                    // Get the index of the derivative functional value
                    indexBasis = getIGABasis()->indexDerivativeBasisFunction(derivDegree, k, l,
                            counterBasis);

                    // Compute iteratively the x-coordinate of the point
                    _baseVectors[counterBaseVector * noCoordinates + 0] +=
                            _localBasisFunctionsAndDerivatives[indexBasis]
                                    * ControlPointNet[CPindex]->getX();
                    // Compute iteratively the y-coordinate of the point
                    _baseVectors[counterBaseVector * noCoordinates + 1] +=
                            _localBasisFunctionsAndDerivatives[indexBasis]
                                    * ControlPointNet[CPindex]->getY();
                    // Compute iteratively the z-coordinate of the point
                    _baseVectors[counterBaseVector * noCoordinates + 2] +=
                            _localBasisFunctionsAndDerivatives[indexBasis]
                                    * ControlPointNet[CPindex]->getZ();

                    // Update the base vector counter
                    counterBaseVector++;
                }
            }

            // Update basis function's counter
            counterBasis++;
        }
    }
}

void IGAPatchSurface::computeBaseVectors(double* _baseVectors, double _u, int _spanU, double _v,
        int _spanV) const {

    // Get the polynomial degree of the basis in each direction
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);

    // Derivative order of the basis functions needed for the computation of the base vectors
    int derivDegree = 1;

    // Compute the local basis functions and their derivatives
    double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1) * (derivDegree + 2)
            * noLocalBasisFunctions / 2];
    IGABasis->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
            derivDegree, _u, _spanU, _v, _spanV);

    // Compute the base vectors at the given surface parameters
    computeBaseVectors(_baseVectors, localBasisFunctionsAndDerivatives, _spanU, _spanV);
}

int IGAPatchSurface::indexDerivativeBaseVector(int _derivDegree, int _uDerivIndex, int _vDerivIndex,
        int _componentIndex, int _baseVecIndex) {
    /*
     * Returns the correct index when sorting the base vectors and the their derivatives in an 1D pointer array with the rule:
     * (_derivDegree - _vDerivIndex) * (_derivDegree - _vDerivIndex + 1) * noCoordinates * noBaseVec / 2 +
     *  + _uDerivIndex * noCoordinates * noBaseVec + _componentIndex * noBaseVec + _baseVecIndex
     */

    // Read input
    if (_uDerivIndex + _vDerivIndex > _derivDegree) {
        ERROR_OUT() << "in IGAPatchSurface::indexDerivativeBaseVector" << endl;
        ERROR_OUT() << "It has been requested the " << _uDerivIndex
                << "-th partial derivative w.r.t. u and the " << _vDerivIndex
                << "-th partial derivative w.r.t. v of the base vectors but " << endl;
        ERROR_OUT() << "the maximum absolute derivative selected is of " << _derivDegree
                << "-th order" << endl;
        exit(-1);
    }

    // Number of tangent base vectors
    int noBaseVec = 2;

    // Number of coordinates
    int noCoordinates = 3;

    // Compute the index of the functional value
    return (_derivDegree - _vDerivIndex) * (_derivDegree - _vDerivIndex + 1) * noCoordinates
            * noBaseVec / 2 + _uDerivIndex * noCoordinates * noBaseVec + _componentIndex * noBaseVec
            + _baseVecIndex;
}

void IGAPatchSurface::computeBaseVectorsAndDerivatives(double* _baseVectorsAndDerivatives,
        double*_localBasisFunctionsAndDerivatives, int _derivDegree, int _uKnotSpanIndex,
        int _vKnotSpanIndex) {
    /*
     * Returns the base vectors and their j-th derivative w.r.t. v-parametric line and i-th derivative w.r.t. u-parametric line partial derivatives.
     * On the input arrays the following rules must be fulfiled:
     *
     * · The array _baseVectorsAndDerivatives = new double[(_derivDegree + 1) * (_derivDegree + 2) * noCoord * noBaseVec / 2] will be filled up
     *   with the components of the base vectors and their derivatives up to absolute order _derivDegree meaning that 0 <= i + j <= _derivDegree
     *
     * · The array _localBasisFunctionsAndDerivatives = new double[(_derivDegree + 2) * (_derivDegree + 3) * noBasisFcts / 2] contains the IGA basis
     *   functions and their derivatives up to absolute order _derivDegree + 1 since the base vectors contain by definition the first derivatives of
     *   IGA basis functions
     */

    // Read input
    assert(_baseVectorsAndDerivatives != NULL);
    assert(_localBasisFunctionsAndDerivatives != NULL);

    // The polynomial degrees of the patch
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // Initialize the index related to the basis functions
    int indexBasis = 0;

    // Initialize the index related to the CP's
    int indexCP = 0;

    // Initialize the index related to the derivatives of the base vectors
    int indexBaseVct = 0;

    // The derivative order of the basis functions
    int derivDegreeBasis = _derivDegree + 1;

    // Initialize factor by which to multiply the derivative of the basis function
    double factor = 0.0;

    // The number of coordinates for the base vectors
    int noCoordinates = 3;

    // The number of the base vectors to be returned
    int noBaseVct = 2;

    // Initialize the output array
    int counterBaseVec = 0;
    for (int i = 0; i <= _derivDegree; i++) {
        for (int j = 0; j <= _derivDegree - i; j++) {
            for (int k = 0; k < noCoordinates; k++) {
                for (int l = 0; l < noBaseVct; l++) {
                    _baseVectorsAndDerivatives[counterBaseVec] = 0.0;
                    counterBaseVec++;
                }
            }
        }
    }

    // Initialize the counter of the basis functions
    int counterBasis = 0;

    // Loop over all the non-zero contributions in v-direction
    for (int vBasis = 0; vBasis <= qDegree; vBasis++) {
        // Loop over all the non-zero contributions in u-direction
        for (int uBasis = 0; uBasis <= pDegree; uBasis++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            indexCP = (_vKnotSpanIndex - qDegree + vBasis) * uNoControlPoints
                    + (_uKnotSpanIndex - pDegree + uBasis);

            // Loop over all the derivatives in u-direction
            for (int i = 0; i <= _derivDegree; i++) {
                // Loop over all the derivatives in v-direction
                for (int j = 0; j <= _derivDegree - i; j++) {
                    // Loop over all the coordinates
                    for (int k = 0; k < noCoordinates; k++) {
                        // Loop over all the base vectors
                        for (int l = 0; l < noBaseVct; l++) {
                            // Compute the index related to the derivatives of the base vectors
                            indexBaseVct = indexDerivativeBaseVector(_derivDegree, i, j, k, l);

                            // · Case l = 0 :: --> i + 1 - l = i + 1 && j + l --> j giving the derivative of the base vector g1 = dX/du
                            // · Case l = 1 :: --> i + 1 - l = i && j + l --> j + 1 giving the derivative of the base vector g2 = dX/dv

                            // Compute the index of the basis functions and their derivatives
                            indexBasis = IGABasis->indexDerivativeBasisFunction(derivDegreeBasis,
                                    i + 1 - l, j + l, counterBasis);

                            // Factor by which to multiply the derivative of the basis function
                            if (k == 0)
                                factor = ControlPointNet[indexCP]->getX();
                            else if (k == 1)
                                factor = ControlPointNet[indexCP]->getY();
                            else
                                factor = ControlPointNet[indexCP]->getZ();

                            // Add the contribution from each basis function in the interval
                            _baseVectorsAndDerivatives[indexBaseVct] +=
                                    _localBasisFunctionsAndDerivatives[indexBasis] * factor;
                        }
                    }
                }
            }
            // Update basis function's counter
            counterBasis++;
        }
    }
}

void IGAPatchSurface::computeSurfaceNormalVectorAndDerivatives(double* _surfNormalVctAndDervs, double* _baseVctsAndDerivs, int _derivDegreeBaseVec){
    /*
     * Returns the surface normal base vector and its first order parametric derivatives given the surface base vectors and
     * their first order parametric derivatives. The output array is filled up the components of the surface normal base
     * vector and its derivatives as follows:
     *
     * · The array _surfNormalVctAndDervs = double[9] with
     *   _surfNormalVctAndDervs = [surfNormalVct[0] surfNormalVct[1] surfNormalVct[2] …
     *                             dSurfNormalVctdU[0] dSurfNormalVctdU[1] dSurfNormalVctdU[2] …
     *                             dSurfNormalVctdV[0] dSurfNormalVctdV[1] dSurfNormalVctdV[2]]
     */

    // Initialize auxiliary variables
    int indexBaseVct;

    // Number of Cartesian coodinates
    int noCoord = 3;

    // Get the base vector along u
    double baseVctU[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 0, 0, i, 0);
        baseVctU[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Get the base vector along v
    double baseVctV[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 0, 0, i, 1);
        baseVctV[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Get the derivative of base vector along u with respect to u parametric direction
    double dBaseVctUdU[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 1, 0, i, 0);
        dBaseVctUdU[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Get the derivative of base vector along v with respect to u parametric direction
    double dBaseVctVdU[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 1, 0, i, 1);
        dBaseVctVdU[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Get the derivative of base vector along u with respect to v parametric direction
    double dBaseVctUdV[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 0, 1, i, 0);
        dBaseVctUdV[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Get the derivative of base vector along v with respect to v parametric direction
    double dBaseVctVdV[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 0, 1, i, 1);
        dBaseVctVdV[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Compute the not normalized surface normal vector
    double surfNormalVctTilde[noCoord];
    MathLibrary::crossProduct(surfNormalVctTilde, baseVctU, baseVctV);

    // Compute the differential area
    double dA = MathLibrary::vector2norm(surfNormalVctTilde,noCoord);

    // Compute the surface normal vector
    for(int i = 0; i < noCoord; i++){
        _surfNormalVctAndDervs[i] = surfNormalVctTilde[i]/dA;
    }

    // Compute the derivative of the surface normal vector with respect to the u parametric coordinate
    double dBaseVctUdUTimesBaseVctV[noCoord];
    double baseVctUTimesdBaseVctVdU[noCoord];
    MathLibrary::crossProduct(dBaseVctUdUTimesBaseVctV, dBaseVctUdU, baseVctV);
    MathLibrary::crossProduct(baseVctUTimesdBaseVctVdU, baseVctU, dBaseVctVdU);
    // A1TimesA2Komma1 = cross(dGCovariant(:,1),GCovariant(:,2)) + cross(GCovariant(:,1),dGCovariant(:,3));
    MathLibrary::computeDenseVectorAddition(dBaseVctUdUTimesBaseVctV, baseVctUTimesdBaseVctVdU,1.0,noCoord);
    // A3TildeDotA1TimesA2Komma1 = A3Tilde'*A1TimesA2Komma1
    double A3TildeDotA1TimesA2Komma1 = MathLibrary::computeDenseDotProduct(noCoord,surfNormalVctTilde,dBaseVctUdUTimesBaseVctV);
    // dA3Komma1 = (1/dA)*A1TimesA2Komma1 - (1/dA^3)*(A3Tilde'*A1TimesA2Komma1)*A3Tilde;
    for(int i = 0; i < noCoord; i++){
        _surfNormalVctAndDervs[1*noCoord + i] = (1/dA)*dBaseVctUdUTimesBaseVctV[i] - (1/(dA*dA*dA))*A3TildeDotA1TimesA2Komma1*surfNormalVctTilde[i];
    }

    // Compute the derivative of the surface normal vector with respect to the u parametric coordinate
    double dBaseVctUdVTimesBaseVctV[noCoord];
    double baseVctUTimesdBaseVctVdV[noCoord];
    MathLibrary::crossProduct(dBaseVctUdVTimesBaseVctV, dBaseVctUdV, baseVctV);
    MathLibrary::crossProduct(baseVctUTimesdBaseVctVdV, baseVctU, dBaseVctVdV);
    // A1TimesA2Komma2 = cross(dGCovariant(:,3),GCovariant(:,2)) + cross(GCovariant(:,1),dGCovariant(:,2));
    MathLibrary::computeDenseVectorAddition(dBaseVctUdVTimesBaseVctV, baseVctUTimesdBaseVctVdV,1.0,noCoord);
    double A3TildeDotA1TimesA2Komma2 = MathLibrary::computeDenseDotProduct(noCoord,surfNormalVctTilde,dBaseVctUdVTimesBaseVctV);
    // dA3Komma2 = (1/dA)*A1TimesA2Komma2 - (1/dA^3)*(A3Tilde'*A1TimesA2Komma2)*A3Tilde;
    for(int i = 0; i < noCoord; i++){
        _surfNormalVctAndDervs[2*noCoord + i] = (1/dA)*dBaseVctUdVTimesBaseVctV[i] - (1/(dA*dA*dA))*A3TildeDotA1TimesA2Komma2*surfNormalVctTilde[i];
    }
}

void IGAPatchSurface::computeCovariantMetricTensor(double* _covariantMetricTensor, double* _baseVctsAndDerivs, int _derivDegreeBaseVec){
    /*
     * Returns the coefficients of the covariant metric tensor.
     *
     * · The array _covariantMetricTensor = double[4] is sorted as
     *   _surfNormalVctAndDervs = [baseVctU*baseVctU baseVctU*baseVctV baseVctV*baseVctU baseVctV*baseVctV]
     */

    // Initialize auxiliary variables
    int indexBaseVct;

    // Number of Cartesian coodinates
    int noCoord = 3;

    // Get the base vector along u
    double baseVctU[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 0, 0, i, 0);
        baseVctU[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Get the base vector along v
    double baseVctV[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 0, 0, i, 1);
        baseVctV[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Compute the entries of the covariant metric tensor
    _covariantMetricTensor[0] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,baseVctU,baseVctU);
    _covariantMetricTensor[1] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,baseVctU,baseVctV);
    _covariantMetricTensor[2] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,baseVctV,baseVctU);
    _covariantMetricTensor[3] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,baseVctV,baseVctV);
}

void IGAPatchSurface::computeContravariantBaseVectors(double* _contravariantBaseVcts, double* _covariantMetricTensor, double* _baseVctsAndDerivs, int _derivDegreeBaseVec) {
    /*
     * Returns the Cartesian components of the contravariant base vectors in an array of constant size 6
     *
     * · The array _contravariantBaseVcts = double[4] is sorted as
     *   _contravariantBaseVcts = [baseVctContraU[0] baseVctContraU[1] baseVctContraU[2] baseVctContraV[0] baseVctContraV[1] baseVctContraV[2]]
     */

    // Initialize flag for the solvability of the linear equation system
    bool isLinearSystemSolvable;

    // Initialize auxiliary variables
    int indexBaseVct;

    // Number of Cartesian coodinates
    int noCoord = 3;

    // Number of parametric coordinates
    int noParameticCoord = 2;

    // Initialize the components of the base vectors for each Cartesian direction seperately
    double baseVctComponents[noParameticCoord];

    // Get the base vector along u
    double baseVctU[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 0, 0, i, 0);
        baseVctU[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Get the base vector along v
    double baseVctV[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 0, 0, i, 1);
        baseVctV[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Loop over all Cartesian directions
    for(int i = 0; i < noCoord; i++){
        // Get the right-hand side of the linear equation system
        baseVctComponents[0] = baseVctU[i];
        baseVctComponents[1] = baseVctV[i];

        // Solve for the Cartesian coordinates of the contravariant base vectors
        isLinearSystemSolvable = EMPIRE::MathLibrary::solve2x2LinearSystem(_covariantMetricTensor, baseVctComponents, EMPIRE::MathLibrary::EPS);

        // Assign the values to the return array
        _contravariantBaseVcts[i] = baseVctComponents[0];
        _contravariantBaseVcts[noCoord + i] = baseVctComponents[1];
    }
}

void IGAPatchSurface::computeContravariantCurvatureTensor(double* _contravariantCurvatureTensor, double* _surfNormalVctAndDervs, double* _baseVctsAndDerivs, int _derivDegreeBaseVec){
    /*
     * Returns the components of the curvature tensor in the contravariant basis in an array of constant size 4
     *
     * · The array _contravariantCurvatureTensor = double[6] is sorted as
     *   _contravariantCurvatureTensor = [dbaseVctContraUdU dbaseVctContraUdV dbaseVctContraVdU dbaseVctContraVdV]*_surfNormalVctAndDervs
     */

    // Initialize auxiliary variables
    int indexBaseVct;

    // Number of Cartesian coodinates
    int noCoord = 3;

    // Get the derivative of base vector along u with respect to u parametric direction
    double dBaseVctUdU[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 1, 0, i, 0);
        dBaseVctUdU[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Get the derivative of base vector along u with respect to v parametric direction
    double dBaseVctUdV[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 0, 1, i, 0);
        dBaseVctUdV[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Get the derivative of base vector along v with respect to u parametric direction
    double dBaseVctVdU[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 1, 0, i, 1);
        dBaseVctVdU[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Get the derivative of base vector along v with respect to v parametric direction
    double dBaseVctVdV[noCoord];
    for(int i = 0; i < noCoord; i++){
        indexBaseVct = this->indexDerivativeBaseVector(_derivDegreeBaseVec, 0, 1, i, 1);
        dBaseVctVdV[i] = _baseVctsAndDerivs[indexBaseVct];
    }

    // Compute the components of the curvature tensor in the contravariant basis
    _contravariantCurvatureTensor[0] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,dBaseVctUdU,_surfNormalVctAndDervs);
    _contravariantCurvatureTensor[1] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,dBaseVctUdV,_surfNormalVctAndDervs);
    _contravariantCurvatureTensor[2] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,dBaseVctVdU,_surfNormalVctAndDervs);
    _contravariantCurvatureTensor[3] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,dBaseVctVdV,_surfNormalVctAndDervs);
}

bool IGAPatchSurface::computePointProjectionOnPatch(double& _u, double& _v, double* _P,
        bool& _flagConverge, int _maxIt, double _orthoTol, double _distTol) {

    /*
     * Returns the projection of a point _P on the NURBS patch given an initial guess for the surface parameters _u, _v via references:
     * _P = double[3]
     * Return value is a bool flag on the convergence of the Newton-Raphson iterations.
     *
     * Function layout :
     *
     * 1. Read input and initialize the data
     *
     * 2. Loop over all the Newton-Raphson iterations
     *    2i. Update the iteration counter
     *   2ii. Find the span of the given surface parameters
     *  2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
     *   2iv. Compute the Cartesian components of the point on the surface
     *    2v. Compute the distance vector between the vector to be projected and the estimated one
     *   2vi. Compute the 2-norm of the distance vector
     *  2vii. Compute the base vectors and their derivatives
     * 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
     *   2ix. Compute the cosine of the angle with respect to u-parametric line
     *    2x. Compute the cosine of the angle with respect to v-parametric line
     *   2xi. Check the orthogonality condition and if it is fulfilled break the loop
     *  2xii. Compute the entries of the Jacobian matrix
     * 2xiii. Compute the entries of the right-hand side vector
     *  2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved
     *   2xv. Update the surface parameters u += du and v += dv
     *  2xvi. Check and modify the surface parameters if they stay out of their knot spans
     *
     * 3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
     *
     * 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
     */

    // 1. Read input and initialize the data
    // Read input
    assert(_P != NULL);

    // Initialize the flag to true
    bool flagNewtonRaphson = true;
    _flagConverge = true;
	
	// Flag for solvability of linear equation system
	bool isSystemSolved = true;

    // Flag for solving linear equation system
    bool isSystemSolvable = true;
    bool conditionFirstColumnZero = false;
    bool conditionSecondColumnZero = false;
    bool conditionFirstRowZero = false;
    bool conditionSecondRowZero = false;

    const double epsJ = 1e-6;
    const double epsDuv = 1e-10;
    //bool fixU = false;
    //bool fixV = false;

    // Initialize number of spatial dimensions
    int noSpatialDimensions = 3;

    // Save the location of the point to be projected
    double point[noSpatialDimensions];
    for (int i = 0; i < noSpatialDimensions; i++)
        point[i] = _P[i];

    // Initialize the distance vector
    double distanceVector[3];

    // Initialize the indices of the base vectors and their derivatives
    int indexGu = 0;
    int indexGv = 0;
    int indexDGuDu = 0;
    int indexDGvDv = 0;
    int indexDGuDv = 0;
    //int indexDGvDu = indexDGuDv;

    // Initialize the base vectors and their derivatives
    double Gu[3];
    double Gv[3];
    double DGuDu[3];
    double DGvDv[3];
    double DGuDv[3];
    //double* DGvDu=DGuDv;

    // Initialize the dot products
    double distanceVector2norm = 0.0;
    double GuXdistanceVector = 0.0;
    double squareGu2norm = 0.0;
    double Gu2norm = 0.0;
    double GvXdistanceVector = 0.0;
    double squareGv2norm = 0.0;
    double Gv2norm = 0.0;

    // Initialize Jacobian matrix
    double dR[4];

    // Initialize right-hand side solution vector
    double R[2];

    // Initialize flag on the solvability of the 2x2 equation system
    bool flagLinearSystem = 1;

    // Initialize the cosines w.r.t. each parametric line
    double cosu = 0.0;
    double cosv = 0.0;

    // Initialize the knot span indices
    int uKnotSpan = 0;
    int vKnotSpan = 0;

    // The NURBS polynomial degrees
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // The lengths of the knot vectors to the NURBS patch
    //int lengthUKnotVct = IGABasis->getUBSplineBasis1D()->getNoKnots();
    //int lengthVKnotVct = IGABasis->getVBSplineBasis1D()->getNoKnots();

    // Local number of basis functions
    int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

    // Initialize the Newton-Raphson iteration counter
    int counter = 0;

    // Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
    int derivDegreeBasis = 2;

    // The number of the base vectors
    int noBaseVcts = 2;

    // Initialize the array of the IGA basis functions and their derivatives
    double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1) * (derivDegreeBasis + 2)
            * noLocalBasisFcts / 2];

    // Number of derivatives for the base vectors
    int derivDegreeBaseVcts = derivDegreeBasis - 1;

    // Initialize the array of the base vectors and their derivatives
    double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1) * (derivDegreeBaseVcts + 2)
            * noSpatialDimensions * noBaseVcts / 2];

    // 2. Loop over all the Newton-Raphson iterations
    while (counter <= _maxIt) {
        // 2i. Update the iteration counter
        counter++;

        // 2ii. Find the span of the given surface parameters
        uKnotSpan = IGABasis->getUBSplineBasis1D()->findKnotSpan(_u);
        vKnotSpan = IGABasis->getVBSplineBasis1D()->findKnotSpan(_v);

        // 2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
        IGABasis->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs, derivDegreeBasis, _u,
                uKnotSpan, _v, vKnotSpan);

        computeCartesianCoordinates(_P, basisFctsAndDerivs, derivDegreeBasis, uKnotSpan, vKnotSpan);

        // 2v. Compute the distance vector between the vector to be projected and the estimated one
        for (int i = 0; i < noSpatialDimensions; i++)
            distanceVector[i] = _P[i] - point[i];

        // 2vi. Compute the 2-norm of the distance vector
        distanceVector2norm = EMPIRE::MathLibrary::vector2norm(distanceVector, noSpatialDimensions);

        if (distanceVector2norm < _distTol)
            break;

        // 2vii. Compute the base vectors and their derivatives
        computeBaseVectorsAndDerivatives(baseVecAndDerivs, basisFctsAndDerivs, derivDegreeBaseVcts,
                uKnotSpan, vKnotSpan);

        // 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
        for (int i = 0; i < noSpatialDimensions; i++) {
            // On the base vector Gu = dR/du
            indexGu = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 0);
            Gu[i] = baseVecAndDerivs[indexGu];

            // On the base vector Gv = dR/dv
            indexGv = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 1);
            Gv[i] = baseVecAndDerivs[indexGv];

            // On the derivative of the base vector Gu w.r.t. u namely dGu/du = d^2R/du^2
            indexDGuDu = indexDerivativeBaseVector(derivDegreeBaseVcts, 1, 0, i, 0);
            DGuDu[i] = baseVecAndDerivs[indexDGuDu];

            // On the derivative of the base vector Gv w.r.t. u namely dGv/dv = d^2R/dv^2
            indexDGvDv = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 1);
            DGvDv[i] = baseVecAndDerivs[indexDGvDv];

            // On the mixed derivative of the base vectors namely d^2Gu/dv = d^2Gv/du = d^2R/dudv
            indexDGuDv = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 0);
            DGuDv[i] = baseVecAndDerivs[indexDGuDv];
        }

        // 2ix. Compute the cosine of the angle with respect to u-parametric line
        GuXdistanceVector = MathLibrary::computeDenseDotProduct(noSpatialDimensions, Gu, distanceVector);
        Gu2norm = MathLibrary::vector2norm(Gu, noSpatialDimensions);
        squareGu2norm = Gu2norm * Gu2norm;
        cosu = fabs(GuXdistanceVector) / Gu2norm / distanceVector2norm;

        // 2x. Compute the cosine of the angle with respect to v-parametric line
        GvXdistanceVector = MathLibrary::computeDenseDotProduct(noSpatialDimensions, Gv, distanceVector);
        Gv2norm = MathLibrary::vector2norm(Gv,noSpatialDimensions);
        squareGv2norm = Gv2norm*Gv2norm;
        cosv = fabs(GvXdistanceVector) / Gv2norm / distanceVector2norm;

        // 2xi. Check the orthogonality condition and if it is fulfilled break the loop
        if (cosu <= _orthoTol && cosv <= _orthoTol)
            break;

        // 2xii. Compute the entries of the Jacobian matrix
        dR[0] = squareGu2norm
                + EMPIRE::MathLibrary::computeDenseDotProduct(noSpatialDimensions, DGuDu, distanceVector);
        dR[1] = EMPIRE::MathLibrary::computeDenseDotProduct(noSpatialDimensions, Gu, Gv)
                + EMPIRE::MathLibrary::computeDenseDotProduct(noSpatialDimensions, DGuDv, distanceVector);
        dR[2] = dR[1];
        dR[3] = squareGv2norm
                + EMPIRE::MathLibrary::computeDenseDotProduct(noSpatialDimensions, DGvDv, distanceVector);

        // 2xiii. Compute the entries of the right-hand side vector
        R[0] = -EMPIRE::MathLibrary::computeDenseDotProduct(noSpatialDimensions, Gu, distanceVector);
        R[1] = -EMPIRE::MathLibrary::computeDenseDotProduct(noSpatialDimensions, Gv, distanceVector);


        conditionFirstRowZero = false;
        if ((fabs(dR[0]) < epsJ && fabs(dR[1]) < epsJ) ){
            conditionFirstRowZero = true;
        }

        conditionSecondRowZero = false;
        if (fabs(dR[2]) < epsJ && fabs(dR[3]) < epsJ){
            conditionSecondRowZero = true;
        }

        conditionFirstColumnZero = false;
        if ((fabs(dR[0]) < epsJ && fabs(dR[2]) < epsJ) ){
            conditionFirstColumnZero = true;
        }

        conditionSecondColumnZero = false;
        if ((fabs(dR[0]) < epsJ && fabs(dR[2]) < epsJ) ){
            conditionSecondColumnZero = true;
        }

        // Check if the system is solvable by checking the condition of the diagonal entries
        if (conditionFirstRowZero || conditionSecondRowZero || conditionFirstColumnZero || conditionSecondColumnZero){
            isSystemSolvable = false;
        }

        if (isSystemSolvable){
            // 2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved

            // Solve the equation system
            flagLinearSystem = EMPIRE::MathLibrary::solve2x2LinearSystem(dR, R, EMPIRE::MathLibrary::EPS );

            // Check if the equation system has been successfully solved
            if (!flagLinearSystem) {
                WARNING_OUT() << "Error in IGAPatchSurface::computePointProjectionOnPatch" << endl;
                WARNING_OUT()
                        << "The 2x2 equation system to find the updates of the surface parameters"
                        << endl;
                WARNING_OUT()
                        << "for the orthogonal projection of a point on the NURBS patch has been"
                        << endl;
                WARNING_OUT() << "detected not solvable up to tolerance" << EMPIRE::MathLibrary::EPS << endl;
                // exit(-1);
				isSystemSolved = false;
				break;
            }
        } else {
            if (conditionFirstRowZero){
                R[1] = R[1] / dR[3];
                R[0] = 0.0;
            } else if(conditionSecondRowZero){
                R[0] = R[0] / dR[0];
                R[1] = 0.0;
            } else if(conditionFirstColumnZero){
                R[1] = (R[0] + R[1])/(dR[1]+dR[3]);
                R[0] = 0.0;
            } else if(conditionSecondColumnZero){
                R[0] = (R[0] + R[1])/(dR[0]+dR[2]);
                R[1] = 0.0;
            }
        }

//        if (fabs(dR[0]) < epsJ || fixU) {
//            std::cout << "fabs(dR[0]) < epsJ || fixU dR: " << dR[0] <<" " << dR[1] <<" " << dR[2] << " " << dR[3] <<std::endl;
//            R[0] = 0.0;
//            R[1] = R[1] / dR[3];
//            fixU = false;
//            fixV = true;
//        } else if (fabs(dR[3]) < epsJ || fixV) {
//            std::cout << "fabs(dR[3]) < epsJ || fixV dR: " << dR[0] <<" " << dR[1] <<" " << dR[2] << " " << dR[3] <<std::endl;
//            // According to Fabien that must be R[0] / dR[0];
//            // R[0] = R[0] / dR[1];  // According to Chenshen
//            R[0] = R[0] / dR[0]; // According to Fabien
//            R[1] = 0.0;
//            fixU = true;
//            fixV = false;
//        } else {
//            std::cout << "Solving the linear system !"<<std::endl;
//            // 2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved

//            // Solve the equation system
//            flagLinearSystem = EMPIRE::MathLibrary::solve2x2LinearSystem(dR, R, EMPIRE::MathLibrary::EPS );

//            // Check if the equation system has been successfully solved
//            if (!flagLinearSystem) {
//                ERROR_OUT() << "Error in IGAPatchSurface::computePointProjectionOnPatch" << endl;
//                ERROR_OUT()
//                        << "The 2x2 equation system to find the updates of the surface parameters"
//                        << endl;
//                ERROR_OUT()
//                        << "for the orthogonal projection of a point on the NURBS patch has been"
//                        << endl;
//                ERROR_OUT() << "detected not solvable up to tolerance" << EMPIRE::MathLibrary::EPS << endl;
//                exit(-1);
//            }
//        }

        // 2xv. Update the surface parameters u += du and v += dv
        _u += R[0];
        _v += R[1];

        // 2xvi. Check and modify the surface parameters if they stay out of their knot spans
    	IGABasis->getUBSplineBasis1D()->clampKnot(_u);
    	IGABasis->getVBSplineBasis1D()->clampKnot(_v);
    }
////     3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
    if (counter > _maxIt) {
        flagNewtonRaphson = false;
        if (R[0] * R[0] + R[1] * R[1] < epsDuv)
            _flagConverge = true;
        else
            _flagConverge = false;
    } else {
        flagNewtonRaphson = true;
    }
    if (!isSystemSolved) {
		flagNewtonRaphson = false;
	}
    
    // 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
    // Clear the memory on the heap
    delete[] basisFctsAndDerivs;
    delete[] baseVecAndDerivs;

    // Return the flag
    return flagNewtonRaphson;
}

bool IGAPatchSurface::computePointProjectionOnPatch(double& _u, double& _v, double* _P, int _maxIt, double _orthoTol, double _distTol) {
    /*
     *  Returns the projection of a point _P on the NURBS patch given an initial guess for the surface parameters _u, _v via references:
     * _P = double[3]. Its return value is a bool flag on the convergence of the Newton-Raphson iterations. Makes use of the previosuly
     * implemented function IGAPatchSurface::computePointProjectionOnPatch
     */

    // Initialize the boolean on the convergence of the Newton iterations
    bool tmp = false;

    // Compute the closest point projection using the Newton-Rapshon algorithm
    return computePointProjectionOnPatch(_u, _v, _P, tmp, _maxIt, _orthoTol, _distTol);
}

bool IGAPatchSurface::computeForcedPointProjectionOnPatch(double& _u, double& _v, double* _P, int _relMaxIt, double _relOrthoTol, double _relDistTol) {

    /*
     * Forces the point projection on patch
     */

    // Initialize the projection flag
    bool isProjected = false;

    // Reset the Initialize new variables for so that the initial guess remains unchanged
    double initU;
    double initV;


    // Force projection Project and return the patch parameters of the projection
    int tryCtr = 0;
    while (!isProjected) {

        initU = _u;
        initV = _v;

        // Relax the orthogonality tolerance for projection at each try by 10.0
        isProjected = computePointProjectionOnPatch(initU, initV, _P, _relMaxIt, TOL_ORTHOGONALITY*pow(10.0,tryCtr*2), TOL_DISTANCE*pow(10.0,tryCtr));

        // If projection is successful then update the input variables
        if (isProjected) {
            _u = initU;
            _v = initV;
        }
        tryCtr++;

        // If maximum relaxed tolerance is reached then raise error
        if (TOL_ORTHOGONALITY*pow(10.0,tryCtr*2) > _relOrthoTol || TOL_DISTANCE*pow(10.0,tryCtr) > _relDistTol)
            ERROR_OUT("IGAPatchSurface::computeForcedPointProjectionOnPatch: Relaxed tolerance exceeded! Point not projected!");
    }

    return isProjected;
}

bool IGAPatchSurface::computePointProjectionOnTrimmingCurve(double& _projectedUTilde,
                                                            double* _P, int _patchBLIndex, int _patchBLTrCurveIndex){

    /*
     * This function receives a cartesian point for projection onto the trimming curve on this patch
     * The point is firstly projected onto the patch. After the patch parameters of the point are found, it is projected onto the curve
     * in the patch parameter space.
     * The initial guess for the projection onto the patch is performed making use of the linearization of the trimming curve.
     * The minimum cartesian distance between the point and the linearization vertices is taken as an initial guess.
     */

    // Initialize coordinates
    const int noCoordParam = 2;

    // Initialize variables and flags
    bool isProjectedOnPatch = false;
    bool isProjectedOnCurve = false;
    double UV[noCoordParam];

    // Find the initial guess for projection onto the patch for the first point
    findInitialGuess4PointProjectionOnTrimmingCurve(_projectedUTilde, UV[0], UV[1], _P, _patchBLIndex, _patchBLTrCurveIndex);

    // Force projection and return the patch parameters of the projection
    isProjectedOnPatch = computeForcedPointProjectionOnPatch(UV[0], UV[1], _P);

    if (isProjectedOnPatch) {
        isProjectedOnCurve = Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).computePointProjectionOn2DCurve(_projectedUTilde, UV, noCoordParam);
        if (!isProjectedOnCurve)
            WARNING_OUT("In \"computePointProjectionOnTrimmingCurve\": Point was not projected on curve and skipped!");
    } else {
        WARNING_OUT("In \"computePointProjectionOnTrimmingCurve\": Point was not projected on patch and skipped!");
        isProjectedOnCurve = false;
    }

    return isProjectedOnCurve;
}

bool IGAPatchSurface::computePointProjectionOnTrimmingCurve(double& _projectedUTilde,
                                                            double* _P, IGAPatchCurve* _curve) {

    /*
     * This function receives a cartesian point for projection onto the trimming curve on this patch
     * The point is firstly projected onto the patch. After the patch parameters of the point are found, it is projected onto the curve
     * in the patch parameter space.
     * The initial guess for the projection onto the patch is performed making use of the linearization of the trimming curve.
     * The minimum cartesian distance between the point and the linearization vertices is taken as an initial guess.
     */

    // Initialize coordinates
    const int noCoordParam = 2;

    // Initialize variables and flags
    bool isProjectedOnPatch = false;
    bool isProjectedOnCurve = false;
    double UV[noCoordParam];

    // Find the initial guess for projection onto the patch for the first point
    findInitialGuess4PointProjectionOnTrimmingCurve(_projectedUTilde, UV[0], UV[1], _P, _curve);

    // Force projection and return the patch parameters of the projection
    isProjectedOnPatch = computeForcedPointProjectionOnPatch(UV[0], UV[1], _P);

    if (isProjectedOnPatch) {
        isProjectedOnCurve = _curve->computePointProjectionOn2DCurve(_projectedUTilde, UV, noCoordParam);
        if (!isProjectedOnCurve)
            WARNING_OUT("In \"computePointProjectionOnTrimmingCurve\": Point was not projected on curve and skipped!");
    } else {
        WARNING_OUT("In \"computePointProjectionOnTrimmingCurve\": Point was not projected on patch and skipped!");
        isProjectedOnCurve = false;
    }

    return isProjectedOnCurve;
}

bool IGAPatchSurface::solvePointProjectionOnPatchBoundaryBisection(
        double& _u, double& _v, double& _lambda, double& _distance, double* _P1,
        double* _P2, int _maxIt, double _tol) {
    // Check input
       assert(_P1 != NULL);
       assert(_P2 != NULL);

       // Initialize variables
       int numCoord = 3;
       bool isIn;
       bool isConverged = false;
       double UV1[2] = { _u, _v };
       double UV[2] = { _u, _v };
       double P[numCoord];
       double P1[numCoord];
       double P2[numCoord];
       double Q[numCoord],Q1[numCoord];
       double P1P[numCoord];
       double P1P2[numCoord];
       double QP[numCoord];
       int counter = 0;

       // Initialize the end points of the segment
       for (int i = 0; i < numCoord; i++) {
           P1[i] = _P1[i];
           P2[i] = _P2[i];
       }

       // Loop over all fixed point iterations
       while (counter < _maxIt) {
           // Update the iteration counter
           counter++;

           // Compute P out of P1 and P2, the vector P1P and update the Q as P
           for (int i = 0; i < numCoord; i++) {
               P[i] = 0.5 * (P1[i] + P2[i]);
               P1P[i] = P[i] - P1[i];
               Q[i] = P[i];
           }

           // Check convergence based on whether the segment P1P is small enough
           if (MathLibrary::vector2norm(P1P,numCoord) <= _tol) {
               isConverged = true;
               break;
           }

           // Try to project Q on the patch
           isIn = computePointProjectionOnPatch(UV[0], UV[1], Q);

           // Update the P1 and P2 vertices accordingly
           if (isIn) {
               for (int i = 0; i < numCoord; i++) {
                   P1[i] = P[i];
                   Q1[i] = Q[i];
               }
               UV1[0] = UV[0];
               UV1[1] = UV[1];
           } else {
               for (int i = 0; i < numCoord; i++)
                   P2[i] = P[i];
               UV[0] = UV1[0];
               UV[1] = UV1[1];
           }
       }

       // Compute the output variables if the algorithm converged
       if (isConverged) {
           // Compute the distance QP and the ratio on the line being inside the patch
           for (int i = 0; i < numCoord; i++) {
               P1P[i] = P1[i] - _P1[i];
               P1P2[i] = _P2[i] - _P1[i];
               QP[i] = P1[i] - Q1[i];
           }
           _u = UV1[0];
           _v = UV1[1];
           _distance = MathLibrary::vector2norm(QP,numCoord);
           _lambda =  MathLibrary::computeDenseDotProduct(numCoord, P1P2, P1P)
                   / MathLibrary::computeDenseDotProduct(numCoord, P1P2, P1P2);
       }

       // Return the convergence flag
       return isConverged;
}

char IGAPatchSurface::computePointProjectionOnPatchBoundaryBisection(double& _u, double& _v, double& _lambda,
        double& _distance, double* _P1, double* _P2, int _maxIt, double _tol){
    DEBUG_OUT()<<"\t======================================================"<<endl;
    DEBUG_OUT() << "in IGAPatchSurface::computePointProjectionOnPatchBoundaryBisection"<<endl;
    DEBUG_OUT()<<"\tPROJECT line on boundary using Bisection for line"<<endl
            <<"\t\t(("<<_P1[0]<<" , "<<_P1[1]<<" , "<<_P1[2]<<");"
            <<"("<<_P2[0]<<" , "<<_P2[1]<<" , "<<_P2[2]<<")) "<<endl
            <<"\t\twith initial guess  projection of P1 : (u,v)=("<<_u<<" , "<<_v<<")"<<endl;

    double distance = _distance;
    double lambda = _lambda;
    bool isConverged;
    double u = _u;
    double v = _v;

    // Compute point projection from the line to the NURBS patch boundary
    isConverged = solvePointProjectionOnPatchBoundaryBisection(u, v, lambda, distance, _P1, _P2, _maxIt, _tol);

    if(isConverged) {
        char edge = getEdge(u, v);
        DEBUG_OUT()<<"\t-------------------------------------------------------"<<endl;
        DEBUG_OUT()<<"\tAlgorithm has CONVERGED on edge "<<int(edge)<<" and distance to patch is "<<distance<<endl
                <<"\t\t and ratio P1P/P1P2 is " << _lambda << endl
                <<"\t\t and parametric value are (u,v)("<<u<<" , "<<v<<")"<<endl;
        DEBUG_OUT()<<"\t======================================================"<<endl;

        // Fix the lambda to 0 or 1 if outside the [0,1] domain
        if(lambda - 1.0 > 0) lambda = 1.0;
        if(lambda < 0)       lambda = 0.0;

        // Fix the parameter running on the boundary
        if (edge == EDGE_U0)
            u = getIGABasis()->getUBSplineBasis1D()->getFirstKnot();
        else if (edge == EDGE_UN)
            u = getIGABasis()->getUBSplineBasis1D()->getLastKnot();
        else if (edge == EDGE_V0)
            v = getIGABasis()->getVBSplineBasis1D()->getFirstKnot();
        else if (edge == EDGE_VN)
            v = getIGABasis()->getVBSplineBasis1D()->getLastKnot();

        // Set the output variables
        _lambda = lambda;
        _distance = distance;
        _u = u;
        _v = v;
        return edge;

    }
    DEBUG_OUT() << "in IGAPatchSurface::computePointProjectionOnPatchBoundaryBisection"<<endl;
    DEBUG_OUT()<<"\tAlgorithm has NOT CONVERGED. Relax bisection parameters in XML input file!"<<endl;
    return false;
}

bool IGAPatchSurface::solvePointProjectionOnPatchBoundaryNewtonRaphson(
        double& _t, double& _lambda, double& _distance, double* _P1,
        double* _P2, int _edge, int _maxIt, double _tol) {

    assert(_P1 != NULL);
    assert(_P2 != NULL);

    // Initialize the convergence flag to false.
    // The flag is not updated if the convergence was not met
    bool isConverged = false;

    // Initialize number of spatial dimensions
    int noCoord = 3;

    // Initialize the distance vector
    double P1Q[3];
    // Normal to the patch = Gu x Gv
    double normalSurface[3];
    // Line to be projected vector
    double P1P2[3];
    for (int i = 0; i < noCoord; i++)
        P1P2[i] = _P2[i] - _P1[i];
    // The normal of the plane made up by the vectors (GuxGv) and Q
    double n[3];

    // Initialize the base vectors and their derivatives
    double Gu[3];
    double Gv[3];
    double DGuDu[3];
    double DGvDv[3];
    double DGuDv[3];

    // Initialize Jacobian matrix
    double dR;

    // Initialize right-hand side solution vector
    double R=0,R_previous1=0,R_previous2=0;

    // Initialize the knot span indices
    int uKnotSpan = 0;
    int vKnotSpan = 0;

    // The NURBS polynomial degrees
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // Local number of basis functions
    int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

    // Initialize the Newton-Raphson iteration counter
    int counter = 0;

    // Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
    int derivDegreeBasis = 2;

    // The number of the base vectors
    int noBaseVcts = 2;

    // Initialize the array of the NURBS basis functions and their derivatives
    double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1)
            * (derivDegreeBasis + 2) * noLocalBasisFcts / 2];

    // Number of derivatives for the base vectors
    int derivDegreeBaseVcts = derivDegreeBasis - 1;

    // Initialize the array of the base vectors and their derivatives
    double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1)
            * (derivDegreeBaseVcts + 2) * noCoord * noBaseVcts / 2];

    int direction;
    double u;
    double v;
    switch (_edge) {
    case 0:
        u = _t;
        v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
        direction = 0;
        break;
    case 1:
        u = _t;
        v = IGABasis->getVBSplineBasis1D()->getLastKnot();
        direction = 0;
        break;
    case 2:
        u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
        v = _t;
        direction = 1;
        break;
    case 3:
        u = IGABasis->getUBSplineBasis1D()->getLastKnot();
        v = _t;
        direction = 1;
        break;
    }

    IGABasis->getUBSplineBasis1D()->clampKnot(u);
    IGABasis->getVBSplineBasis1D()->clampKnot(v);

    double Q[3];
    // 2. Loop over all the Newton-Raphson iterations
    while (counter < _maxIt) {
        // 2i. Update the iteration counter
        counter++;
        // 2ii. Find the span of the given surface parameters
        uKnotSpan = IGABasis->getUBSplineBasis1D()->findKnotSpan(u);
        vKnotSpan = IGABasis->getVBSplineBasis1D()->findKnotSpan(v);
        // 2iii. Compute the NURBS basis functions and their derivatives at current (_u,_v) pair of surface parameters
        IGABasis->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs,
                derivDegreeBasis, u, uKnotSpan, v, vKnotSpan);
        // 2iv. Compute the base vectors at current (_u,_v) pair of surface parameters
        computeBaseVectorsAndDerivatives(baseVecAndDerivs, basisFctsAndDerivs,
                derivDegreeBaseVcts, uKnotSpan, vKnotSpan);
        // 2v. Compute the Cartesian components of the point on the surface
        computeCartesianCoordinates(Q, basisFctsAndDerivs, derivDegreeBasis,
                uKnotSpan, vKnotSpan);
        // 2vi. Compute the distance vector between the vector to be projected and the estimated one
        for (int i = 0; i < noCoord; i++) {
            P1Q[i] = Q[i] - _P1[i];
            Gu[i] = baseVecAndDerivs[indexDerivativeBaseVector(
                    derivDegreeBaseVcts, 0, 0, i, 0)];
            Gv[i] = baseVecAndDerivs[indexDerivativeBaseVector(
                    derivDegreeBaseVcts, 0, 0, i, 1)];
            DGuDu[i] = baseVecAndDerivs[indexDerivativeBaseVector(
                    derivDegreeBaseVcts, 1, 0, i, 0)];
            DGvDv[i] = baseVecAndDerivs[indexDerivativeBaseVector(
                    derivDegreeBaseVcts, 0, 1, i, 1)];
            DGuDv[i] = baseVecAndDerivs[indexDerivativeBaseVector(
                    derivDegreeBaseVcts, 0, 1, i, 0)];
        }
        // Update normal and compute the residual
        MathLibrary::crossProduct(normalSurface, Gu, Gv);
        MathLibrary::crossProduct(n, P1P2, normalSurface);
        R_previous2 = R_previous1;
        R_previous1 = R;
        R = MathLibrary::computeDenseDotProduct(noCoord, P1Q, n);
        // 2vii. Compute the stopping criteria
        // If the residual is low enough stop and set the flag
        if(fabs(R) < _tol) {
            isConverged = true;
            break;
        }
        // If oscillating stop
        if(R==R_previous2) {
            break;
        }

        // Declare temporary cross product required for derivation of normal n
        double product1[3];
        double product2[3];
        double product3[3];
        // Compute derivatives of the normal n into product1 product2 and product3 and then sums them correctly
        if (direction == 0) {
            MathLibrary::crossProduct(product1, DGuDu, Gv);
            MathLibrary::crossProduct(product2, Gu, DGuDv);
            dR = MathLibrary::computeDenseDotProduct(noCoord, Gu, n);
        } else {
            MathLibrary::crossProduct(product1, DGuDv, Gv);
            MathLibrary::crossProduct(product2, Gu, DGvDv);
            dR = MathLibrary::computeDenseDotProduct(noCoord, Gv, n);
        }
        for (int i = 0; i < noCoord; i++)
            product1[i] += product2[i];
        MathLibrary::crossProduct(product3, P1P2, product1);
        dR += MathLibrary::computeDenseDotProduct(noCoord, P1Q, product3);

        // Apply the variation dw on parameter and check if going outside of knot span then break
        double dw=-R/dR;
        if (direction == 0) {
            u += dw;
        } else {
            v += dw;
        }
        // 2vii. Check and modify the surface parameters if they stay out of their knot spans
        IGABasis->getUBSplineBasis1D()->clampKnot(u);
        IGABasis->getVBSplineBasis1D()->clampKnot(v);
    }

    if (direction == 0)
        _t = u;
    else
        _t = v;

    // Compute point on the line and thus get the ratio of P1P/P1P2
    double normNormalSurface = MathLibrary::vector2norm(normalSurface,noCoord);
    double unitNormalSurface[3];
    double P1P[3];
    for (int i = 0; i < noCoord; i++)
        unitNormalSurface[i] = normalSurface[i]/normNormalSurface;
    // Project P1Q onto the normal of the patch
    double h10 = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, P1Q, unitNormalSurface);
    for (int i = 0; i < noCoord; i++) {
        P1P[i] = P1Q[i];
        P1P[i] -= h10 * unitNormalSurface[i];
    }

    // Compute the collinearity of P1P and P1P2
    double unitP1P[3], unitP1P2[3];
    double normP1P = MathLibrary::vector2norm(P1P,noCoord);
    double normP1P2 = MathLibrary::vector2norm(P1P2,noCoord);
    for (int i = 0; i < noCoord; i++) {
        unitP1P[i] = P1P[i] / normP1P;
        unitP1P2[i] = P1P2[i] / normP1P2;
    }
    double collinearity = fabs(EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, unitP1P, unitP1P2));

    // Update the flag depending on the collinearity.
    // This check is performed here since the point P is reconstructed and not returned to the upper level
    isConverged = isConverged && (1 - collinearity <= _tol );

    // Compute the output variables if the algorithm validly converged
    if (isConverged) {
        // Compute output parameters _lambda and _distance
        _lambda = normP1P / normP1P2;
        double QP[3];
        for (int i = 0; i < noCoord; i++) {
            QP[i] = _P1[i] + _lambda * (_P2[i] - _P1[i]) - Q[i];
        }
        _distance = MathLibrary::vector2norm(QP,noCoord);
    }

    // 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
    // Clear the memory on the heap
    delete[] basisFctsAndDerivs;
    delete[] baseVecAndDerivs;

    // Return the flag
    return isConverged;
}

bool IGAPatchSurface::solvePointProjectionOnPatchBoundaryNewtonRaphsonClosestDistance(
        double& _t, double& _lambda, double& _distance, double* _P1,
        double* _P2, int _edge, int _maxIt, double _tol) {

    // Check input
    assert(_P1 != NULL);
    assert(_P2 != NULL);
    assert(_edge >= 0 && _edge <= 4);

    // Initialize variables
    double u,v;
    double lambda = 0.0;
    bool isNRConverged = false;
    bool isOnU;
    double distance;

    // Initialize number of spatial dimensions
    int dim = 3;

    // Line to be projected vector
    double P1P2[3];
    for (int i = 0; i < dim; i++)
        P1P2[i] = _P2[i] - _P1[i];

    // Initialize the base vectors and their derivatives
    double G[dim];
    double GSum[dim];
    double DGSum[dim];
    double DG[dim];

    int knotSpanIndexFixed;

    // Initialize the knot span indices
    int uKnotSpan = 0;
    int vKnotSpan = 0;

    // The NURBS polynomial degrees
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // Local number of basis functions
    int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

    // Initialize the Newton-Raphson iteration counter
    int counter = 0;

    // Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
    int derivDegreeBasis = 2;

    // The number of the base vectors
    int noBaseVcts = 2;

    // Initialize the array of the NURBS basis functions and their derivatives
    double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1) * (derivDegreeBasis + 2) * noLocalBasisFcts / 2];

    // Number of derivatives for the base vectors
    int derivDegreeBaseVcts = derivDegreeBasis - 1;

    // Initialize the array of the base vectors and their derivatives
    double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1)
            * (derivDegreeBaseVcts + 2) * dim * noBaseVcts / 2];

    // Get the fixed and the running parameters on the edge
    switch (_edge) {
    case 0:
        u = _t;
        v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
        isOnU = true;
        break;
    case 1:
        u = _t;
        v = IGABasis->getVBSplineBasis1D()->getLastKnot();
        isOnU = true;
        break;
    case 2:
        u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
        v = _t;
        isOnU = false;
        break;
    case 3:
        u = IGABasis->getUBSplineBasis1D()->getLastKnot();
        v = _t;
        isOnU = false;
        break;
    }
    IGABasis->getUBSplineBasis1D()->clampKnot(u);
    IGABasis->getVBSplineBasis1D()->clampKnot(v);
    if (isOnU) {
        knotSpanIndexFixed = IGABasis->getVBSplineBasis1D()->findKnotSpan(v);
    } else {
        knotSpanIndexFixed = IGABasis->getUBSplineBasis1D()->findKnotSpan(u);
    }

    // The point to be projected on the line
    double P[3];

    // The projected point on the patch surface
    double Q[3];
    double PQ[3];
    double R[2];
    double J[4];
    J[3] = EMPIRE::MathLibrary::computeDenseDotProduct(dim, P1P2, P1P2);
    bool isLinearEquationSolved = false;

    // 2. Loop over the Newton-Raphson iterations
    while (counter < _maxIt) {
        // 2i. Update the iteration counter
        counter++;

        // 2ii. Compute the NURBS basis functions and their derivatives at current (u,v) pair of surface parameters
        if (isOnU) {
            uKnotSpan = IGABasis->getUBSplineBasis1D()->findKnotSpan(u);
            vKnotSpan = knotSpanIndexFixed;
        } else {
            uKnotSpan = knotSpanIndexFixed;
            vKnotSpan = IGABasis->getVBSplineBasis1D()->findKnotSpan(v);
        }
        IGABasis->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs,
                derivDegreeBasis, u, uKnotSpan, v, vKnotSpan);

        // 2iv. Compute the base vectors at current (u,v) pair of surface parameters
        computeBaseVectorsAndDerivatives(baseVecAndDerivs, basisFctsAndDerivs,
                derivDegreeBaseVcts, uKnotSpan, vKnotSpan);

        // 2v. Compute the Cartesian coordinates of the point on the surface
        computeCartesianCoordinates(Q, basisFctsAndDerivs, derivDegreeBasis,
                uKnotSpan, vKnotSpan);

        // Compute the base vector and its derivative along the running parameter
//        for (int i = 0 ; i < dim; i++) {
//            GSum[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 0)] +
//                    baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 1)];
//            if (isOnU) {
//                G[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 0)];
//                DGSum[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 1, 0, i, 0)] +
//                        baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 0)];
//            } else {
//                G[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 1)];
//                DGSum[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 1)] +
//                        baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 0)];
//            }
//        }
        for (int i = 0 ; i < dim; i++) {
            if (isOnU) {
                G[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 0)];
                DG[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 1, 0, i, 0)];
            } else {
                G[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 1)];
                DG[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 0)];
            }
        }

        // Update the P and compute the distance between P and Q
        for (int i = 0; i < dim; i++) {
            P[i] = (1.0 - lambda) * _P1[i] + lambda * _P2[i];
            PQ[i] = Q[i] - P[i];
        }

        // Compute the residual vector
        R[0] = -EMPIRE::MathLibrary::computeDenseDotProduct(dim, G, PQ);
        R[1] = -EMPIRE::MathLibrary::computeDenseDotProduct(dim, P1P2, PQ);
        distance = EMPIRE::MathLibrary::vector2norm(PQ, dim);

        // Check convergence
        if ((distance <=  _tol || EMPIRE::MathLibrary::vector2norm(R, 2) <=  _tol) && counter > 1){
            isNRConverged = true;
            break;
        }

        // Compute the Jacobian
        J[0] = EMPIRE::MathLibrary::computeDenseDotProduct(dim, PQ, DG) + EMPIRE::MathLibrary::computeDenseDotProduct(dim, G, G);
        J[1] = EMPIRE::MathLibrary::computeDenseDotProduct(dim, P1P2, G);
        J[2] = EMPIRE::MathLibrary::computeDenseDotProduct(dim, G, P1P2);

        // Solve the 2x2 system
        isLinearEquationSolved = EMPIRE::MathLibrary::solve2x2LinearSystem(J, R);

        // Update the parametric coordinates
        if (isOnU)
            u += R[0];
        else
            v += R[0];
        lambda += R[1];

        // Clamp the running curve parameter to the knot limits
        if (isOnU) {
            IGABasis->getUBSplineBasis1D()->clampKnot(u);
        } else {
            IGABasis->getVBSplineBasis1D()->clampKnot(v);
        }

        // Clamp the running line parameter to the 0-1 limits
        if (lambda < 0)
            lambda = 0.0;

        if (lambda > 1.0)
            lambda = 1.0;
    }

    // Update output variables
    if(isOnU)
        _t = u;
    else
        _t = v;
    _lambda = lambda;
    _distance = distance;

    // Delete pointers
    delete[] basisFctsAndDerivs;
    delete[] baseVecAndDerivs;

    // Return convergence flag
    return isNRConverged;
}

char IGAPatchSurface::computePointProjectionOnPatchBoundaryNewtonRhapson(double& _u, double& _v, double& _lambda,
        double& _distance, double* _P1, double* _P2, int _maxIt, double _tol) {
    DEBUG_OUT()<<"\t======================================================"<<endl;
    DEBUG_OUT() << "in IGAPatchSurface::computePointProjectionOnPatchBoundaryNewtonRhapson"<<endl;
    DEBUG_OUT()<<"\tPROJECT line on boundary using Newton-Rhapson for line"<<endl
            <<"\t\t(("<<_P1[0]<<" , "<<_P1[1]<<" , "<<_P1[2]<<");"
            <<"("<<_P2[0]<<" , "<<_P2[1]<<" , "<<_P2[2]<<")) "<<endl
            <<"\t\twith initial guess is projection of P1 : (u,v)=("<<_u<<" , "<<_v<<")"<<endl;
    double u1 = _u;
    double v1 = _v;
    double t;
    double maxDistance = _distance;
    double distance = numeric_limits<double>::max();
    _distance = numeric_limits<double>::max();
    double lambda = numeric_limits<double>::max();
    bool isConverged = false;

    char edgeOut=0, _edge=0;
    // Loop over all the edges of the NURBS patch (for tensor product surfaces there are 4 edges)
    for (int edge = 0; edge < 4; edge++) {
        // Store intermediate results
        double u[3];
        double v[3];
        bool hasConverged=false;

        // Do the test for every edge for each extremity of the boundary patch
        for(int point=0;point<3;point++) {
            // Find the fixed and the running parameter on the patch boundary
            if (edge == 0 || edge == 1)
                if(point==0)
                    t=IGABasis->getUBSplineBasis1D()->getFirstKnot();
                else
                    t=IGABasis->getUBSplineBasis1D()->getLastKnot();
            else
                if(point==0)
                    t=IGABasis->getVBSplineBasis1D()->getFirstKnot();
                else
                    t=IGABasis->getVBSplineBasis1D()->getLastKnot();

            // If extremities of edge have not converged then try with initial guess
            if(point==2 && !hasConverged) {
                if (edge == 0 || edge == 1)
                    t = u1;
                else
                    t = v1;
            }
            // Compute point projection from the line to the NURBS patch boundary
            isConverged = solvePointProjectionOnPatchBoundaryNewtonRaphson(t, lambda, distance, _P1, _P2, edge, _maxIt, _tol);
//            isConverged = solvePointProjectionOnPatchBoundaryNewtonRaphsonClosestDistance(t, lambda, distance, _P1, _P2, edge, _maxIt, _tol);

            // Fix possible numerical error
            if(lambda > 1.0 && lambda-1.0 <= 1e-3) {
                lambda = 1.0;
                DEBUG_OUT("In IGAPatchSurface::computePointProjectionOnPatchBoundaryNewtonRhapson, line parameter clamped to 1 !");
            }
            if(lambda < 0 && lambda >= -1e-3) {
                lambda = 0.0;
                DEBUG_OUT("In IGAPatchSurface::computePointProjectionOnPatchBoundaryNewtonRhapson, line parameter clamped to 0 !");
            }

            // Update the flag depending on the validity and quality of the projection compared to the previous found results
            isConverged = isConverged && lambda >= 0.0 && lambda <= 1.0 && distance <= _distance && distance <= maxDistance;

            if (isConverged) {
                switch (edge) {
                case 0:
                    edgeOut=EDGE_V0;
                    IGABasis->getUBSplineBasis1D()->clampKnot(t);
                    u[point] = t;
                    v[point] = IGABasis->getVBSplineBasis1D()->getFirstKnot();
                    break;
                case 1:
                    edgeOut=EDGE_VN;
                    IGABasis->getUBSplineBasis1D()->clampKnot(t);
                    u[point] = t;
                    v[point] = IGABasis->getVBSplineBasis1D()->getLastKnot();
                    break;
                case 2:
                    edgeOut=EDGE_U0;
                    IGABasis->getVBSplineBasis1D()->clampKnot(t);
                    u[point] = IGABasis->getUBSplineBasis1D()->getFirstKnot();
                    v[point] = t;
                    break;
                case 3:
                    edgeOut=EDGE_UN;
                    IGABasis->getVBSplineBasis1D()->clampKnot(t);
                    u[point] = IGABasis->getUBSplineBasis1D()->getLastKnot();
                    v[point] = t;
                    break;
                }
                // If the point is the same as the initial guess from input function
                bool validPoint1=(u[point]!=u1 || v[point]!=v1); //Different from entry point then true else false
                bool validPoint2=(point==1)?(u[0]==u1 && v[0]==v1 && u[1]==u1 && v[1]==v1):false; //Both are the same then true else false
                // If it is not a point to take into account, continue
                if(!(validPoint1 || validPoint2) && point<2) continue;

                _u = u[point];
                _v = v[point];
                _distance = distance;
                _lambda = lambda;
                _edge = _edge | edgeOut;
                DEBUG_OUT()<<"\t-------------------------------------------------------"<<endl;
                DEBUG_OUT()<<"\tAlgorithm has CONVERGED for initial guess point "<<point<<" on edge "<<int(edgeOut)<<endl
                          <<"\t\t and new distance found to patch is "<<_distance<<endl
                         <<"\t\t and ratio P1P/P1P2 is "<<_lambda<<endl
                        <<"\t\t and parametric value are (u,v)("<<_u<<" , "<<_v<<")"<<endl;
            }
        }
    }
    if(!_edge) {
        _u = u1;
        _v = v1;
        DEBUG_OUT() << "in IGAPatchSurface::computePointProjectionOnPatchBoundaryNewtonRhapson"<<endl;
        DEBUG_OUT()<<"\tAlgorithm has NOT CONVERGED. Relax newtonRaphsonBoundary parameters in XML input file!"<<endl;
    }
    DEBUG_OUT()<<"\t======================================================"<<endl;
    return _edge;
}

void IGAPatchSurface::findInitialGuess4PointProjection(double& _u, double& _v, double* _P,
        int _uDiv, int _vDiv) {

    assert(_P != NULL);

    const int noSpatialDimensions = 3;

    // Number of Basis function for U and V
    int nU = IGABasis->getUBSplineBasis1D()->computeNoBasisFunctions();
    int nV = IGABasis->getVBSplineBasis1D()->computeNoBasisFunctions();

    // The NURBS polynomial degrees for U and V
    int pU = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int pV = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // Knot vector for U and V
    double *knotU = IGABasis->getUBSplineBasis1D()->getKnotVector();
    double *knotV = IGABasis->getVBSplineBasis1D()->getKnotVector();

    double coords[3];
    double minDis = numeric_limits<double>::max();
    double Dis;

    double u0 = knotU[pU];
    double v0 = knotV[pV];
    double uEnd = knotU[nU];
    double vEnd = knotV[nV];
    double du = (uEnd - u0) / _uDiv;
    double dv = (vEnd - v0) / _vDiv;
    double uv[2];

    for (int i = 1; i < _uDiv - 1; i++)
        for (int j = 1; j < _vDiv - 1; j++) {
            uv[0] = u0 + du * i;
            uv[1] = v0 + dv * j;
            computeCartesianCoordinates(coords, uv);

            for (int k = 0; k < noSpatialDimensions; k++)
                coords[k] -= _P[k];
            Dis = MathLibrary::vector2norm(coords, noSpatialDimensions);
            if (Dis < minDis) {
                minDis = Dis;
                _u = uv[0];
                _v = uv[1];
            }

        }
}

bool IGAPatchSurface::findInitialGuess4PointProjectionOnTrimmingCurve(double& _uTilde, double& _u, double& _v,
                                                                      double* _P, int _patchBLIndex, int _patchBLTrCurveIndex){
    /*
     * This function computes an initial guess for projection of a cartesian point onto a trimming curve belonging to a patch.
     * The algorithm advances on the linearization vertices of the trimming curve and returns
     * the patch parameters of the vertex which has the closest cartesian distance to the point to be projected.
     */

    const int noCoordParam = 2;
    const int noCoord = 3;

    double theVertex[noCoordParam];
    double xyzCoord[noCoord];
    double distance = 1e3;
    double tmpDistance;
    int xyzCtr = 0;
    bool isFound;
    double minNumVertices = 51;

    // Number of existing linearziation vertices
    int numVertices = Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getPolyline()->size()/noCoordParam;

    // Candidates for the initial guess
    std::vector<double> candidatesUTilde;
    std::vector<double> candidatesUV;

    // If number of linearization vertices are less than 10 then create equidistant points on the curve parameter space to get an initial guess
    if (numVertices<minNumVertices){
        double u0 = Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getIGABasis()->getFirstKnot();
        double u1 = Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getIGABasis()->getLastKnot();
        double du = (u1-u0)/(minNumVertices-1);
        for(int i=0;i<minNumVertices;i++) {
            double knot = u0 + i*du;
            double parametricCoordinates[2] = {0};
            candidatesUTilde.push_back(knot);
            Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).computeCartesianCoordinates(parametricCoordinates,knot);
            candidatesUV.push_back(parametricCoordinates[0]);
            candidatesUV.push_back(parametricCoordinates[1]);
        }
    } // else use the linearization vertices for getting an initial guess
    else {
        candidatesUTilde = *Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getPolylineKnots();
        candidatesUV = *Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getPolyline();
    }

    // For each vertex
    for(int vertexCtr=0; vertexCtr<candidatesUTilde.size(); vertexCtr++){

        // Get the vertex
        theVertex[0] = candidatesUV[vertexCtr*noCoordParam];
        theVertex[1] = candidatesUV[vertexCtr*noCoordParam+1];

        // Compute the physical coordinates and the distance
        computeCartesianCoordinates(xyzCoord, theVertex);
        for(xyzCtr=0; xyzCtr < noCoord; xyzCtr++)
            xyzCoord[xyzCtr] = xyzCoord[xyzCtr] - _P[xyzCtr];

        // Compute the length of the distance vector
        tmpDistance = MathLibrary::vector2norm(xyzCoord,noCoord);

        // Update the distance and the patch parameters
        if (tmpDistance<distance) {
            distance = tmpDistance;
            _uTilde = candidatesUTilde[vertexCtr];
            _u = theVertex[0];
            _v = theVertex[1];
            isFound = true;
        }
    }


    if (!isFound)
        ERROR_BLOCK_OUT("IGAPatchSurface","findInitialGuess4PointProjectionOnTrimmingCurve","Initial guess for projection cannot be found!");

}

bool IGAPatchSurface::findInitialGuess4PointProjectionOnTrimmingCurve(double& _uTilde, double& _u, double& _v,
                                                                      double* _P, IGAPatchCurve* _curve){
    /*
     * This function computes an initial guess for projection of a cartesian point onto a trimming curve belonging to a patch.
     * The algorithm advances on the linearization vertices of the trimming curve and returns
     * the patch parameters of the vertex which has the closest cartesian distance to the point to be projected.
     */

    const int noCoordParam = 2;
    const int noCoord = 3;

    double theVertex[noCoordParam];
    double xyzCoord[noCoord];
    double distance = 1e3;
    double tmpDistance;
    int xyzCtr = 0;
    bool isFound;
    double minNumVertices = 51;

    // Number of existing linearziation vertices
    int numVertices = _curve->getPolyline()->size()/noCoordParam;

    // Candidates for the initial guess
    std::vector<double> candidatesUTilde;
    std::vector<double> candidatesUV;

    // If number of linearization vertices are less than 10 then create equidistant points on the curve parameter space to get an initial guess
    if (numVertices<minNumVertices){
        double u0 = _curve->getIGABasis()->getFirstKnot();
        double u1 = _curve->getIGABasis()->getLastKnot();
        double du = (u1-u0)/(minNumVertices-1);
        for(int i=0;i<minNumVertices;i++) {
            double knot = u0 + i*du;
            double parametricCoordinates[2] = {0};
            candidatesUTilde.push_back(knot);
            _curve->computeCartesianCoordinates(parametricCoordinates,knot);
            candidatesUV.push_back(parametricCoordinates[0]);
            candidatesUV.push_back(parametricCoordinates[1]);
        }
    } // else use the linearization vertices for getting an initial guess
    else {
        candidatesUTilde = *_curve->getPolylineKnots();
        candidatesUV = *_curve->getPolyline();
    }

    // For each vertex
    for(int vertexCtr=0; vertexCtr<candidatesUTilde.size(); vertexCtr++){

        // Get the vertex
        theVertex[0] = candidatesUV[vertexCtr*noCoordParam];
        theVertex[1] = candidatesUV[vertexCtr*noCoordParam+1];

        // Compute the physical coordinates and the distance
        computeCartesianCoordinates(xyzCoord, theVertex);
        for(xyzCtr=0; xyzCtr < noCoord; xyzCtr++)
            xyzCoord[xyzCtr] = xyzCoord[xyzCtr] - _P[xyzCtr];

        // Compute the length of the distance vector
        tmpDistance = MathLibrary::vector2norm(xyzCoord,noCoord);

        // Update the distance and the patch parameters
        if (tmpDistance<distance) {
            distance = tmpDistance;
            _uTilde = candidatesUTilde[vertexCtr];
            _u = theVertex[0];
            _v = theVertex[1];
            isFound = true;
        }
    }


    if (!isFound)
        ERROR_BLOCK_OUT("IGAPatchSurface","findInitialGuess4PointProjectionOnTrimmingCurve","Initial guess for projection cannot be found!");

}

void IGAPatchSurface::computeCartesianCoordinatesAndNormalVector(double* _coords, double* _normal,
        double _u, double _v)const {

    // The knot spans
    int spanU = IGABasis->getUBSplineBasis1D()->findKnotSpan(_u);
    int spanV = IGABasis->getVBSplineBasis1D()->findKnotSpan(_v);

    // Compute the Cartesian coordinates of (_u,_v)
    computeCartesianCoordinates(_coords, _u, spanU, _v, spanV);

    // Compute the base vectors
    double baseVec[6];
    computeBaseVectors(baseVec, _u, spanU, _v, spanV);

    // Compute the cross product of the surface base vectors to get the surface normal
    _normal[0] = baseVec[1] * baseVec[5] - baseVec[2] * baseVec[4];
    _normal[1] = baseVec[2] * baseVec[3] - baseVec[0] * baseVec[5];
    _normal[2] = baseVec[0] * baseVec[4] - baseVec[1] * baseVec[3];
}
void IGAPatchSurface::computeCartesianCoordinatesAndNormalVector(double* _coords, double* _normal,
        double _u, double _v, int _spanU, int _spanV) const {
    // Compute the Cartesian coordinates of (_u,_v)
    computeCartesianCoordinates(_coords, _u, _spanU, _v, _spanV);

    // Compute the base vectors
    double baseVec[6];
    computeBaseVectors(baseVec, _u, _spanU, _v, _spanV);

    // Compute the cross product of the surface base vectors to get the surface normal
    MathLibrary::computeVectorCrossProduct(&baseVec[0], &baseVec[3],_normal);
}

void IGAPatchSurface::computeKnotIntersectionsWithTrimmingCurve(std::vector<double>& _uTilde,
                                                                int _patchBLIndex, int _patchBLTrCurveIndex){

    /*
     * 1. Get the knots and the trimming curve to perform the intersection
     * 2. Compute U and V intersections separately. At the end of the process _uTilde is sorted in ascending order
     * 3. Given the curve parameters, compute UV parameters on the patch and the global Cartesian coordinates
     */

    // Define an equality tolerance which is actually the tolerance of the bisection algorithm that is used here
    double tol = 1e-6;

    // Initialize parameter line direction
    int dirU = 0;
    int dirV = 1;

    // Compute knot intersections
    double knot;
    // Compute U knot intersections in V direction
    for (int knotCtr = 0; knotCtr < IGABasis->getUBSplineBasis1D()->getNoKnots(); knotCtr++){
        if (knotCtr == 0) knot = IGABasis->getUBSplineBasis1D()->getKnotVector()[knotCtr];
        else if (knot == IGABasis->getUBSplineBasis1D()->getKnotVector()[knotCtr]) continue;
        else knot = IGABasis->getUBSplineBasis1D()->getKnotVector()[knotCtr];

        Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).computeIntersectionsWithKnotBisection(_uTilde, dirV, knot);
    }

    // Compute V knot intersections in U direction
    for (int knotCtr = 0; knotCtr < IGABasis->getVBSplineBasis1D()->getNoKnots(); knotCtr++){
        if (knotCtr == 0) knot = IGABasis->getVBSplineBasis1D()->getKnotVector()[knotCtr];
        else if (knot == IGABasis->getVBSplineBasis1D()->getKnotVector()[knotCtr]) continue;
        else knot = IGABasis->getVBSplineBasis1D()->getKnotVector()[knotCtr];

        Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).computeIntersectionsWithKnotBisection(_uTilde, dirU, knot);
    }

    // Check if intersections exist and the first/last knot is included
    if (!_uTilde.empty()) {
        EMPIRE::MathLibrary::sortRemoveDuplicates(_uTilde);
        // Check for the first knot
        if (fabs (*_uTilde.begin() - Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getIGABasis()->getFirstKnot())  > tol) {
            _uTilde.insert(_uTilde.begin(), Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getIGABasis()->getFirstKnot());
            DEBUG_OUT("In \"IGAPatchSurface::computeKnotIntersectionsWithTrimmingCurve\": Intersections exist but the first knot is not found!");
            DEBUG_OUT("Adding the first knot.");
        }
        // Check for the last knot
        if (fabs (*(_uTilde.end()-1) - Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getIGABasis()->getLastKnot()) > tol) {
            _uTilde.push_back(Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getIGABasis()->getLastKnot());
            DEBUG_OUT("In \"IGAPatchSurface::computeKnotIntersectionsWithTrimmingCurve\": Intersections exist but the last knot is not found!");
            DEBUG_OUT("Adding the last knot.");
        }
    }
    // Check if the trimming curve is not intersected by any knot then add the first and the last knots
    else if (_uTilde.empty()) {
        DEBUG_OUT("In \"IGAPatchSurface::computeKnotIntersectionsWithTrimmingCurve\": No intersections found!");
        DEBUG_OUT("Adding the first and the last knot.");
        _uTilde.push_back(Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getIGABasis()->getFirstKnot());
        _uTilde.push_back(Trimming.getLoop(_patchBLIndex).getIGACurve(_patchBLTrCurveIndex).getIGABasis()->getLastKnot());
    }
}

void IGAPatchSurface::computeKnotIntersectionsWithTrimmingCurve(std::vector<double>& _uTilde,
                                                                IGAPatchCurve* _theCurve){

    /*
     * 1. Get the knots and the trimming curve to perform the intersection
     * 2. Compute U and V intersections separately. At the end of the process _uTilde is sorted in ascending order
     * 3. Given the curve parameters, compute UV parameters on the patch and the global Cartesian coordinates
     */

    // Define an equality tolerance which is actually the tolerance of the bisection algorithm that is used here
    double tol = 1e-6;

    // Initialize parameter line direction
    int dirU = 0;
    int dirV = 1;

    // Compute knot intersections
    double knot;
    // Compute U knot intersections in V direction
    for (int knotCtr = 0; knotCtr < IGABasis->getUBSplineBasis1D()->getNoKnots(); knotCtr++){
        if (knotCtr == 0) knot = IGABasis->getUBSplineBasis1D()->getKnotVector()[knotCtr];
        else if (knot == IGABasis->getUBSplineBasis1D()->getKnotVector()[knotCtr]) continue;
        else knot = IGABasis->getUBSplineBasis1D()->getKnotVector()[knotCtr];

        _theCurve->computeIntersectionsWithKnotBisection(_uTilde, dirV, knot);
    }

    // Compute V knot intersections in U direction
    for (int knotCtr = 0; knotCtr < IGABasis->getVBSplineBasis1D()->getNoKnots(); knotCtr++){
        if (knotCtr == 0) knot = IGABasis->getVBSplineBasis1D()->getKnotVector()[knotCtr];
        else if (knot == IGABasis->getVBSplineBasis1D()->getKnotVector()[knotCtr]) continue;
        else knot = IGABasis->getVBSplineBasis1D()->getKnotVector()[knotCtr];

        _theCurve->computeIntersectionsWithKnotBisection(_uTilde, dirU, knot);
    }

    // Check if intersections exist and the first/last knot is included
    if (!_uTilde.empty()) {
        EMPIRE::MathLibrary::sortRemoveDuplicates(_uTilde);
        // Check for the first knot
        if (fabs (*_uTilde.begin() - _theCurve->getIGABasis()->getFirstKnot())  > tol) {
            _uTilde.insert(_uTilde.begin(), _theCurve->getIGABasis()->getFirstKnot());
            DEBUG_OUT("In \"IGAPatchSurface::computeKnotIntersectionsWithTrimmingCurve\": Intersections exist but the first knot is not found!");
            DEBUG_OUT("Adding the first knot.");
        }
        // Check for the last knot
        if (fabs (*(_uTilde.end()-1) - _theCurve->getIGABasis()->getLastKnot()) > tol) {
            _uTilde.push_back(_theCurve->getIGABasis()->getLastKnot());
            DEBUG_OUT("In \"IGAPatchSurface::computeKnotIntersectionsWithTrimmingCurve\": Intersections exist but the last knot is not found!");
            DEBUG_OUT("Adding the last knot.");
        }
    }
    // Check if the trimming curve is not intersected by any knot then add the first and the last knots
    else if (_uTilde.empty()) {
        DEBUG_OUT("In \"IGAPatchSurface::computeKnotIntersectionsWithTrimmingCurve\": No intersections found!");
        DEBUG_OUT("Adding the first and the last knot.");
        _uTilde.push_back(_theCurve->getIGABasis()->getFirstKnot());
        _uTilde.push_back(_theCurve->getIGABasis()->getLastKnot());
    }
}

char IGAPatchSurface::getEdge(const double _u, const double _v, const double _tolerance) {
	char edge = 0;
	if(fabs(_u - getIGABasis(0)->getFirstKnot()) < _tolerance) {
		edge = edge | EDGE_U0;
	}
	else if(fabs(_u - getIGABasis(0)->getLastKnot()) < _tolerance) {
		edge = edge | EDGE_UN;
	}
	if(fabs(_v - getIGABasis(1)->getFirstKnot()) < _tolerance) {
		edge = edge | EDGE_V0;
	}
	else if(fabs(_v - getIGABasis(1)->getLastKnot()) < _tolerance) {
		edge = edge | EDGE_VN;
	}
	return edge;
}
std::vector<std::pair<double,double> > IGAPatchSurface::getCorner(const char _edgeIn, const char _edgeOut, const bool _isCounterclockwise) {
	std::vector<std::pair<double,double> > corners;
	// If share one common edge means point does not exist
	// Or that the corner is already present
	if(_edgeIn & _edgeOut)
		return corners;
	// If the input is correct
	if(_edgeIn == 0 || _edgeOut == 0) {
		WARNING_BLOCK_OUT("IGAPatchSurface","getCorner","One of the input edge is not defined!");
		return corners;
	}
	double u0 = getIGABasis()->getUBSplineBasis1D()->getFirstKnot();
	double uN = getIGABasis()->getUBSplineBasis1D()->getLastKnot();
	double v0 = getIGABasis()->getVBSplineBasis1D()->getFirstKnot();
	double vN = getIGABasis()->getVBSplineBasis1D()->getLastKnot();
	const char corner = _edgeIn | _edgeOut;
	if(corner ==(EDGE_U0 | EDGE_V0)) {
		corners.push_back(make_pair(u0,v0));
		return corners;
	}
	if(corner ==(EDGE_U0 | EDGE_VN)) {
		corners.push_back(make_pair(u0,vN));
		return corners;
	}
	if(corner ==(EDGE_UN | EDGE_V0)) {
		corners.push_back(make_pair(uN,v0));
		return corners;
	}
	if(corner ==(EDGE_UN | EDGE_VN)) {
		corners.push_back(make_pair(uN,vN));
		return corners;
	}
	if(_edgeIn == EDGE_U0 && _edgeOut == EDGE_UN && _isCounterclockwise) {
		corners.push_back(make_pair(uN,vN));
		corners.push_back(make_pair(u0,vN));
		return corners;
	}
	if(_edgeIn == EDGE_U0 && _edgeOut == EDGE_UN && !_isCounterclockwise) {
		corners.push_back(make_pair(uN,v0));
		corners.push_back(make_pair(u0,v0));
		return corners;
	}
	if(_edgeIn == EDGE_UN && _edgeOut == EDGE_U0 && _isCounterclockwise) {
		corners.push_back(make_pair(u0,v0));
		corners.push_back(make_pair(uN,v0));
		return corners;
	}
	if(_edgeIn == EDGE_UN && _edgeOut == EDGE_U0 && !_isCounterclockwise) {
		corners.push_back(make_pair(u0,vN));
		corners.push_back(make_pair(uN,vN));
		return corners;
	}
	if(_edgeIn == EDGE_V0 && _edgeOut == EDGE_VN && _isCounterclockwise) {
		corners.push_back(make_pair(u0,vN));
		corners.push_back(make_pair(u0,v0));
		return corners;
	}
	if(_edgeIn == EDGE_V0 && _edgeOut == EDGE_VN && !_isCounterclockwise) {
		corners.push_back(make_pair(uN,vN));
		corners.push_back(make_pair(uN,v0));
		return corners;
	}
	if(_edgeIn == EDGE_VN && _edgeOut == EDGE_V0 && _isCounterclockwise) {
		corners.push_back(make_pair(uN,v0));
		corners.push_back(make_pair(uN,vN));
		return corners;
	}
	if(_edgeIn == EDGE_VN && _edgeOut == EDGE_V0 && !_isCounterclockwise) {
		corners.push_back(make_pair(u0,v0));
		corners.push_back(make_pair(u0,vN));
		return corners;
	}
	if(_edgeIn == EDGE_VN && _edgeOut == EDGE_V0 && _isCounterclockwise) {
		corners.push_back(make_pair(uN,v0));
		corners.push_back(make_pair(uN,vN));
		return corners;
	}
	if(_edgeIn == EDGE_VN && _edgeOut == EDGE_V0 && !_isCounterclockwise) {
		corners.push_back(make_pair(u0,v0));
		corners.push_back(make_pair(u0,vN));
		return corners;
	}
	ERROR_OUT()<<"No corner found to add in polygon"<<endl;
	ERROR_OUT()<<"Edge going IN the patch is ["<<int(_edgeIn)<<"] and Edge going OUT is [" << int(_edgeOut) << "] with direction "
			<<(_isCounterclockwise?"counterclockwise":"clockwise")<<endl;
	exit(EXIT_FAILURE);
	return corners;
}

Message &operator<<(Message &message, const IGAPatchSurface &mesh) {
//  message << "\t" << "IGA Patch name: " << mesh.name << endl;

    message << "\t\tpDegree:  " << mesh.getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
            << endl;
    message << "\t\tqDegree:  " << mesh.getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
            << endl;

    message << "\t\tKnots Vector U: \t";
    for (int i = 0; i < mesh.getIGABasis()->getUBSplineBasis1D()->getNoKnots(); i++)
        message << mesh.getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] << "  ";
    message << endl;

    message << "\t\tKnots Vector V: \t";
    for (int i = 0; i < mesh.getIGABasis()->getVBSplineBasis1D()->getNoKnots(); i++)
        message << mesh.getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] << "  ";
    message << endl;

    message << "\t\t" << "number of control points U: " << mesh.getUNoControlPoints() << endl;
    message << "\t\t" << "number of control points V: " << mesh.getVNoControlPoints() << endl;

    message << "\t\tControl Points Net: " << endl;
    int count = 0;
    for (int j = 0; j < mesh.getVNoControlPoints(); j++) {
        message << "\t\t";
        for (int i = 0; i < mesh.getUNoControlPoints(); i++) {
            message << mesh.getControlPointNet()[count]->getX() << ", "
                    << mesh.getControlPointNet()[count]->getY() << ", "
                    << mesh.getControlPointNet()[count]->getZ() << "\t";
            count++;
        }
        message << endl;
    }
//    //Printf patch to be written in C++ unit style
//    count=0;
//    for (int j = 0; j < mesh.getVNoControlPoints(); j++) {
//        message << "\t\t";
//        for (int i = 0; i < mesh.getUNoControlPoints(); i++) {
//            message << "controlPoints[4*"<<(i+j* mesh.getUNoControlPoints())<<"+0]="<<mesh.getControlPointNet()[count]->getX() << ";" <<endl
//            		<< "controlPoints[4*"<<(i+j* mesh.getUNoControlPoints())<<"+1]="<<mesh.getControlPointNet()[count]->getY() << ";" <<endl
//            		<< "controlPoints[4*"<<(i+j* mesh.getUNoControlPoints())<<"+2]="<<mesh.getControlPointNet()[count]->getZ() << ";" <<endl
//            		<< "controlPoints[4*"<<(i+j* mesh.getUNoControlPoints())<<"+3]="<<mesh.getControlPointNet()[count]->getW() << ";" <<endl;
//            count++;
//        }
//        message << endl;
//    }
//    for(int j = 0; j < mesh.getVNoControlPoints()*mesh.getUNoControlPoints(); j++) {
//        message << "dofIndexNet["<<j<<"]="<<mesh.getControlPointNet()[j]->getDofIndex() << ";" <<endl;
//    }
    if(mesh.isTrimmed()) {
		message << mesh.getTrimming();
	}
    message << "\t" << "---------------------------------End Patch" << endl;
    return message;
}

}/* namespace EMPIRE */


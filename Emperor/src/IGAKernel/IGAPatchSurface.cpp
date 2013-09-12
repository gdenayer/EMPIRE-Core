// Inclusion of standard libraries
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// Inclusion of user defined libraries
#include "IGAPatchSurface.h"
#include "IGAMath.h"
#include "DataField.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

const int IGAPatchSurface::MAX_NUM_ITERATIONS = 20;

const int IGAPatchSurface::REGULAR_NUM_ITERATIONS = 6;

const double IGAPatchSurface::EPS_ORTHOGONALITY_CONDITION = 1e-9;

const double IGAPatchSurface::EPS_DISTANCE = 1e-9;

IGAPatchSurface::IGAPatchSurface(std::string _name, int _IDBasis, int _pDegree,
		int _uNoKnots, double* _uKnotVector, int _qDegree, int _vNoKnots,
		double* _vKnotVector, int _uNoControlPoints, int _vNoControlPoints,
		IGAControlPoint** _controlPointNet) :
		AbstractMesh(_name, EMPIRE_Mesh_IGAMesh), uNoControlPoints(
				_uNoControlPoints), vNoControlPoints(_vNoControlPoints) {

	// Read input
	bool ucondition = _uNoControlPoints != _uNoKnots - _pDegree - 1;
	bool vcondition = _vNoControlPoints != _vNoKnots - _qDegree - 1;

	if (ucondition || vcondition) {
		cout << endl;
		cout << endl;
		cout << "Error in IGAPatchSurface::IGAPatchSurface" << endl;
		cout
				<< "Number of Control Points, number of knots and polynomial degree do not match!"
				<< endl;
		cout << endl;
		cout << endl;
		exit(-1);
	}

	// Figure out whether the patch has a B-Spline or a NURBS underlying basis
	int isNurbs = 0;
	int counter = 0;
	for (int i = 0; i < uNoControlPoints; i++) {
		for (int j = 0; j < vNoControlPoints; j++) {
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
		IGABasis = new BSplineBasis2D(_IDBasis, _pDegree, _uNoKnots,
				_uKnotVector, _qDegree, _vNoKnots, _vKnotVector);
	} else {
		double* controlPointWeights = new double[uNoControlPoints
				* vNoControlPoints];
		for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
			controlPointWeights[i] = _controlPointNet[i]->getW();
		IGABasis = new NurbsBasis2D(_IDBasis, _pDegree, _uNoKnots, _uKnotVector,
				_qDegree, _vNoKnots, _vKnotVector, _uNoControlPoints,
				_vNoControlPoints, controlPointWeights);
	}

	// On the Control Point net
	assert(_controlPointNet!=NULL);
	ControlPointNet = _controlPointNet;
}

IGAPatchSurface::IGAPatchSurface(std::string _name, int _IDBasis, int _pDegree,
		int _uNoKnots, double* _uKnotVector, int _qDegree, int _vNoKnots,
		double* _vKnotVector, int _uNoControlPoints, int _vNoControlPoints,
		double* _controlPointNet) :
		AbstractMesh(_name, EMPIRE_Mesh_IGAMesh), uNoControlPoints(
				_uNoControlPoints), vNoControlPoints(_vNoControlPoints) {

	assert(_uKnotVector!=NULL);
	assert(_vKnotVector!=NULL);
	assert(_controlPointNet!=NULL);

	// Read input
	bool ucondition = _uNoControlPoints != _uNoKnots - _pDegree - 1;
	bool vcondition = _vNoControlPoints != _vNoKnots - _qDegree - 1;

	if (ucondition || vcondition) {
		cout << endl;
		cout << endl;
		cout << "Error in IGAPatchSurface::IGAPatchSurface" << endl;
		cout
				<< "Number of Control Points, number of knots and polynomial degree do not match!"
				<< endl;
		cout << endl;
		cout << endl;
		exit(-1);
	}

	ControlPointNet = new IGAControlPoint*[uNoControlPoints * vNoControlPoints];

	// Figure out whether the patch has a B-Spline or a NURBS underlying basis
	int isNurbs = 0;
	int counter = 0;
	for (int counter = 0; counter < uNoControlPoints * vNoControlPoints;
			counter++) {
		ControlPointNet[counter] = new IGAControlPoint(counter,
				_controlPointNet[counter * 4],
				_controlPointNet[counter * 4 + 1],
				_controlPointNet[counter * 4 + 2],
				_controlPointNet[counter * 4 + 3]);
		if (_controlPointNet[counter * 4 + 3] != 1.0)
			isNurbs = 1;
	}

	// Create the NURBS or the B-Spline underlying basis
	if (!isNurbs) {
		IGABasis = new BSplineBasis2D(_IDBasis, _pDegree, _uNoKnots,
				_uKnotVector, _qDegree, _vNoKnots, _vKnotVector);
	} else {
		double* controlPointWeights = new double[uNoControlPoints
				* vNoControlPoints];
		for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
			controlPointWeights[i] = ControlPointNet[i]->getW();
		IGABasis = new NurbsBasis2D(_IDBasis, _pDegree, _uNoKnots, _uKnotVector,
				_qDegree, _vNoKnots, _vKnotVector, _uNoControlPoints,
				_vNoControlPoints, controlPointWeights);
	}

}

IGAPatchSurface::~IGAPatchSurface() {
//	delete IGABasis;

}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates,
		double _uPrm, int _uKnotSpanIndex, double _vPrm, int _vKnotSpanIndex) {
	/*
	 *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
	 *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
	 */

	// Read input
	assert(_cartesianCoordinates!=NULL);

	// Initialize the coordinates of the point
	for (int i = 0; i < 3; i++)
		_cartesianCoordinates[i] = 0;

	// Compute the local basis functions in the vicinity of the point
	int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
	int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
	int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);
	double* localBasisFunctions = new double[noLocalBasisFunctions];

	IGABasis->computeLocalBasisFunctions(localBasisFunctions, _uPrm,
			_uKnotSpanIndex, _vPrm, _vKnotSpanIndex);

	// Initialize the Control Point index
	int CPindex = 0;

	// Initialize a basis functions counter
	int counter_basis = 0;

	// Loop over all the non-zero contributions
	for (int j = 0; j <= qDegree; j++) {
		for (int i = 0; i <= pDegree; i++) {

			// Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[i*m+j]
			CPindex = (_uKnotSpanIndex - pDegree + i) * vNoControlPoints
					+ (_vKnotSpanIndex - qDegree + j);

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
		double* _localBasisFunctions, int _uKnotSpanIndex,
		int _vKnotSpanIndex) {
	/*
	 *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
	 *  It is also expected that the local basis functions have been precomputed outside the scope of this function and are given as arguments.
	 *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
	 */

	// Read input
	assert(_cartesianCoordinates!=NULL);
	assert(_localBasisFunctions!=NULL);

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

			// Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[i*m+j]
			CPindex = (_uKnotSpanIndex - pDegree + i) * vNoControlPoints
					+ (_vKnotSpanIndex - qDegree + j);

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
		double* _localCoordinates) {
	int _uKnotSpanIndex = IGABasis->getUBSplineBasis1D()->findKnotSpan(
			_localCoordinates[0]);
	int _vKnotSpanIndex = IGABasis->getVBSplineBasis1D()->findKnotSpan(
			_localCoordinates[1]);
	IGAPatchSurface::computeCartesianCoordinates(_cartesianCoordinates,
			_localCoordinates[0], _uKnotSpanIndex, _localCoordinates[1],
			_vKnotSpanIndex);
}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates,
		double* _localBasisFctsAndDerivs, int _derivDegree, int _uKnotSpanIndex,
		int _vKnotSpanIndex) {
	/*
	 *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
	 *  It is also expected that the local basis functions have been precomputed outside the scope of this function and are given as arguments.
	 *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
	 */

	// Read input
	assert(_cartesianCoordinates!=NULL);
	assert(_localBasisFctsAndDerivs!=NULL);

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

			// Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[i * m + j]
			CPindex = (_uKnotSpanIndex - pDegree + i) * vNoControlPoints
					+ (_vKnotSpanIndex - qDegree + j);

			// Update the basis function index
			indexBasis = IGABasis->indexDerivativeBasisFunction(_derivDegree,
					derivIndex, derivIndex, counter_basis);

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

void IGAPatchSurface::computeBaseVectors(double* _baseVectors,
		double* _localBasisFunctionsAndDerivatives, int _uKnotSpanIndex,
		int _vKnotSpanIndex) {
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
	assert(_baseVectors!=NULL);
	assert(_localBasisFunctionsAndDerivatives!=NULL);

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
	int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

	// Initialize the Control Point index
	int CPindex = 0;

	// Initialize a basis functions counter
	int counterBasis = 0;

	// Initialize the derivative counter
	int counterBaseVector = 0;

	// Loop over all the non-zero contributions
	for (int j = 0; j <= qDegree; j++) {
		for (int i = 0; i <= pDegree; i++) {

			// Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[i*m+j]
			CPindex = (_uKnotSpanIndex - pDegree + i) * vNoControlPoints
					+ (_vKnotSpanIndex - qDegree + j);

			// Reset the base vector's counter
			counterBaseVector = 0;

			// Loop over all the first derivatives of the basis functions and sum up the contributions
			for (int l = 0; l <= derivDegree; l++) {
				for (int k = 0; k <= derivDegree - l; k++) {

					// Skip the case when k=l=0 i.e. the case where the basis functions themselves will be issued
					if (k == l)
						continue;

					// Get the index of the derivative functional value
					indexBasis = getIGABasis()->indexDerivativeBasisFunction(
							derivDegree, k, l, counterBasis);

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

int IGAPatchSurface::indexDerivativeBaseVector(int _derivDegree,
		int _uDerivIndex, int _vDerivIndex, int _componentIndex,
		int _baseVecIndex) {
	/*
	 * Returns the correct index when sorting the base vectors and the their derivatives in an 1D pointer array with the rule:
	 * (_derivDegree - _vDerivIndex) * (_derivDegree - _vDerivIndex + 1) * noCoordinates * noBaseVec / 2 +
	 *  + _uDerivIndex * noCoordinates * noBaseVec + _componentIndex * noBaseVec + _baseVecIndex
	 */

	// Read input
	if (_uDerivIndex + _vDerivIndex > _derivDegree) {
		cout << endl;
		cout << endl;
		cout << "Error in IGAPatchSurface::indexDerivativeBaseVector" << endl;
		cout << "It has been requested the " << _uDerivIndex
				<< "-th partial derivative w.r.t. u and" << endl;
		cout << "the " << _vDerivIndex
				<< "-th partial derivative w.r.t. v of the base vectors but "
				<< endl;
		cout << "the maximum absolute derivative selected is of "
				<< _derivDegree << "-th order" << endl;
		cout << endl;
		cout << endl;
	}

	// Number of tangent base vectors
	int noBaseVec = 2;

	// Number of coordinates
	int noCoordinates = 3;

	// Compute the index of the functional value
	return (_derivDegree - _vDerivIndex) * (_derivDegree - _vDerivIndex + 1)
			* noCoordinates * noBaseVec / 2
			+ _uDerivIndex * noCoordinates * noBaseVec
			+ _componentIndex * noBaseVec + _baseVecIndex;
}

void IGAPatchSurface::computeBaseVectorsAndDerivatives(
		double* _baseVectorsAndDerivatives,
		double*_localBasisFunctionsAndDerivatives, int _derivDegree,
		int _uKnotSpanIndex, int _vKnotSpanIndex) {
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
	assert(_baseVectorsAndDerivatives!=NULL);
	assert(_localBasisFunctionsAndDerivatives!=NULL);

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

			// Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[i*m+j]
			indexCP = (_uKnotSpanIndex - pDegree + uBasis) * vNoControlPoints
					+ (_vKnotSpanIndex - qDegree + vBasis);

			// Loop over all the derivatives in u-direction
			for (int i = 0; i <= _derivDegree; i++) {
				// Loop over all the derivatives in v-direction
				for (int j = 0; j <= _derivDegree - i; j++) {
					// Loop over all the coordinates
					for (int k = 0; k < noCoordinates; k++) {
						// Loop over all the base vectors
						for (int l = 0; l < noBaseVct; l++) {
							// Compute the index related to the derivatives of the base vectors
							indexBaseVct = indexDerivativeBaseVector(
									_derivDegree, i, j, k, l);

							// · Case l = 0 :: --> i + 1 - l = i + 1 && j + l --> j giving the derivative of the base vector g1 = dX/du
							// · Case l = 1 :: --> i + 1 - l = i && j + l --> j + 1 giving the derivative of the base vector g2 = dX/dv

							// Compute the index of the basis functions and their derivatives
							indexBasis = IGABasis->indexDerivativeBasisFunction(
									derivDegreeBasis, i + 1 - l, j + l,
									counterBasis);

							// Factor by which to multiply the derivative of the basis function
							if (k == 0)
								factor = ControlPointNet[indexCP]->getX();
							else if (k == 1)
								factor = ControlPointNet[indexCP]->getY();
							else
								factor = ControlPointNet[indexCP]->getZ();

							// Add the contribution from each basis function in the interval
							_baseVectorsAndDerivatives[indexBaseVct] +=
									_localBasisFunctionsAndDerivatives[indexBasis]
											* factor;
						}
					}
				}
			}
			// Update basis function's counter
			counterBasis++;
		}
	}
}

bool IGAPatchSurface::computePointProjectionOnPatch(double& _u, double& _v,	double* _P, bool& _flagConverge) {
	/*
	 * Returns the projection of a point _P on the NURBS patch given an initial guess for the surface parameters _u, _v via references:
	 * _P = double[3]
	 * Return value is a bool flag on the convergence of the Newton-Rapson iterations.
	 *
	 * Function layout :
	 *
	 * 1. Read input and initialize the data
	 *
	 * 2. Loop over all the Newton-Rapson iterations
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
	assert(_P!=NULL);

//	cout << "initial u: " << _u << endl;
//	cout << "initial v: " << _v << endl;
//	cout << "Cartesian Coordinate: " << _P[0] << ", " << _P[1] << ", " << _P[2]
//			<< endl;

	// Initialize the flag to true
	bool flagNewtonRapson = true;
	double du,dv;

	// Initialize number of spatial dimensions
	int noSpatialDimensions = 3;

	// Save the location of the point to be projected
	double point[noSpatialDimensions];
	for (int i = 0; i < noSpatialDimensions; i++)
		point[i] = _P[i];

	// Initialize the distance vector
	double distanceVector[3];

	// Initialize the indices of the base vectors and their derivatives
	int indexUBaseVec = 0;
	int indexVBaseVec = 0;
	int indexdUBaseVecdu = 0;
	int indexdVBaseVecdv = 0;
	int indexdBaseVecMix = 0;

	// Initialize the base vectors and their derivatives
	double uBaseVec[3];
	double vBaseVec[3];
	double uBaseVecdu[3];
	double vBaseVecdv[3];
	double dBaseVecMix[3];

	// Initialize the dot products
	double distanceVector2norm = 0.0;
	double uBaseVecXdistanceVector = 0.0;
	double squareUBaseVec2norm = 0.0;
	double uBaseVec2norm = 0.0;
	double vBaseVecXdistanceVector = 0.0;
	double squareVBaseVec2norm = 0.0;
	double vBaseVec2norm = 0.0;

	// Initialize Jacobian matrix
	double jacobianMatrix[4];

	// Initialize right-hand side solution vector
	double rightSideVct[2];

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
	int lengthUKnotVct = IGABasis->getUBSplineBasis1D()->getNoKnots();
	int lengthVKnotVct = IGABasis->getVBSplineBasis1D()->getNoKnots();

	// Local number of basis functions
	int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

	// Initialize the Newton-Rapson iteration counter
	int counter = 0;

	// Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Rapson iterations)
	int derivDegreeBasis = 2;

	// The number of the base vectors
	int noBaseVcts = 2;

	// Initialize the array of the IGA basis functions and their derivatives
	double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1)
			* (derivDegreeBasis + 2) * noLocalBasisFcts / 2];

	// Number of derivatives for the base vectors
	int derivDegreeBaseVcts = derivDegreeBasis - 1;

	// Initialize the array of the base vectors and their derivatives
	double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1)
			* (derivDegreeBaseVcts + 2) * noSpatialDimensions * noBaseVcts / 2];

	// 2. Loop over all the Newton-Rapson iterations
	while (counter <= MAX_NUM_ITERATIONS) {

		// 2i. Update the iteration counter
		counter++;

		// 2ii. Find the span of the given surface parameters
		uKnotSpan = IGABasis->getUBSplineBasis1D()->findKnotSpan(_u);
		vKnotSpan = IGABasis->getVBSplineBasis1D()->findKnotSpan(_v);

		// 2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
		IGABasis->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs,
				derivDegreeBasis, _u, uKnotSpan, _v, vKnotSpan);

		// 2iv. Compute the Cartesian components of the point on the surface
		computeCartesianCoordinates(_P, basisFctsAndDerivs, derivDegreeBasis,
				uKnotSpan, vKnotSpan);

		// 2v. Compute the distance vector between the vector to be projected and the estimated one
		for (int i = 0; i < noSpatialDimensions; i++)
			distanceVector[i] = _P[i] - point[i];

		// 2vi. Compute the 2-norm of the distance vector
		distanceVector2norm = square2normVector(noSpatialDimensions,
				distanceVector);
		distanceVector2norm = sqrt(distanceVector2norm);

		if (distanceVector2norm < EPS_DISTANCE)
			break;

		// 2vii. Compute the base vectors and their derivatives
		computeBaseVectorsAndDerivatives(baseVecAndDerivs, basisFctsAndDerivs,
				derivDegreeBaseVcts, uKnotSpan, vKnotSpan);

		// 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
		for (int i = 0; i < noSpatialDimensions; i++) {
			// On the base vector Gu = dR/du
			indexUBaseVec = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0,
					i, 0);
			uBaseVec[i] = baseVecAndDerivs[indexUBaseVec];

			// On the base vector Gv = dR/dv
			indexVBaseVec = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0,
					i, 1);
			vBaseVec[i] = baseVecAndDerivs[indexVBaseVec];

			// On the derivative of the base vector Gu w.r.t. u namely dGu/du = d^2R/du^2
			indexdUBaseVecdu = indexDerivativeBaseVector(derivDegreeBaseVcts, 1,
					0, i, 0);
			uBaseVecdu[i] = baseVecAndDerivs[indexdUBaseVecdu];

			// On the derivative of the base vector Gv w.r.t. u namely dGv/dv = d^2R/dv^2
			indexdVBaseVecdv = indexDerivativeBaseVector(derivDegreeBaseVcts, 0,
					1, i, 1);
			vBaseVecdv[i] = baseVecAndDerivs[indexdVBaseVecdv];

			// On the mixed derivative of the base vectors namely d^2Gu/dv = d^2Gv/du = d^2R/dudv
			indexdBaseVecMix = indexDerivativeBaseVector(derivDegreeBaseVcts, 0,
					1, i, 0);
			dBaseVecMix[i] = baseVecAndDerivs[indexdBaseVecMix];
		}

		// 2ix. Compute the cosine of the angle with respect to u-parametric line
		uBaseVecXdistanceVector = dotProduct(noSpatialDimensions, uBaseVec,
				distanceVector);
		squareUBaseVec2norm = square2normVector(noSpatialDimensions, uBaseVec);
		uBaseVec2norm = sqrt(squareUBaseVec2norm);
		cosu = fabs(uBaseVecXdistanceVector) / uBaseVec2norm
				/ distanceVector2norm;

		// 2x. Compute the cosine of the angle with respect to v-parametric line
		vBaseVecXdistanceVector = dotProduct(noSpatialDimensions, vBaseVec,
				distanceVector);
		squareVBaseVec2norm = square2normVector(noSpatialDimensions, vBaseVec);
		vBaseVec2norm = sqrt(squareVBaseVec2norm);
		cosv = fabs(vBaseVecXdistanceVector) / vBaseVec2norm
				/ distanceVector2norm;

		// 2xi. Check the orthogonality condition and if it is fulfilled break the loop

		if (cosu <= EPS_ORTHOGONALITY_CONDITION
				&& cosv <= EPS_ORTHOGONALITY_CONDITION)
			break;


		// 2xii. Compute the entries of the Jacobian matrix
		jacobianMatrix[0] = squareUBaseVec2norm
				+ dotProduct(noSpatialDimensions, uBaseVecdu, distanceVector);
		jacobianMatrix[1] = dotProduct(noSpatialDimensions, uBaseVec, vBaseVec)
				+ dotProduct(noSpatialDimensions, dBaseVecMix, distanceVector);
		jacobianMatrix[2] = jacobianMatrix[1];
		jacobianMatrix[3] = squareVBaseVec2norm
				+ dotProduct(noSpatialDimensions, vBaseVecdv, distanceVector);

		// 2xiii. Compute the entries of the right-hand side vector
		rightSideVct[0] = -dotProduct(noSpatialDimensions, uBaseVec,
				distanceVector);
		rightSideVct[1] = -dotProduct(noSpatialDimensions, vBaseVec,
				distanceVector);

		// 2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved


		// Solve the equation system
		flagLinearSystem = solve2x2linearSystem(rightSideVct, jacobianMatrix);


		// Check if the equation system has been successfully solved
		if (!flagLinearSystem) {
			cout << endl;
			cout << endl;
			cout << "Error in IGAPatchSurface::computePointProjectionOnPatch"
					<< endl;
			cout
					<< "The 2x2 equation system to find the updates of the surface parameters"
					<< endl;
			cout
					<< "for the orthogonal projection of a point on the NURBS patch has been"
					<< endl;
			cout << "detected not solvable up to tolerance" << EPS << endl;
			cout << endl;
			cout << endl;
			exit(-1);
		}

		du = _u;
		dv = _v;


		// 2xv. Update the surface parameters u += du and v += dv
		_u += rightSideVct[0];
		_v += rightSideVct[1];

		// 2xvi. Check and modify the surface parameters if they stay out of their knot spans
		if (_u > IGABasis->getUBSplineBasis1D()->getKnotVector()[lengthUKnotVct - 1])
			_u = IGABasis->getUBSplineBasis1D()->getKnotVector()[lengthUKnotVct - 1];
		if (_u < IGABasis->getUBSplineBasis1D()->getKnotVector()[0])
			_u = IGABasis->getUBSplineBasis1D()->getKnotVector()[0];
		if (_v > IGABasis->getVBSplineBasis1D()->getKnotVector()[lengthVKnotVct	- 1])
			_v = IGABasis->getVBSplineBasis1D()->getKnotVector()[lengthVKnotVct	- 1];
		if (_v < IGABasis->getVBSplineBasis1D()->getKnotVector()[0])
			_v = IGABasis->getVBSplineBasis1D()->getKnotVector()[0];

		du -= _u;
		dv -= _v;
	}


	// 3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
	if (counter > MAX_NUM_ITERATIONS) {
		flagNewtonRapson = 0;
		if (du * du + dv * dv < 1e-10){
			cout << "Warning in function IGAPatchSurface::computePointProjectionOnPatch " << endl;
			cout << "Newton-Rapson iterations converged to a non-projection point" << endl;
			cout << "(U, V) : (" << _u << " , " << _v << ")" << endl;

		} else {
			cout << endl;
			cout << "Error in function IGAPatchSurface::computePointProjectionOnPatch" << endl;
			cout << "Newton-Rapson iterations did not converge up to iteration number:  "<< MAX_NUM_ITERATIONS << endl;
			cout << "Point: (" << _P[0] << "," << _P[1] << "," << _P[2] << ")" << endl;
			cout << endl;
//			exit(-1);
		}
	} else if (counter > REGULAR_NUM_ITERATIONS) {
		cout << endl;
		cout << "Warning in function IGAPatchSurface::computePointProjectionOnPatch" << endl;
		cout << "Number of iterations for the Newton-Rapson algorithm to converge:  " << counter - 1 << endl;
		cout << endl;
	}

	// 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
	// Clear the memory on the heap
	delete[] basisFctsAndDerivs;
	delete[] baseVecAndDerivs;


	// Return the flag
	return flagNewtonRapson;
}

bool IGAPatchSurface::computePointProjectionOnPatch(double& _u, double& _v,	double* _P){
	bool tmp = false;
	return computePointProjectionOnPatch(_u,_v,_P,tmp);
}

void IGAPatchSurface::findNearestKnotIntersection(double& _u, double& _v,
		double* _P) {

	assert(_P != NULL);

//	IGABasis->getUBSplineBasis1D()->printKnotVector();
//	IGABasis->getVBSplineBasis1D()->printKnotVector();

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

	for (int i = pU; i < nU; i++)
		for (int j = pV; j < nV; j++) {
			computeCartesianCoordinates(coords, knotU[i], i, knotV[j], j);

			for (int k = 0; k < noSpatialDimensions; k++)
				coords[k] -= _P[k];

			Dis = square2normVector(noSpatialDimensions, coords);
			Dis = sqrt(Dis);

			if (Dis < minDis) {

				minDis = Dis;
				_u = knotU[i];
				_v = knotV[j];
			}

		}

}

void IGAPatchSurface::printControlPointNet() {

	cout << endl;
	cout << "---------------------------------------------" << endl;
	cout << "Debugging information in class IGAPatchSurface" << endl;
	cout << endl;
	cout << "IGAPatchSurface::uNoControlPoints = " << getUNoControlPoints()
			<< endl;
	cout << "IGAPatchSurface::vNoControlPoints = " << getVNoControlPoints()
			<< endl;
	cout << endl;
	for (int i = 0; i < getUNoControlPoints() * getVNoControlPoints(); i++)
		cout << "CP " << ControlPointNet[i]->getId() << ": ( "
				<< ControlPointNet[i]->getX() << " , "
				<< ControlPointNet[i]->getY() << " , "
				<< ControlPointNet[i]->getZ() << " ) --Weight--> "
				<< ControlPointNet[i]->getW() << endl;
	;
	cout << "_____________________________________________" << endl;
	cout << endl;
}

void IGAPatchSurface::computeBoundingBox() {
	if (boundingBox.isComputed)
		return;
	boundingBox.xmin = ControlPointNet[0]->getX();
	boundingBox.xmin = ControlPointNet[0]->getX();
	boundingBox.xmin = ControlPointNet[0]->getY();
	boundingBox.xmin = ControlPointNet[0]->getY();
	boundingBox.xmin = ControlPointNet[0]->getZ();
	boundingBox.xmin = ControlPointNet[0]->getZ();
	for (int i = 1; i < getNoControlPoints(); i++) {
		double x = ControlPointNet[i]->getX();
		double y = ControlPointNet[i]->getY();
		double z = ControlPointNet[i]->getZ();
		if (x < boundingBox.xmin)
			boundingBox.xmin = x;
		else if (x > boundingBox.xmax)
			boundingBox.xmax = x;
		if (y < boundingBox.ymin)
			boundingBox.ymin = y;
		else if (y > boundingBox.ymax)
			boundingBox.ymax = y;
		if (z < boundingBox.zmin)
			boundingBox.zmin = z;
		else if (z > boundingBox.zmax)
			boundingBox.zmax = z;
	}
	boundingBox.isComputed = true;
}

void IGAPatchSurface::addDataField(string dataFieldName,
		EMPIRE_DataField_location location,
		EMPIRE_DataField_dimension dimension,
		EMPIRE_DataField_typeOfQuantity typeOfQuantity) {

	int numLocatiions = -1;
	if (location == EMPIRE_DataField_atNode)
		numLocatiions = getNoControlPoints();
	else
		assert(false);

	assert(nameToDataFieldMap.find(dataFieldName) == nameToDataFieldMap.end());
	DataField *dataField = new DataField(dataFieldName, location, numLocatiions,
			dimension, typeOfQuantity);
	nameToDataFieldMap.insert(
			pair<string, DataField*>(dataFieldName, dataField));

}

void IGAPatchSurface::printSelf(){

	cout << "\t\tpDegree:  "
			<< IGABasis->getUBSplineBasis1D()->getPolynomialDegree()
			<< endl;
	cout << "\t\tqDegree:  "
			<< IGABasis->getVBSplineBasis1D()->getPolynomialDegree()
			<< endl;

	cout << "\t\tKnots Vector U: \t";
	for (int i = 0; i < IGABasis->getUBSplineBasis1D()->getNoKnots();
			i++)
		cout << IGABasis->getUBSplineBasis1D()->getKnotVector()[i]
				<< "  ";
	cout << endl;

	cout << "\t\tKnots Vector V: \t";
	for (int i = 0; i < IGABasis->getVBSplineBasis1D()->getNoKnots();
			i++)
		cout << IGABasis->getVBSplineBasis1D()->getKnotVector()[i]
				<< "  ";
	cout << endl;

	cout << "\t\t" << "number of control points U: "
			<< uNoControlPoints << endl;
	cout << "\t\t" << "number of control points V: "
			<< vNoControlPoints << endl;

	cout << "\t\tControl Points Net: " << endl;
	int count = 0;
	for (int i = 0; i < uNoControlPoints; i++) {
		cout << "\t\t";
		for (int j = 0; j < vNoControlPoints; j++) {
			cout << ControlPointNet[count] -> getId() << ":"
					<< ControlPointNet[count]->getX() << ", "
					<< ControlPointNet[count]->getY() << ", "
					<< ControlPointNet[count]->getZ() << "  "
					<<  ControlPointNet[count]->getW() << "\t";
			count++;
		}
		cout << endl;
	}

	cout << "\t" << "---------------------------------" << endl;
}

Message &operator<<(Message &message, IGAPatchSurface &mesh) {
//	message << "\t" << "IGA Patch name: " << mesh.name << endl;

	message << "\t\tpDegree:  "
			<< mesh.getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
			<< endl;
	message << "\t\tqDegree:  "
			<< mesh.getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
			<< endl;

	message << "\t\tKnots Vector U: \t";
	for (int i = 0; i < mesh.getIGABasis()->getUBSplineBasis1D()->getNoKnots();
			i++)
		message << mesh.getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i]
				<< "  ";
	message << endl;

	message << "\t\tKnots Vector V: \t";
	for (int i = 0; i < mesh.getIGABasis()->getVBSplineBasis1D()->getNoKnots();
			i++)
		message << mesh.getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i]
				<< "  ";
	message << endl;

	message << "\t\t" << "number of control points U: "
			<< mesh.getUNoControlPoints() << endl;
	message << "\t\t" << "number of control points V: "
			<< mesh.getVNoControlPoints() << endl;

	message << "\t\tControl Points Net: " << endl;
	int count = 0;
	for (int i = 0; i < mesh.getUNoControlPoints(); i++) {
		cout << "\t\t";
		for (int j = 0; j < mesh.getVNoControlPoints(); j++) {
			message << mesh.getControlPointNet()[count]->getX() << ", "
					<< mesh.getControlPointNet()[count]->getY() << ", "
					<< mesh.getControlPointNet()[count]->getZ() << "\t";
			count++;
		}
		message << endl;
	}

	message() << "\t" << "---------------------------------" << endl;
	return message;
}

}/* namespace EMPIRE */


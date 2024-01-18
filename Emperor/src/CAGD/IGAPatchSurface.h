/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Ragnar Björnsson Munich
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
/***********************************************************************************************//**
 * \file IGAPatchSurface.h
 * This file holds the class IGAPatchSurface.h
 * \date 28/5/2013
 **************************************************************************************************/

#ifndef IGAPatchSurface_H_
#define IGAPatchSurface_H_

// Inclusion of user defined libraries
#include "AbstractMesh.h"
#include "BoundingBox.h"
#include "NurbsBasis2D.h"
#include "IGAControlPoint.h"
#include "IGAPatchSurfaceTrimming.h"
#include <limits>
#include <set>
#include <vector>

namespace EMPIRE {
class DataField;
class Message;

/********//**
 * \brief class IGAPatchSurface is a specialization of the class AbstractMesh used for IGA Mesh with two parameters like shell elements
 ***********/

class IGAPatchSurface {

protected:
    /// The basis functions of the 2D NURBS patch
    BSplineBasis2D* IGABasis;

    /// Number of Control Points in u-direction
    int uNoControlPoints;

    /// Number of Control Points in v-direction
    int vNoControlPoints;

    /// The set of the Control Points of the patch
    IGAControlPoint** ControlPointNet;

    /// The class holding the trimming information
    IGAPatchSurfaceTrimming Trimming;

    /// The bounding box of the patch
    AABB boundingBox;

    /// The constructor and the destructor and the copy constructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _IDBasis The id of the underlying basis to the IGA 2D patch
     * \param[in] _pDegree The polynomial degree of the IGA 2D patch in the u-direction
     * \param[in] _uNoKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 2D patch in the u-direction
     * \param[in] _qDegree The polynomial degree of the IGA 2D patch in the v-direction
     * \param[in] _vNoKnots The number of knots for the knot vector in the v-direction
     * \param[in] _vKnotVector The underlying knot vector of the IGA 2D patch in the v-direction
     * \param[in] _uNoControlPoints The number of the Control Points for the 2D NURBS patch in the u-direction
     * \param[in] _vNoControlPoints The number of the Control Points for the 2D NURBS patch in the v-direction
     * \param[in] _controlPointNet The set of the Control Points related to the 2D NURBS patch
     * \author Andreas Apostolatos
     ***********/
    IGAPatchSurface(int, int, int, double*, int, int, double*, int, int, IGAControlPoint**);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    ~IGAPatchSurface();

    void computeBoundingBox();
    /// Trimming related functions
public:
    /***********************************************************************************************
     * \brief Setup information about the loop soon to be received
     * \param[in] inner 0 for outter and 1 for inner
     * \param[in] numCurves Number of curves to be received for this loop
     * \author Fabien Pean
     ***********/
    void addTrimLoop(int inner, int numCurves);
    /***********************************************************************************************
     * \brief Add a Nurbs curve for the current loop and its attached information
     * \param[in] _direction The direction of the curve if is following standard or not
     * \param[in] _pDegree The polynomial degree of the IGA 1D curve in the u-direction
     * \param[in] _uNoKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 1D curve in the u-direction
     * \param[in] _uNoControlPoints The number of the Control Points for the 1D NURBS patch in the u-direction
     * \param[in] _controlPointNet The set of the Control Points related to the 1D NURBS patch
     * \author Fabien Pean
     ***********/
    void addTrimCurve(int _direction, int _pDegree, int _uNoKnots, double* _uKnotVector,
                      int _uNoControlPoints, double* _controlPointNet);

    /***********************************************************************************************
     * \brief Linearize all the trimming loops and curves of the given mesh and patch
     * \author Fabien Pean
     ***********/
    inline void linearizeTrimming() {
        Trimming.linearizeLoops();
    }

    /// Basis related functions
public:
    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters are known
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the Patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _uPrm The parameter on the u-coordinate line
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vPrm The parameter on the v-coordinate line
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeCartesianCoordinates(double*, double, int, double, int)const;

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters and the local basis functions are given
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _localBasisFunctions The local basis functions
     * \compute the knot span Index inside the function. Convenient but in-efficient.
     * \author Chenshen Wu
     ***********/
    void computeCartesianCoordinates(double*, double*)const;

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters and the local basis functions are given
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _localBasisFunctions The local basis functions
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeCartesianCoordinates(double*, double*, int, int)const;

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters and the local basis functions and their derivatives are given
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _localBasisFctsAndDerivs The local basis functions and their derivatives
     * \param[in] _derivDegree The derivative degree up to which the basis functions have been computed contained ion array _localBasisFctsAndDerivs
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeCartesianCoordinates(double*, double*, int, int, int)const;

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a trimming curve that lies on the patch
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the patch whose curve parameters are _uPrm
     * \param[in] _localBasisFctsAndDerivs The local basis functions and their derivatives
     * \param[in] _uPrm The curve parameter
     * \param[in] _patchBLIndex The boundary loop index in which the trimming curve exists
     * \param[in] _patchBLTrCurveIndex The index of trimming curve in the boundary loop
     * \author Altug Emiroglu
     ***********/
    void computeCartesianCoordinates(double* _cartesianCoordinates, double _uPrm, int _patchBLIndex, int _patchBLTrCurveIndex) const;

    void computeCartesianCoordinates(double* _cartesianCoordinates, double _uPrm, IGAPatchCurve* _curve) const;

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of the base vectors at a given pair of surface parameters given the basis functions and their derivatives
     * \param[in/out] _baseVectors The Cartesian coordinates of the base vectors on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _localBasisFunctionsAndDerivatives The local basis functions and their derivatives
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeBaseVectors(double* _baseVectors, double* _localBasisFunctionsAndDerivatives,
            int _uKnotSpanIndex, int _vKnotSpanIndex)const;

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of the base vectors at a given pair of surface parameters given the basis functions and their derivatives
     * \param[in/out] _baseVectors The Cartesian coordinates of the base vectors on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _uPrm The parameter on the u-coordinate line
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vPrm The parameter on the v-coordinate line
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Chenshen Wu
     ***********/
    void computeBaseVectors(double* _baseVectors, double _uPrm, int _uKnotSpanIndex, double _vPrm,
            int _vKnotSpanIndex)const;

    /***********************************************************************************************
     * \brief Returns the index of the i-th partial derivative w.r.t. to u , j-th partial derivative w.r.t v of the _componentIndex-th component to the l-th base vector
     * \param[out] The index of the _uDerivIndex-th partial derivative w.r.t. to u , j-th partial derivative w.r.t v of the k-th component to the l-th base vector
     * \param[in] _derivDegree The absolute order of the partial derivatives to the base vectors
     * \param[in] _uDerivIndex The order of the partial derivative w.r.t. u-parametric coordinate, _uDerivIndex = 0 , … , _derivDegree
     * \param[in] _vDerivIndex The order of the partial derivative w.r.t. v-parametric coordinate, _uDerivIndex = 0 , … , _derivDegree - _uDerivIndex
     * \param[in] _componentIndex The component of the base vector, _componentIndex = 1, 2, 3
     * \param[in] _baseVecIndex The base vector for which to compute the derivatives _baseVecIndex = 1, 2
     * \author Andreas Apostolatos
     ***********/
    int indexDerivativeBaseVector(int, int, int, int, int);

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of the base vectors and their derivatives at a given pair of surface parameters given the basis functions and their derivatives
     * \param[in/out] _baseVectorsAndDerivatives The Cartesian coordinates of the point on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _localBasisFunctionsAndDerivatives The local basis functions and their derivatives
     * \param[in] _derivDegree The derivative order with respect to both u-,v- directions
     * \param[in] _uPrm The parameter on the u-coordinate line
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vPrm The parameter on the v-coordinate line
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeBaseVectorsAndDerivatives(double*, double*, int, int, int);

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of the surface normal vector and its first derivatives at a given pair of surface parameters given the surface base vectors
     * \param[in/out] _surfNormalVctAndDervs The surface normal vector and its derivatives, should be static array of size 9
     * \param[in] _baseVctsAndDerivs The base vectors and their first order parametric derivatives
     * \param[in] _derivDegreeBaseVec The derivative order up to which the base vectors have been computed
     * \author Andreas Apostolatos
     ***********/
    void computeSurfaceNormalVectorAndDerivatives(double*, double*, int);

    /***********************************************************************************************
     * \brief Returns the coefficients of the covariant metric tensor
     * \param[in/out] _covariantMetricTensor The coefficients of the covariant metric tensor in an array of constant size 4
     * \param[in] _baseVctsAndDerivs The base vectors and their first order parametric derivatives
     * \param[in] _derivDegreeBaseVec The derivative order up to which the base vectors have been computed
     * \author Andreas Apostolatos
     ***********/
    void computeCovariantMetricTensor(double*, double*, int);

    /***********************************************************************************************
     * \brief Returns the contravariant base vectors of the surface
     * \param[in/out] _contravariantBaseVcts The Cartesian components of the contravariant base vectors stored in an array of constant size 6
     * \param[in] _covariantMetricTensor The coefficients of the covariant metric tensor stored in an array of constant sizes 4
     * \param[in] _baseVctsAndDerivs The base vectors and their first order parametric derivatives
     * \param[in] _derivDegreeBaseVec The derivative order up to which the base vectors have been computed
     * \author Andreas Apostolatos
     ***********/
    void computeContravariantBaseVectors(double*, double*, double*, int);

    /***********************************************************************************************
     * \brief Returns the coefficients of the curvature tensor in the contravariant basis
     * \param[in/out] _contravariantCurvatureTensor The components of the curvature tensor in the contravariant basis
     * \param[in] _surfNormalVctAndDervs
     * \param[in] _baseVctsAndDerivs The base vectors and their first order parametric derivatives
     * \param[in] _derivDegreeBaseVec The derivative order up to which the base vectors have been computed
     * \author Andreas Apostolatos
     ***********/
    void computeContravariantCurvatureTensor(double*, double*, double*, int);

    /// Projection related functions
public:
    /***********************************************************************************************
     * \brief Computes the orthogonal projection of point of the 3D Euclidean space onto the NURBS pacth
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _P Given the Cartesian components of the point to be projected on the NURBS patch it is returned the Cartesian components of its orthogonal projection
     * \param[in/out] _flagConverge Flag indicating whether the Newton iterations have converged true/false
     * \param[in]	  _maxIt The number of iteration to do in the scheme
     * \param[in]	  _tol The tolerance for which the scheme stops
     * \return The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \author Andreas Apostolatos
     ***********/
    bool computePointProjectionOnPatch(double&, double&, double*, bool&, const int _maxIt=MAX_NUM_ITERATIONS, const double _tol=TOL_ORTHOGONALITY, const double _distTol=TOL_DISTANCE);

    /***********************************************************************************************
     * \brief Computes the orthogonal projection of point of the 3D Euclidean space onto the NURBS pacth (overloaded)
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _P Given the Cartesian components of the point to be projected on the NURBS patch it is returned the Cartesian components of its orthogonal projection
     * \param[in]	  _maxIt The number of iteration to do in the scheme
     * \param[in]	  _tol The tolerance for which the scheme stops
     * \return The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \author Andreas Apostolatos
     ***********/
    bool computePointProjectionOnPatch(double&, double&, double*, const int _maxIt=MAX_NUM_ITERATIONS, const double _tol=TOL_ORTHOGONALITY, const double _distTol=TOL_DISTANCE);

    /***********************************************************************************************
     * \brief Computes the forced orthogonal projection of point of the 3D Euclidean space onto the NURBS pacth
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _P Given the Cartesian components of the point to be projected on the NURBS patch it is returned the Cartesian components of its orthogonal projection
     * \param[in]	  _relMaxIt The relaxed number of iteration to do in the scheme
     * \param[in]	  _relTol The relaxed tolerance for which the scheme stops
     * \return The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \author Altug Emiroglu
     ***********/
    bool computeForcedPointProjectionOnPatch(double& _u, double& _v, double* _P, const int _relMaxIt=REL_MAX_NUM_ITERATIONS, const double _relTol=REL_TOL_ORTHOGONALITY, const double _distTol=REL_TOL_DISTANCE);

    /***********************************************************************************************
     * \brief Computes the orthogonal projection of point of the 3D Euclidean space onto the trimming curve
     * \param[in/out] _projectedUTilde Curve parameter of the projected point on the trimming curve
     * \param[in]     _coordsXYZ Cartesian coordinates of the points to be projected
     * \param[in]	  _patchBLIndex The index of the boundary loop which contains the trimming curve
     * \param[in]	  _patchBLTrCurveIndex The index of the trimming curve in the boundary loop
     * \return The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    bool computePointProjectionOnTrimmingCurve(double& _projectedUTilde, double* _P, int _patchBLIndex, int _patchBLTrCurveIndex);

    /***********************************************************************************************
     * \brief Computes the orthogonal projection of point of the 3D Euclidean space onto the trimming curve
     * \param[in/out] _projectedUTilde Curve parameter of the projected point on the trimming curve
     * \param[in]     _coordsXYZ Cartesian coordinates of the points to be projected
     * \param[in]	  _curve The curve on which the point should be projected
     * \return The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    bool computePointProjectionOnTrimmingCurve(double& _projectedUTilde, double* _P, IGAPatchCurve* _curve);

    /***********************************************************************************************
     * \brief Solves the Newton-Raphson problem for projecting a line segment on a patch boundary by finding the location where P1Q vector lies in the plane (P1P2 - n)
     * \param[in/out] _t The running parameter on the given NURBS patch boundary
     * \param[in/out] _lambda The ratio between the line segment that is projected on the NURBS patch to the complete line segment
     * \param[in/out] _distance The orthogonal distance from the NURBS surface to the line segment
     * \param[in/out] _distanceActual The actual ortogonal distance (only used for validity checks)
     * \param[in] _P1 The first point of the line segment
     * \param[in] _P2 The second point of the line segment
     * \param[in] _edge (0,1,2,3) --> (uRunsvStart,uRunsvEnd,uStartvRuns,uEndvRuns)
     * \param[in] _maxIt The number of iteration to do in the scheme
     * \param[in] _tol The tolerance for which the scheme stops
     * \return The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \author Fabien Pean
     ***********/
    bool solvePointProjectionOnPatchBoundaryNewtonRaphson(double& _t, double& _lambda,
            double& _distance, double* _P1, double* _P2, int _edge,
            const int _maxIt=MAX_NUM_ITERATIONS, const double _tol=TOL_ORTHOGONALITY);

    /***********************************************************************************************
     * \brief Solves the Newton-Raphson problem for computing the closest point projection of a straight line over a patch boundary by orthogonality
     * \param[in/out] _t The running parameter on the given NURBS patch boundary
     * \param[in/out] _lambda The ratio between the line segment that is projected on the NURBS patch to the complete line segment
     * \param[in/out] _distance The orthogonal distance from the NURBS surface to the line segment
     * \param[in/out] _distanceActual The actual ortogonal distance (only used for validity checks)
     * \param[in] _P1 The first point of the line segment
     * \param[in] _P2 The second point of the line segment
     * \param[in] _edge (0,1,2,3) --> (uRunsvStart,uRunsvEnd,uStartvRuns,uEndvRuns)
     * \param[in] _maxIt The number of iteration to do in the scheme
     * \param[in] _tol The tolerance for which the scheme stops
     * \return The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    bool solvePointProjectionOnPatchBoundaryNewtonRaphsonClosestDistance(double& _t, double& _lambda,
            double & _distance, double* _P1, double* _P2, int _indexEdge,
            const int _maxIt=MAX_NUM_ITERATIONS, double _tol=TOL_ORTHOGONALITY);

    /***********************************************************************************************
     * \brief Solves the bisection problem for computing the closest point projection of a straight line over a patch boundary
     * \param[in/out] _t The running parameter on the given NURBS patch boundary
     * \param[in/out] _lambda The ratio between the line segment that is projected on the NURBS patch to the complete line segment
     * \param[in/out] _distance The orthogonal distance from the NURBS surface to the line segment
     * \param[in/out] _distanceActual The actual distance of the point found on the patch and the finite element edge
     * \param[in] _P1 The first point of the line segment
     * \param[in] _P2 The second point of the line segment
     * \param[in] _maxIt The number of iteration to do in the scheme
     * \param[in] _tol The tolerance for which the scheme stops
     * \return The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \author Fabien Pean
     ***********/
    bool solvePointProjectionOnPatchBoundaryBisection(double& _u,double& _v, double& _lambda,
            double& _distance, double* _P1, double* _P2,
            const int _maxIt=MAX_NUM_ITERATIONS, const double _tol=TOL_ORTHOGONALITY);

    /***********************************************************************************************
     * \brief Returns the point on the given NURBS patch boundary which defines an orthogonal projection from the given line to the NURBS boundary
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _lambda The ratio between the line segment that is projected on the NURBS patch to the complete line segment
     * \param[in/out] _distance The orthogonal distance from the NURBS surface to the line segment
     * \param[in/out] _distActual The actual orthogonal distance
     * \param[in] _P1 The first point of the line segment
     * \param[in] _P2 The second point of the line segment
     * \param[in] _maxIt The number of iteration to do in the scheme
     * \param[in] _tol The tolerance for which the scheme stops
     * \return The id of the edge it is crossing, indicating also if scheme has converged (no edge = no convergence)
     * \author Fabien Pean
     ***********/
    char computePointProjectionOnPatchBoundaryNewtonRhapson(double& _u, double& _v, double& _lambda,
            double& _distance, double* _P1, double* _P2,
            const int _maxIt=MAX_NUM_ITERATIONS, const double _tol=TOL_ORTHOGONALITY);

    /***********************************************************************************************
     * \brief Returns the point on the given NURBS patch boundary which defines an orthogonal projection from the given line to the NURBS boundary using the bisection method
     * \param[out] The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _ratio The ratio between the line segment that is projected on the NURBS patch to the complete line segment
     * \param[in/out] _distance The orthogonal distance from the NURBS surface to the line segment
     * \param[in/out] _distanceActual The actual distance of the found point on the patch and the finite element edge
     * \param[in] _P1 The first point of the line segment
     * \param[in] _P2 The second point of the line segment
     * \param[in] _maxIt The number of iteration to do in the scheme
     * \param[in] _tol The tolerance for which the scheme stops
     * \return The id of the edge it is crossing, indicating also if scheme has converged (no edge = no convergence)
     * \author Fabien Pean
     ***********/
    char computePointProjectionOnPatchBoundaryBisection(double& _u, double& _v, double& _ratio,
            double& _distance, double* _P1, double* _P2,
            const int _maxIt=MAX_NUM_ITERATIONS, const double _tol=TOL_ORTHOGONALITY);

    /***********************************************************************************************
     * \brief Find the nearest knot intersection on the patch as an initial guess for the projection
     * \param[in/out] _u Given is the u-surface parameter of the nearest knot intersection.
     * \param[in/out] _v Given is the v-surface parameter of the nearest knot intersection.
     * \param[in] _P Given the Cartesian components of the point to be projected on the NURBS patch
     * \param[in] _uDiv Given the number of division in u-direction
     * \param[in] _vDiv Given the number of division in u-direction
     * \author Chenshen Wu
     ***********/
    void findInitialGuess4PointProjection(double& _u, double& _v, double* _P, int uDiv = 5,
            int _vDiv = 5);

    /***********************************************************************************************
     * \brief Find the nearest vertex from the linearization of the trimming curve as an initial guess for the projection
     * \param[in/out] _u Given is the u-surface parameter of the nearest linearization vertex.
     * \param[in/out] _v Given is the v-surface parameter of the nearest linearization vertex.
     * \param[in] _P Given the Cartesian components of the point to be projected on the trimming curve
     * \param[in] _patchBLIndex The index of the boundary loop which contains the trimming curve
     * \param[in] _patchBLTrCurveIndex The index of the trimming curve in the boundary loop
     * \author Altug Emiroglu
     ***********/
    bool findInitialGuess4PointProjectionOnTrimmingCurve(double& _uTilde, double& _u, double& _v, double* _P, int _patchBLIndex, int _patchBLTrCurveIndex);

    /***********************************************************************************************
     * \brief Find the nearest vertex from the linearization of the trimming curve as an initial guess for the projection
     * \param[in/out] _u Given is the u-surface parameter of the nearest linearization vertex.
     * \param[in/out] _v Given is the v-surface parameter of the nearest linearization vertex.
     * \param[in] _P Given the Cartesian components of the point to be projected on the trimming curve
     * \param[in] _curve The trimming curve to project the point onto
     * \author Altug Emiroglu
     ***********/
    bool findInitialGuess4PointProjectionOnTrimmingCurve(double& _uTilde, double& _u, double& _v,
                                                                          double* _P, IGAPatchCurve* _curve);

    /***********************************************************************************************
     * \brief Find the nearest knot intersection on the patch as an initial guess for the projection
     * \param[in/out] _coords The Cartesian coordinates of the point on the patch
     * \param[in/out] _normal The normal to the patch vector
     * \param[in] _u Given is the u-surface parameter
     * \param[in] _v Given is the v-surface parameter
     * \author Chenshen Wu
     ***********/
    void computeCartesianCoordinatesAndNormalVector(double* _coords, double* _normal, double _u,
            double _v)const;
    void computeCartesianCoordinatesAndNormalVector(double* _coords, double* _normal,
            double _u, double _v, int _spanU, int _spanV)const;

    /// Intersection related functions
public:
    /***********************************************************************************************
     * \brief Find the knot intersections of a trimming curve and return the coordinates
     * in the curve parameter space
     * \param[in/out] _uTilde The curve parameters of the intersections
     * \param[in] _patchBLIndex Given is the u-surface parameter
     * \param[in] _patchBLTrCurveIndex Given is the v-surface parameter
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    void computeKnotIntersectionsWithTrimmingCurve(std::vector<double>& _uTilde,
                                                   int _patchBLIndex, int _patchBLTrCurveIndex);

    /***********************************************************************************************
     * \brief Find the knot intersections of the given trimming curve and return the coordinates
     * in the curve parameter space
     * \param[in/out] _uTilde The curve parameters of the intersections
     * \param[in] _theCurve The curve on which the knot intersections to be computed
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    void computeKnotIntersectionsWithTrimmingCurve(std::vector<double>& _uTilde,
                                                   IGAPatchCurve* _theCurve);


    /// Postprocessing functions
public:
    /***********************************************************************************************
     * \brief Returns the approximate value of a field given its values on the Control Points at the specified parametric location
     * \param[in] _u Given is the u-surface parameter
     * \param[in] _v Given is the v-surface parameter
     * \param[in] _valuesOnCP The scalar values on the Control Points
     * \param[out] The approximate value of the scalar on (_u,_v)
     * \author Chenshen Wu
     ***********/
    double computePostprocessingScalarValue(double _u, double _v, double* _valuesOnCP);

    /// Get and set functions
public:
    char getEdge(const double _u, const double _v, const double _tolerance = 1e-3);
    /***********************************************************************************************
     * \brief Get the corners between edge 1 and edge 2 following the orientation
     * \param[in]	_edgeIn		The id of the edge going in the patch
     * \param[in]	_edgeOut	The id of the edge going out of the patch
     * \param[in]	_isCounterclockwise		The orientation of the polygon
     * \return The list of corner in between edgeIn and edgeOut
     * \author Fabien Pean
     ***********/
    std::vector<std::pair<double,double> > getCorner(const char _edge1, const char _edge2, const bool _isCounterclockwise);

    /***********************************************************************************************
     * \brief Get the underlying IsoGeometric basis of the patch
     * \author Andreas Apostolatos
     ***********/
    inline BSplineBasis2D* getIGABasis() const {
        return IGABasis;
    }
    inline BSplineBasis1D* getIGABasis(bool _u0v1) const {
        if(_u0v1==0)	return IGABasis->getUBSplineBasis1D();
        else 			return IGABasis->getVBSplineBasis1D();
    }
    inline BSplineBasis1D* operator[](const char _uOrV) const {
        if(_uOrV=='u')		return IGABasis->getUBSplineBasis1D();
        else if(_uOrV=='v') return IGABasis->getVBSplineBasis1D();
        assert(0);			return NULL;
    }

    /***********************************************************************************************
     * \brief Get the number of the Control Points of the patch in u-direction
     * \author Andreas Apostolatos
     ***********/
    inline int getUNoControlPoints() const {
        return uNoControlPoints;
    }

    /***********************************************************************************************
     * \brief Get the number of the Control Points of the patch in v-direction
     * \author Andreas Apostolatos
     ***********/
    inline int getVNoControlPoints() const {
        return vNoControlPoints;
    }

    /***********************************************************************************************
     * \brief Get the number of the Control Points of the patch
     * \author Chenshen Wu
     ***********/
    inline int getNoControlPoints() const {
        return uNoControlPoints * vNoControlPoints;
    }

    /***********************************************************************************************
     * \brief Get the Control Points of the patch
     * \author Andreas Apostolatos
     ***********/
    inline IGAControlPoint** getControlPointNet() const {
        return ControlPointNet;
    }
    inline IGAControlPoint* operator[](int i) const {
        return ControlPointNet[i];
    }

    /***********************************************************************************************
     * \brief Find know span on u direction
     * \author Chenshen Wu
     ***********/
    inline int findSpanU(double _u) const {
        return getIGABasis()->getUBSplineBasis1D()->findKnotSpan(_u);
    }

    /***********************************************************************************************
     * \brief Find know span on v direction
     * \author Chenshen Wu
     ***********/
    inline int findSpanV(double _v) const {
        return getIGABasis()->getVBSplineBasis1D()->findKnotSpan(_v);
    }

    /***********************************************************************************************
     * \brief Get Trimming class
     * \author Fabien Pean
     ***********/
    inline IGAPatchSurfaceTrimming& getTrimming() {
        return Trimming;
    }
    inline const IGAPatchSurfaceTrimming& getTrimming() const {
        return Trimming;
    }
    /***********************************************************************************************
     * \brief Check if patch is trimmed
     * \author Fabien Pean
     ***********/
    inline bool isTrimmed() const {
        return Trimming.isTrimmed();
    }

    inline double getBoundingBox(int id) {
        return boundingBox[id];
    }

    inline const BoundingBox& getBoundingBox() {
        return boundingBox;
    }

    /// The maximum number of Newton-Raphson iterations for the computation of the orthogonal projection of point on the NURBS patch
    static int MAX_NUM_ITERATIONS;
    /// The maximum relaxed number of Newton-Raphson iterations for the computation of the forced orthogonal projection of point on the NURBS patch
    static int REL_MAX_NUM_ITERATIONS;
    /// The orthogonality tolerance for the Newton-Raphson iterations for the computation of the orthogonal projection of point on the NURBS patch
    static double TOL_ORTHOGONALITY;
    /// The maximum relaxed orthogonality tolerance for the Newton-Raphson iterations for the computation of the forced orthogonal projection of point on the NURBS patch
    static double REL_TOL_ORTHOGONALITY;
    /// The distance tolerance for the Newton-Raphson iterations for the computation of the orthogonal projection of point on the NURBS patch
    static double TOL_DISTANCE;
    /// The maximum relaxed distance tolerance for the Newton-Raphson iterations for the computation of the forced orthogonal projection of point on the NURBS patch
    static double REL_TOL_DISTANCE;
    /// The tolerance for the ratio parameter in line projection on boundary
    static double TOL_PARAMETER_ON_STRAIGHT_EDGE;
    /// The tolerance for the parameter in the patch parameter space
    static double TOL_PARAMETER_ON_PATCH_PARAMETER_SPACE;

    static const char EDGE_U0;
    static const char EDGE_UN;
    static const char EDGE_V0;
    static const char EDGE_VN;
    static const char EDGES[4];

};

/***********************************************************************************************
 * \brief Allows for nice debug output
 * \author Fabien Pean, Chenshen Wu
 ***********/
Message &operator<<(Message &message, const IGAPatchSurface &mesh);

}/* namespace EMPIRE */

#endif /* IGAPatchSurface_H_ */


/******************************************************************************//**
 * \file IGAMath.h
 * The header file of math functions for mortar isogeometric mapping and for the basis functions
 * \date 1/29/2013
 *********************************************************************************/

#ifndef IGAMATH_H_
#define IGAMATH_H_

namespace EMPIRE {

/// The binomial coefficients stored in a vector format. The precomputed values are up to k=50
extern const double binomialCoefficients[2500];

/// The tolerance for the linear equation system solvers
extern const double EPS;

/***********************************************************************************************
 * \brief Computes the index function for the binomial coefficients (_i;_j)
 * \param[out] The index function for the binomial coefficients (_i;_j)
 * \param[in] _i The integer on the nominator
 * \param[in] _j The integer on the denominator
 * \author Andreas Apostolatos
 ***********/
int indexBinomialCoefficients(int, int);

/***********************************************************************************************
 * \brief Compute the square of the Euclidean distance of two points in n-D space.
 * \param[out] The square of the Euclidean distance between two points in the n-D space
 * \param[in] _length The length of the n-dimensional space
 * \param[in] _Pi The first point
 * \param[in] _Pj The second point
 * \author Andreas Apostolatos
 ***********/
double squareEuclideanDistance(int, double*, double*);

/***********************************************************************************************
 * \brief Compute the square of the 2-norm of a vector in the n-D space
 * \param[out] The square of the 2-norm of a vector in the n-D space
 * \param[in] _length The dimensinality of the n-D space
 * \param[in] _vector A vector in the nD space
 * \author Andreas Apostolatos
 ***********/
double square2normVector(int, double*);

/***********************************************************************************************
 * \brief Compute the dot product between two vectors in the n-D space
 * \param[out] The dot product between two vectors in the n-D space
 * \param[in] _length The dimensinality of the n-D space
 * \param[in] _vecI The 1st vector
 * \param[in] _vecJ The 2nd vector
 * \author Andreas Apostolatos
 ***********/
double dotProduct(int, double*, double*);

/***********************************************************************************************
 * \brief Solve a 2x2 linear equation system
 * \param[out] Flag on whether the linear system is solvable up to tolerance EPS or not
 * \param[in/out] _b The right-hand side vector where the solution is also stored, _b = double[2]
 * \param[in] _A The 2x2 matrix stored in a vector format namely A[i][j] = V[2 * i + j]
 * \author Andreas Apostolatos
 ***********/
bool solve2x2linearSystem(double*, double*);

}/* namespace EMPIRE */

#endif /* IGAMATH_H_ */

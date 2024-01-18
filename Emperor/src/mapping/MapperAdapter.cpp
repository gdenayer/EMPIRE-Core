/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Munich
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
#include "MapperAdapter.h"
#include "MortarMapper.h"
#include "IGAMortarMapper.h"
#include "IGABarycentricMapper.h"
#include "NearestNeighborMapper.h"
#include "BarycentricInterpolationMapper.h"
#include "NearestElementMapper.h"
#include "CurveSurfaceMapper.h"
#include "AbstractMesh.h"
#include "FEMesh.h"
#include "IGAMesh.h"
#include "SectionMesh.h"
#include "IGAPatchSurface.h" 	
#include "DataField.h"
#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include "AuxiliaryParameters.h"
#include "AbstractMapper.h"

using namespace std;

namespace EMPIRE {

MapperAdapter::MapperAdapter(std::string _name, AbstractMesh *_meshA, AbstractMesh *_meshB) :
        name(_name), meshA(_meshA), meshB(_meshB) {
    mapperImpl = NULL;
}

MapperAdapter::~MapperAdapter() {
    delete mapperImpl;
}

void MapperAdapter::initMortarMapper(bool oppositeSurfaceNormal, bool dual,
        bool enforceConsistency) {

    MortarMapper::mklSetNumThreads = AuxiliaryParameters::mklSetNumThreads;
    MortarMapper::mapperSetNumThreads = AuxiliaryParameters::mapperSetNumThreads;

    FEMesh *a = NULL;
    FEMesh *b = NULL;

    assert(meshA->type == EMPIRE_Mesh_FEMesh || meshA->type == EMPIRE_Mesh_SectionMesh);
    assert(meshB->type == EMPIRE_Mesh_FEMesh || meshB->type == EMPIRE_Mesh_SectionMesh);

    FEMesh *feMeshA = dynamic_cast<FEMesh *>(meshA);
    FEMesh *feMeshB = dynamic_cast<FEMesh *>(meshB);

    if (feMeshA->triangulate() == NULL) {
        a = feMeshA;
    } else {
        a = feMeshA->triangulate();
    }

    if (feMeshB->triangulate() == NULL) {
        b = feMeshB;
    } else {
        b = feMeshB->triangulate();
    }

    assert(mapperImpl == NULL);

    mapperImpl = new MortarMapper(a->numNodes, a->numElems, a->numNodesPerElem, a->nodes,
            a->nodeIDs, a->elems, b->numNodes, b->numElems, b->numNodesPerElem, b->nodes,
            b->nodeIDs, b->elems, oppositeSurfaceNormal, dual, enforceConsistency);
    MortarMapper* mapper = dynamic_cast<MortarMapper*>(mapperImpl);
    mapper->writeMode = this->writeMode;
    mapper->buildCouplingMatrices();
}

void MapperAdapter::initIGAMortarMapper(bool _enforceConsistency, double _tolConsistency,
                                        double _maxProjectionDistance, int _noInitialGuess, double _maxProjectionDistanceOnDifferentPatches,
                                        int _noIterationsNewton, double _tolProjectionNewtonRaphson,
                                        int _noIterationsNewtonRaphsonBoundary, double _tolProjectionNewtonRaphsonBoundary,
                                        int _noIterationsBisection, double _tolProjectionBisection,
                                        bool _isAutomaticNoGPTriangle, int _noGPTriangle, bool _isAutomaticNoGPQuadrilateral, int _noGPQuadrilateral,
                                        bool _isWeakCurveDirichletConditions, bool _isAutomaticPenaltyParametersWeakCurveDirichletConditions, bool _isPrimPrescribedWeakCurveDirichletConditions, bool _isSecBendingPrescribedWeakCurveDirichletConditions, bool _isSecTwistingPrescribedWeakCurveDirichletConditions, double _alphaPrimWeakCurveDirichletConditions, double _alphaSecBendingWeakCurveDirichletConditions, double _alphaSecTwistingWeakCurveDirichletConditions,
                                        bool _isWeakSurfaceDirichletConditions, bool _isAutomaticPenaltyParametersWeakSurfaceDirichletConditions, bool _isPrimPrescribedWeakSurfaceDirichletConditions, double _alphaPrimWeakSurfaceDirichletConditions,
                                        bool _isWeakPatchContinuityConditions, bool _isAutomaticPenaltyParametersWeakContinuityConditions, bool _isPrimCoupledWeakContinuityConditions, bool _isSecBendingCoupledWeakContinuityConditions, bool _isSecTwistingCoupledWeakContinuityConditions, double _alphaPrimWeakContinuityConditions, double _alphaSecBendingWeakContinuityConditions, double _alphaSecTwistingWeakContinuityConditions,
                                        bool _isStrongCurveDirichletConditions,
                                        bool _isErrorComputation, bool _isDomainError, bool _isCurveError, bool _isInterfaceError) {

    assert((meshA->type == EMPIRE_Mesh_FEMesh && meshB->type == EMPIRE_Mesh_IGAMesh) ||
           (meshB->type == EMPIRE_Mesh_FEMesh && meshA->type == EMPIRE_Mesh_IGAMesh));

    mapperImpl = new IGAMortarMapper(name, meshA, meshB);
    IGAMortarMapper* mapper = dynamic_cast<IGAMortarMapper*>(mapperImpl);
    mapper->writeMode = this->writeMode;
    mapper->setParametersConsistency(_enforceConsistency, _tolConsistency);
    mapper->setParametersProjection(_maxProjectionDistance, _noInitialGuess, _maxProjectionDistanceOnDifferentPatches);
    mapper->setParametersNewtonRaphson(_noIterationsNewton, _tolProjectionNewtonRaphson);
    mapper->setParametersNewtonRaphsonBoundary(_noIterationsNewtonRaphsonBoundary, _tolProjectionNewtonRaphsonBoundary);
    mapper->setParametersBisection(_noIterationsBisection, _tolProjectionBisection);
    mapper->setParametersIntegration(_isAutomaticNoGPTriangle, _noGPTriangle, _isAutomaticNoGPQuadrilateral, _noGPQuadrilateral);
    mapper->setParametersWeakCurveDirichletConditions(_isWeakCurveDirichletConditions, _isAutomaticPenaltyParametersWeakCurveDirichletConditions,
                                                      _isPrimPrescribedWeakCurveDirichletConditions, _isSecBendingPrescribedWeakCurveDirichletConditions,
                                                      _isSecTwistingPrescribedWeakCurveDirichletConditions, _alphaPrimWeakCurveDirichletConditions,
                                                      _alphaSecBendingWeakCurveDirichletConditions, _alphaSecTwistingWeakCurveDirichletConditions);
    mapper->setParametersWeakSurfaceDirichletConditions(_isWeakSurfaceDirichletConditions, _isAutomaticPenaltyParametersWeakSurfaceDirichletConditions,
                                                        _isPrimPrescribedWeakSurfaceDirichletConditions, _alphaPrimWeakSurfaceDirichletConditions);
    mapper->setParametersWeakPatchContinuityConditions(_isWeakPatchContinuityConditions, _isAutomaticPenaltyParametersWeakContinuityConditions,
                                                       _isPrimCoupledWeakContinuityConditions, _isSecBendingCoupledWeakContinuityConditions,
                                                       _isSecTwistingCoupledWeakContinuityConditions, _alphaPrimWeakContinuityConditions,
                                                       _alphaSecBendingWeakContinuityConditions, _alphaSecTwistingWeakContinuityConditions);
    mapper->setParametersStrongCurveDirichletConditions(_isStrongCurveDirichletConditions);
    mapper->setParametersErrorComputation(_isErrorComputation, _isDomainError, _isCurveError, _isInterfaceError);
    mapper->initialize();
    mapper->buildCouplingMatrices();
}

void MapperAdapter::initIGABarycentricMapper(double _maxProjectionDistance, int _numRefinementForIntialGuess,
    double _maxDistanceForProjectedPointsOnDifferentPatches, int _newtonRaphsonMaxIt, double _newtonRaphsonTol) {
    bool meshAIGA = (meshA->type == EMPIRE_Mesh_IGAMesh);
    bool meshBIGA = (meshB->type == EMPIRE_Mesh_IGAMesh);
    if (meshAIGA && !meshBIGA) {
        assert(meshB->type == EMPIRE_Mesh_FEMesh || meshB->type == EMPIRE_Mesh_SectionMesh);
        mapperImpl = new IGABarycentricMapper(name, dynamic_cast<IGAMesh *>(meshA),
                dynamic_cast<FEMesh *>(meshB), true);
    } else if (!meshAIGA && meshBIGA) {
        assert(meshA->type == EMPIRE_Mesh_FEMesh || meshA->type == EMPIRE_Mesh_SectionMesh);
        mapperImpl = new IGABarycentricMapper(name, dynamic_cast<IGAMesh *>(meshB),
                dynamic_cast<FEMesh *>(meshA), false);
    } else {
        ERROR_OUT() << "Error in MapperAdapter::initIGABarycentricMapper" << endl;
        ERROR_OUT() << "Wrong type of mesh! Put a NURBS mesh and a FE mesh!" << endl;
        exit(-1);
    }
    IGABarycentricMapper* mapper = dynamic_cast<IGABarycentricMapper*>(mapperImpl);
    mapper->writeMode = this->writeMode;
    mapper->setParametersProjection(_maxProjectionDistance, _numRefinementForIntialGuess,
            _maxDistanceForProjectedPointsOnDifferentPatches);
    mapper->setParametersNewtonRaphson(_newtonRaphsonMaxIt, _newtonRaphsonTol);
    mapper->buildCouplingMatrices();
}

void MapperAdapter::initNearestNeighborMapper() {
    assert(meshA->type == EMPIRE_Mesh_FEMesh || meshA->type == EMPIRE_Mesh_SectionMesh);
    assert(meshB->type == EMPIRE_Mesh_FEMesh || meshB->type == EMPIRE_Mesh_SectionMesh);

    FEMesh *a = dynamic_cast<FEMesh *>(meshA);
    FEMesh *b = dynamic_cast<FEMesh *>(meshB);
    assert(mapperImpl == NULL);

    mapperImpl = new NearestNeighborMapper(a->numNodes, a->nodes, b->numNodes, b->nodes);
    NearestNeighborMapper* mapper = dynamic_cast<NearestNeighborMapper*>(mapperImpl);
    mapper->writeMode = this->writeMode;
    mapper->buildCouplingMatrices();
}

void MapperAdapter::initBarycentricInterpolationMapper() {
    assert(meshA->type == EMPIRE_Mesh_FEMesh || meshA->type == EMPIRE_Mesh_SectionMesh);
    assert(meshB->type == EMPIRE_Mesh_FEMesh || meshB->type == EMPIRE_Mesh_SectionMesh);

    FEMesh *a = dynamic_cast<FEMesh *>(meshA);
    FEMesh *b = dynamic_cast<FEMesh *>(meshB);
    assert(mapperImpl == NULL);

    mapperImpl = new BarycentricInterpolationMapper(a->numNodes, a->nodes, b->numNodes, b->nodes);
    BarycentricInterpolationMapper* mapper = dynamic_cast<BarycentricInterpolationMapper*>(mapperImpl);
    mapper->writeMode = this->writeMode;
    mapper->buildCouplingMatrices();
}

void MapperAdapter::initNearestElementMapper() {
    assert(meshA->type == EMPIRE_Mesh_FEMesh || meshA->type == EMPIRE_Mesh_SectionMesh);
    assert(meshB->type == EMPIRE_Mesh_FEMesh || meshB->type == EMPIRE_Mesh_SectionMesh);

    FEMesh *a = dynamic_cast<FEMesh *>(meshA);
    FEMesh *b = dynamic_cast<FEMesh *>(meshB);
    assert(mapperImpl == NULL);
    NearestElementMapper::mapperSetNumThreads = AuxiliaryParameters::mapperSetNumThreads;

    mapperImpl = new NearestElementMapper(a->numNodes, a->numElems, a->numNodesPerElem, a->nodes,
            a->nodeIDs, a->elems, b->numNodes, b->numElems, b->numNodesPerElem, b->nodes,
            b->nodeIDs, b->elems);
    NearestElementMapper* mapper = dynamic_cast<NearestElementMapper*>(mapperImpl);
    mapper->writeMode = this->writeMode;
    mapper->buildCouplingMatrices();
}

void MapperAdapter::initCurveSurfaceMapper(EMPIRE_CurveSurfaceMapper_type type) {
    assert(meshA->type == EMPIRE_Mesh_FEMesh || meshA->type == EMPIRE_Mesh_SectionMesh);
    assert(meshB->type == EMPIRE_Mesh_SectionMesh);

    FEMesh *a = dynamic_cast<FEMesh *>(meshA);
    SectionMesh *b = dynamic_cast<SectionMesh *>(meshB);
    assert(mapperImpl == NULL);
    mapperImpl = new CurveSurfaceMapper(type, a->numNodes, a->numElems, a->nodes, a->nodeIDs,
            a->elems, b->numNodes, b->nodes, b->getNumSections(), b->getNumRootSectionNodes(),
            b->getNumNormalSectionNodes(), b->getNumTipSectionNodes(), b->getRotationGlobal2Root(),
            b->getTranslationGlobal2Root());
}

void MapperAdapter::consistentMapping(const DataField *fieldA, DataField *fieldB) {
    assert(mapperImpl != NULL);

    // 1. Do the consistent mapping
    if (dynamic_cast<CurveSurfaceMapper *>(mapperImpl) != NULL) { // CurveSurfaceMappers map DOFs together
        assert(fieldA->dimension == EMPIRE_DataField_doubleVector);
        assert(fieldB->dimension == EMPIRE_DataField_vector);
        mapperImpl->consistentMapping(fieldA->data, fieldB->data);
    }
    else if(dynamic_cast<IGAMortarMapper *>(mapperImpl) != NULL && (dynamic_cast<IGAMortarMapper *>(mapperImpl)->getIsExpanded())) {
        mapperImpl->consistentMapping(fieldA->data, fieldB->data);
    }
    else if(dynamic_cast<IGABarycentricMapper *>(mapperImpl) != NULL) {
        mapperImpl->consistentMapping(fieldA->data, fieldB->data);
    }
    else { // CurveSurfaceMappers map DOFs seperately
        int numNodesA, numNodesB;
        if (meshA->type == EMPIRE_Mesh_FEMesh || meshA->type == EMPIRE_Mesh_SectionMesh)
            numNodesA = dynamic_cast<FEMesh *>(meshA)->numNodes;
        else if (meshA->type == EMPIRE_Mesh_IGAMesh)
            numNodesA = dynamic_cast<IGAMesh *>(meshA)->getNumNodes();
        else {
            ERROR_OUT() << "Error in MapperAdapter::consistentMapping" << endl;
            ERROR_OUT() << "Wrong type of meshA!" << endl;
            exit(-1);
        }
        if (meshB->type == EMPIRE_Mesh_FEMesh || meshB->type == EMPIRE_Mesh_SectionMesh)
            numNodesB = dynamic_cast<FEMesh *>(meshB)->numNodes;
        else if (meshB->type == EMPIRE_Mesh_IGAMesh)
            numNodesB = dynamic_cast<IGAMesh *>(meshB)->getNumNodes();
        else {
            ERROR_OUT() << "Error in MapperAdapter::consistentMapping" << endl;
            ERROR_OUT() << "Wrong type of meshB!" << endl;
            exit(-1);
        }

        assert(fieldA->dimension == fieldB->dimension);
        assert(fieldA->numLocations == numNodesA);
        assert(fieldB->numLocations == numNodesB);

        double *fieldADOFi = new double[fieldA->numLocations];
        double *fieldBDOFi = new double[fieldB->numLocations];
        int numDOFs = fieldA->dimension;
        for (int i = 0; i < numDOFs; i++) {
            for (int j = 0; j < fieldA->numLocations; j++)
                fieldADOFi[j] = fieldA->data[j * numDOFs + i];
            mapperImpl->consistentMapping(fieldADOFi, fieldBDOFi);
            for (int j = 0; j < fieldB->numLocations; j++)
                fieldB->data[j * numDOFs + i] = fieldBDOFi[j];
        }
        delete[] fieldADOFi;
        delete[] fieldBDOFi;
    }

    // 2. Compute the mapping error
    if (dynamic_cast<IGAMortarMapper *>(mapperImpl) != NULL && dynamic_cast<IGAMortarMapper *>(mapperImpl)->getIsErrorComputation())
        mapperImpl->computeErrorsConsistentMapping(fieldA->data, fieldB->data);
}

void MapperAdapter::conservativeMapping(const DataField *fieldB, DataField *fieldA) {
    assert(mapperImpl != NULL);
    if (dynamic_cast<CurveSurfaceMapper *>(mapperImpl) != NULL) { // CurveSurfaceMappers map DOFs together
        assert(fieldA->dimension == EMPIRE_DataField_doubleVector);
        assert(fieldB->dimension == EMPIRE_DataField_vector);
        mapperImpl->conservativeMapping(fieldB->data, fieldA->data);
    } 
    else if(dynamic_cast<IGAMortarMapper *>(mapperImpl) != NULL && (dynamic_cast<IGAMortarMapper *>(mapperImpl)->getIsExpanded())) {
        mapperImpl->conservativeMapping(fieldB->data, fieldA->data);
    }
    else if(dynamic_cast<IGABarycentricMapper *>(mapperImpl) != NULL) {
        mapperImpl->conservativeMapping(fieldB->data, fieldA->data);
    }
    else { // CurveSurfaceMappers map DOFs seperately

        int numNodesA, numNodesB;
        if (meshA->type == EMPIRE_Mesh_FEMesh || meshA->type == EMPIRE_Mesh_SectionMesh)
            numNodesA = dynamic_cast<FEMesh *>(meshA)->numNodes;
        else if (meshA->type == EMPIRE_Mesh_IGAMesh)
            numNodesA = dynamic_cast<IGAMesh *>(meshA)->getNumNodes();
        else {
            ERROR_OUT() << "Error in MapperAdapter::conservativeMapping" << endl;
            ERROR_OUT() << "Wrong type of meshA!" << endl;
            exit(-1);
        }
        if (meshB->type == EMPIRE_Mesh_FEMesh || meshB->type == EMPIRE_Mesh_SectionMesh)
            numNodesB = dynamic_cast<FEMesh *>(meshB)->numNodes;
        else if (meshB->type == EMPIRE_Mesh_IGAMesh)
            numNodesB = dynamic_cast<IGAMesh *>(meshB)->getNumNodes();
        else {
            ERROR_OUT() << "Error in MapperAdapter::conservativeMapping" << endl;
            ERROR_OUT() << "Wrong type of meshB!" << endl;
            exit(-1);
        }

        assert(fieldA->dimension == fieldB->dimension);
        assert(fieldA->numLocations == numNodesA);
        assert(fieldB->numLocations == numNodesB);

        double *fieldADOFi = new double[fieldA->numLocations];
        double *fieldBDOFi = new double[fieldB->numLocations];
        int numDOFs = fieldA->dimension;
        for (int i = 0; i < numDOFs; i++) {
            for (int j = 0; j < fieldB->numLocations; j++)
                fieldBDOFi[j] = fieldB->data[j * numDOFs + i];
            mapperImpl->conservativeMapping(fieldBDOFi, fieldADOFi);
            for (int j = 0; j < fieldA->numLocations; j++)
                fieldA->data[j * numDOFs + i] = fieldADOFi[j];
        }
        delete[] fieldADOFi;
        delete[] fieldBDOFi;
    }
}

} /* namespace EMPIRE */

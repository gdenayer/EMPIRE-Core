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
#include "Empire.h"
#include "ClientMetaDatabase.h"
#include "ClientCommunication.h"
#include <assert.h>
#include <string.h>
#include <iostream>
#include "EMPIRE_API.h"
using namespace std;

namespace EMPIRE {

Empire::Empire() {

}

Empire::~Empire() {
}

void Empire::initEnvironment(char *inputFileName) {
    ClientMetaDatabase::init(inputFileName);
}

string Empire::getUserDefinedText(string elementName) {
    return ClientMetaDatabase::getSingleton()->getUserDefinedText(elementName);
}

void Empire::connect() {
    ClientCommunication::getSingleton()->connect();
}

void Empire::disconnect() {
    ClientCommunication::getSingleton()->disconnect();
    //ClientMetaDatabase::getSingleton()->free();
    //ClientCommunication::getSingleton()->free();
}

void Empire::sendMesh(int numNodes, int numElems, double *nodes, int *nodeIDs, int *numNodesPerElem,
        int *elems) {
    const int BUFFER_SIZE = 2;
    const int DIMENSION = 3;
    int meshInfo[BUFFER_SIZE] = { numNodes, numElems };
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE, meshInfo);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(numNodes * DIMENSION, nodes);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(numNodes, nodeIDs);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(numElems, numNodesPerElem);
    int count = 0;
    for (int i = 0; i < numElems; i++)
        count += numNodesPerElem[i];
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(count, elems);
}

void Empire::sendSectionMesh(int numNodes, int numElems, double *nodes, int *nodeIDs,
        int *numNodesPerElem, int *elems, int numSections, int numRootSectionNodes,
        int numNormalSectionNodes, int numTipSectionNodes, double *rotationGlobal2Root,
        double *translationGlobal2Root) {
    sendMesh(numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elems);

    const int SECTION_INFO_SIZE = 4;
    int sectionInfo[4]; // number of sections, number of root section nodes, number of normal section nodes, number of tip section nodes
    sectionInfo[0] = numSections;
    sectionInfo[1] = numRootSectionNodes;
    sectionInfo[2] = numNormalSectionNodes;
    sectionInfo[3] = numTipSectionNodes;
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(SECTION_INFO_SIZE, sectionInfo);

    ClientCommunication::getSingleton()->sendToServerBlocking<double>(9, rotationGlobal2Root);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(3, translationGlobal2Root);
}

void Empire::recvMesh(int *numNodes, int *numElems, double **nodes, int **nodeIDs,
        int **numNodesPerElem, int **elems) {
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>(1, numNodes);
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>(1, numElems);
    const int DIMENSION = 3;

    *nodes = new double[(*numNodes) * DIMENSION];
    *nodeIDs = new int[*numNodes];
    *numNodesPerElem = new int[*numElems];

    ClientCommunication::getSingleton()->receiveFromServerBlocking<double>((*numNodes) * DIMENSION,
            *nodes);
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>((*numNodes), *nodeIDs);
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>((*numElems),
            *numNodesPerElem);
    int count = 0;
    for (int i = 0; i < *numElems; i++)
        count += (*numNodesPerElem)[i];
    *elems = new int[count];
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>(count, *elems);
}

void Empire::sendIGAMesh(int _numPatches, int _numNodes) {
    const int BUFFER_SIZE = 2;
    int meshInfo[BUFFER_SIZE] = { _numPatches, _numNodes };
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE, meshInfo);
}

void Empire::sendIGAPatch(int _pDegree, int _uNumKnots, double* _uKnotVector, int _qDegree,
        int _vNumKnots, double* _vKnotVector, int _uNumControlPoints, int _vNumControlPoints,
        double* _cpNet, int* _nodeNet) {
    const int BUFFER_SIZE = 6;
    int meshInfo[BUFFER_SIZE] = { _pDegree, _uNumKnots, _qDegree, _vNumKnots, _uNumControlPoints,
            _vNumControlPoints };
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE, meshInfo);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_uNumKnots, _uKnotVector);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_vNumKnots, _vKnotVector);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(
            _uNumControlPoints * _vNumControlPoints * 4, _cpNet);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(
            _uNumControlPoints * _vNumControlPoints, _nodeNet);
}

void Empire::sendIGATrimmingInfo(int _isTrimmed, int _numLoops) {
    const int BUFFER_SIZE = 2;
    int trimInfo[BUFFER_SIZE] = { _isTrimmed, _numLoops };
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE, trimInfo);

}

void Empire::sendIGATrimmingLoopInfo(int _inner, int _numCurves) {
    const int BUFFER_SIZE = 2;
    int trimInfo[BUFFER_SIZE] = { _inner, _numCurves };
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE, trimInfo);
}

void Empire::sendIGATrimmingCurve(int _direction, int _pDegree, int _uNumKnots,
        double* _uKnotVector, int _uNumControlPoints, double* _cpNet) {
    const int BUFFER_SIZE = 4;
    int trimInfo[BUFFER_SIZE] = { _direction, _pDegree, _uNumKnots, _uNumControlPoints };
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE, trimInfo);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_uNumKnots, _uKnotVector);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_uNumControlPoints * 4,
            _cpNet);
}

void Empire::sendIGANumDirichletConditions(int _numDirichletConditions){
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1,&_numDirichletConditions);
}

void Empire::sendIGADirichletConditionInfo(int _patchIndex, int _patchBLIndex, int _patchBLTrCurveIndex,
                                        int _isGPProvided){
    const int BUFFER_SIZE = 3;
    int conditionInfo[BUFFER_SIZE] = { _patchIndex, _patchBLIndex, _patchBLTrCurveIndex };
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE,conditionInfo);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1,&_isGPProvided);

}

void Empire::sendIGADirichletConditionData(int _trCurveNumGP,
                                        double *_trCurveGPs, double *_trCurveGPWeights,
                                        double *_trCurveGPTangents,
                                        double *_trCurveGPJacobianProducts){

    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1,&_trCurveNumGP);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_trCurveNumGP*2,_trCurveGPs);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_trCurveNumGP,_trCurveGPWeights);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_trCurveNumGP*3,_trCurveGPTangents);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_trCurveNumGP,_trCurveGPJacobianProducts);

}

void Empire::sendIGANumPatchConnections(int _numConnections){
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1,&_numConnections);
}

void Empire::sendIGAPatchConnectionInfo(int _masterPatchIndex, int _masterPatchBLIndex, int _masterPatchBLTrCurveIndex,
                                        int _slavePatchIndex, int _slavePatchBLIndex, int _slavePatchBLTrCurveIndex,
                                        int _isGPProvided){
    const int BUFFER_SIZE = 6;
    int couplingInfo[BUFFER_SIZE] = { _masterPatchIndex, _masterPatchBLIndex, _masterPatchBLTrCurveIndex,
                                      _slavePatchIndex, _slavePatchBLIndex, _slavePatchBLTrCurveIndex };
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE,couplingInfo);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1,&_isGPProvided);

}

void Empire::sendIGAPatchConnectionData(int _trCurveNumGP,
                                        double *_trCurveMasterGPs, double *_trCurveSlaveGPs, double *_trCurveGPWeights,
                                        double *_trCurveMasterGPTangents, double *_trCurveSlaveGPTangents,
                                        double *_trCurveGPJacobianProducts){

    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1,&_trCurveNumGP);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_trCurveNumGP*2,_trCurveMasterGPs);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_trCurveNumGP*2,_trCurveSlaveGPs);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_trCurveNumGP,_trCurveGPWeights);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_trCurveNumGP*3,_trCurveMasterGPTangents);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_trCurveNumGP*3,_trCurveSlaveGPTangents);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_trCurveNumGP,_trCurveGPJacobianProducts);

}

void Empire::sendIGAPatchCouplingGaussPointsTest(int numPatches, int* numBRepsPerPatch, int totalNumGP, int totalNumBRePs,
        int* allID_slave, int* allNumElemsPerBRep, int* allNumGPsPerElem,
        double* allGPOfBRep_master, double* allGPOfBRep_slave, double* allGPOfBRep_weight,
        double* allTangents_master, double* allTangents_slave, double* allMappings) {
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(numPatches, numBRepsPerPatch);

    int CouplingInfo[2] = {totalNumGP, totalNumBRePs};
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(2, CouplingInfo);

    ClientCommunication::getSingleton()->sendToServerBlocking<int>(totalNumBRePs, allID_slave);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(totalNumBRePs, allNumElemsPerBRep);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(totalNumBRePs, allNumGPsPerElem);

    ClientCommunication::getSingleton()->sendToServerBlocking<double>(2*totalNumGP, allGPOfBRep_master);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(2*totalNumGP, allGPOfBRep_slave);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(totalNumGP, allGPOfBRep_weight);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(3*totalNumGP, allTangents_master);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(3*totalNumGP, allTangents_slave);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(totalNumGP, allMappings);
}

void Empire::sendIGADirichletDofs(int numberOfClampedDofs, int* clampedDofs, int clampedDirections) {
    int dirichletBCInfo[3]={numberOfClampedDofs, clampedDirections};
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(2, dirichletBCInfo);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(numberOfClampedDofs, clampedDofs);
}

void Empire::sendDataField(int sizeOfArray, double *dataField) {
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1, &sizeOfArray);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(sizeOfArray, dataField);
}

void Empire::recvDataField(int sizeOfArray, double *dataField) {
    int sizeOfArrayRecv = 0;
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>(1, &sizeOfArrayRecv);
    assert(sizeOfArray == sizeOfArrayRecv);
    ClientCommunication::getSingleton()->receiveFromServerBlocking<double>(sizeOfArray, dataField);
}

void Empire::sendSignal_double(char *name, int sizeOfArray, double *signal) {
    char nameBuffer[EMPIRE_API_NAME_STRING_LENGTH];
    strcpy(nameBuffer, name);
    ClientCommunication::getSingleton()->sendToServerBlocking<char>(EMPIRE_API_NAME_STRING_LENGTH,
            nameBuffer);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1, &sizeOfArray);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(sizeOfArray, signal);
}

void Empire::recvSignal_double(char *name, int sizeOfArray, double *signal) {
    char nameRecv[EMPIRE_API_NAME_STRING_LENGTH];
    ClientCommunication::getSingleton()->receiveFromServerBlocking<char>(
            EMPIRE_API_NAME_STRING_LENGTH, nameRecv);
    if (strcmp(name, nameRecv) != 0) {
        cout << "Error: signal names are not matching: " << name << " and, " << nameRecv << endl;
        assert(false);
    }
    int sizeOfArrayRecv = 0;
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>(1, &sizeOfArrayRecv);
    assert(sizeOfArray == sizeOfArrayRecv);
    ClientCommunication::getSingleton()->receiveFromServerBlocking<double>(sizeOfArray, signal);
}

void Empire::sendConvergenceSignal(int signal) {
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1, &signal);
}

int Empire::recvConvergenceSignal() {
    int signal;
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>(1, &signal);
    return signal;
}

void Empire::printDataField(char *name, int sizeOfArray, double *dataField) {
    std::cout << name << std::endl;
    for (int i = 0; i < sizeOfArray; i++)
        std::cout << dataField[i] << '\t';
    std::cout << std::endl;
}
}/* namespace EMPIRE */


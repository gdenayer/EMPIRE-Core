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
#include <assert.h>
#include <mpi.h>
// for the function 'gettimeofday' (microsecond resolution)
#include <sys/time.h>
#include <sstream>

#include "Emperor.h"
#include "ServerCommunication.h"
#include "ClientCode.h"
#include "DataOutput.h"
#include "Connection.h"
#include "ConnectionIO.h"
#include "MetaDatabase.h"
#include "EMPEROR_Enum.h"
#include "CopyFilter.h"
#include "SetFilter.h"
#include "MappingFilter.h"
#include "LocationFilter.h"
#include "ScalingFilter.h"
#include "WeakCouplingFilter.h"
#include "DataFieldIntegrationFilter.h"
#include "AdditionFilter.h"
#include "MapperAdapter.h"
#include "NearestNeighborMapper.h"
#include "BarycentricInterpolationMapper.h"
#include "FEMesh.h"
#include "AbstractMesh.h"
#include "Message.h"
#include "DataField.h"
#include "MapperAdapter.h"
#include "Aitken.h"
#include "ConstantRelaxation.h"
#include "IJCSA.h"
#include "GMRES.h"
#include "AbstractCouplingLogic.h"
#include "LinearExtrapolator.h"
#include "CouplingLogicSequence.h"
#include "IterativeCouplingLoop.h"
#include "TimeStepLoop.h"
#include "OptimizationLoop.h"
#include "PseudoCodeOutput.h"
#include "ConvergenceChecker.h"
#include "AuxiliaryFunctions.h"
#include "AuxiliaryParameters.h"
#include "Residual.h"

using namespace std;

namespace EMPIRE {

Emperor::Emperor() {
    globalCouplingLogic = NULL;
}

Emperor::~Emperor() {
    for (map<string, ClientCode*>::iterator it = nameToClientCodeMap.begin();
            it != nameToClientCodeMap.end(); it++) {
        delete it->second;
    }
    for (map<string, DataOutput*>::iterator it = nameToDataOutputMap.begin();
            it != nameToDataOutputMap.end(); it++) {
        delete it->second;
    }
    for (map<string, Connection*>::iterator it = nameToConnetionMap.begin();
            it != nameToConnetionMap.end(); it++) {
        delete it->second;
    }
    for (map<string, MapperAdapter*>::iterator it = nameToMapperMap.begin();
            it != nameToMapperMap.end(); it++) {
        delete it->second;
    }
    for (map<string, AbstractCouplingAlgorithm*>::iterator it = nameToCouplingAlgorithmMap.begin();
            it != nameToCouplingAlgorithmMap.end(); it++) {
        delete it->second;
    }
    for (map<string, AbstractExtrapolator*>::iterator it = nameToExtrapolatorMap.begin();
            it != nameToExtrapolatorMap.end(); it++) {
        delete it->second;
    }
    for (int i = 0; i < couplingLogicVec.size(); i++)
        delete couplingLogicVec[i];
}

void Emperor::initEnvironment(int *argc, char ***argv) {
    /// Check for command line arguments
    if (*argc != 2) {
        ERROR_BLOCK_OUT("Emperor", "initEnvironment", "Please provide a valid input file.");
    } else {
        /// Generate MetaDatabase in order give it to ServerCommunication::getSingleton()
        MetaDatabase::init((*argv)[1]);
        ServerCommunication::init(argc, argv);
    }
    ASCIIART_BLOCK();
    HEADING_OUT(1, "Emperor", "Emperor version " + AuxiliaryParameters::gitTAG + " started!",
            infoOut);

    PseudoCodeOutput *pcOutput = new PseudoCodeOutput(MetaDatabase::getSingleton(),
            "pseudocode.txt");
    pcOutput->writePseudoCode();
    delete pcOutput;
}

void Emperor::startServerListening() {
    ServerCommunication::getSingleton()->startListening();
}

void Emperor::startServerCoupling() {
    connectAllClients();
    // Start time stamps
    // for the function 'gettimeofday'
    struct timeval highrestimeStart, highrestimeEnd;
    stringstream timeMessage;

    // Capture the start time
    gettimeofday(&highrestimeStart, NULL);
    HEADING_OUT(4, "Emperor", "initClientCodes", infoOut);
    initClientCodes();
    // Capture the end time
    gettimeofday(&highrestimeEnd, NULL);
    double duration = AuxiliaryFunctions::highresDiffTime(highrestimeStart, highrestimeEnd);
    timeMessage.str(""); /// delete old message
    timeMessage << "It took " << duration << " seconds for initClientCodes";
    INDENT_OUT(1, timeMessage.str(), infoOut);

    //
    gettimeofday(&highrestimeStart, NULL);
    HEADING_OUT(4, "Emperor", "initDataOutputs", infoOut);
    initDataOutputs();
    gettimeofday(&highrestimeEnd, NULL);
    duration = AuxiliaryFunctions::highresDiffTime(highrestimeStart, highrestimeEnd);
    timeMessage.str(""); /// delete old message
    timeMessage << "It took " << duration << " seconds for initDataOutputs";
    INDENT_OUT(1, timeMessage.str(), infoOut);

    //
    gettimeofday(&highrestimeStart, NULL);
    HEADING_OUT(4, "Emperor", "initMappers", infoOut);
    initMappers();
    gettimeofday(&highrestimeEnd, NULL);
    duration = AuxiliaryFunctions::highresDiffTime(highrestimeStart, highrestimeEnd);
    timeMessage.str("");
    timeMessage << "It took " << duration << " seconds for initMappers";
    INDENT_OUT(1, timeMessage.str(), infoOut);

    //
    gettimeofday(&highrestimeStart, NULL);
    HEADING_OUT(4, "Emperor", "initCouplingAlgorithms", infoOut);
    initCouplingAlgorithms();
    gettimeofday(&highrestimeEnd, NULL);
    duration = AuxiliaryFunctions::highresDiffTime(highrestimeStart, highrestimeEnd);
    timeMessage.str("");
    timeMessage << "It took " << duration << " seconds for initCouplingAlgorithms";
    INDENT_OUT(1, timeMessage.str(), infoOut);

    //
    gettimeofday(&highrestimeStart, NULL);
    HEADING_OUT(4, "Emperor", "initExtrapolators", infoOut);
    initExtrapolators();
    gettimeofday(&highrestimeEnd, NULL);
    duration = AuxiliaryFunctions::highresDiffTime(highrestimeStart, highrestimeEnd);
    timeMessage.str("");
    timeMessage << "It took " << duration << " seconds for initExtrapolators";
    INDENT_OUT(1, timeMessage.str(), infoOut);

    //
    gettimeofday(&highrestimeStart, NULL);
    HEADING_OUT(4, "Emperor", "initConnections", infoOut);
    initConnections();
    gettimeofday(&highrestimeEnd, NULL);
    duration = AuxiliaryFunctions::highresDiffTime(highrestimeStart, highrestimeEnd);
    timeMessage.str("");
    timeMessage << "It took " << duration << " seconds for initConnections";
    INDENT_OUT(1, timeMessage.str(), infoOut);

    //
    gettimeofday(&highrestimeStart, NULL);
    HEADING_OUT(4, "Emperor", "initGlobalCouplingLogic", infoOut);
    initGlobalCouplingLogic();
    gettimeofday(&highrestimeEnd, NULL);
    duration = AuxiliaryFunctions::highresDiffTime(highrestimeStart, highrestimeEnd);
    timeMessage.str("");
    timeMessage << "It took " << duration << " seconds for initGlobalCouplingLogic";
    INDENT_OUT(1, timeMessage.str(), infoOut);

    //
    gettimeofday(&highrestimeStart, NULL);
    HEADING_OUT(4, "Emperor", "doCoSimulation", infoOut);
    doCoSimulation();
    gettimeofday(&highrestimeEnd, NULL);
    duration = AuxiliaryFunctions::highresDiffTime(highrestimeStart, highrestimeEnd);
    timeMessage.str("");
    timeMessage << "It took " << duration << " seconds for doCoSimulation";
    INDENT_OUT(1, timeMessage.str(), infoOut);

    disconnectAllClients();
}

void Emperor::connectAllClients() {
    while (1) {
        if (ServerCommunication::getSingleton()->allClientsConnected()) {
            INFO_OUT() << "All clients successfully connected!" << endl;
            break;
        }
    }
}

void Emperor::disconnectAllClients() {
    ServerCommunication::getSingleton()->disconnectAllClients();
    HEADING_OUT(1, "Emperor", "Emperor ends!", infoOut);
}

void Emperor::initClientCodes() {
    //set<string> clientNames;
    //ServerCommunication::getSingleton()->getClientNames(&clientNames);
    const vector<structClientCode> &settingClientCodesVec =
            MetaDatabase::getSingleton()->settingClientCodeVec;
    for (int i = 0; i < settingClientCodesVec.size(); i++) {
        const structClientCode &settingClientCode = settingClientCodesVec[i];
        string name = settingClientCode.name;
        //assert(clientNames.find(name)!=clientNames.end());
        const vector<structClientCode::structMesh> &settingMeshes = settingClientCode.meshes;
        const vector<structClientCode::structSignal> &settingSignals = settingClientCode.signals;
        const vector<std::string> &initialDataFields = settingClientCode.initialDataFields;

        ClientCode *clientCode = new ClientCode(name);
        nameToClientCodeMap.insert(pair<string, ClientCode*>(name, clientCode));
        clientCode->setServerCommunication(ServerCommunication::getSingleton());
        for (int j = 0; j < settingMeshes.size(); j++) {
            const structClientCode::structMesh &settingMesh = settingMeshes[j];
            if (settingMesh.type == EMPIRE_Mesh_FEMesh) {
                clientCode->recvFEMesh(settingMesh.name, settingMesh.triangulateAll);
            } else if (settingMesh.type == EMPIRE_Mesh_IGAMesh) {
                clientCode->recvIGAMesh(settingMesh.name);
            } else if (settingMesh.type == EMPIRE_Mesh_SectionMesh) {
                clientCode->recvSectionMesh(settingMesh.name, settingMesh.triangulateAll);
            } else if (settingMesh.type == EMPIRE_Mesh_copyFEMesh) {
                ClientCode *clientToCopyFrom = nameToClientCodeMap[settingMesh.clientNameToCopyFrom];
                AbstractMesh *meshToCopyFrom = clientToCopyFrom->getMeshByName(
                        settingMesh.meshNameToCopyFrom);
                clientCode->copyMesh(settingMesh.name, meshToCopyFrom);
                if (settingMesh.sendMeshToClient) {
                    clientCode->sendMesh(settingMesh.name);
                } else {
                    INFO_OUT() << "FEMesh copied but not sent" << endl;
                }
            } else if (settingMesh.type == EMPIRE_Mesh_copyIGAMesh) {
                ClientCode *clientToCopyFrom = nameToClientCodeMap[settingMesh.clientNameToCopyFrom];
                AbstractMesh *meshToCopyFrom = clientToCopyFrom->getMeshByName(
                        settingMesh.meshNameToCopyFrom);
                clientCode->copyMesh(settingMesh.name, meshToCopyFrom);
                if (settingMesh.sendMeshToClient) {
                    assert(false);
                    clientCode->sendMesh(settingMesh.name);
                } else {
                    INFO_OUT() << "IGAMesh copied but not sent" << endl;
                }
            } else {
                assert(false);
            }
            const vector<structClientCode::structMesh::structDataField> &settingDataFields =
                    settingMesh.dataFields;
            AbstractMesh *mesh = clientCode->getMeshByName(settingMesh.name);
            for (int k = 0; k < settingDataFields.size(); k++) {
                const structClientCode::structMesh::structDataField &settingDataField =
                        settingDataFields[k];
                mesh->addDataField(settingDataFields[k].name, settingDataFields[k].location,
                        settingDataFields[k].dimension, settingDataFields[k].typeOfQuantity);
            }

            // Receiving the initial values of the data fields specified. // Aditya
            for (int k = 0; k < initialDataFields.size(); k++) {
                clientCode->recvDataField(settingMesh.name, initialDataFields.at(k));
            }
        }


        for (int j = 0; j < settingSignals.size(); j++) {
            const structClientCode::structSignal &settingSignal = settingSignals[j];
            clientCode->addSignal(settingSignal.name, settingSignal.size3D[0],
                    settingSignal.size3D[1], settingSignal.size3D[2]);
        }
//        nameToClientCodeMap.insert(pair<string, ClientCode*>(name, clientCode));

    }
}

void Emperor::initDataOutputs() {
    const vector<structDataOutput> &settingDataOutputVec =
            MetaDatabase::getSingleton()->settingDataOutputVec;
    int numDataOutputs = settingDataOutputVec.size();
    for (int i = 0; i < numDataOutputs; i++) {
        DataOutput *dataOutput = new DataOutput(settingDataOutputVec[i], nameToClientCodeMap);
        nameToDataOutputMap.insert(
                pair<string, DataOutput*>(settingDataOutputVec[i].name, dataOutput));
    }
}

void Emperor::initMappers() {
    const vector<structMapper> &settingMapperVec = MetaDatabase::getSingleton()->settingMapperVec;
    int numMappers = settingMapperVec.size();
    for (int i = 0; i < numMappers; i++) {
        const structMapper &settingMapper = settingMapperVec[i];
        string name = settingMapper.name;
        const structMeshRef &meshRefA = settingMapper.meshRefA;
        const structMeshRef &meshRefB = settingMapper.meshRefB;

        assert(nameToClientCodeMap.find(meshRefA.clientCodeName) != nameToClientCodeMap.end());
        assert(nameToClientCodeMap.find(meshRefB.clientCodeName) != nameToClientCodeMap.end());
        AbstractMesh *meshA = nameToClientCodeMap.at(meshRefA.clientCodeName)->getMeshByName(
                meshRefA.meshName);
        AbstractMesh *meshB = nameToClientCodeMap.at(meshRefB.clientCodeName)->getMeshByName(
                meshRefB.meshName);

        MapperAdapter *mapper = new MapperAdapter(name, meshA, meshB);
        mapper->setWriteMode(settingMapper.writeMode);
        if (settingMapper.type == EMPIRE_MortarMapper) {
            mapper->initMortarMapper(settingMapper.mortarMapper.oppositeSurfaceNormal,
                    settingMapper.mortarMapper.dual, settingMapper.mortarMapper.enforceConsistency);
        } else if (settingMapper.type == EMPIRE_NearestNeighborMapper) {
            mapper->initNearestNeighborMapper();
        } else if (settingMapper.type == EMPIRE_BarycentricInterpolationMapper) {
            mapper->initBarycentricInterpolationMapper();
        } else if (settingMapper.type == EMPIRE_NearestElementMapper) {
            mapper->initNearestElementMapper();
        } else if (settingMapper.type == EMPIRE_IGAMortarMapper) {
            mapper->initIGAMortarMapper(
                    settingMapper.IGAMortarMapper.propConsistency.enforceConsistency,
                    settingMapper.IGAMortarMapper.propConsistency.tolConsistency,
                    settingMapper.IGAMortarMapper.propProjection.maxProjectionDistance,
                    settingMapper.IGAMortarMapper.propProjection.noInitialGuess,
                    settingMapper.IGAMortarMapper.propProjection.maxProjectionDistanceOnDifferentPatches,
                    settingMapper.IGAMortarMapper.propNewtonRaphson.noIterations,
                    settingMapper.IGAMortarMapper.propNewtonRaphson.tolProjection,
                    settingMapper.IGAMortarMapper.propNewtonRaphsonBoundary.noIterations,
                    settingMapper.IGAMortarMapper.propNewtonRaphsonBoundary.tolProjection,
                    settingMapper.IGAMortarMapper.propBisection.noIterations,
                    settingMapper.IGAMortarMapper.propBisection.tolProjection,
                    settingMapper.IGAMortarMapper.propIntegration.isAutomaticNoGPTriangle,
                    settingMapper.IGAMortarMapper.propIntegration.noGPTriangle,
                    settingMapper.IGAMortarMapper.propIntegration.isAutomaticNoGPQuadrilateral,
                    settingMapper.IGAMortarMapper.propIntegration.noGPQuadrilateral,
                    settingMapper.IGAMortarMapper.propWeakCurveDirichletConditions.isWeakCurveDirichletConditions,
                    settingMapper.IGAMortarMapper.propWeakCurveDirichletConditions.isAutomaticPenaltyParameters,
                    settingMapper.IGAMortarMapper.propWeakCurveDirichletConditions.isPrimPrescribed,
                    settingMapper.IGAMortarMapper.propWeakCurveDirichletConditions.isSecBendingPrescribed,
                    settingMapper.IGAMortarMapper.propWeakCurveDirichletConditions.isSecTwistingPrescribed,
                    settingMapper.IGAMortarMapper.propWeakCurveDirichletConditions.alphaPrim,
                    settingMapper.IGAMortarMapper.propWeakCurveDirichletConditions.alphaSecBending,
                    settingMapper.IGAMortarMapper.propWeakCurveDirichletConditions.alphaSecTwisting,
                    settingMapper.IGAMortarMapper.propWeakSurfaceDirichletConditions.isWeakSurfaceDirichletConditions,
                    settingMapper.IGAMortarMapper.propWeakSurfaceDirichletConditions.isAutomaticPenaltyParameters,
                    settingMapper.IGAMortarMapper.propWeakSurfaceDirichletConditions.isPrimPrescribed,
                    settingMapper.IGAMortarMapper.propWeakSurfaceDirichletConditions.alphaPrim,
                    settingMapper.IGAMortarMapper.propWeakPatchContinuityConditions.isWeakPatchContinuityConditions,
                    settingMapper.IGAMortarMapper.propWeakPatchContinuityConditions.isAutomaticPenaltyParameters,
                    settingMapper.IGAMortarMapper.propWeakPatchContinuityConditions.isPrimCoupled,
                    settingMapper.IGAMortarMapper.propWeakPatchContinuityConditions.isSecBendingCoupled,
                    settingMapper.IGAMortarMapper.propWeakPatchContinuityConditions.isSecTwistingCoupled,
                    settingMapper.IGAMortarMapper.propWeakPatchContinuityConditions.alphaPrim,
                    settingMapper.IGAMortarMapper.propWeakPatchContinuityConditions.alphaSecBending,
                    settingMapper.IGAMortarMapper.propWeakPatchContinuityConditions.alphaSecTwisting,
                    settingMapper.IGAMortarMapper.propStrongCurveDirichletConditions.isStrongCurveDirichletConditions,
                    settingMapper.IGAMortarMapper.propErrorComputation.isErrorComputation,
                    settingMapper.IGAMortarMapper.propErrorComputation.isDomainError,
                    settingMapper.IGAMortarMapper.propErrorComputation.isCurveError,
                    settingMapper.IGAMortarMapper.propErrorComputation.isInterfaceError);
        } else if (settingMapper.type == EMPIRE_IGABarycentricMapper) {
            mapper->initIGABarycentricMapper(
                    settingMapper.IGABarycentricMapper.propProjection.maxProjectionDistance,
                    settingMapper.IGABarycentricMapper.propProjection.noInitialGuess,
                    settingMapper.IGABarycentricMapper.propProjection.maxProjectionDistanceOnDifferentPatches,
                    settingMapper.IGABarycentricMapper.propNewtonRaphson.noIterations,
                    settingMapper.IGABarycentricMapper.propNewtonRaphson.tolProjection);
        } else if (settingMapper.type == EMPIRE_CurveSurfaceMapper) {
            mapper->initCurveSurfaceMapper(settingMapper.curveSurfaceMapper.type);
        } else {
            assert(false);
        }
        nameToMapperMap.insert(pair<string, MapperAdapter*>(name, mapper));
    }
}

void Emperor::initCouplingAlgorithms() {
    const std::vector<structCouplingAlgorithm> &settingCouplingAlgorithmVec =
            MetaDatabase::getSingleton()->settingCouplingAlgorithmVec;
    int numCoupAlgs = settingCouplingAlgorithmVec.size();
    for (int i = 0; i < numCoupAlgs; i++) {
        const structCouplingAlgorithm &settingCouplingAlgorithm = settingCouplingAlgorithmVec[i];
        string name = settingCouplingAlgorithm.name;

        // call constructor
        AbstractCouplingAlgorithm *couplingAlgorithm = NULL;
        if (settingCouplingAlgorithm.type == EMPIRE_Aitken) {
            double initialRelaxationFactor = settingCouplingAlgorithm.aitken.initialRelaxationFactor;
            couplingAlgorithm = new Aitken(name, initialRelaxationFactor);
        } else if (settingCouplingAlgorithm.type == EMPIRE_ConstantRelaxation) {
            double relaxationFactor = settingCouplingAlgorithm.constantRelaxation.relaxationFactor;
            couplingAlgorithm = new ConstantRelaxation(name, relaxationFactor);
        } else if (settingCouplingAlgorithm.type == EMPIRE_IJCSA) {
            couplingAlgorithm = new IJCSA(name);
            // add interfaceJacobianConsts
            for (int j = 0; j < settingCouplingAlgorithm.interfaceJacobians.size(); j++) {
                const structCouplingAlgorithm::structInterfaceJacobian &settingInterfaceJacobian =
                        settingCouplingAlgorithm.interfaceJacobians[j];
                if (IJCSA* couplingAlgorithmIJCSA = dynamic_cast<IJCSA*>(couplingAlgorithm)) {
                    if (settingCouplingAlgorithm.interfaceJacobians[j].isConstant == true) {
                        couplingAlgorithmIJCSA->addInterfaceJacobianEntry(
                                settingInterfaceJacobian.indexRow,
                                settingInterfaceJacobian.indexColumn,
                                settingInterfaceJacobian.value);
                    }
                    if (settingCouplingAlgorithm.interfaceJacobians[j].isAutoDiff == true) {
                        couplingAlgorithmIJCSA->addInterfaceJacobianEntry(
                                settingInterfaceJacobian.indexRow,
                                settingInterfaceJacobian.indexColumn,
                                constructConnectionIO(settingInterfaceJacobian.functionInput),
                                constructConnectionIO(settingInterfaceJacobian.functionOutput),
                                settingInterfaceJacobian.coefficient);
                    }

                } else {
                    assert(false);
                }

            }
        }else if(settingCouplingAlgorithm.type == EMPIRE_GMRES){

        	int maxOuterItter = settingCouplingAlgorithm.gmres.maxOuterItter;
        	int maxInnerItter = settingCouplingAlgorithm.gmres.maxInnerItter;
        	double residual = settingCouplingAlgorithm.gmres.residualTolerance;
        	couplingAlgorithm = new GMRES(name, maxOuterItter, maxInnerItter, residual);
        	int j;
        	for(j=0; j<settingCouplingAlgorithm.gmres.inputs.size(); j++){

        		if (GMRES* couplingAlgorithmGMRES = dynamic_cast<GMRES*>(couplingAlgorithm)){
        			couplingAlgorithmGMRES -> addInputConnection( constructConnectionIO(settingCouplingAlgorithm.gmres.inputs[j]));
        			couplingAlgorithmGMRES -> addOutputConnection( constructConnectionIO(settingCouplingAlgorithm.gmres.outputs[j]));
        		}
        	}

        }else { // Nothing is found
        	assert(false);
        }
        // add residuals
        for (int j = 0; j < settingCouplingAlgorithm.residuals.size(); j++) {
            const structResidual &settingResidual = settingCouplingAlgorithm.residuals[j];
            Residual *residual = new Residual(settingResidual.index);
            for (int k = 0; k < settingResidual.components.size(); k++) {
                const structResidual::structComponent &settingComponent =
                        settingResidual.components[k];
                residual->addComponent(settingComponent.coefficient, settingComponent.timeToUpdate,
                        constructConnectionIO(settingComponent.connectionIO));
            }
            residual->init();
            couplingAlgorithm->addResidual(residual, settingResidual.index);
        }
        // add outputs
        for (int j = 0; j < settingCouplingAlgorithm.outputs.size(); j++) {
            const structCouplingAlgorithm::structOutput &settingOutput =
                    settingCouplingAlgorithm.outputs[j];
            ConnectionIO *io = constructConnectionIO(settingOutput.connectionIO);
            couplingAlgorithm->addOutput(io, settingOutput.index);
        }
        // init all coupling algorithms
        couplingAlgorithm->init();

        nameToCouplingAlgorithmMap.insert(pair<string, AbstractCouplingAlgorithm*>(name, couplingAlgorithm));
    }
}

void Emperor::initExtrapolators() {
    const std::vector<structExtrapolator> &settingExtrapolatorVec =
            MetaDatabase::getSingleton()->settingExtrapolatorVec;
    for (int i = 0; i < settingExtrapolatorVec.size(); i++) {
        const structExtrapolator &settingExtrapolator = settingExtrapolatorVec[i];
        string name = settingExtrapolator.name;
        AbstractExtrapolator *extrapolator = NULL;
        if (settingExtrapolator.type == EMPIRE_LinearExtrapolator) {
            extrapolator = new LinearExtrapolator(name);
        } else {
            assert(false);
        }
        for (int j = 0; j < settingExtrapolator.connectionIOs.size(); j++) {
            extrapolator->addConnectionIO(
                    (constructConnectionIO(settingExtrapolator.connectionIOs[j])));
        }
        extrapolator->init();
        nameToExtrapolatorMap.insert(pair<string, AbstractExtrapolator*>(name, extrapolator));
    }
}

void Emperor::initConnections() {
    const std::vector<structConnection> &settingConnectionVec =
            MetaDatabase::getSingleton()->settingConnectionVec;
    int numConnections = settingConnectionVec.size();
    for (int i = 0; i < numConnections; i++) {
        Connection *connection;
        const structConnection &settingConnection = settingConnectionVec[i];
        string name = settingConnection.name;
        connection = new Connection(name);
        for (int j = 0; j < settingConnection.inputs.size(); j++) {
            const structConnectionIO &settingConnectionIO = settingConnection.inputs[j];
            connection->addInput(constructConnectionIO(settingConnectionIO));

        }
        for (int j = 0; j < settingConnection.outputs.size(); j++) {
            const structConnectionIO &settingConnectionIO = settingConnection.outputs[j];
            connection->addOutput(constructConnectionIO(settingConnectionIO));
        }

        nameToConnetionMap.insert(pair<string, Connection*>(name, connection));

        const vector<structFilter> &filterSequence = settingConnection.filterSequence;
        for (int j = 0; j < filterSequence.size(); j++) {
            const structFilter &settingFilter = filterSequence[j];
            AbstractFilter *filter;
            if (settingFilter.type == EMPIRE_MappingFilter) {
                string mapperName = settingFilter.mappingFilter.mapperName;
                assert(nameToMapperMap.find(mapperName) != nameToMapperMap.end());
                MapperAdapter *mapper = nameToMapperMap.at(mapperName);
                filter = new MappingFilter(mapper);
            } else if (settingFilter.type == EMPIRE_LocationFilter) {
                filter = new LocationFilter();
            } else if (settingFilter.type == EMPIRE_ScalingFilter) {
                filter = new ScalingFilter(settingFilter.scalingFilter.factor);
            } else if (settingFilter.type == EMPIRE_WeakCouplingFilter) {
                filter = new WeakCouplingFilter(settingFilter.weakCouplingFilter.beta);
            } else if (settingFilter.type == EMPIRE_SetFilter) {
                filter = new SetFilter(settingFilter.setFilter.value);
            } else if (settingFilter.type == EMPIRE_CopyFilter) {
                filter = new CopyFilter(settingFilter.copyFilter.signalOffset);
            } else if (settingFilter.type == EMPIRE_DataFieldIntegrationFilter) {
                const structMeshRef &meshRef = settingFilter.dataFieldIntegrationFilter.meshRef;
                assert(
                        nameToClientCodeMap.find(meshRef.clientCodeName)
                                != nameToClientCodeMap.end());
                AbstractMesh *mesh = nameToClientCodeMap.at(meshRef.clientCodeName)->getMeshByName(
                        meshRef.meshName);
                filter = new DataFieldIntegrationFilter(mesh);
            } else if (settingFilter.type == EMPIRE_AdditionFilter) {
                double a = settingFilter.additionFilter.a;
                double b = settingFilter.additionFilter.b;
                filter = new AdditionFilter(a, b);
            } else {
                assert(false);
            }
            for (int k = 0; k < settingFilter.inputs.size(); k++) {
                const structConnectionIO &settingConnectionIO = settingFilter.inputs[k];
                filter->addInput(constructConnectionIO(settingConnectionIO));
            }
            for (int k = 0; k < settingFilter.outputs.size(); k++) {
                const structConnectionIO &settingConnectionIO = settingFilter.outputs[k];
                filter->addOutput(constructConnectionIO(settingConnectionIO));
            }
            filter->init(); // initialize something after the inputs and outputs are set
            connection->addFilter(filter);
        }
    }
}

void Emperor::initGlobalCouplingLogic() {
    globalCouplingLogic = parseStructCouplingLogic(
            MetaDatabase::getSingleton()->settingGlobalCouplingLogic);
}

AbstractCouplingLogic *Emperor::parseStructCouplingLogic(
        structCouplingLogic &settingCouplingLogic) {
    vector<structCouplingLogic> &settingCouplingLogicSequence = settingCouplingLogic.sequence;
    AbstractCouplingLogic *couplingLogic;
    if (settingCouplingLogic.type == EMPIRE_connection) {
        // only a reference to an existing connection
        string connectionName = settingCouplingLogic.connectionRef.connectionName;
        assert(nameToConnetionMap.find(connectionName) != nameToConnetionMap.end());
        couplingLogic = nameToConnetionMap.at(connectionName);
    } else if (settingCouplingLogic.type == EMPIRE_CouplingLogicSequence) {
        couplingLogic = new CouplingLogicSequence();
    } else if (settingCouplingLogic.type == EMPIRE_IterativeCouplingLoop) {
        structCouplingLogic::structIterativeCouplingLoop &settingIterCoupLoop =
                settingCouplingLogic.iterativeCouplingLoop;
        structCouplingLogic::structIterativeCouplingLoop::structConvergenceChecker &settingConvgChecker =
                settingIterCoupLoop.convergenceChecker;
        vector < string > &convergenceObservers = settingIterCoupLoop.convergenceObservers;
        vector < string > &couplingAlgorithmRefs = settingIterCoupLoop.couplingAlgorithmRefs;

        IterativeCouplingLoop *iterativeCouplingLoop = new IterativeCouplingLoop();
        { // convergence checker
            double maxNumOfIterations = settingConvgChecker.maxNumOfIterations;
            ConvergenceChecker *checker = new ConvergenceChecker(maxNumOfIterations);
            for (int i = 0; i < settingConvgChecker.checkResiduals.size(); i++) {
                structCouplingLogic::structIterativeCouplingLoop::structConvergenceChecker::structCheckResidual & settingCheckResidual =
                        settingConvgChecker.checkResiduals[i];
                assert(
                        nameToCouplingAlgorithmMap.find(
                                settingCheckResidual.residualRef.couplingAlgorithmName)
                                != nameToCouplingAlgorithmMap.end());
                AbstractCouplingAlgorithm *coupAlgRef =
                        nameToCouplingAlgorithmMap[settingCheckResidual.residualRef.couplingAlgorithmName];
                checker->addCheckResidual(settingCheckResidual.absoluteTolerance,
                        settingCheckResidual.relativeTolerance, coupAlgRef,
                        settingCheckResidual.residualRef.index, settingCheckResidual.isAbsolute);
            }

            iterativeCouplingLoop->setConvergenceChecker(checker);
        }
        { // add convergence observers
            int tmpSize = convergenceObservers.size();
            for (int i = 0; i < tmpSize; i++) {
                assert(
                        nameToClientCodeMap.find(convergenceObservers[i])
                                != nameToClientCodeMap.end());
                ClientCode *clientCodeTmp = nameToClientCodeMap.at(convergenceObservers[i]);
                iterativeCouplingLoop->addConvergenceObserver(clientCodeTmp);
            }
        }
        { // set coupling algorithm
            int tmpSize = couplingAlgorithmRefs.size();
            for (int i = 0; i < tmpSize; i++) {
                assert(
                        nameToCouplingAlgorithmMap.find(
                                settingIterCoupLoop.couplingAlgorithmRefs[i])
                                != nameToCouplingAlgorithmMap.end());
                AbstractCouplingAlgorithm *couplingAlgorithmTmp = nameToCouplingAlgorithmMap.at(
                        settingIterCoupLoop.couplingAlgorithmRefs[i]);
                iterativeCouplingLoop->addCouplingAlgorithm(couplingAlgorithmTmp);
            }
        }
        { // add dataOutputs
            const vector<string> &dataOutputRefs = settingIterCoupLoop.dataOutputRefs;
            for (int i = 0; i < dataOutputRefs.size(); i++) {
                assert(nameToDataOutputMap.find(dataOutputRefs[i]) != nameToDataOutputMap.end());
                DataOutput *dataOutput = nameToDataOutputMap.at(dataOutputRefs[i]);
                iterativeCouplingLoop->addDataOutput(dataOutput);
            }
        }
        couplingLogic = iterativeCouplingLoop;
    } else if (settingCouplingLogic.type == EMPIRE_TimeStepLoop) {
        structCouplingLogic::structTimeStepLoop &settingTimeStepLoop =
                settingCouplingLogic.timeStepLoop;
        int numTimeSteps = settingTimeStepLoop.numTimeSteps;
        TimeStepLoop *timeStepLoop = new TimeStepLoop(numTimeSteps);
        { // set extrapolater
            if (settingTimeStepLoop.extrapolatorRef.first) {
                string extrapolatorName = settingTimeStepLoop.extrapolatorRef.second;
                assert(nameToExtrapolatorMap.find(extrapolatorName) != nameToExtrapolatorMap.end());
                AbstractExtrapolator *extrapolator = nameToExtrapolatorMap.at(extrapolatorName);
                timeStepLoop->setExtrapolator(extrapolator);
            }
        }
        { // add dataOutputs
            const vector<string> &dataOutputRefs = settingTimeStepLoop.dataOutputRefs;
            for (int i = 0; i < dataOutputRefs.size(); i++) {
                assert(nameToDataOutputMap.find(dataOutputRefs[i]) != nameToDataOutputMap.end());
                DataOutput *dataOutput = nameToDataOutputMap.at(dataOutputRefs[i]);
                timeStepLoop->addDataOutput(dataOutput);
            }
        }
        couplingLogic = timeStepLoop;
    } else if (settingCouplingLogic.type == EMPIRE_OptimizationLoop) {
        structCouplingLogic::structOptimizationLoop &settingOptLoop =
                settingCouplingLogic.optimizationLoop;
        OptimizationLoop *optimizationLoop = new OptimizationLoop(
                settingOptLoop.maxNumOfIterations);
        { // convergence signal sender
            assert(
                    nameToClientCodeMap.find(settingOptLoop.convergenceSignalSender)
                            != nameToClientCodeMap.end());
            ClientCode *clientCodeTmp = nameToClientCodeMap.at(
                    settingOptLoop.convergenceSignalSender);
            optimizationLoop->setConvergenceSignalSender(clientCodeTmp);
        }

        { // convergence signal receivers
            for (int i = 0; i < settingOptLoop.convergenceSignalReceivers.size(); i++) {
                assert(
                        nameToClientCodeMap.find(settingOptLoop.convergenceSignalReceivers[i])
                                != nameToClientCodeMap.end());
                ClientCode *clientCodeTmp = nameToClientCodeMap.at(
                        settingOptLoop.convergenceSignalReceivers[i]);
                optimizationLoop->addConvergenceSignalReceiver(clientCodeTmp);
            }
        }
        { // add dataOutputs
            const vector<string> &dataOutputRefs = settingOptLoop.dataOutputRefs;
            for (int i = 0; i < dataOutputRefs.size(); i++) {
                assert(nameToDataOutputMap.find(dataOutputRefs[i]) != nameToDataOutputMap.end());
                DataOutput *dataOutput = nameToDataOutputMap.at(dataOutputRefs[i]);
                optimizationLoop->addDataOutput(dataOutput);
            }
        }
        couplingLogic = optimizationLoop;
    } else {
        assert(false);
    }

    if (settingCouplingLogic.type != EMPIRE_connection)
        couplingLogicVec.push_back(couplingLogic); // if it is not a connection, delete it in the destructor

    int size = settingCouplingLogicSequence.size();
    for (int i = 0; i < size; i++) {
        couplingLogic->addCouplingLogic(parseStructCouplingLogic(settingCouplingLogicSequence[i]));
    }
    return couplingLogic;
}

void Emperor::doCoSimulation() {
    assert(globalCouplingLogic->size() != 0);
    globalCouplingLogic->doCoupling();
}

ConnectionIO *Emperor::constructConnectionIO(const structConnectionIO &settingConnectionIO) {
    ConnectionIO *io;
    if (settingConnectionIO.type == EMPIRE_ConnectionIO_DataField) {
        io = new ConnectionIO(nameToClientCodeMap, settingConnectionIO.dataFieldRef.clientCodeName,
                settingConnectionIO.dataFieldRef.meshName,
                settingConnectionIO.dataFieldRef.dataFieldName);
    } else if (settingConnectionIO.type == EMPIRE_ConnectionIO_Signal) {
        io = new ConnectionIO(nameToClientCodeMap, settingConnectionIO.signalRef.clientCodeName,
                settingConnectionIO.signalRef.signalName);
    } else {
        assert(false);
    }
    return io;
}

} /* namespace EMPIRE */

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
#include "Aitken.h"
#include "DataField.h"
#include "ConnectionIO.h"
#include "Signal.h"
#include "Residual.h"
#include "MathLibrary.h"

#include <assert.h>
#include <math.h>
#include <sstream>
#include <string.h>
//#include <fstream>
//#include <iomanip>

using namespace std;

namespace EMPIRE {

Aitken::Aitken(std::string _name, double _initRelaxationFactor) :
        AbstractCouplingAlgorithm(_name), INIT_RELAXATION_FACTOR(_initRelaxationFactor) {
    debugMe = false;
    globalResidual = NULL;
    globalResidualOld = NULL;
    tmpVec = NULL;
    
    string filename = "outputEmperor_" + _name + ".dat";
    file.open (filename.c_str());
}

Aitken::~Aitken() {
    delete[] globalResidual;
    delete[] globalResidualOld;
    delete[] tmpVec;
    
    file.close();
}

void Aitken::calcNewValue() {
	/// reset if new time step is started
	if (newTimeStep)
	{
		startNewTimeStep();
	}
//    /// compute the current residuals
//    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++) {
//        it->second->computeCurrentResidual();
//    }
    /// assemble global residual vector
    int oldSize =0;
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++) {
        Residual *residual = it->second;
    	for (int i=0; i<it->second->size; i++) {
    		globalResidual[i+oldSize]=residual->residualVector[i];
    	}
    	oldSize+=it->second->size;;
    }
    /// compute Aitken update
	if (newTimeStep)
	{
	    stringstream toOutput;
	    toOutput << scientific;
	    toOutput << "Initial Aitken relaxation factor: " << relaxationFactor;
	    INDENT_OUT(1, toOutput.str(), infoOut);
		newTimeStep=false;
		cout.unsetf(ios_base::floatfield);
	}else{
	    computeRelaxationFactor();
	}

    /// apply the new output
    assert(outputs.size() == residuals.size());
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++) {
        Residual *residual = it->second;
        assert(outputs.find(it->first) != outputs.end());
        CouplingAlgorithmOutput *output = outputs.find(it->first)->second;

        assert(residual->size == output->size);
        double *newOuput = new double[residual->size];
        // U_i_n+1 = U_i_n + alpha R_i_n
        for (int i=0; i<residual->size; i++) {
            newOuput[i] = output->outputCopyAtIterationBeginning[i]+ relaxationFactor*residual->residualVector[i] ;
        }

        output->overwrite(newOuput);
        delete[] newOuput;
    }

    /// save old values
    MathLibrary::copyDenseVector(globalResidualOld,globalResidual,globalResidualSize);

    relaxationFactorOld=relaxationFactor;
    
    file << currentTimeStep << "\t" << currentIteration << "\t" << relaxationFactor << "\t" << endl;


}

void Aitken::computeRelaxationFactor() {
	MathLibrary::copyDenseVector(tmpVec,globalResidual,globalResidualSize);
	MathLibrary::computeDenseVectorMultiplicationScalar(tmpVec,-1.0,globalResidualSize);
	MathLibrary::computeDenseVectorAddition(tmpVec,globalResidualOld,1.0,globalResidualSize);
	/// tmpVec holds now globalResidualOld - globalResidual
	double denominator = MathLibrary::vector2norm(tmpVec,globalResidualSize);
	denominator*=denominator;
	//double numerator = MathLibrary::computeDenseDotProduct(globalResidualOld, tmpVec, globalResidualSize);
	// Edit Aditya
	double numerator = MathLibrary::computeDenseDotProduct(globalResidualSize,globalResidualOld,tmpVec);
    if (denominator>1e-30){
        relaxationFactor = relaxationFactorOld * (numerator/denominator);
    }else{
        stringstream toOutput;
        toOutput << scientific;
        toOutput << "denominator = " << denominator << " <= 1.0e-30 close to zero !!!";
        WARNING_BLOCK_OUT("Aitken", "computeRelaxationFactor()", toOutput.str());
    }
    stringstream toOutput;
    toOutput << scientific;
    toOutput << "Aitken relaxation factor: " << relaxationFactor;
    INDENT_OUT(1, toOutput.str(), infoOut);
    cout.unsetf(ios_base::floatfield);
}

void Aitken::startNewTimeStep() {
	relaxationFactor   =INIT_RELAXATION_FACTOR;
	relaxationFactorOld=INIT_RELAXATION_FACTOR;
	memset(globalResidual   ,0,sizeof(double)*globalResidualSize);
	memset(globalResidualOld,0,sizeof(double)*globalResidualSize);
}

void Aitken::init() {
    // determine global residual vector size
    globalResidualSize =0;
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++) {
    	globalResidualSize +=it->second->size;
    }
    globalResidual    = new double [globalResidualSize];
    globalResidualOld = new double [globalResidualSize];
    tmpVec            = new double [globalResidualSize];
    startNewTimeStep();
}

void Aitken::calcCurrentResidual() {
    /// compute the current residuals
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++) {
        it->second->computeCurrentResidual();
    }
}






//// Writing out the ouputCopyAtIterationBeginning
//std::ofstream file;
//char fileName[22];
//int n = sprintf(fileName,"outputCopyAtIterationBeginning_%d_%d",currentTimeStep,currentIteration);
//file.open (fileName);
//for (int i = 0; i < output->size; i++) {
//    file << output->outputCopyAtIterationBeginning[i+0]<<" "<<output->outputCopyAtIterationBeginning[i+1]<<" "<<output->outputCopyAtIterationBeginning[i+2]<<"\n";
//}
//
//
//// Writing out the ouputCopyAtIterationBeginning
//std::ofstream testfile;
//char testfileName[22];
//int in = sprintf(testfileName,"residualVector_%d_%d",currentTimeStep,currentIteration);
//testfile.open (testfileName);
//for (int i = 0; i < output->size; i++) {
//	testfile << residual->residualVector[i+0]<<" "<<residual->residualVector[i+1]<<" "<<residual->residualVector[i+2]<<"\n";
//}


} /* namespace EMPIRE */

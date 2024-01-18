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
#ifndef WEAKIGAPATCHCONTINUITYCONDITION_H_
#define WEAKIGAPATCHCONTINUITYCONDITION_H_

#include <string>
#include <map>
#include <vector>
#include "EMPEROR_Enum.h"
#include "AbstractCondition.h"

namespace EMPIRE {

class Message;
class IGAPatchCurve;
class IGAPatchSurface;

/********//**
 * \brief Class WeakIGAPatchContinuityCondition is the class that holds and applies the weak continuity conditions for IGA
 ***********/
class WeakIGAPatchContinuityCondition: public AbstractCondition {

protected:

    /// The master curve
    IGAPatchCurve* masterCurve;

    /// The slave curve
    IGAPatchCurve* slaveCurve;

    /// If the curves belong to a trimming loop of a patch
    bool isTrimmingCurve;

    /// If the GP data for the coupling condition is initialized
    bool isGPDataInitialized;

    /// Number of GPs on the coupled trimming curve
    int trCurveNumGP;

    /// GP weights on the coupled trimming curve
    double* trCurveGPWeights;

    /// Jacobian products at the GPs on the coupled trimming curve
    double* trCurveGPJacobianProducts;

    /// The index of the master patch in the sending order
    int masterPatchIndex;

    /// The index of the master patch boundary loop in the sending order
    int masterPatchBLIndex;

    /// The index of the master patch boundary loop trimming curve in the sending order
    int masterPatchBLTrCurveIndex;

    /// GP parametric coordinates on the trimming curve of the master patch
    double* trCurveMasterGPs;

    /// Trimming curve tangents in the parametric space of the master patch
    double* trCurveMasterGPTangents;

    /// The index of the slave patch in the sending order
    int slavePatchIndex;

    /// The index of the slave patch boundary loop in the sending order
    int slavePatchBLIndex;

    /// The index of the slave patch boundary loop trimming curve in the sending order
    int slavePatchBLTrCurveIndex;

    /// GP parametric coordinates on the trimming curve of the slave patch
    double* trCurveSlaveGPs;

    /// Trimming curve tangents in the parametric space of the slave patch
    double* trCurveSlaveGPTangents;

public:

    /***********************************************************************************************
     * brief Constructor, initializing the class
     * \param[in] _ID ID of the condition
     * \param[in] _masterPatchIndex The index of the patch in the EMPIRE data structure
     * \param[in] _pMaster The polynomial degree of the IGA 1D curve in the u-direction
     * \param[in] _uNoKnotsMaster The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVectorMaster The underlying knot vector of the IGA 1D curve in the u-direction
     * \param[in] _uNoControlPointsMaster The number of the Control Points for the 1D NURBS patch in the u-direction
     * \param[in] _controlPointNetMaster The set of the Control Points related to the 1D NURBS patch
     * \param[in] _slavePatchIndex The index of the patch in the EMPIRE data structure
     * \param[in] _pDegreeSlave The polynomial degree of the IGA 1D curve in the u-direction
     * \param[in] _uNoKnotsSlave The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVectorSlave The underlying knot vector of the IGA 1D curve in the u-direction
     * \param[in] _uNoControlPointsSlave The number of the Control Points for the 1D NURBS patch in the u-direction
     * \param[in] _controlPointNetSlave The set of the Control Points related to the 1D NURBS patch
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    WeakIGAPatchContinuityCondition(int _ID,
                                    int _masterPatchIndex, int _pMaster, int _uNoKnotsMaster, double* _uKnotVectorMaster, int _uNoControlPointsMaster, double* _controlPointNetMaster,
                                    int _slavePatchIndex,  int _pSlave, int _uNoKnotsSlave, double* _uKnotVectorSlave, int _uNoControlPointsSlave, double* _controlPointNetSlave);

    /***********************************************************************************************
     * brief Constructor, initializing the class
     * \param[in] _ID ID of the condition
     * \param[in] _masterPatchIndex The index of the patch in the EMPIRE data structure
     * \param[in] _masterCurve The master curve
     * \param[in] _slavePatchIndex The index of the patch in the EMPIRE data structure
     * \param[in] _slaveCurve The slave curve
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    WeakIGAPatchContinuityCondition(int _ID,
                                    int _masterPatchIndex, IGAPatchCurve* _masterCurve,
                                    int _slavePatchIndex, IGAPatchCurve* _slaveCurve);

    /***********************************************************************************************
     * \brief Constructor, initializing the class
     * \param[in] _ID ID of the condition
     * \param[in] _masterPatchIndex The index of the master patch in the EMPIRE data structure
     * \param[in] _masterPatchBLIndex The index of the master patch boundary loop in the EMPIRE data structure
     * \param[in] _masterPatchBLTrCurveIndex The index of the master patch trimming curve in the current boundary loop in the EMPIRE data structure
     * \param[in] _slavePatchIndex The index of the slave patch in the EMPIRE data structure
     * \param[in] _slavePatchBLIndex The index of the slave patch boundary loop in the EMPIRE data structure
     * \param[in] _slavePatchBLTrCurveIndex The index of the slave patch trimming curve in the current boundary loop in the EMPIRE data structure
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    WeakIGAPatchContinuityCondition(int _ID,
                    int _masterPatchIndex, int _masterPatchBLIndex, int _masterPatchBLTrCurveIndex,
                    int _slavePatchIndex, int _slavePatchBLIndex, int _slavePatchBLTrCurveIndex);

    /***********************************************************************************************
     * \brief Set the GP data for the coupling condition
     * \param[in] _trCurveNumGP The total number of GPs on the trimming curve
     * \param[in] _trCurveMasterGPs The parametric coordinates of the GPs on the master trimming curve in the master patch parameter space
     * \param[in] _trCurveSlaveGPs The parametric coordinates of the GPs on the slave trimming curve in the slave patch parameter space
     * \param[in] _trCurveGPWeights The GP weights
     * \param[in] _trCurveMasterGPTangents The tangent to the trimming curve vector in the parameter space of the master patch
     * \param[in] _trCurveSlaveGPTangents The tangent to the trimming curve vector in the parameter space of the slave patch
     * \param[in] _trCurveGPJacobianProducts The Jacobian products for the transformation of the integrals
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void addWeakContinuityConditionGPData(int _trCurveNumGP,
                                          double* _trCurveMasterGPs, double* _trCurveSlaveGPs, double* _trCurveGPWeights,
                                          double* _trCurveMasterGPTangents, double* _trCurveSlaveGPTangents,
                                          double* _trCurveGPJacobianProducts);
    
    /***********************************************************************************************
     * \brief Set the GP data for the coupling condition
     * \param[in] _masterPatch The master patch
     * \param[in] _slavePatch The slave patch
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void createGPData(const std::vector<IGAPatchSurface*>& _surfacePatches);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    ~WeakIGAPatchContinuityCondition();

    /***********************************************************************************************
     * \brief get masterPatchIndex
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getMasterPatchIndex() const {
        return masterPatchIndex;
    }

    /***********************************************************************************************
     * \brief get slavePatchIndex
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getSlavePatchIndex() const {
        return slavePatchIndex;
    }

    /***********************************************************************************************
     * \brief get masterPatchBLIndex
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getMasterPatchBLIndex() const {
        return masterPatchBLIndex;
    }

    /***********************************************************************************************
     * \brief get slavePatchBLIndex
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getSlavePatchBLIndex() const {
        return slavePatchBLIndex;
    }

    /***********************************************************************************************
     * \brief get masterPatchBLTrCurveIndex
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getMasterPatchBLTrCurveIndex() const {
        return masterPatchBLTrCurveIndex;
    }

    /***********************************************************************************************
     * \brief get slavePatchBLTrCurveIndex
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getSlavePatchBLTrCurveIndex() const {
        return slavePatchBLTrCurveIndex;
    }

    /***********************************************************************************************
     * \brief get trCurveNumGP
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getTrCurveNumGP() const {
        return trCurveNumGP;
    }

    /***********************************************************************************************
     * \brief get trCurveMasterGPs
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getTrCurveMasterGPs() const {
        return trCurveMasterGPs;
    }

    /***********************************************************************************************
     * \brief get trCurveSlaveGPs
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getTrCurveSlaveGPs() const {
        return trCurveSlaveGPs;
    }

    /***********************************************************************************************
     * \brief get trCurveGPWeights
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getTrCurveGPWeights() const {
        return trCurveGPWeights;
    }

    /***********************************************************************************************
     * \brief get trCurveSlaveGPTangents
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getTrCurveMasterGPTangents() const {
        return trCurveMasterGPTangents;
    }

    /***********************************************************************************************
     * \brief get trCurveSlaveGPTangents
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getTrCurveSlaveGPTangents() const {
        return trCurveSlaveGPTangents;
    }

    /***********************************************************************************************
     * \brief get trCurveGPJacobianProducts
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getTrCurveGPJacobianProducts() const {
        return trCurveGPJacobianProducts;
    }

};

} /* namespace EMPIRE */

#endif /* WEAKIGAPATCHCONTINUITYCONDITION_H_ */

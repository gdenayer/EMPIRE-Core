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
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include "KinematicMotion.h"
#include <iostream>
#include <string>
#include <math.h>

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the class KinematicMotion.
 ***********/
class TestKinematicMotion: public CppUnit::TestFixture {
private:

public:
    void setUp() {
    }
    void tearDown() {
    }

    /***********************************************************************************************
     * \brief Test constructor
     ***********/
    void testConstructor() {
        const double TOL = 1E-3;
        KinematicMotion *km = new KinematicMotion();
        double point[] = { 1.0, 2.0, 3.0 };
        km->move(point);
        CPPUNIT_ASSERT(fabs(point[0] - 1.0) < TOL);
        CPPUNIT_ASSERT(fabs(point[1] - 2.0) < TOL);
        CPPUNIT_ASSERT(fabs(point[2] - 3.0) < TOL);
        delete km;
    }

    /***********************************************************************************************
     * \brief Test motion and coordinates transformation
     ***********/
    void testMotionAndCoorTransformation() {
        // test coordinate transformation, and move point (the inverse is also tested)
        // define a local coordinate system A, assume O is the global coordinate system
        const double TOL = 1E-3;

        KinematicMotion *km_O_A = new KinematicMotion();

        double translationVec[] = { 2.0, 2.0, 0.0 };
        double zAxis[] = { 0.0, 0.0, 1.0 };
        km_O_A->addRotation(zAxis, true, M_PI / 4.0);
        km_O_A->addTranslation(translationVec);
        KinematicMotion *km_A_O = km_O_A->newInverse();

        { // x_O -> x_A
            double x_O[] = { 1.0, 1.0, 0.0 };
            km_A_O->move(x_O); // x_O -> x_A
            double *ref_x_A = x_O;
            CPPUNIT_ASSERT(fabs(ref_x_A[0] - (-1.41421)) < TOL);
            CPPUNIT_ASSERT(fabs(ref_x_A[1] - 0.0) < TOL);
            CPPUNIT_ASSERT(fabs(ref_x_A[2] - 0.0) < TOL);
        }
        { // x_A -> x_O
            double x_A[] = { -1.41421, 0.0, 0.0 };
            km_O_A->move(x_A); // x_A -> x_O
            double *ref_x_O = x_A;
            CPPUNIT_ASSERT(fabs(ref_x_O[0] - 1.0) < TOL);
            CPPUNIT_ASSERT(fabs(ref_x_O[1] - 1.0) < TOL);
            CPPUNIT_ASSERT(fabs(ref_x_O[2] - 0.0) < TOL);
        }
        { // move x_O in O
            double x_O[] = { 1.0, 1.0, 0.0 };
            km_O_A->move(x_O);
            CPPUNIT_ASSERT(fabs(x_O[0] - 2.0) < TOL);
            CPPUNIT_ASSERT(fabs(x_O[1] - 3.41421) < TOL);
            CPPUNIT_ASSERT(fabs(x_O[2] - 0.0) < TOL);
        }
        { // move x_A in A
            double x_A[] = { -1.41421, 0.0, 0.0 };
            km_A_O->move(x_A);
            CPPUNIT_ASSERT(fabs(x_A[0] - (-3.82843)) < TOL);
            CPPUNIT_ASSERT(fabs(x_A[1] - 1.0) < TOL);
            CPPUNIT_ASSERT(fabs(x_A[2] - 0.0) < TOL);
        }

        delete km_O_A;
        delete km_A_O;
    }

    /***********************************************************************************************
     * \brief Test adding motion
     ***********/
    void testAddMotion() {
        // test add transformation
        const double TOL = 1E-3;

        double point[] = { 1.0, 0.0, 0.0 };
        KinematicMotion *km1 = new KinematicMotion();
        double zAxis[] = { 0.0, 0.0, 1.0 };
        km1->addRotation(zAxis, true, M_PI / 3.0);
        double translationVec1[] = { 0.0, -0.866, 0.0 };
        km1->addTranslation(translationVec1);

        KinematicMotion *km2 = new KinematicMotion();
        km2->addRotation(zAxis, true, M_PI / 3.0);
        double translationVec2[] = { 0.0, -0.433, 0.0 };
        km2->addTranslation(translationVec2);

        km1->addKinematicMotion(km2);
        km1->move(point);

        CPPUNIT_ASSERT(fabs(point[0] - 0.25) < TOL);
        CPPUNIT_ASSERT(fabs(point[1] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(point[2] - 0.0) < TOL);

        delete km1, km2;

    }
    /***********************************************************************************************
     * \brief Test adding rotation with different axises (non-commutative)
     ***********/
    void testAddRotation() {
        const double TOL = 1E-3;

        double point[] = { 1.0, 0.0, 0.0 };

        KinematicMotion *km = new KinematicMotion();

        double zAxis[] = { 0.0, 0.0, 1.0 };
        km->addRotation(zAxis, true, M_PI / 2.0);

        double xAxis[] = { 1.0, 0.0, 0.0 };
        km->addRotation(xAxis, true, M_PI / 2.0);

        km->move(point);
        //KinematicMotion::printVector(point);
        CPPUNIT_ASSERT(fabs(point[0] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(point[1] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(point[2] - 1.0) < TOL);

        delete km;
    }

    /***********************************************************************************************
     * \brief Test computation of rotation vector according to rotation matrix
     ***********/
    void testRotationVector() {
        const double TOL = 1E-10;
        { // a general angle
            KinematicMotion rot1;
            KinematicMotion rot2;

            double xAxis[] = { 1.0, 0.0, 0.0 };
            rot1.addRotation(xAxis, true, M_PI / 6.0);

            double zAxis[] = { 0.0, 0.0, 1.0 };
            rot1.addRotation(zAxis, true, M_PI / 3.0);

            //rot1.print();

            double rotationVector[3];
            rot1.getRotationVector(rotationVector);

            double angle = sqrt(
                    rotationVector[0] * rotationVector[0] + rotationVector[1] * rotationVector[1]
                            + rotationVector[2] * rotationVector[2]);
            if (angle < 1E-20) // 1E-20 is a physical "small" length, or a physical small angle
                angle = 0.0;
            for (int j = 0; j < 3; j++)
                rotationVector[j] /= angle;

            rot2.addRotation(rotationVector, true, angle);

            //rot2.print();

            for (int i = 0; i < 9; i++) {
                CPPUNIT_ASSERT(
                        fabs(rot1.getRotationMatrix()[i] - rot2.getRotationMatrix()[i]) < TOL);
            }
        }
        { // angle 0 and pi
            KinematicMotion rot;
            double rotVec[3];

            double xAxis[] = { 1.0, 0.0, 0.0 };
            rot.addRotation(xAxis, true, 1E-20);
            rot.getRotationVector(rotVec);

            for (int i = 0; i < 3; i++) {
                CPPUNIT_ASSERT(rotVec[i] == 0.0);
            }

            rot.resetRotation();
            rot.addRotation(xAxis, true, M_PI - 1E-20);
            rot.getRotationVector(rotVec);
            CPPUNIT_ASSERT(fabs(rotVec[0] - M_PI) < TOL);
            CPPUNIT_ASSERT(fabs(rotVec[1] - 0.0) < TOL);
            CPPUNIT_ASSERT(fabs(rotVec[2] - 0.0) < TOL);

            double tmpAxis1[] = { 1.0, 1.0, 0.0 };
            rot.resetRotation();
            rot.addRotation(tmpAxis1, false, M_PI - 1E-20);
            rot.getRotationVector(rotVec);
            CPPUNIT_ASSERT(fabs(rotVec[0] - sqrt(2.0)/2.0*M_PI) < TOL);
            CPPUNIT_ASSERT(fabs(rotVec[1] - sqrt(2.0)/2.0*M_PI) < TOL);
            CPPUNIT_ASSERT(fabs(rotVec[2] - 0.0) < TOL);

            double tmpAxis2[] = { 1.0, 1.0, 1.0 };
            rot.resetRotation();
            rot.addRotation(tmpAxis2, false, M_PI - 1E-20);
            rot.getRotationVector(rotVec);
            CPPUNIT_ASSERT(fabs(rotVec[0] - sqrt(3.0)/3.0*M_PI) < TOL);
            CPPUNIT_ASSERT(fabs(rotVec[1] - sqrt(3.0)/3.0*M_PI) < TOL);
            CPPUNIT_ASSERT(fabs(rotVec[2] - sqrt(3.0)/3.0*M_PI) < TOL);
        }
    }

    CPPUNIT_TEST_SUITE (TestKinematicMotion);
    CPPUNIT_TEST (testConstructor);
    CPPUNIT_TEST (testMotionAndCoorTransformation);
    CPPUNIT_TEST (testAddMotion);
    CPPUNIT_TEST (testAddRotation);
    CPPUNIT_TEST (testRotationVector);CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestKinematicMotion);

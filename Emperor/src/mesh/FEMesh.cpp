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
#include <iostream>
#include <map>
#include "FEMesh.h"
#include "Message.h"
#include "DataField.h"
#include "TriangulatorAdaptor.h"

namespace EMPIRE {

using namespace std;

FEMesh::FEMesh(std::string _name, int _numNodes, int _numElems, bool _triangulateAll) :
        AbstractMesh(_name), numNodes(_numNodes), numElems(_numElems) {
    type = EMPIRE_Mesh_FEMesh;
    boundingBox.isComputed(false);
    nodes = new double[numNodes * 3];
    nodeIDs = new int[numNodes];
    numNodesPerElem = new int[numElems];
    elems = NULL;

    elemIDs = new int[numElems];
    for (int i = 0; i < numElems; i++)
        elemIDs[i] = i + 1; // set element ID by hand instead of receiving from the client, element ID starts from 1.

    tobeTriangulated = false;
    triangulateAll = _triangulateAll;
    triangulatedMesh = NULL;
}

FEMesh::~FEMesh() {
    delete[] nodes;
    delete[] nodeIDs;
    delete[] numNodesPerElem;
    delete[] elems;
    delete[] elemIDs;
    if (triangulatedMesh != NULL)
        delete triangulatedMesh;
}

void FEMesh::initElems() {
    assert(elems==NULL);
    elemsArraySize = 0;

    for (int i = 0; i < numElems; i++) {
        int numNodesThisElem = numNodesPerElem[i];
        elemsArraySize += numNodesThisElem;
    }
    elems = new int[elemsArraySize];
    for (int i = 0; i < numElems; i++) {
        int numNodesThisElem = numNodesPerElem[i];
        if (numNodesThisElem > 4) {
            tobeTriangulated = true;
            break;
        }
    }
    if (triangulateAll)
        tobeTriangulated = true;
}

void FEMesh::addDataField(string dataFieldName, EMPIRE_DataField_location location,
        EMPIRE_DataField_dimension dimension, EMPIRE_DataField_typeOfQuantity typeOfQuantity) {
    int numLocatiions = -1;
    if (location == EMPIRE_DataField_atNode)
        numLocatiions = numNodes;
    if (location == EMPIRE_DataField_atElemCentroid)
        numLocatiions = numElems;
    assert(nameToDataFieldMap.find(dataFieldName) == nameToDataFieldMap.end());
    DataField *dataField = new DataField(dataFieldName, location, numLocatiions, dimension,
            typeOfQuantity);
    nameToDataFieldMap.insert(pair<string, DataField*>(dataFieldName, dataField));
}

FEMesh *FEMesh::triangulate() {
    if (!tobeTriangulated)
        return NULL;
    if (triangulatedMesh != NULL)
        return triangulatedMesh;

    map<int, int> *nodeIDToNodePosMap = new map<int, int>;
    for (int i = 0; i < numNodes; i++)
        nodeIDToNodePosMap->insert(nodeIDToNodePosMap->end(), pair<int, int>(nodeIDs[i], i));

    vector<int> *elemsTri = new vector<int>;
    vector<int> *numNodesPerElemTri = new vector<int>;
    int count = 0;
	bool isCompletelyTriangulated;
    for (int i = 0; i < numElems; i++) {
        int numNodesThisElem = numNodesPerElem[i];
        if (numNodesThisElem == 3) { // put triangles in the new mesh
            numNodesPerElemTri->push_back(numNodesThisElem);
            for (int j = 0; j < numNodesThisElem; j++)
                elemsTri->push_back(elems[count + j]);
        } else if ((numNodesThisElem == 4) && (!triangulateAll)) { // if not triangulateAll, put quads in the new mesh
                numNodesPerElemTri->push_back(numNodesThisElem);
                for (int j = 0; j < numNodesThisElem; j++)
                    elemsTri->push_back(elems[count + j]);
        } else {
            // use third party triangulation algorithm
            TriangulatorAdaptor *triangulator = new TriangulatorAdaptor();
            for (int j = 0; j < numNodesThisElem; j++) {
                int nodeID = elems[count + j];
                int nodePos = nodeIDToNodePosMap->at(nodeID);
                double x = nodes[nodePos * 3 + 0];
                double y = nodes[nodePos * 3 + 1];
                double z = nodes[nodePos * 3 + 2];
                triangulator->addPoint(x, y, z);
            }
            int numTriangles = numNodesThisElem - 2;
            int triangleIndexes[numTriangles * 3];
            isCompletelyTriangulated = triangulator->triangulate(triangleIndexes);
            for (int j = 0; j < numTriangles; j++) {
                numNodesPerElemTri->push_back(3);
                for (int k = 0; k < 3; k++) {
                    int nodeID = elems[count + triangleIndexes[j * 3 + k]];
                    elemsTri->push_back(nodeID);
                }
            }
            delete triangulator;
        }
        count += numNodesThisElem;
    }
	assert(isCompletelyTriangulated==true);
    delete nodeIDToNodePosMap;
    assert(count == elemsArraySize);

    triangulatedMesh = new FEMesh(name + "_triangulated", numNodes, numNodesPerElemTri->size());
    for (int i = 0; i < numNodes; i++) {
        triangulatedMesh->nodeIDs[i] = nodeIDs[i];
    }
    for (int i = 0; i < numNodes * 3; i++) {
        triangulatedMesh->nodes[i] = nodes[i];
    }
    for (int i = 0; i < numNodesPerElemTri->size(); i++) {
        triangulatedMesh->numNodesPerElem[i] = numNodesPerElemTri->at(i);
    }
    triangulatedMesh->initElems();
    for (int i = 0; i < triangulatedMesh->elemsArraySize; i++) {
        triangulatedMesh->elems[i] = elemsTri->at(i);
    }

    delete numNodesPerElemTri;
    delete elemsTri;
    assert(triangulatedMesh->tobeTriangulated == false);
    return triangulatedMesh;
}

void FEMesh::computeBoundingBox() {
    if (boundingBox.isComputed())
        return;
    boundingBox[0] = nodes[0 * 3 + 0];
    boundingBox[1] = nodes[0 * 3 + 0];
    boundingBox[2] = nodes[0 * 3 + 1];
    boundingBox[3] = nodes[0 * 3 + 1];
    boundingBox[4] = nodes[0 * 3 + 2];
    boundingBox[5] = nodes[0 * 3 + 2];
    for (int i = 1; i < numNodes; i++) {
        double x = nodes[i * 3 + 0];
        double y = nodes[i * 3 + 1];
        double z = nodes[i * 3 + 2];
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

void FEMesh::validateMesh() {
    int count1 = 0;
    int count2 = 0;
    // Check that two nodes index are not identical in a same element
    for(int elem=0; elem < numElems; elem++) {
    	for(int node1=0; node1 < numNodesPerElem[elem]-1; node1++){
        	for(int node2=node1+1; node2 < numNodesPerElem[elem]; node2++){
        		bool identical=elems[count1+node1]==elems[count1+node2];
        		if(identical)
        			ERROR_BLOCK_OUT("FEMesh","validateMesh","Mesh not valid. Duplicated node in element");
        		//assert(!identical);
        	}
    	}
    	count1+=numNodesPerElem[elem];
    }
    // Check that two elements are not identical
    // !WARNING! Check only element together in the same order
	/* count1=0;
    for(int elem1=0;elem1<numElems-1;elem1++) {
    	count2=count1+numNodesPerElem[elem1];
    	for(int elem2=elem1+1;elem2<numElems;elem2++) {
    		if(numNodesPerElem[elem1]==numNodesPerElem[elem2]) {
    			bool different=true;
    			for(int k=0;k<numNodesPerElem[elem1];k++) {
            		different = (elems[count1+k]!=elems[count2+k]);
            		if(different) break;
    			}
    			if(!different)
        			ERROR_BLOCK_OUT("FEMesh","validateMesh","Mesh not valid. Two strictly identical element present");
    			//assert(different);
    		}
    		count2 += numNodesPerElem[elem2];
    	}
    	count1 += numNodesPerElem[elem1];
    } */ // commented by Tianyang, since it is an O(N^2) operation.
    //TODO: move this function to the place where logically it requires element identity check, no place else should call this function.
}

void revertSurfaceNormalOfFEMesh(FEMesh *mesh) {
    int count = 0;
    for (int i = 0; i < mesh->numElems; i++) {
        int numNodesThisElem = mesh->numNodesPerElem[i];
        int tmpElem[numNodesThisElem];
        int *elems = &(mesh->elems[count]);
        for (int j = 0; j < numNodesThisElem; j++)
            tmpElem[j] = elems[numNodesThisElem - j - 1];
        for (int j = 0; j < numNodesThisElem; j++)
            elems[j] = tmpElem[j];
        count += numNodesThisElem;
    }
}

Message &operator<<(Message &message, FEMesh &mesh) {
	// First validate the mesh
	mesh.validateMesh();
	// Display the mesh information
    message << "\t+" << "FEMesh name: " << mesh.name << endl;
    message << "\t\t+" << "no. of nodes: " << mesh.numNodes << endl;
    message << "\t\t+" << "no. of elements: " << mesh.numElems << endl;
    message << "\t\t+" << "Nodes:" << endl;
    for (int i = 0; i < mesh.numNodes; i++) {
        message << "\t\t+" << "\t" << mesh.nodeIDs[i];
        for (int j = 0; j < 3; j++) {
            message << "\t" << mesh.nodes[i * 3 + j];
        }
        message << endl;
    }
    message << "\t\t+" << "Elements:" << endl;
    int count = 0;
    for (int i = 0; i < mesh.numElems; i++) {
        message << "\t\t+" << '\t' << i + 1;
        for (int j = 0; j < mesh.numNodesPerElem[i]; j++) {
            message << "\t" << mesh.elems[count + j];
        }
        count += mesh.numNodesPerElem[i];
        message << endl;
    }
    message() << "\t+" << "---------------------------------" << endl;

    return message;
}

} /* namespace EMPIRE */

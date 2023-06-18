/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   JointProbability.cpp
 * Copyright (C) 2018
 * Jenny Farmer jfarmer6@uncc.edu
 * Donald Jacobs djacobs1@uncc.edu
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in 
 * the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with 
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "JointProbability.h"

JointProbability::JointProbability(vector <Variable> variables, int nSamples, int pdfSize) {
    this->variables = variables;
    this->pdfSize = pdfSize;
    nVariables = variables.size();
    this->nSamples = nSamples;
    
    int gridPoints = pdfSize + 1;    
    gridSize = pdfSize;    
    
    gridPointsCourse.reserve(gridPoints);
    double dx = 1.0 / gridSize;
    double sum = 0;
    for (int i = 0; i < gridSize; i++) {
        gridPointsCourse.push_back(sum);
        sum += dx;
    }      
    gridPointsCourse.push_back(1.0);
      
    gridPointsFine.reserve(pdfSize);
    dx = 1.0 / (pdfSize - 1);
    sum = 0;
    for (int i = 0; i < (pdfSize - 1); i++) {
        gridPointsFine.push_back(sum);
        sum += dx;
    }
    gridPointsFine.push_back(1.0);
       
    if (nVariables > 1) {
        int v1 = nVariables - 1;
        int nLoops = pow(gridSize, v1);
        vector < vector < int > > grid(nLoops, vector<int> (v1));
        grids = grid;
        for (int v = 0; v < v1; v++) {
            int repeats = pow(gridSize, v);
            int loopIndex = 0;
            while (loopIndex < nLoops) {
                for (int g = 0; g < gridSize; g++) {
                    for (int i = 0; i < repeats; i++) {
                        grids[loopIndex++][v] = g;
                    }
                }
            }
        }        
    }    
}

JointProbability::~JointProbability() {
    delete [] jointPDF;
}
vector <double> JointProbability::getJP() {
    vector <double> VjointPDF;
    VjointPDF.reserve(matrixSize);
    for (int i = 0; i < matrixSize; i++) {
        VjointPDF.push_back(jointPDF[i]);
    }
    return VjointPDF;
}
        
vector <double> JointProbability::getRange() {
    vector <double> range;
    range.reserve(nVariables * pdfSize);
    for (int v = 0; v < nVariables; v++) {
        for (int i = 0; i < pdfSize; i++) {
            if (nVariables > 2) {                
                range.push_back(variables[v].meanSampleGrid[i]);
            } else {
                range.push_back(variables[v].xPDF[i]);
            }
        }
    }
    return range;
}


void JointProbability::calculate() {
    int nVariables = variables.size();
    matrixSize = pow(pdfSize, nVariables);
    jointPDF = new double[matrixSize];
    vector <double> pdf(pdfSize);
    
    for (int v = 0; v < nVariables; v++) {
        variables[v].calculateMarginals(gridPointsCourse, gridPointsFine);    
    }         
     
    int * repeats = new int[nVariables];
    repeats[nVariables - 1] = 1;
    for(int v = (nVariables - 2); v >= 0; v--) {
        repeats[v] = repeats[v + 1] * pdfSize;
    }
        
    pdf = variables[0].pdf;
    int repeat = repeats[0];
    for (int i = 0; i < pdfSize; i++) {
        for (int j = 0; j < repeat; j++) {
            jointPDF[i * repeat + j] = pdf[i];            
        }        
    }
          
    int iMatrix;
    
    for (int v = 1; v < nVariables; v++) {
        int nLoops = pow(gridSize, v);
        int nRepeats = repeats[v];
        iMatrix = 0;
        for (int iloop = 0; iloop < nLoops; iloop++) {
            int idx0 = grids[iloop][v - 1];
            vector <int> rows = variables[0].indices[idx0];
            for (int r = (v - 2); r >= 0; r--) {
                int idx = grids[iloop][r];
                rows = rowsIntersect(rows, variables[r + 1].indices[idx]);
            }
            pdf = variables[v].calculatePDF(rows);
            for (int i = 0; i < pdfSize; i++) {
                for (int iRepeats = 0; iRepeats < nRepeats; iRepeats++) {
                    jointPDF[iMatrix++] *= pdf[i];   
                }
            }            
        }
    }
           
    delete [] repeats;
}                    

vector <int> JointProbability::rowsIntersect(vector <int> v1, vector <int> v2) { 
    vector<int>::iterator it;
    int maxSize = v1.size() + v2.size();
    vector <int> v(maxSize);
    it = set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), v.begin()); 
    v.resize(it - v.begin());
    return v;
}


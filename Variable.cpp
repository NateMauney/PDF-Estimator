/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   Variable.cpp
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

#include "Variable.h"

Variable::Variable(InputParameters &input, vector <double> sample, string name, bool writeMarginal) {
    this->input = input;
    this->sample = sample;
    this->name = name;       
    this->writeMarginal = writeMarginal;
}

Variable::~Variable() {
}


void Variable::calculateMarginals(vector <double> gridPointsCourse, vector <double> gridPointsFine) {   
    
    InputData data = InputData(input);
    data.out.debug = input.debug;    
    data.setData(sample);
    data.processData();
    
    MinimizeScore minimumPDF = MinimizeScore();
    minimumPDF.out.debug = input.debug;
    minimumPDF.minimize(input, data); 
    WriteResults write;
    if (writeMarginal) {
        string file = "MarginalPDF_" + name;
        write.writeSolution(input, data, minimumPDF, file);   
    } else {
        write.createSolution(input, data, minimumPDF);
    }
        
    vector <double> samplePoints = interpolateGrid(write.CDF, write.x, gridPointsCourse);
    input.estimatePoints = true;
    input.estimatedPoints = samplePoints;    
    write.createSolution(input, data, minimumPDF);
    xGrid = write.xPoints;
    
    
    vector <double> pdfPoints = interpolateGrid(write.CDF, write.x, gridPointsFine);
    input.estimatedPoints = pdfPoints;    
    write.createSolution(input, data, minimumPDF);
    xPDF = write.xPoints;
    pdf = write.PDFPoints;   
    pdfSize = xPDF.size();
    
    calculateIndices(samplePoints);
    
    int nSamplePoints = samplePoints.size();
    meanSampleGrid.reserve(nSamplePoints - 1);
    meanSampleGrid.push_back((samplePoints[0] + samplePoints[1]) / 2);
    for (int i = 1; i < (nSamplePoints - 1); i++) {
        meanSampleGrid.push_back((samplePoints[i] + samplePoints[i + 1]) / 2);
    }    
    
}

vector <double> Variable::calculatePDF(vector <int> index) {
    int n = index.size();
    if (n < 5) {
        vector <double> pzero(pdfSize);
        return pzero;
    }
    vector <double> subsample;
    subsample.reserve(n);
    for (int i = 0; i < n; i++) {
        subsample.push_back(sample[index[i]]);        
    }        
    
    InputData data = InputData(input);
    data.out.debug = input.debug;    
    data.setData(subsample);
    data.processData();
    
    MinimizeScore minimumPDF = MinimizeScore();
    minimumPDF.out.debug = input.debug;
    minimumPDF.minimize(input, data);
    WriteResults write;
    write.createSolution(input, data, minimumPDF);         
    
    return write.PDFPoints;    
}

vector < vector < double > > Variable::interpolatePDF(int gridIndex1, int nTimes, int pdfIndex, vector <double> pdf0, vector <double> pdf1) {
    
    vector < vector < double > > interpolatePDFs(nTimes, vector <double> (pdfSize)); 
    int gridIndex0 = gridIndex1 - 1;
    
    for (int i = 0; i < nTimes; i++) {
        for (int j = 0; j < pdfSize; j++) {
            double value = pdf0[j] + (xPDF[pdfIndex] - meanSampleGrid[gridIndex0]) * 
                (pdf1[j] - pdf0[j]) / (meanSampleGrid[gridIndex1] - meanSampleGrid[gridIndex0]);
            if (value < 0) {
                value = 0;
            }
            interpolatePDFs[i][j] = value;
        }
        pdfIndex++;
    }
    
    return(interpolatePDFs);
}

vector <double> Variable::interpolateGrid(vector <double> x, vector <double> y, vector <double> gridPoints) {

    vector <double> interp;
    int n = x.size();    
    int nGrids = gridPoints.size();    
    interp.reserve(nGrids);
    
    interp.push_back(y[0]);
    x[0] = 0;
    x[n-1] = 1.0;
    int i = 1;   
    
    for (int gridIndex = 1; gridIndex < (nGrids - 1); gridIndex++) {
        while (x[i] < gridPoints[gridIndex]) {
            if (++i == (n - 1)){
                break;
            } 
        }
        double p = y[i - 1] + (gridPoints[gridIndex] - x[i - 1]) * (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        interp.push_back(p);
    }
    
    interp.push_back(y[n-1]);
    return(interp);
}

void Variable::calculateIndices(vector <double> samplePoints) {
    int n = sample.size();
    int g = xGrid.size();
    vector <int> index;
    vector <int> v;
    for (int i = 0; i < n; i++) {
        for (int j = 1; j < g; j++) {
            if (sample[i] < samplePoints[j]) {
                index.push_back(j - 1);
                break;
            }
            if (j == (g - 1)) {
                index.push_back(-1);
            }
        }
    }
    for (int j = 0; j < g; j++) {
        for (int i = 0; i < n; i++) {
            if (index[i] == j) {
                v.push_back(i);
            }
        }        
        indices.push_back(v);
        v.clear();
    }
}

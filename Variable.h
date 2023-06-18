/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   Variable.h
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

#ifndef VARIABLE_H
#define VARIABLE_H

#include "InputParameters.h"
#include "ScoreQZ.h"
#include "InputData.h"
#include "MinimizeScore.h"
#include "WriteResults.h"

class Variable {
public:
    
    Variable(InputParameters& input, vector <double> sample, string name, bool writeMarginal);
    virtual ~Variable();
    
    void calculateMarginals(vector <double> gridPointsCourse, vector <double> gridPointsFine);
    vector <double> calculatePDF(vector <int> index);
    vector < vector < double > > interpolatePDF(int gridIndex0, int nTimes, int pdfIndex, vector <double> pdf0, vector <double> pdf1);
   
    string name;
        
    vector <double> xPDF;
    vector <double> pdf;
        
    vector <double> meanSampleGrid;
    
    vector < vector < int > > indices;
        
    OutputControl out;
private:
    bool writeMarginal;
    vector <double> sample;
    InputParameters input;
    
    vector <double> xGrid;
    int pdfSize;
    
    vector <double> interpolateGrid(vector <double> x, vector <double> y, vector <double> gridPoints);
    void calculateIndices(vector <double> samplePoints);
    
};

#endif /* VARIABLE_H */


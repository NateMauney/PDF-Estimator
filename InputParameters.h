/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   InputParameters.h
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

#ifndef INPUTPARAMETERS_HPP
#define	INPUTPARAMETERS_HPP

#include <string>
#include <unistd.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include "OutputControl.h"
using namespace std;

class InputParameters {
 
public:     
    
    string debugOpt;
    bool   debug;
   
    string inputPath;
    string inputFile;
    string outputFile;
    string outputPath;
    bool    writeFile;
    bool    writeHeader;
    bool    writeFailed;
    string  writeOpt;
    string  headerOpt;
    bool    writeQQ;
    bool    writeSQR;
    string  qqFile;
    string  sqrFile;
    bool    adaptive;
        
    float   lowerBound;
    float   upperBound;
    bool    lowerBoundSpecified;
    bool    upperBoundSpecified;
    
    string  scoreType;
    double  SURDMinimum;
    double  SURDTarget;
    double  SURDMaximum;
    int     initPartitionSize;
    int     startSolutionNumber;
    int     integrationPoints;
    int     maxLagrange;
    int     minLagrange;
    int     nLagrangeAdd;
    double  outlierCutoff;
    bool    smooth;
    
    double  fractionLagrangeAdd;
    double  initSigma;
    double  finalSigma;
    double  decayFactor;
    int     loopMax;
    
    void setEstimationPoints(vector <double> x);    
    vector <double> estimatedPoints;
    bool estimatePoints;
        
    InputParameters();
//    InputParameters(const InputParameters& orig);
    virtual ~InputParameters();
    bool userInput(int argc, char** argv);    
    
    OutputControl out;
    
private:    
    void printUsage();
    
};

#endif	/* INPUTPARAMETERS_HPP */


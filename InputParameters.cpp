/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   InputParameters.cpp
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

#include "InputParameters.h"

InputParameters::InputParameters() {  
    
    debug = false;
   
    inputPath = "";
    inputFile = ".txt";
    outputPath = "";
    writeFile = true;
    writeHeader = true;
    writeQQ = false;
    writeSQR = false;
    writeFailed = true;
    qqFile = "";
    sqrFile = "";
    adaptive = true;
    
    lowerBoundSpecified = false;
    upperBoundSpecified = false;
    
    scoreType = "QZ";
    SURDMinimum = 5;
    SURDTarget  = 70;
    SURDMaximum = 100;
    initPartitionSize = 1025;
    startSolutionNumber = 0;
    integrationPoints = -1;
    maxLagrange = 200;//2 for power
    minLagrange = 1;
    nLagrangeAdd = 5;
    outlierCutoff = 7.0;
    smooth = true;
    
    fractionLagrangeAdd = 0.1;
    initSigma = 0.1;
    finalSigma =0.001;
    decayFactor = sqrt(2.0);
    loopMax = 100;                  // updated from 30; September 2020
    
    estimatePoints = false;
    
 }

//InputParameters::InputParameters(const InputParameters& orig) {
//}

InputParameters::~InputParameters() {
}


void InputParameters::setEstimationPoints(vector<double> x) {
    estimatedPoints.resize(x.size());
    estimatedPoints =  x;
    estimatePoints = true;
}

bool InputParameters::userInput(int argc, char**  argv){
    
    int c;
    bool inputEntered = false;
    
    while ((c = getopt(argc, argv, "f:o:w:h:q:r:l:u:s:p:n:m:z:a:b:t:c:d:e:x:g:v:")) != -1)
    switch (c){    
         case 'g':
            debugOpt = optarg;
            if (debugOpt == "on") {
                debug = true;
                out.debug = true;
                out.print("debug = on");
            }
            break;
        case 'f':
            inputFile = optarg;
            out.print("Input data file name = " + inputFile);
            inputEntered = true;
            break; 
        case 'o':
            outputFile = optarg;
            out.print("Output data file name = " + outputFile);
            break; 
        case 'a':
            inputPath = optarg;
            out.print("Input data path = " + inputPath);
            break; 
        case 'b':
            outputPath = optarg;
            out.print("Output data path = " + outputPath);
            break; 
        case 'w':
            writeOpt = optarg;
            if (writeOpt == "off") {
                writeFile = false;
                out.print("Write File = off");
            }
            break;
        case 'x':
            writeOpt = optarg;
            if (writeOpt == "off") {
                writeFailed = false;
                out.print("Write Failed Solutions = off");
            }
            break;        
        case 'h':
            headerOpt = optarg;
            if (headerOpt == "off") {
                writeHeader = false;
                out.print("Write File Header = off");
            }
            break;
        case 'q':
            qqFile = optarg;
            writeQQ = true;
            out.print("Write QQ File = " + qqFile);
            break;
        case 'r':
            sqrFile = optarg;
            writeSQR = true;
            out.print("Write SQR File = " + sqrFile);
            break;
        case 'l': 
            lowerBound = atof(optarg);
            out.print("lower bound = ", lowerBound);
            lowerBoundSpecified = true;
            break; 
        case 'u':
            upperBound = atof(optarg);
            out.print("upper bound = ", upperBound);
            upperBoundSpecified = true;
            break;  
        case 'v':
            scoreType = optarg;
            out.print("Scoring method = " + scoreType);
            break;       
        case 's':                                                               
            SURDTarget = atof(optarg);
            if (SURDTarget < 1) {
                out.print("WARNING: coverage must be between 1 and 100; setting to 1");
                SURDTarget = 1;
            } else if (SURDTarget > 100) {
                out.print("WARNING:  coverage must be between 1 and 100; setting to 100");
                SURDTarget = 100;
            } else {
                out.print("coverage = ", SURDTarget);
            }
            break;
        case 'd':                                                               
            SURDMinimum = atof(optarg);
            if (SURDMinimum < 1) {
                out.print("WARNING: coverage must be between 1 and 100; setting to 1");
                SURDMinimum = 1;
            } else if (SURDMinimum > 100) {
                out.print("WARNING:  coverage must be between 1 and 100; setting to 100");
                SURDMinimum = 100;
            } else {
                out.print("maximum coverage = ", SURDMinimum);
            }
            break;
        case 'e':                                                               
            SURDMaximum = atof(optarg);
            if (SURDMaximum < 1) {
                out.print("WARNING: coverage must be between 1 and 100; setting to 1");
                SURDMaximum = 1;
            } else if (SURDMaximum > 100) {
                out.print("WARNING:  coverage must be between 1 and 100; setting to 100");
                SURDMaximum = 100;
            } else {
                out.print("minimum coverage = ", SURDMaximum);
            }
            break;
        case 'p':                                                               
            integrationPoints = atoi(optarg);
            out.print("integration points = ", integrationPoints);
            break;
        case 'n':                                                               
            maxLagrange = atoi(optarg);
            out.print("maximum Lagrange = ", maxLagrange);
            break;       
        case 'm':                                                               
            minLagrange = atoi(optarg);
            out.print("minimum Lagrange = ", minLagrange);
            break;                
        default:
            out.print("Invalid parameter flag: ", c);
            printUsage();            
            return false;
    }
    if (!inputEntered) {
        out.print("Input file name required");
        printUsage();
        return false;
    }
   
    if (SURDMaximum < SURDMinimum) {
        out.print("maximum coverage cannot be less than minimum coverage values");
        return false;
    }
    if (SURDMaximum < SURDTarget) {
        out.print("maximum coverage cannot be less than target coverage values");
        return false;
    }
    if (SURDMinimum > SURDTarget) {
        out.print("minimum coverage cannot be greater than target coverage values");
        return false;
    }
    return true;
}

  


void InputParameters::printUsage() {
    out.print("Usage:");
    out.print("getpdf -f <filename> [-option <argument>]");
    
    out.print("Options:");
    out.print(" -f    input filename (REQUIRED)");
    out.print(" -o    main output filename");
    out.print( " -w    write main output file [on/off]");
    out.print( " -h    include header info in main output file [on/off]");
    out.print( " -q    QQ filename");
    out.print( " -r    SQR filename");
    out.print( " -l    lower bound");
    out.print( " -u    upper bound");
    out.print( " -s    score threshold percentage [1-100]");
    out.print( " -p    minimum number of integration points");
    out.print( " -n    maximum number of Lagrange multipliers");
    out.print( " -m    minimum number of Lagrange multipliers");
    out.print( " -y    penalty flag [on/off]");
    out.print( " -g    debug [on/off]");
}


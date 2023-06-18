/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   WriteResults.cpp
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

#include "WriteResults.h"

WriteResults::WriteResults() {
}


WriteResults::~WriteResults() {
}

void WriteResults::writeSolution(const InputParameters& input, const InputData& data, MinimizeScore& solution, 
                                 string fileNameAdd) {
    writeSolution(input, data, solution, ScoreQZ(), 0, fileNameAdd);
}

void WriteResults::writeSolution(const InputParameters& input, const InputData& data, MinimizeScore& solution, 
                                 Score score, bool failed, string fileNameAdd) {
        
    createSolution(input, data, solution);
    
    ofstream outFile;
    if (input.writeFile) {
        ostringstream surdString; 
        surdString << input.SURDTarget;  
        string filename;
        if (!(fileNameAdd == "")) {
            filename = fileNameAdd + "_" + input.inputFile;
        } else {
            if (input.outputFile == "") {
                if (failed) {
                    filename = "FAILED_PDF_" + surdString.str() + "_" + input.inputFile;
                } else {
                    filename = "PDF_" + surdString.str() + "_" + input.inputFile;
                }
            } else {
                filename = input.outputFile;
            }
        }
        filename = input.outputPath + filename;
        outFile.open((filename).c_str());
        if(!outFile.is_open()){
            out.print("Failed to open data file " + filename);
            return;
        }
    }
    
    if ((input.writeFile) && (input.writeHeader)) {
        if (failed) {
            outFile << "***** FAILED SOLUTION **********\n\n";
        }
#ifdef clock
        float time = solution.duration;
        outFile << "#   calculation time: " << time << " seconds\n";
        if (time >= 60) {
            outFile << "#                    (" << time/60 << " minutes)\n";
        }
        outFile << "#\n";
#endif          
        outFile << "# PARAMETERS\n";            
        outFile << "#\n";            
        outFile << "#   minimum confidence:           "  << input.SURDTarget << "%\n";    
        
        outFile << "#   minimum data value:           "  << data.minimumRaw << "\n";  
        outFile << "#   maximum data value:           "  << data.maximumRaw << "\n";   
            
        if (input.lowerBoundSpecified) {
            outFile << "#   requested left boundary:      " << input.lowerBound << "\n";                
        } else {
            outFile << "#   requested left boundary:      infinite\n" ;
        }
          
        if (input.upperBoundSpecified) {
            outFile << "#   requested right boundary:      " << input.upperBound << "\n";                
        } else {
            outFile << "#   requested right boundary:     infinite\n" ;
        }
                      
        outFile << "#   right outliers removed:       "  << data.nRightOutliers << "\n";  
        outFile << "#   left outliers removed:        "  << data.nLeftOutliers << "\n";        
        outFile << "#   recalculated left boundary:   "  << data.minimumCalc << "\n";                    
        outFile << "#   recalculated right boundary:  "  << data.maximumCalc << "\n";     
        outFile << "#   number of sample points used: "  << data.N << "\n";                                 
        outFile << "#   starting partition:           "  << input.initPartitionSize << "\n";  
        outFile << "#   max lagrange :                "  << input.maxLagrange << "\n";         
        outFile << "#   attempts per step size:       "  << input.loopMax << "\n";                                        
        outFile << "#   integration points:           "  << data.nPoints << "\n";      
        outFile << "#   integration points adjusted:  "  << data.nPointsAdjust << "\n";                
        outFile << "#\n";            
            
        outFile << "# SCORES\n";                    
        outFile << "#\n";               
//        outFile << "#   variance :         " << solution.bestScore << "\n";  
//        outFile << "#   QZ score:     " << score.getLikelihood() << "\n";
//        outFile << "#   SURD threshold:       " << score.getConfidence(score.getLikelihood()) << "%\n";
        outFile << "#\n";           
        outFile << "#\n";           
          
    }
    
    
    if ((input.writeFile) && (input.writeHeader)) {  
        outFile << "# LAGRANGE MULTIPLIERS (" << solution.mode << ")\n";            
        outFile << "#\n";            
        for (int j = 0; j < solution.mode; j++) {
            outFile << "#   " << L[j] << "\n";     
        }
        outFile << "#\n";            
        outFile << "# PLOT POINTS\n";            
        outFile << "#\n";                      
        
     }
    
    if (input.estimatePoints) {
        outFile << "#  x                     PDF(x)\n";   
        for (unsigned int k = 0; k < PDFPoints.size(); k++) {
            if (input.writeFile) {
                outFile << xPoints[k] << "   " <<  PDFPoints[k] << "\n";
            }
        }
    } else {
        outFile << "#  x                     PDF(x)                     CDF(x) \n";   
        for (unsigned int k = 0; k < PDF.size(); k++) {
            if (input.writeFile) {
                outFile << x[k] << "   " <<  PDF[k] << "   " << CDF[k] << "\n";
            }
        }
    }
    
    if (input.writeFile) {
        outFile.close();
    }
    
}


void WriteResults::createSolution(const InputParameters& input, const InputData& data, MinimizeScore& solution) {
       
    double max = data.maximumCalc;
    double min = data.minimumCalc;  
      
    double dzSize;
    double * dz;
    double * dzPoint;
    int nSize = input.estimatedPoints.size();
    dzPoint = new double[(int) nSize];
    if (input.estimatePoints) {        
        dzPoint[0] = (input.estimatedPoints[0] - min);
        dzPoint[nSize - 1] = (max - input.estimatedPoints[nSize - 1]);
        for (int k = 1; k < (nSize - 1); k++) {
            dzPoint[k] = input.estimatedPoints[k + 1] - input.estimatedPoints[k];
        }
    } 
    if (input.adaptive) {
        double * dr;
        dr = data.dz; 
        double dzUniform = (max - min);
        dzSize = (data.nPointsAdjust) * 2 - 2;
        dz = new double[(int) dzSize];
        for (int i = 0; i < dzSize; i++) {
            dz[i] = dr[i] * dzUniform;
        }
    } else {      
        dzSize = data.nPoints;     
        double dzUniform = (max - min) * 1.0/dzSize;     
        dzSize++;    
        dz = new double[(int) dzSize];      
        for (int i = 0; i < dzSize; i++) {
            dz[i] = dzUniform;
        }
    }
    
        
    vector <double> termsT;;         
    vector <double> termsP; 
    double p = 0;
    double termsSum = 0;     
    double z = 0;
    double q = min;
    vector <double> lagrange = solution.getLagrange();
    for (int k = 0; k < dzSize; k++) {
        z = (2*q - max - min) / (max - min);           
        termsT.clear();
        termsT.push_back(1.0);
        termsT.push_back(z);        
        p = lagrange[0];        
        for (int t = 1; t < solution.mode; t++) {                
            p += termsT[t] * lagrange[t];
            termsT.push_back(2 * z * termsT[t] - termsT[t-1]);
        }  
        p = exp(p);
        p /= (max - min)/2;        
        termsP.push_back(p);
        termsSum += p * dz[k];
        q += dz[k];          
    }    
    double lambdaZero =  -log(termsSum);    
    double normConstant = termsSum;
    
    L.push_back(lambdaZero);
    for (int j = 1; j < solution.mode; j++) {   
        L.push_back(lagrange[j]);
//        out.print("Lagrange   ", L[j]);
    }
    
    
    if (input.estimatePoints) { 
        xPoints.clear();
        PDFPoints.clear();
        int nEstimate = input.estimatedPoints.size();
        for (int k = 0; k < nEstimate; k++) {
            q = input.estimatedPoints[k];
            if ((q < min) || (q > max)) {
                p = 0;
            } else {
                z = (2*q - max - min) / (max - min);    
                termsT.clear();
                termsT.push_back(1.0);
                termsT.push_back(z);        
                p = lagrange[0];
                for (int t = 1; t < solution.mode; t++) {                
                    p += termsT[t] * lagrange[t];
                    termsT.push_back(2 * z * termsT[t] - termsT[t-1]);
                }  
                p = exp(p);
                p /= (max - min)/2;   
                p /= normConstant;
            }
            xPoints.push_back(q);
            PDFPoints.push_back(p);
        }
    } 
         
    termsSum = 0;
    q = min;
    
    x.clear();
    CDF.clear();
    PDF.clear();
    
    for (int k = 0; k < dzSize; k++) {
        double pk = termsP[k] / normConstant;
        termsSum += pk * dz[k];
        CDF.push_back(termsSum);  
        x.push_back(q);
        PDF.push_back(pk);                
        q += dz[k];
    }             
        
    delete [] dz;
    delete [] dzPoint;
}



void WriteResults::writeColumn(string filename, vector <double> r, int length) {   
    ofstream outFile;
    outFile.open(filename.c_str());
    for (int i = 0; i < length; i++) {
        outFile << r[i] <<  "\n";
    }
    outFile.close(); 
}

void WriteResults::writeColumn(string filename, vector <int> r, int length) {   
    ofstream outFile;
    outFile.open(filename.c_str());
    for (int i = 0; i < length; i++) {
        outFile << r[i] <<  "\n";
    }
    outFile.close(); 
}

void WriteResults::writeColumn(string filename, double r[], int length) {
    ofstream outFile;
    outFile.open(filename.c_str());
    for (int i = 0; i < length; i++) {
        outFile << r[i] <<  "\n";
    }
    outFile.close();
}

void WriteResults::writeColumn(string filename, int r[], int length) {
    ofstream outFile;
    outFile.open(filename.c_str());
    for (int i = 0; i < length; i++) {
        outFile << r[i] <<  "\n";
    }
    outFile.close();
}

void WriteResults::writeQQ(string filename, double r[], int length, bool sqr) {
    ofstream outFile;
    outFile.open(filename.c_str());
    if(!outFile.is_open()){
        out.print("Failed to open data file " + filename);
        return;
    }
   
    for (int i = 0; i < length; i++) {
        double position = (i + 1) * 1.0 / (length + 1);
        if (sqr) {
            SQR.push_back(sqrt(length + 2) * (r[i] - position));
            outFile << setprecision(12) << r[i] << " " << position << " " << SQR[i] <<  "\n";
        } else {
            outFile << setprecision(12) << position << " " << r[i] <<  "\n";
        }
    }    
    outFile.close();
}

void WriteResults::createQQ(double r[], int length) {
   
    for (int i = 0; i < length; i++) {
        double position = (i + 1) * 1.0 / (length + 1);
        SQR.push_back(sqrt(length + 2) * (r[i] - position));
    }    
}
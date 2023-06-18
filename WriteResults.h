/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   WriteResults.h
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

#ifndef WRITERESULTS_HPP
#define	WRITERESULTS_HPP

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>      
#include <time.h>
#include "MinimizeScore.h"
#include "OutputControl.h"

//#define clock

using namespace std;

class WriteResults {
public:    
    WriteResults();
    WriteResults(const WriteResults& orig);
    virtual ~WriteResults();
   
    void writeSolution(const InputParameters& input, const InputData& data, MinimizeScore& solution, 
                       Score score, bool failed, string fileNameAdd);
    void writeSolution(const InputParameters& input, const InputData& data, MinimizeScore& solution, 
                       string fileNameAdd);
    void createSolution(const InputParameters& input, const InputData& data, MinimizeScore& solution);
    
    void writeColumn(string filename, double r[], int length);
    void writeColumn(string filename, int r[], int length);
    void writeColumn(string filename, vector <double> r, int length);
    void writeColumn(string filename, vector <int> r, int length);
    
    void writeQQ(string filename, double r[], int length, bool sqr);
    void createQQ(double r[], int length);
    
    vector <double> x;
    vector <double> PDF;   
    vector <double> CDF;   
    vector <double> xPoints;
    vector <double> PDFPoints;  
    vector <double> SQR;
    vector <double> L;
    vector <double> R;
    
    OutputControl out;

};

#endif	/* WRITERESULTS_HPP */


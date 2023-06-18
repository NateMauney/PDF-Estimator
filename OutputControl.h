/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   OutputControl.h
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

#ifndef OUTPUTCONTROL_H
#define	OUTPUTCONTROL_H

#define outputCommandLine
//#define outputMatlab
//#define outputR


#include <string>

#ifdef outputMatlab
#include "cppmex/mexMatlabEngine.hpp"
#include "MatlabDataArray/ArrayFactory.hpp"
#endif

#ifdef outputR
#include "R_ext/Print.h"
#endif

#ifdef outputCommandLine
#include <iostream>
#endif

using namespace std;
class OutputControl {
public:
    OutputControl();
    OutputControl(const OutputControl& orig);
    virtual ~OutputControl();    
    bool debug;
    
#ifdef outputMatlab    
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr; 
    
    void displayError(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr, std::string errorMessage);
    void displayOnMATLAB(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr, std::string message);  
    void setPtr(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr); 
#endif
    
    
void print(string output);
void print(string output, int value);
void print(string output, double value);
void error(string output);
void error(string output, int value);
void error(string output, double value);       

};

#endif	/* OUTPUTCONTROL_H */



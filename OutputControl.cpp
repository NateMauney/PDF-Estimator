/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   OutputControl.cpp
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


#include "OutputControl.h"



OutputControl::OutputControl() {
    debug = false;
}

OutputControl::OutputControl(const OutputControl& orig) {
    debug = false;
}

OutputControl::~OutputControl() {
}

#ifdef outputMatlab
void OutputControl::displayError(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr, std::string errorMessage) {    
    matlab::data::ArrayFactory factory;
    matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"), 0, std::vector<matlab::data::Array>({
    factory.createScalar(errorMessage) }));
}

void OutputControl::displayOnMATLAB(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr, std::string message) {
    if (debug) {
        matlab::data::ArrayFactory factory;
        matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("fprintf"),0, std::vector<matlab::data::Array>
            ({ factory.createScalar(message)}));
        matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("fprintf"),0, std::vector<matlab::data::Array>
            ({ factory.createScalar("\n")}));
    }
}


void OutputControl::setPtr(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr) {
    this->matlabPtr = matlabPtr; 
}

 void OutputControl::print(string output) { 
     displayOnMATLAB(matlabPtr, output); 
 }
 
 void OutputControl::error(string output) { 
     displayError(matlabPtr, output); 
 }

void OutputControl::print(string output, int value) {
    ostringstream strOut;
    strOut << output << ": " << value;                                      
    print(strOut.str());
} 

void OutputControl::print(string output, double value) {
    ostringstream strOut;
    strOut << output << ": " << value;                                      
    print(strOut.str());
}

void OutputControl::error(string output, int value)  {
    ostringstream strOut;
    strOut << output << ": " << value;                                      
    error(strOut.str());
}

void OutputControl::error(string output, double value)  {
    ostringstream strOut;
    strOut << output << ": " << value;                                      
    error(strOut.str());
}
#endif


#ifdef outputR


void OutputControl::error(string output) {
    REprintf("%s\n", output.c_str());
}

void OutputControl::error(string output, int value) {
    REprintf("%s: %d\n", output.c_str(), value);
}

void OutputControl::error(string output, double value) {
    REprintf("%s: %f\n", output.c_str(), value);
}


void OutputControl::print(string output) {
    if (debug) {
        Rprintf("%s\n", output.c_str());
    }
}

void OutputControl::print(string output, int value) {
    if (debug) {
        Rprintf("%s: %d\n", output.c_str(), value);
    }
}

void OutputControl::print(string output, double value) {
    if (debug) {
        Rprintf("%s: %f\n", output.c_str(), value);
    }
}

#endif

#ifdef outputCommandLine

void OutputControl::print(string output) {
    if (debug) {
        cout << output << "\n";
    }
}

void OutputControl::print(string output, int value) {
    if (debug) {
        cout << output << ": " << value << "\n";
    }
}

void OutputControl::print(string output, double value) {
    if (debug) {
        cout << output << ": " << value << "\n";
    }
}

void OutputControl::error(string output) {
    cout << output << "\n";
}

void OutputControl::error(string output, int value) {
   cout << output << ": " << value << "\n";
}

void OutputControl::error(string output, double value) {
    cout << output << ": " << value << "\n";
}

#endif
    
   
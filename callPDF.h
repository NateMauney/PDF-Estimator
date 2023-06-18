/* 
 * File:   callPDF.hpp
 * Author: jenny
 *
 * Created on September 18, 2018, 1:20 PM
 */

#include <vector>
#include "InputParameters.h"
#include "InputData.h"
#include "ScoreQZ.h"
#include "MinimizeScore.h"
#include "WriteResults.h"


#ifndef CALLPDF_HPP
#define	CALLPDF_HPP

class callPDF {
public:
    callPDF();
    callPDF(const callPDF& orig);
    virtual ~callPDF();
    void makeCall(double * sampleData, int sampleLength, double * estimationPoints, int estimationLength, int isSpecifyPoints, double low, double high, int isLow, int isHigh, double target, int points, int lagrangeMin, int lagrangeMax, int outlierCutoff, int debug, int smooth);
    
    vector <double> Vcdf;
    vector <double> Vpdf;
    vector <double> Vx; 
    vector <double> VpdfPoints;
    vector <double> Vsqr;
    vector <double> Vlagrange;
    vector <double> Vr;
    
    double N;
    
    bool solutionFailed;
    
    double solutionThreshold;
        
private:
    OutputControl out;
};

#endif	/* CALLPDF_HPP */


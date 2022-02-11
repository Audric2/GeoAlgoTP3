#ifndef OTSU_HISTOGRAME_HPP
#define OTSU_HISTOGRAME_HPP

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double otsu(std::vector<int> histogramme){
    int size = histogramme.size();
    int sum1 = 0,sum0 = 0,w1 = 0,w0 = 0;
    for(int i = 0;i<size;i++){
        w1 += i * histogramme[i];
        sum1 += histogramme[i];
    }
    	
    double maxVar = 0;
    int maxVarindex = 0;
	for(int t=0;t<size;t++){
        w0 += t * histogramme[t];
		w1 -= t * histogramme[t];
		sum0 += histogramme[t];
		sum1 -= histogramme[t];
		if(histogramme[t]!=0 && sum0 != 0 && sum1 != 0){
		    // Calcul de Wikipedia 
            const double moy0 = w0/sum0;
            const double moy1 = w1/sum1;
            const double w0s = (double)sum0/(double)(sum0+sum1);
            const double w1s = (double)sum1/(double)(sum0+sum1);
            const double varianceInter = w0s * w1s * pow(moy0-moy1,2);
            // On recupere le maximum si il a change
            if(varianceInter > maxVar){
                maxVar = varianceInter;
                maxVarindex = t;
            }
        }
    }
    return maxVarindex;
}

#endif


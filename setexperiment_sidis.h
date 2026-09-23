#ifndef _SET_EXPERIMENT_SIDIS__
#define _SET_EXPERIMENT_SIDIS__

#include <sidis.h>
#include <common_sidis.h>

#ifndef  __CONSTANTS_H_
#include "constants.h"
#endif

#ifndef  __HADRON_H_
#include "hadron.h"
#endif

#ifndef  __PARAMETERS_H_
#include "parameters.h"
#endif

//void SetExperiment_SIDIS(int Experiment);
int SetExperiment_SIDIS(const char * argv);

int SetBeam_SIDIS(const char * argv);
int SetTarget_SIDIS(const char * argv);
int SetDependence_SIDIS(const char * argv);
int SetHadron_SIDIS(const char * arv);


#endif // #ifndef _SET_EXPERIMENT_SIDIS__

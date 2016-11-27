#ifndef PARCOURS_LOGADD_H
#define PARCOURS_LOGADD_H

#define logUnderflowThreshold 7.5
#define posteriorMatchThreshold 0.01

#define LOG_ZERO (double )-std::numeric_limits<double>::infinity()

double LogAdd(double x, double y);

#endif

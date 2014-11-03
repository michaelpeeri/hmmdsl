//-------------------------------------------------------------------------------------
// Header for asa121.cpp
//
// Obtained from:  http://people.sc.fsu.edu/~burkardt/cpp_src/asa121/asa121.C
//  On April 15th, 2010
//
// This file (and the website it comes from) don't contain any information about
// its licensing. It is probably in the public domain.
//
// The author is John Burkardt, closely based on the publication by Schneider.
// He (Burkardt) is working at Virginia Tech. The website has been updated recently,
// but there is no contact information.
//
// Accuracy: According to the paper, should be at least 6 digits.
//
//-------------------------------------------------------------------------------------
void timestamp ( void );
double trigam ( double x, int *ifault );
void trigamma_values ( int *n_data, double *x, double *fx );

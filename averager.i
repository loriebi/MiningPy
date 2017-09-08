/* File : averager.i */
%module averager
%include "carrays.i"
%include "std_vector.i"

%template(vecFl) std::vector<float>;
%template(vecDb) std::vector<double>;
%template(vecUint) std::vector<unsigned int>;

%{
#include "Averager.h"
%}


/* %array_class(float, floatArray); */

/* Let's just grab the original header file here */

%include "Averager.h"


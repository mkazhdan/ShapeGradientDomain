#ifndef PRE_PROCESSOR_INCLUDED
#define PRE_PROCESSOR_INCLUDED

#if _WIN32 || _WIN64
#define NOMINMAX
#endif // _WIN32 || _WIN64

#define NEW_CODE				// Experimental code

//#define EIGEN_USE_MKL_ALL		// Accelerate Eigen performance by using the MKL


#endif // PRE_PROCESSOR_INCLUDED
// Demonstrates how to use time.h to return cputime

#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif
	double cputime_();
#ifdef __cplusplus
}
#endif

double cputime_() 
{
	return  (double) clock() / (double) CLOCKS_PER_SEC;
}

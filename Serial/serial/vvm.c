#ifdef __cplusplus
extern "C" {
#endif
    void vvm_(  int *len, double *va, double *vb, double *ma);
#ifdef __cplusplus
    }
#endif

// Computes the tensor product of two vectors 

void  vvm_( int *len, double *va, double *vb, double *ma){

	int i, j;

	int alength = *len;

	for (i=0; i<alength; i++) {
		for (j=0; j<alength; j++) {
			*(ma+(alength*i)+j) = *(va+i) * *(vb+j);
		}
	}

}

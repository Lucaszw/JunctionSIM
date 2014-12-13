/* File: cubic.h. -- #defines, typedefs, and prototypes for cubic_io.c &
   cubic_extrema.c. -- by Mike J. Courtney
*/
/* status defines */
//#define SUCCESS 1
#define FAILURE 0
#define FAILURE1    -1
#define FAILURE2    -2
/* flag defines */
#define TRUE    1
#define FALSE   0
/* typedefs */
typedef int BOOL;
typedef unsigned char BYTE;
/* structures */
struct point {
  struct point *next;
  float x;
  float y;
};
/* prototypes */
BOOL ComputeSecDerivs (unsigned int, float *, float *, float *);
BOOL FindCubicExtrema (unsigned int, float *,float *, struct point *);
BOOL FindQuadRoots (float, float, float, float *, float *);
void ComputeY (unsigned int, float *, float, float *, float *, float *);

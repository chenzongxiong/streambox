#ifndef _ALBERTA_DEMO_GRAPHICS_H_
#define _ALBERTA_DEMO_GRAPHICS_H_

#include <alberta/alberta.h>

#ifdef __cplusplus
extern "C" {
#elif 0
} /* otherwise some editors attempt automatic indentation based on above '{' */
#endif

/* function for displaying mesh, discrete solution, and/or estimate
 * defined in graphics.c
 */
extern void graphics(MESH *mesh, DOF_REAL_VEC *u_h, REAL (*get_est)(EL *el),
		     REAL (*u)(const REAL_D x), REAL time);

extern void graphics_d(MESH *mesh, DOF_REAL_VEC_D *u_h, DOF_REAL_VEC *p_h,
		       REAL (*get_est)(EL *el),
		       const REAL *(*u)(const REAL_D val, REAL_D x), REAL time);
  
#ifdef __cplusplus
}
#endif

#endif

#if HAVE_CONFIG_H
# include "config.h"
#endif

//#include <stdio.h>
#include <alberta/alberta.h>

#include "graphics.h" /* just to keep the proto-type consistent */

/*---8<---------------------------------------------------------------------*/
/*---   graphical output of mesh, discrete solution, estimate, and error ---*/
/*--------------------------------------------------------------------->8---*/

#if !ALBERTA_USE_GRAPHICS

void graphics(MESH *mesh, DOF_REAL_VEC *u_h, REAL (*get_est)(EL *el),
	      REAL (*u)(const REAL_D x), REAL time)
{
  FUNCNAME("graphics");
  MSG("Graphics disabled by compile-time switch !ALBERTA_USE_GRAPHICS");
}


void graphics_d(MESH *mesh, DOF_REAL_VEC_D *u_h, DOF_REAL_VEC *p_h,
		REAL (*get_est)(EL *el),
		const REAL *(*u)(const REAL_D val, REAL_D x), REAL time)
{
  FUNCNAME("graphics_d");
  MSG("Graphics disabled by compile-time switch !ALBERTA_USE_GRAPHICS");
}

#else

# if !HAVE_LIBGLTOOLS
/*---8<---------------------------------------------------------------------*/
/*---   simple GL graphics ...                                           ---*/
/*---   nothing will be done in 3d                                       ---*/
/*--------------------------------------------------------------------->8---*/

void graphics(MESH *mesh, DOF_REAL_VEC *u_h, REAL (*get_est)(EL *el),
	      REAL (*u)(const REAL_D x), REAL time)
{
  /*printf("*** in routine graphics ***\n");*/
  FUNCNAME("graphics");
  static bool first = true;
  static GRAPH_WINDOW win_est, win_val, win_mesh, win_err;
  static DOF_REAL_VEC   *u_diff = NULL;
  int    refine = 0;

  if (first) {
    int size[4] = { 600,600,600,600 };
    char   geom[128] = "500x500+0+0";
//    GET_PARAMETER(1, "graphic windows", "%d %d %d %d", size, size+1, size+2, size+3);

    if (size[0] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[0], size[0], 0, 0);
      win_val = graph_open_window("ALBERTA values", geom, NULL, mesh);
    }
    if (size[1] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[0], size[0], size[0], 0);
      win_est = graph_open_window("ALBERTA estimate", geom, NULL, mesh);
    }
    if (size[2] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[2], size[2], size[0]+size[1], 0);
      win_mesh = graph_open_window("ALBERTA mesh", geom, NULL, mesh);
    }
    if (size[3] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[3], size[3], 0, size[0]);
      win_err = graph_open_window("ALBERTA error (u_h-u)", geom, NULL, mesh);
    }

    first = false;
  }

  if (mesh && win_mesh) {
    graph_clear_window(win_mesh, rgb_white);
    graph_mesh(win_mesh, mesh, rgb_black, 2);
  }

  if (u_h && win_val) {
    graph_clear_window(win_val, rgb_white);
//    refine = u_h->fe_space->bas_fcts->degree;
    /*if (refine<2) */ refine=0;
    graph_drv(win_val, u_h, 0.0, 0.0, refine);
  }

  if (get_est && win_est) {
    graph_clear_window(win_est, rgb_white);
    graph_mesh(win_est, mesh, rgb_blue, 0);
    graph_el_est(win_est, mesh, get_est, 0.0, 0.0);
  }

  if(u && u_h && win_err) {
    if(!u_diff)
      u_diff = get_dof_real_vec("u_h - u", u_h->fe_space);

    interpol(u, u_diff);
    dof_xpay(-1.0, u_h, u_diff);
    graph_drv(win_err, u_diff, 0.0, 0.0, refine);
  }

  WAIT;

  return;
}

void graphics_d(MESH *mesh, DOF_REAL_VEC_D *u_h, DOF_REAL_VEC *p_h,
		REAL (*get_est)(EL *el),
		const REAL *(*u)(const REAL_D val, REAL_D x), REAL time)
{
  FUNCNAME("graphics_d");
  static bool first = true;
  static GRAPH_WINDOW win_est, win_val, win_mesh, win_err, win_pressure;
  static DOF_REAL_D_VEC *u_diff = NULL;
  int refine = 0;

  if (first) {
    int size[6] = { 0, };
    char   geom[128] = "500x500+0+0";
    GET_PARAMETER(1, "graphic windows", "%d %d %d %d %d %d",
		  size, size+1, size+2, size+3, size+4, size+5);

    if (size[0] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[0], size[0], 0, 0);
      win_val = graph_open_window("ALBERTA values", geom, NULL, mesh);
    }

    if (size[1] > 0) {
      WARNING("Plotting of vector-valued "
	      "functions as flow-field not supported.\n");
    }

    if (size[2] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[2], size[2], size[0]+size[1], 0);
      win_pressure = graph_open_window("ALBERTA pressure", geom, NULL, mesh);
    }

    if (size[3] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[3], size[3], 0, size[0]);
      win_est = graph_open_window("ALBERTA estimate", geom, NULL, mesh);
    }

    if (size[4] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[4], size[4], size[3], size[0]);
      win_mesh = graph_open_window("ALBERTA mesh", geom, NULL, mesh);
    }

    if (size[5] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[5], size[5], size[3]+size[4], size[0]);
      win_err = graph_open_window("ALBERTA error ||u_h-u||", geom, NULL, mesh);
    }

    first = false;
  }

  if (mesh && win_mesh) {
    graph_clear_window(win_mesh, rgb_white);
    graph_mesh(win_mesh, mesh, rgb_black, 2);
  }

  if (u_h && win_val) {
    graph_clear_window(win_val, rgb_white);
    refine = u_h->fe_space->bas_fcts->degree;
    if (refine<2) refine=0;
    graph_drv_d(win_val, u_h, 0.0, 0.0, refine);
  }

  if (p_h && win_pressure) {
    graph_clear_window(win_pressure, rgb_white);
    refine = p_h->fe_space->bas_fcts->degree;
    if (refine < 2) refine=0;
    graph_drv(win_pressure, p_h, 0.0, 0.0, refine);
  }

  if (get_est && win_est) {
    graph_clear_window(win_est, rgb_white);
    graph_mesh(win_est, mesh, rgb_blue, 0);
    graph_el_est(win_est, mesh, get_est, 0.0, 0.0);
  }

  if(u && u_h && win_err) {
    if(!u_diff)
      u_diff = get_dof_real_d_vec("u_h - u", u_h->fe_space);

    interpol_d(u, u_diff);
    dof_xpay_d(-1.0, u_h, u_diff);
    graph_drv_d(win_err, u_diff, 0.0, 0.0, refine);
  }

  WAIT;
}

# else /* gltools version follows below */

void graphics(MESH *mesh, DOF_REAL_VEC *u_h, REAL (*get_est)(EL *el),
	      REAL (*u)(const REAL_D x) , REAL time)
{
  FUNCNAME("graphics");
  static bool first = true;
  static GLTOOLS_WINDOW  win_est=NULL, win_val=NULL, win_mesh=NULL, win_err=NULL;
  static DOF_REAL_VEC   *u_diff = NULL;
  static REAL  min, max;

  if (first) {
    int size[4] = { 0, };
    char   geom[128] = "500x500+0+0";
    GET_PARAMETER(1, "graphic windows", "%d %d %d %d", size, size+1, size+2, size+3);

    GET_PARAMETER(1, "gltools range", "%e %e", &min, &max);

    if (size[0] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[0], size[0], 0, 0);
      win_val = open_gltools_window("ALBERTA values", geom, NULL, mesh, true);
    }

    if (size[1] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[1], size[1], size[0], 0);
      win_est = open_gltools_window("ALBERTA estimate", geom, NULL, mesh, true);
    }

    if (size[2] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[2], size[2], size[0]+size[1], 0);
      win_mesh = open_gltools_window("ALBERTA mesh", geom, NULL, mesh, true);
    }

    if (size[3] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[3], size[3], 0, size[0]);
      win_err=open_gltools_window("ALBERTA error (u_h-u)", geom, NULL, mesh, true);
    }
    first = false;
  }


  if (mesh && win_mesh) {
    gltools_mesh(win_mesh, mesh, 0, time);
  }


  if (u_h && win_val) {
    gltools_drv(win_val, u_h, min, max, time);
  }

  if (get_est && win_est ) {
    gltools_est(win_est, mesh, get_est, 0.0, -1.0, time);
  }

  if(u && u_h && win_err) {
    if(!u_diff)
      u_diff = get_dof_real_vec("u_h - u", u_h->fe_space);

    interpol(u, u_diff);
    dof_xpay(-1.0, u_h, u_diff);
    gltools_drv(win_err, u_diff, 0.0, -1.0, time);
  }
}

void graphics_d(MESH *mesh, DOF_REAL_VEC_D *u_h, DOF_REAL_VEC *p_h,
		REAL (*get_est)(EL *el),
		const REAL *(*u)(const REAL_D val, REAL_D x),
		REAL time)
{
  FUNCNAME("graphics_d");
  static bool first = true;
  static GLTOOLS_WINDOW win_est, win_val, win_vec, win_pressure;
  static GLTOOLS_WINDOW win_mesh, win_err;
  static DOF_REAL_VEC_D *u_diff;
  static REAL  min, max;

  if (first) {
    int size[6] = { 0, };
    char   geom[128] = "500x500+0+0";
    GET_PARAMETER(1, "graphic windows", "%d %d %d %d %d %d",
		  size, size+1, size+2, size+3, size+4, size+5);

    GET_PARAMETER(1, "gltools range", "%e %e", &min, &max);

    if (size[0] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[0], size[0], 0, 0);
      win_val = open_gltools_window("ALBERTA values", geom, NULL, mesh, true);
    }

    if (size[1] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[1], size[1], size[0], 0);
      win_vec = open_gltools_window("ALBERTA values (vectors)",
				    geom, NULL, mesh, true);
    }

    if (size[2] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[2], size[2], size[0]+size[1], 0);
      win_pressure = open_gltools_window("ALBERTA pressure",
					 geom, NULL, mesh, true);
    }

    if (size[3] > 0) {
      sprintf(geom, "%dx%d+%d+%d",
	      size[3], size[3], 0, size[0]);
      win_est = open_gltools_window("ALBERTA estimate", geom, NULL, mesh, true);
    }

    if (size[4] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[4], size[4], size[3], size[0]);
      win_mesh = open_gltools_window("ALBERTA mesh", geom, NULL, mesh, true);
    }

    if (size[5] > 0) {
      sprintf(geom, "%dx%d+%d+%d", size[5], size[5], size[3]+size[4], size[0]);
      win_err=open_gltools_window("ALBERTA error ||u_h-u||",
				  geom, NULL, mesh, true);
    }
    first = false;
  }


  if (mesh && win_mesh) {
    gltools_mesh(win_mesh, mesh, 0, time);
  }

  if (u_h && win_val) {
    gltools_drv_d(win_val, u_h, min, max, time);
  }

  if (u_h && win_vec) {
    gltools_vec(win_vec, u_h, min, max, time);
  }

  if (p_h && win_pressure) {
    gltools_drv(win_pressure, p_h, min, max, time);
  }

  if (get_est && win_est) {
    gltools_est(win_est, mesh, get_est, 0.0, -1.0, time);
  }

  if(u && u_h && win_err) {
    if(!u_diff) {
      u_diff = get_dof_real_vec_d("u_h - u", u_h->fe_space);
    }

    interpol_dow(u, u_diff);
    dof_xpay_dow(-1.0, u_h, u_diff);
    gltools_drv_d(win_err, u_diff, 0.0, -1.0, time);
  }
}

# endif /* HAVE_LIBGLTOOLS */

#endif /* USE_GRAPHICS */

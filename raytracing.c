#include "raytracing.h"

int counter()
{
  static int ctr = 0;
  return ctr++;
}

void init_ray_ribbon(Transmitter * tx, Receiver * rx, Perfectreflector * patcharray,
		     Rayribbon * rb, int num_segments, int num_reflectors,
		     int num_reflections)
{
  if (rb == NULL)
  {
    fprintf(stderr, "Null pointer ray ribbon. Cannot proceed.\n");
    return;
  }
  
  double phi, theta;
  Ribbonnode *rnprev = 0;
  int ctr;
  
  for (phi=0; phi<2*PI; phi += 2*PI/num_segments)
  {
    for (theta=-PI/2; theta<PI/2; theta += PI/num_segments)
    {
      Ribbonnode * rn = init_ribbonnode();
      cblas_dcopy(3, tx->gn.smm.position, 1, rn->current->point, 1);
      double direction[3] = {cos(phi), sin(phi)*cos(theta), sin(phi)*sin(theta)};
      cblas_dcopy(3, direction, 1, rn->current->unit_direction, 1);

      bool hit_destination = process_vertical_chain(rn, patcharray, num_reflectors, num_reflections);

      if (hit_destination)
      {
	if (rb->head == NULL)
	{
	  rb->head = rn;
	}
	else
	{
	  //rn->left = rnprev;
	  if (rnprev != NULL) rnprev->right = rn;
	}
	rnprev = rn;

	// check number of non collinear points; if more than 4, break
      }
      else
      {
	destroy_ray_ribbon_vertical_down(rn);
      }

    }
  }
}

Ribbonnode * init_ribbonnode()
{
  Ribbonnode * rn;
  rn = malloc(sizeof(Ribbonnode));
  rn->current = malloc(sizeof(Halfinfiniteray));
  rn->down = rn->right = 0;
  rn->ctr = counter();
  rn->num_reflections = 0;
  rn->hit_destination_patch = false;
  rn->surface_index = -1;

  return rn;
}

void destroy_ray_ribbon(Rayribbon * rb)
{
  Ribbonnode * rn = 0;
  if (rb != NULL)
  {
    rn = rb->head;
    while (rn != NULL)
    {
      Ribbonnode * rnv = rn->right;
      destroy_ray_ribbon_vertical_down(rn);
      rn = rnv;
    }
  }
}

/* void unlink_ray_ribbon_node(Rayribbon * rb, Ribbonnode * rn) */
/* { */
/*   // make sure rn is the highest node possible in the vertical strip */
/*   /\* while (rn->up != NULL) *\/ */
/*   /\* { *\/ */
/*   /\*   rn = rn->up; *\/ */
/*   /\* } *\/ */
  
/*   if (rn->right == NULL && rn->left == NULL) */
/*   { */
/*     fprintf(stderr, "Both null\n"); */
/*     rb->head = NULL; */
/*     // if node to the left or node to the right is null */
/*     // set head and tail of ribbon to NULL */
/*   } */
/*   else if (rn->right == NULL && rn->left != NULL) */
/*   { */
/*     fprintf(stderr, "Right null\n");     */
/*   } */
/*   else if (rn->left == NULL && rn->right != NULL) */
/*   { */
/*     fprintf(stderr, "Left null\n");         */
/*     rb->head = rn->right; */
/*     rn->right->left == NULL; */
/*   } */
/*   else */
/*   { */
/*     fprintf(stderr, "None null\n");             */
/*     // if nodes to the left and right are not null */
/*     rn->right->left = rn->left; */
/*     rn->left->right = rn->right; */
/*   } */
/*   //destroy_ray_ribbon_vertical_down(rn); */
/* } */

/* void destroy_ray_ribbon_vertical(Ribbonnode * rn) */
/* { */
/*   if (rn == NULL) return; */
/*   while (rn->up != NULL) */
/*   { */
/*     rn = rn->up; */
/*   } */
/*   destroy_ray_ribbon_vertical_down(rn); */
  
/* } */

void destroy_ray_ribbon_vertical_down(Ribbonnode * rn)
{
  if (rn == NULL) return;
  if (rn->current != NULL) free(rn->current);
  Ribbonnode * rndown = rn->down;
  free(rn);
  destroy_ray_ribbon_vertical_down(rndown);
} 

complex double compute_intersection(Halfinfiniteray * hr, Perfectreflector * pr)
{
  double t, sgn;
  double diff[3];
  cblas_dcopy(3, hr->point, 1, diff, 1);
  cblas_daxpy(3, -1, pr->center_point, 1, diff, 1);
  t = -cblas_ddot(3, diff, 1, pr->unit_normal, 1)
    /cblas_ddot(3, hr->unit_direction, 1, pr->unit_normal, 1);
  //fprintf(stderr, "t obtained is %lf\n", t);

  // check if t lies within the bounds of the patch


  if (t < INFINITY)
  {
    // Verify the signs
    cblas_daxpy(3, t, hr->unit_direction, 1, diff, 1);
    
    //print_vector(diff);
    //print_vector(pr->unit_length_normal);
    //print_vector(pr->unit_width_normal);
    
    double lengtht = cblas_ddot(3, diff, 1, pr->unit_length_normal, 1);
    double widtht = cblas_ddot(3, diff, 1, pr->unit_width_normal, 1);
    /* fprintf(stderr, "lengtht and widtht obtained is %lf, %lf. Limits are %lf, %lf.\n", */
    /* 	  lengtht, widtht, pr->length/2, pr->width/2); */
    if (abs(lengtht) > pr->length/2 || abs(widtht) > pr->width/2 || t<0)
    {
      t = INFINITY;
    }
  }

  if (t < INFINITY)
  {
    cblas_dcopy(3, hr->point, 1, hr->end_pt, 1);
    cblas_daxpy(3, t, hr->unit_direction, 1, hr->end_pt, 1);
  }
  
  sgn = cblas_ddot(3, hr->unit_direction, 1, pr->unit_normal, 1);
  return t + I*sgn; // if sgn positive then ray is blocked
}

/* void process_rayribbon(Rayribbon * rb, Perfectreflector * pr, */
/* 		       int num_reflectors, int num_reflections) // last one is the destination patch */
/* { */
/*   Ribbonnode * rn = rb->head; */

/*   while (rn != NULL) */
/*   { */
/*     fprintf(stderr, "Processing %d\n", rn->ctr); */
/*     bool hit_destination = process_vertical_chain(rn, pr, num_reflectors, num_reflections); */
/*     Ribbonnode * rnnext = rn->right; */
    
/*     if (!hit_destination) */
/*     { */
/*       fprintf(stderr, "Removing  %d\n", rn->ctr); */
/*       unlink_ray_ribbon_node(rb, rn); */
/*       if (rn->left != NULL) rnnext = rn->left->right;       */
/*       destroy_ray_ribbon_vertical_down(rn); */
/*       fprintf(stderr, "Right is %p\n", (void *) rnnext); */
/*     } */
/*     else */
/*     { */
/*       fprintf(stderr, "Bingo! Keeping  %d\n", rn->ctr); */
/*     } */
/*     rn = rnnext; */
    
/*   } */
/* } */

bool process_vertical_chain(Ribbonnode * rn, Perfectreflector * pr,
			    int num_reflectors, int num_reflections)
{
  // this function computes whether a ray can hit the destination after a max num_reflections
  int ctr=0, ctrindex=-1;
  double tmin = INFINITY, sgn;
  
  while(ctr<num_reflectors)
  {
    complex double dbl = compute_intersection(rn->current, pr+ctr);
    if (creal(dbl) < tmin)
    {
      tmin = creal(dbl);
      sgn = cimag(dbl);
      // printf("Hit reflector %d\n", ctr);
      ctrindex = ctr;
    }
    ctr++;    
  }

  if (sgn>0 || tmin>1e4) return false;
  if (ctrindex == num_reflectors-1)
  {
    rn->hit_destination_patch = true;
    return true;
  }

  if (rn->num_reflections > num_reflections) return false;
  
  // only case remaining is if there is intersection with reflector and number of reflections is small
  Ribbonnode * rn_next = malloc(sizeof(Ribbonnode));
  rn_next->right = rn_next->down = 0;
  rn_next->current = calloc(1, sizeof(Halfinfiniteray));
  rn_next->hit_destination_patch = false;
  rn_next->num_reflections = rn->num_reflections + 1;

  // update ray directions for rn_next

  // first update starting point
  /* for(ctr=0; ctr<3; ctr++) */
  /* { */
  /*   rn_next->current->point[ctr] = rn->current->point[ctr]; */
  /*   rn_next->current->unit_direction[ctr] = rn->current->unit_direction[ctr]; */
  /* } */

  cblas_dcopy(3, rn->current->point, 1, rn_next->current->point, 1);
  cblas_daxpy(3, tmin, rn->current->unit_direction, 1, rn_next->current->point, 1);
  rn_next->surface_index = ctrindex;
  //fprintf(stderr, "Next ribbonnode, tmin: %lf\n", tmin);
  //print_vector(rn_next->current->point);

  // next update direction
  cblas_dcopy(3, rn->current->unit_direction, 1, rn_next->current->unit_direction, 1);
  double factor = -2*cblas_ddot(3, rn->current->unit_direction, 1, (pr+ctrindex)->unit_normal, 1);
  cblas_daxpy(3, factor, (pr+ctrindex)->unit_normal, 1, rn_next->current->unit_direction, 1);

  // update pointers 
  rn->down = rn_next;

  return process_vertical_chain(rn_next, pr, num_reflectors, num_reflections);
}

void print_vector(const double * db)
{
  int ctr;
  for(ctr=0;ctr<3;ctr++)
  {
    fprintf(stderr, "%lf ",*(db+ctr)); 
  }
  fprintf(stderr, "\n");
}

void print_rayribbon(Rayribbon * rb)
{
  Ribbonnode * rn = rb->head;
  int ctr=0, ctr1;
  while (rn != NULL)
  {
    ctr1=0;
    Ribbonnode * rn_next = rn;

    while(rn_next != NULL)
    {
      if (rn_next->hit_destination_patch)
      {
	print_ribbonnode(rn_next);
	//fprintf(stderr, "Bingo! ");
	//fprintf(stderr, "Node %d, %d: ", ctr, ctr1);
	//print_vector(rn_next->current->unit_direction);
      }
      rn_next = rn_next->down;
      ctr1++;
    }
    rn = rn->right;
    ctr++;
  }
}

void print_vertical_strip(Ribbonnode * rn)
{
  int ctr=0;
  while (rn != NULL)
  {
    fprintf(stderr, "Level %d\n", ctr++);
    print_ribbonnode(rn);
    rn = rn->down;
  }
}

void print_ribbonnode(Ribbonnode * rn)
{
  if (rn == NULL) return;
  fprintf(stderr, "Starting point: ");
  int ctr=0;
  for(ctr=0; ctr<3; ctr++)
  {
    fprintf(stderr, "%lf ", rn->current->point[ctr]);
  }
  
  fprintf(stderr, "Unit direction: ");
  for(ctr=0; ctr<3; ctr++)
  {
    fprintf(stderr, "%lf ", rn->current->unit_direction[ctr]);
  }

  fprintf(stderr, "Ending point: ");
  for(ctr=0; ctr<3; ctr++)
  {
    fprintf(stderr, "%lf ", rn->current->end_pt[ctr]);
  }
  
  fprintf(stderr, "Hit destination: %d, num reflections: %d, Surface index: %d",
	  rn->hit_destination_patch, rn->num_reflections, rn->surface_index);
  fprintf(stderr, "\n");
}

int count_segments(Ribbonnode * rn)
{
  int ctr=0;

  while (rn != NULL)
  {
    rn = rn->down;
    ctr++;
  }
  return ctr;
}

void compute_average_ribbonnode(Ribbonnode * rn, Ribbonnode ** node_array, double * weights)
{
  double thet, phi;
  double thetaav=0, phiav=0;
  int ctr;
  for (ctr=0; ctr<3; ctr++)
  {
    invert_spherical_angles(node_array[ctr]->current->unit_direction, &phi, &thet);
    //fprintf(stderr, " Phi: %lf, Theta: %lf\n", phi, thet);
    //print_vector(node_array[ctr]->current->unit_direction);
    
    thetaav += *(weights+ctr) * thet;
    phiav += *(weights+ctr) * phi;
  }
  
  //fprintf(stderr, " Phi: %lf, Theta: %lf\n", phi, thet);
  
  rn->current->unit_direction[0] = cos(phiav);
  rn->current->unit_direction[1] = sin(phiav)*cos(thetaav);
  rn->current->unit_direction[2] = sin(phiav)*sin(thetaav);
  rn->hit_destination_patch = 0;
  rn->num_reflections = 0;
  rn->ctr = 1;

  cblas_dcopy(3, node_array[1]->current->point, 1, rn->current->point, 1);
}

void invert_spherical_angles(double * unit_vector, double * phi, double * thet)
{
  * thet = atan(unit_vector[2]/unit_vector[1]);
  * phi = acos(unit_vector[0]);
  if ( sin(*phi)*sin(*thet)/unit_vector[2] < 0 ) *phi = 2*PI- (*phi);
}

void compute_averaging_coefficients(double * point, Ribbonnode ** node_array, double * weights)
{
  gsl_matrix * mat = gsl_matrix_alloc(3, 2);
  gsl_permutation * perm = gsl_permutation_alloc(3);
  
  gsl_vector * x = gsl_vector_alloc(2);
  gsl_vector * b = gsl_vector_alloc(3);
  gsl_vector * tau = gsl_vector_alloc(2);
  gsl_vector * residual = gsl_vector_alloc(3);
  
  int c0, c1;
  for (c0=0; c0<3; c0++)
  {
    gsl_vector_set(b, c0, *(point+c0) - node_array[0]->down->current->end_pt[c0]);
    for (c1=0; c1<2; c1++)
    {
      gsl_matrix_set(mat, c0, c1, node_array[c1+1]->down->current->end_pt[c0] - node_array[0]->down->current->end_pt[c0]);
    }
  }

  /* gsl_linalg_LU_decomp (mat, perm, &c0); */
  /* gsl_linalg_LU_solve (mat, perm, b, x); */

  gsl_linalg_QR_decomp(mat, tau);
  gsl_linalg_QR_lssolve(mat, tau, b, x, residual);
  
  *(weights) = 1 - gsl_vector_get(x, 0) - gsl_vector_get(x, 1);
  *(weights + 1) = gsl_vector_get(x, 0);
  *(weights + 2) = gsl_vector_get(x, 1);

  gsl_matrix_free(mat);
  gsl_permutation_free(perm);
  gsl_vector_free(x);
  gsl_vector_free(b);
  gsl_vector_free(tau);
  gsl_vector_free(residual);
}

Ribbonnode * refine_ribbonnode(Ribbonnode ** node_array, double * point, Ribbonnode * rn,
			       Perfectreflector * pr)
{
  Ribbonnode * node_array_mod[3];
  int ctr;
  for(ctr=0; ctr<3; ctr++)
  {
    node_array_mod[ctr] = node_array[ctr];
  }

  return _refine_ribbonnode_helper(node_array_mod, point, rn, false, pr);
}

static Ribbonnode * _refine_ribbonnode_helper(Ribbonnode ** node_array, double * point, Ribbonnode * rn, bool null_node_array,
				       Perfectreflector * pr)
{
  // null node array should be called with false otherwise segmentation fault or memory leak will occur!

  if (null_node_array) destroy_ray_ribbon_vertical_down(node_array[0]); 

  if (rn != NULL && isclose(rn, point))
  {
    return rn;
  }
  
  if (rn != NULL && !isclose(rn, point))
  {
    node_array[0] = rn;
    null_node_array = true;
    rn = init_ribbonnode();
  }
  
  if (rn == NULL )
  {
    rn = init_ribbonnode();
  }
  
  double weights[3];
  compute_averaging_coefficients(point, node_array, weights);
  compute_average_ribbonnode(rn, node_array, weights);
  bool has_hit = process_vertical_chain(rn, pr, 3, 3); //update with num reflectors TODO

  if (has_hit)
  {
    return _refine_ribbonnode_helper(node_array, point, rn, null_node_array, pr);
  }
  else
  {
    return NULL;
  }
}

bool isclose(Ribbonnode * rn, double * point)
{
  double diff[3];
  while (rn->down != NULL)
  {
    rn = rn->down;
  }
  cblas_dcopy(3, rn->current->end_pt, 1, diff, 1);
  cblas_daxpy(3, -1, point, 1, diff, 1);
  if (cblas_dnrm2(3, diff, 1) < 1e-5)  return true;
  return false;
}

long type_vertical_strip(Ribbonnode * rn)
{
  int ctr=0;
  while (rn != NULL)
  {
    ctr = MAX_SURFACES*ctr + rn->surface_index+1;
    rn = rn->down;
  }
  return ctr;
}

Ribbonnode ** vertical_strip_for_points(Ribbonnode ** nodearray, double ** points,
					int num_points, Perfectreflector * pr)
{
  int ctr=0;
  Ribbonnode ** vertical_strips = malloc(num_points*sizeof(Ribbonnode *));
  for (ctr=0; ctr<num_points; ctr++)
  {
    *(vertical_strips + ctr) = refine_ribbonnode(nodearray, points[ctr], NULL, pr);
  }
  return vertical_strips;
}

Path * generate_all_paths(Transmitter * tn, Receiver * rxarray, Perfectreflector * pr)
{
  // make all the transmitter, receiver and reflector arrays null terminated; that way code is simplified
  // look up best practices for documenting
  // look up best practices for testing
}

#include "raytracing.h"

void init_ray_ribbon(Transmitter * tx, Receiver * rx, Perfectreflector * patcharray,
		     Rayribbon * rb, int num_segments, int num_reflectors,
		     int num_reflections)
{
  if (rb == NULL)
  {
    return;
    fprintf(stderr, "Null pointer. Cannot proceed.\n");
  }
  
  double phi, theta;
  Ribbonnode *rnprev = 0;
  int ctr = 0;
  
  for (phi=0; phi<2*PI; phi += 2*PI/num_segments)
  {
    for (theta=-PI/2; theta<PI/2; theta += PI/num_segments)
    {
      Ribbonnode * rn;
      rn = malloc(sizeof(Ribbonnode));
      rn->right = rn->left = rn->up = rn->down = 0;
      rn->hit_destination_patch = false;
      rn->num_reflections = 0;
      
      rn->current = malloc(sizeof(Halfinfiniteray));
      rn->ctr = ctr++;
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
	  rn->left = rnprev;
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
  rb->tail = rnprev;
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

void unlink_ray_ribbon_node(Rayribbon * rb, Ribbonnode * rn)
{
  // make sure rn is the highest node possible in the vertical strip
  /* while (rn->up != NULL) */
  /* { */
  /*   rn = rn->up; */
  /* } */
  
  if (rn->right == NULL && rn->left == NULL)
  {
    fprintf(stderr, "Both null\n");
    rb->head = NULL;
    rb->tail = NULL;
    // if node to the left or node to the right is null
    // set head and tail of ribbon to NULL
  }
  else if (rn->right == NULL && rn->left != NULL)
  {
    fprintf(stderr, "Right null\n");    
    rb->tail = rn->left;
    rb->tail->right = 0;
  }
  else if (rn->left == NULL && rn->right != NULL)
  {
    fprintf(stderr, "Left null\n");        
    rb->head = rn->right;
    rn->right->left == NULL;
  }
  else
  {
    fprintf(stderr, "None null\n");            
    // if nodes to the left and right are not null
    rn->right->left = rn->left;
    rn->left->right = rn->right;
  }
  //destroy_ray_ribbon_vertical_down(rn);
}

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
  rn_next->right = rn_next->left = rn_next->up = rn_next->down = 0;
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
  //fprintf(stderr, "Next ribbonnode, tmin: %lf\n", tmin);
  //print_vector(rn_next->current->point);

  // next update direction
  cblas_dcopy(3, rn->current->unit_direction, 1, rn_next->current->unit_direction, 1);
  double factor = -2*cblas_ddot(3, rn->current->unit_direction, 1, (pr+ctrindex)->unit_normal, 1);
  cblas_daxpy(3, factor, (pr+ctrindex)->unit_normal, 1, rn_next->current->unit_direction, 1);

  // update pointers 
  rn->down = rn_next;
  rn_next->up = rn;

  return process_vertical_chain(rn_next, pr, num_reflectors, num_reflections);
}

void print_vector(const double * db)
{
  int ctr=0;
  for(;ctr<3;ctr++)
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
      if (rn_next->hit_destination_patch || 1)
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
  fprintf(stderr, "Hit destination: %d, num reflections: %d", rn->hit_destination_patch, rn->num_reflections);
  fprintf(stderr, "\n");
}

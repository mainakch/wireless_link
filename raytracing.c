#include "raytracing.h"

static Ribbonnode *_refine_ribbonnode(Ribbonnode ** node_array,
				      const double * point,
				      Ribbonnode * rn,
				      bool null_node_array,
				      const Perfectreflector **pr);

int counter()
{
	static int ctr = 0;
	return ctr++;
}

void init_ray_ribbon(Transmitter * tx, Receiver * rx,
		     const Perfectreflector ** ref_arr,
		     Rayribbon * rb,
		     int num_divs, int num_ref)
{
	if (rb == NULL)
	{
		fprintf(stderr, "Null pointer ray ribbon. Cannot proceed.\n");
		return;
	}
	
	double phi, theta;
	Ribbonnode *rnprev = 0;
	for (phi=0; phi<2*PI; phi += 2*PI/num_divs)
	{
		for (theta=-PI/2; theta<PI/2; theta += PI/num_divs)
		{
			Ribbonnode * rn = init_ribbonnode();
			cblas_dcopy(3, tx->gn.smm.position, 1,
				    rn->current->point, 1);
			double direction[3] = {cos(phi),
					       sin(phi)*cos(theta),
					       sin(phi)*sin(theta)};
			cblas_dcopy(3, direction, 1,
				    rn->current->unit_direction, 1);
			
			bool hit_des = process_vertical_chain(rn,
							      ref_arr,
							      num_ref);
			
			if (hit_des)
			{
				if (rb->head == NULL)
				{
					rb->head = rn;
				}
				else
				{
					//rn->left = rnprev;
					if (rnprev != NULL)
						rnprev->right = rn;
				}
				rnprev = rn;
				
				// check number of non collinear
				// points; if more than 4, break
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

void destroy_ray_ribbon_vertical_down(Ribbonnode * rn)
{
	if (rn == NULL) return;
	if (rn->current != NULL) free(rn->current);
	Ribbonnode * rndown = rn->down;
	free(rn);
	destroy_ray_ribbon_vertical_down(rndown);
} 

complex double compute_intersection(Halfinfiniteray * hr,
				    const Perfectreflector * pr)
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
		
		double lengtht = cblas_ddot(3, diff, 1,
					    pr->unit_length_normal, 1);
		double widtht = cblas_ddot(3, diff, 1,
					   pr->unit_width_normal, 1);

		if (abs(lengtht) > pr->length/2 ||
		    abs(widtht) > pr->width/2 ||
		    t<0)
		{
			t = INFINITY;
		}
	}
	
	if (t < 1e4)
	{
		cblas_dcopy(3, hr->point, 1, hr->end_pt, 1);
		cblas_daxpy(3, t, hr->unit_direction, 1, hr->end_pt, 1);
	}
	
	sgn = cblas_ddot(3, hr->unit_direction, 1, pr->unit_normal, 1);
	return t + I*sgn; // if sgn positive then ray is blocked
}

bool process_vertical_chain(Ribbonnode * rn,
			    const Perfectreflector ** pr,
			    int num_reflections)
{
	// this function computes whether a ray can hit the
	// destination after a max num_reflections
	int ctr=0, ctrindex=-1, num_reflectors=0;
	double tmin = INFINITY, sgn;
	
	const Perfectreflector * prsurf = *(pr + ctr);
	while(prsurf != NULL)
	{
		complex double dbl = compute_intersection(rn->current, prsurf);
		if (creal(dbl) < tmin && ctr != rn->surface_index)
		{
			tmin = creal(dbl);
			sgn = cimag(dbl);
			ctrindex = ctr;
		}
		ctr++;
		prsurf = *(pr+ctr);
	}
	num_reflectors = ctr;
	
	if (sgn>0 || tmin>1e4) return false;
	if (ctrindex == num_reflectors-1)
	{
		rn->hit_destination_patch = true;
		return true;
	}
	
	if (rn->num_reflections > num_reflections) return false;
	
	// only case remaining is if there is intersection with
	// reflector and number of reflections is small
	Ribbonnode * rn_next = malloc(sizeof(Ribbonnode));
	rn_next->right = rn_next->down = 0;
	rn_next->current = calloc(1, sizeof(Halfinfiniteray));
	rn_next->hit_destination_patch = false;
	rn_next->num_reflections = rn->num_reflections + 1;
	
	cblas_dcopy(3, rn->current->point, 1,
		    rn_next->current->point, 1);
	cblas_daxpy(3, tmin, rn->current->unit_direction, 1,
		    rn_next->current->point, 1);
	rn_next->surface_index = ctrindex;
	
	// next update direction
	cblas_dcopy(3, rn->current->unit_direction, 1,
		    rn_next->current->unit_direction, 1);
	const Perfectreflector * prsurface = pr[ctrindex]; 
	double factor = -2*cblas_ddot(3, rn->current->unit_direction,
				      1, prsurface->unit_normal, 1);
	cblas_daxpy(3, factor, prsurface->unit_normal, 1,
		    rn_next->current->unit_direction, 1);
	
	// update pointers 
	rn->down = rn_next;
	return process_vertical_chain(rn_next, pr, num_reflections);
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

void print_rayribbon(const Rayribbon * rb)
{
	Ribbonnode * rn = rb->head;
	int ctr = 0, ctr1;
	while (rn != NULL)
	{
		ctr1 = 0;
		Ribbonnode * rn_next = rn;

		while(rn_next != NULL)
		{
			if (rn_next->hit_destination_patch)
			{
				print_ribbonnode(rn_next);
			}
			rn_next = rn_next->down;
			ctr1++;
		}
		rn = rn->right;
		ctr++;
	}
}

void print_vertical_strip(const Ribbonnode * rn)
{
	int ctr=0;
	while (rn != NULL)
	{
		fprintf(stderr, "Level %d\n", ctr++);
		print_ribbonnode(rn);
		rn = rn->down;
	}
}

void print_ribbonnode(const Ribbonnode * rn)
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
		fprintf(stderr, "%lf ",
			rn->current->unit_direction[ctr]);
	}
	
	fprintf(stderr, "Ending point: ");
	for(ctr=0; ctr<3; ctr++)
	{
		fprintf(stderr, "%lf ", rn->current->end_pt[ctr]);
	}
	
	fprintf(stderr,
		"Hit dest: %d, num reflec: %d, Surf index: %d",
		rn->hit_destination_patch,
		rn->num_reflections, rn->surface_index);
	fprintf(stderr, "\n");
}

int count_segments(const Ribbonnode * rn)
{
	int ctr=0;
	
	while (rn != NULL)
	{
		rn = rn->down;
		ctr++;
	}
	return ctr;
}

void compute_average_ribbonnode(Ribbonnode * rn,
				const Ribbonnode ** node_array,
				double * weights)
{
	double thet, phi;
	double thetaav=0, phiav=0;
	int ctr;
	for (ctr=0; ctr<3; ctr++)
	{
		const Ribbonnode *node = node_array[ctr];
		invert_spherical_angles(node->current->unit_direction,
					&phi, &thet);
		thetaav += *(weights+ctr) * thet;
		phiav += *(weights+ctr) * phi;
	}
	rn->current->unit_direction[0] = cos(phiav);
	rn->current->unit_direction[1] = sin(phiav)*cos(thetaav);
	rn->current->unit_direction[2] = sin(phiav)*sin(thetaav);
	rn->hit_destination_patch = 0;
	rn->num_reflections = 0;
	rn->ctr = 1;
	
	cblas_dcopy(3, node_array[1]->current->point, 1,
		    rn->current->point, 1);
}
                                                                     
void invert_spherical_angles(const double * unit_vector, double * phi,
			     double * thet)
{
	* thet = atan(unit_vector[2]/unit_vector[1]);
	* phi = acos(unit_vector[0]);
	if ( sin(*phi)*sin(*thet)/unit_vector[2] < 0 )
		*phi = 2*PI- (*phi);
}

void compute_averaging_coefficients(const double * point,
				    const Ribbonnode ** node_array,
				    double * weights)
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
		const Ribbonnode *node = node_array[0]->down;
		gsl_vector_set(b, c0, *(point+c0)
			       - node->current->end_pt[c0]);
		for (c1=0; c1<2; c1++)
		{
			const Ribbonnode *nodeset = node_array[c1+1]->down;
			gsl_matrix_set(mat, c0, c1,
				       nodeset->current->end_pt[c0]
				       - node->current->end_pt[c0]);
		}
	}
	
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

Ribbonnode * refine_ribbonnode(Ribbonnode ** node_array, const double * point,
			       Ribbonnode * rn,
			       const Perfectreflector ** pr)
{
	Ribbonnode * node_array_mod[3];
	int ctr;
	for(ctr=0; ctr<3; ctr++)
	{
		node_array_mod[ctr] = node_array[ctr];
	}
	
	return _refine_ribbonnode(node_array_mod, point, rn, false, pr);
}

static Ribbonnode *_refine_ribbonnode(Ribbonnode ** node_array,
				      const double * point,
				      Ribbonnode * rn,
				      bool null_node_array,
				      const Perfectreflector ** pr)
{
	// null node array should be called with false otherwise
	// segmentation fault or memory leak will occur!
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
	const Ribbonnode ** node_array_const = (const Ribbonnode **) node_array;
	compute_averaging_coefficients(point, node_array_const, weights);
	compute_average_ribbonnode(rn, node_array_const, weights);
	bool has_hit = process_vertical_chain(rn, pr, 3);
	
	if (has_hit)
	{
		return _refine_ribbonnode(node_array,
					  point, rn, null_node_array, pr);
	}
	else
	{
		return NULL;
	}
}

bool isclose(const Ribbonnode * rn, const double * point)
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

long type_vertical_strip(const Ribbonnode * rn)
{
	int ctr=0;
	while (rn != NULL)
	{
		ctr = MAX_SURFACES*ctr + rn->surface_index+1;
		rn = rn->down;
	}
	return ctr;
}

Ribbonnode ** vertical_strip_for_points(Ribbonnode ** nodearray,
					const double ** points,
					int num_points,
					const Perfectreflector ** pr)
{
	int ctr=0;
	Ribbonnode ** vertical_strips = malloc(num_points*sizeof(Ribbonnode *));
	for (ctr=0; ctr<num_points; ctr++)
	{
		*(vertical_strips + ctr) = refine_ribbonnode(nodearray,
							     points[ctr],
							     NULL, pr);
	}
	return vertical_strips;
}

/* Path * generate_all_paths(Transmitter * tn, Receiver * rxarray, Perfectreflector * pr) */
/* { */
/*   // make all the transmitter, receiver and reflector arrays null terminated; that way code is simplified */
/*   // look up best practices for documenting */
/*   // look up best practices for testing */
/* } */
 

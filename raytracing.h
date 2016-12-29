#include "models.h"

#define MAX_SURFACES 10

typedef struct
{
  double point[3];
  double unit_direction[3];
  double end_pt[3];

} Halfinfiniteray;

typedef struct Ribbonnode_t
{
  Halfinfiniteray * current;
  struct Ribbonnode_t * right;
  struct Ribbonnode_t * down;
  bool hit_destination_patch;
  int num_reflections;
  int ctr;
  int surface_index;
} Ribbonnode;

typedef struct
{
  Ribbonnode * head;
} Rayribbon;

typedef struct
{
  Rayribbon * rb;
  double doppler;
  double delay;
  double gain;
  double phase;
} Path;

void init_ray_ribbon(Transmitter * tx, Receiver * rx, Perfectreflector ** patcharray,
		     Rayribbon * rb, int num_segments, int num_reflections);
Ribbonnode * init_ribbonnode();
void destroy_ray_ribbon(Rayribbon * rb);
void unlink_ray_ribbon_node(Rayribbon * rb, Ribbonnode * rn);
// void destroy_ray_ribbon_vertical(Ribbonnode * rb);
void destroy_ray_ribbon_vertical_down(Ribbonnode * rb); //assume that each path hits destination patch only once
complex double compute_intersection(Halfinfiniteray * hr, Perfectreflector * pr);
                                                               // if any
/* void process_rayribbon(Rayribbon * rb, Perfectreflector * pr, */
/* 		       int num_reflectors, int max_num_reflections); */
bool process_vertical_chain(Ribbonnode * rn, Perfectreflector ** pr, int num_reflections);
void print_vector(const double * db);
void print_rayribbon(Rayribbon * rb);
void print_vertical_strip(Ribbonnode * rn);
void print_ribbonnode(Ribbonnode * rn);
int count_segments(Ribbonnode * rn);
void compute_average_ribbonnode(Ribbonnode * rn, Ribbonnode ** node_array, double * weights);
void invert_spherical_angles(double * unit_vector, double * phi, double * theta);
void compute_averaging_coefficients(double * point, Ribbonnode ** node_array, double * weights);
Ribbonnode * refine_ribbonnode(Ribbonnode ** node_array, double * point, Ribbonnode * rn,
			       Perfectreflector ** pr);
static Ribbonnode * _refine_ribbonnode_helper(Ribbonnode ** node_array, double * point, Ribbonnode * rn, bool null_node_array,
			       Perfectreflector ** pr);
bool isclose(Ribbonnode * rn, double * point);
long type_vertical_strip(Ribbonnode * rn);
Ribbonnode ** vertical_strip_for_points(Ribbonnode ** nodearray, double ** points,
				       int num_points, Perfectreflector ** pr);
Path * generate_all_paths(Transmitter * tn, Receiver * rxarray, Perfectreflector * pr);
			  

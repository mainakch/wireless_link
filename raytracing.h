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

void init_ray_ribbon(Transmitter * tx, Receiver * rx, const Perfectreflector ** patcharray,
		     Rayribbon * rb, int num_segments, int num_reflections);
Ribbonnode * init_ribbonnode();
void destroy_ray_ribbon(Rayribbon * rb);
void unlink_ray_ribbon_node(Rayribbon * rb, Ribbonnode * rn);
void destroy_ray_ribbon_vertical_down(Ribbonnode * rb); //assume that each path hits destination patch only once

complex double compute_intersection(Halfinfiniteray * hr, const Perfectreflector * pr);
bool process_vertical_chain(Ribbonnode * rn, const Perfectreflector ** pr, int num_reflections);
void print_vector(const double * db);
void print_rayribbon(const Rayribbon * rb);
void print_vertical_strip(const Ribbonnode * rn);
void print_ribbonnode(const Ribbonnode * rn);
int count_segments(const Ribbonnode * rn);
void compute_average_ribbonnode(Ribbonnode * rn, const Ribbonnode ** node_array, double * weights);
void invert_spherical_angles(const double * unit_vector, double * phi, double * theta);
void compute_averaging_coefficients(const double * point, const Ribbonnode ** node_array, double * weights);
Ribbonnode * refine_ribbonnode(Ribbonnode ** node_array, const double * point, Ribbonnode * rn,
			       const Perfectreflector ** pr);
bool isclose(const Ribbonnode * rn, const double * point);
long type_vertical_strip(const Ribbonnode * rn);
Ribbonnode ** vertical_strip_for_points(Ribbonnode ** nodearray, const double ** points,
				       int num_points, const Perfectreflector ** pr);
Path * generate_all_paths(Transmitter * tn, Receiver * rxarray, Perfectreflector * pr);
			  

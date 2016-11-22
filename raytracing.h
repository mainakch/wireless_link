#include "models.h"

typedef struct
{
  double point[3];
  double unit_direction[3];

} Halfinfiniteray;

typedef struct Ribbonnode_t
{
  Halfinfiniteray * current;
  struct Ribbonnode_t * right;
  struct Ribbonnode_t * left;
  struct Ribbonnode_t * up;
  struct Ribbonnode_t * down;
  bool hit_destination_patch;
  int num_reflections;
  int ctr;
} Ribbonnode;

typedef struct
{
  Ribbonnode * head;
  Ribbonnode * tail; 
} Rayribbon;

void init_ray_ribbon(Transmitter * tx, Receiver * rx, Perfectreflector * patcharray,
		     Rayribbon * rb, int num_segments, int num_reflectors,
		     int num_reflections);
void destroy_ray_ribbon(Rayribbon * rb);
void unlink_ray_ribbon_node(Rayribbon * rb, Ribbonnode * rn);
// void destroy_ray_ribbon_vertical(Ribbonnode * rb);
void destroy_ray_ribbon_vertical_down(Ribbonnode * rb); //assume that each path hits destination patch only once
complex double compute_intersection(Halfinfiniteray * hr, Perfectreflector * parr);
                                                               // if any
/* void process_rayribbon(Rayribbon * rb, Perfectreflector * pr, */
/* 		       int num_reflectors, int max_num_reflections); */
bool process_vertical_chain(Ribbonnode * rn, Perfectreflector * pr,
			    int num_reflectors, int num_reflections);
void print_vector(const double * db);
void print_rayribbon(Rayribbon * rb);
void print_ribbonnode(Ribbonnode * rn);

#include "models.h"

#define MAX_SURFACES 10

struct half_infinite_ray
{
        double point[3];
        double unit_direction[3];
        double end_pt[3];

};

struct ribbon_node
{
        struct half_infinite_ray *current;
        struct ribbon_node *right;
        struct ribbon_node *down;
        bool hit_destination_patch;
        int num_reflections;
        int ctr;
        int surface_index;
};

struct ray_ribbon
{
        struct ribbon_node *head;
};

struct path
{
        struct ray_ribbon *rb;
        double doppler;
        double delay;
        double gain;
        double phase;
};

void init_ray_ribbon(struct transmitter *tx, struct receiver *rx,
                     const struct perfect_reflector **patcharray,
                     struct ray_ribbon *rb, int num_segments,
                     int num_reflections);
struct ribbon_node *init_ribbonnode();
void destroy_ray_ribbon(struct ray_ribbon *rb);
void unlink_ray_ribbon_node(struct ray_ribbon *rb, struct ribbon_node *rn);
//assume that each path hits destination patch only once
void destroy_ray_ribbon_vertical_down(struct ribbon_node *rb);

complex double compute_intersection(struct half_infinite_ray *hr,
                                    const struct perfect_reflector *pr);
bool process_vertical_chain(struct ribbon_node *rn,
                            const struct perfect_reflector **pr,
                            int num_reflections);
void print_vector(const double *db);
void print_rayribbon(const struct ray_ribbon *rb);
void print_vertical_strip(const struct ribbon_node *rn);
void print_ribbonnode(const struct ribbon_node *rn);
int count_segments(const struct ribbon_node *rn);
void compute_average_ribbonnode(struct ribbon_node *rn,
                                const struct ribbon_node **node_array,
                                double *weights);
void invert_spherical_angles(const double *unit_vector, double *phi,
                             double *theta);
void compute_averaging_coefficients(const double *point,
                                    const struct ribbon_node **node_array,
                                    double *weights);
struct ribbon_node *refine_ribbonnode(struct ribbon_node **node_array,
                               const double *point,
                               struct ribbon_node *rn,
                               const struct perfect_reflector **pr);
bool isclose(const struct ribbon_node *rn, const double *point);
long type_vertical_strip(const struct ribbon_node *rn);
struct ribbon_node **vertical_strip_for_points(struct ribbon_node **nodearray,
                                        const double **points,
                                       int num_points,
                                        const struct perfect_reflector **pr);
struct path *generate_all_paths(struct transmitter *tn, struct receiver *rxarray,
                          struct perfect_reflector *pr);

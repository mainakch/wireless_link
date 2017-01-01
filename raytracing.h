#include "models.h"

#define MAX_SURFACES 10
#define NUM_TYPES 100

struct half_infinite_ray {
        double point[3];
        double unit_direction[3];
        double end_pt[3];

};

struct ribbon_node {
        struct half_infinite_ray *current;
        struct ribbon_node *down;
        bool hit_destination_patch;
        int num_reflections;
        int ctr;
        int surface_index;
};

struct ray_ribbon {
        struct ribbon_node *head;
};

struct ray_ribbon_array {
        struct ray_ribbon **ribbons;
        int max_len;
        int current_len;
};

struct path {
        struct ray_ribbon *rb;
        double doppler;
        double delay;
        double gain;
        double phase;
};

struct set_of_types {
        long *types;
};


struct set_of_types *init_set_of_types(int number);
void destroy_set_of_types(struct set_of_types *st);
void add_entry(long number, struct set_of_types *st);

void init_ray_ribbon_array(int number, struct ray_ribbon_array *rarr);
void populate_ray_ribbon_array(struct transmitter *tx, struct receiver *rx,
                               const struct perfect_reflector **patcharray,
                               struct ray_ribbon_array *rarr,
                               int num_segments, int num_reflections);
void populate_ray_ribbon_array_long(struct transmitter *tx, struct receiver *rx,
                                    const struct perfect_reflector **ref_arr,
                                    struct ray_ribbon_array *rarr,
                                    int num_ref,
                                    const double phi_start,
                                    const double phi_end,
                                    const double phi_delta,
                                    const double theta_start,
                                    const double theta_end,
                                    const double theta_delta);
void populate_ray_ribbon_array_full(const struct transmitter *tx,
                                    const struct perfect_reflector **ref_arr,
                                    int num_ref, int num_points,
                                    const complex double *angles,
                                    struct ray_ribbon_array *rarr);
struct ray_ribbon *init_ray_ribbon(struct ribbon_node *rn);
struct ribbon_node *init_ribbon_node();
void destroy_ray_ribbon(struct ray_ribbon *rb);
void destroy_ray_ribbon_nodes(struct ray_ribbon *rb);
void destroy_ray_ribbon_array(struct ray_ribbon_array *array);
void destroy_ray_ribbon_array_all_but_first(struct ray_ribbon_array *array);
void destroy_chain_of_ribbon_nodes(struct ribbon_node *rn);

bool add_ray_ribbon(struct ray_ribbon_array *array, struct ray_ribbon *rb);
complex double compute_intersection(struct half_infinite_ray *hr,
                                    const struct perfect_reflector *pr);
bool process_vertical_chain(struct ribbon_node *rn,
                            const struct perfect_reflector **pr,
                            int num_reflections);
void print_vector(const double *db);
void print_ray_ribbon(const struct ray_ribbon *rb);
void print_ray_ribbon_array(const struct ray_ribbon_array *rarr);
void print_vertical_strip(const struct ribbon_node *rn);
void print_ribbon_node(const struct ribbon_node *rn);
int count_segments(const struct ribbon_node *rn);
void compute_average_ribbon_node(struct ribbon_node *rn,
                                const struct ray_ribbon_array *rba,
                                double *weights);
void invert_spherical_angles(const double *unit_vector, double *phi,
                             double *theta);
void compute_averaging_coefficients(const double *point,
                                    const struct ray_ribbon_array *rba,
                                    double *weights);
struct ray_ribbon_array *generate_nearby_ribbons(const struct transmitter *tx,
                                                 const struct
                                                 perfect_reflector **ref_arr,
                                                 int num_ref,
                                                 const struct ray_ribbon *rb);
struct ray_ribbon *refine_ray_ribbon(const struct transmitter *tx,
                                     const struct ray_ribbon *rb,
                                     const struct receiver *rx,
                                     const struct perfect_reflector **pr);
bool is_close_ribbon(const struct ray_ribbon *rb, const double *point);
bool isclose(const struct ribbon_node *rn, const double *point);
long type_ray_ribbon(const struct ray_ribbon *rb);
struct ribbon_node **vertical_strip_for_points(struct ribbon_node **nodearray,
                                        const double **points,
                                       int num_points,
                                        const struct perfect_reflector **pr);
struct ray_ribbon_array *throw_three_dim_ray_ribbon(struct transmitter *tn,
                                                    struct receiver *rx,
                                                    const struct
                                                    perfect_reflector **p,
                                                    int num_ref,
                                                    const double phi_start,
                                                    const double phi_end,
                                                    const double phi_incr,
                                                    const double thet_start,
                                                    const double thet_end,
                                                    const double thet_incr);
struct set_of_types *get_unique_types(const struct ray_ribbon_array *arr);
void print_ray_ribbon_types(const struct ray_ribbon_array *arr);
struct ribbon_node *get_last_ribbon_node(const struct ray_ribbon *rb);

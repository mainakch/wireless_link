#include "models.h"

#define NUM_TYPES 100
#define _RAYTRACING_DEBUG 0

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
        double delay;
        double doppler;
        // phases not normalized by 2 PI
        double integrated_doppler_phase;
        double gain;
        double reflection_phase;
        const struct transmitter *start_gn;
        const struct receiver *end_gn;
};

struct ray_ribbon_array {
        struct ray_ribbon **ribbons;
        int max_len;
        int current_len;
};

struct receiver_ray_ribbon {
        // linked list of receiver ray ribbons
        struct ray_ribbon *ribbon;
        struct signal_buffer *signal;
        double delay;
        double doppler;
        // phases not normalized by 2 PI
        double integrated_doppler_phase;
        double gain;
        double reflection_phase;
        bool active;
        const struct transmitter *start_gn;
};

struct signal_buffer {
        // linked list of signal buffer
        struct signal_buffer *next;
        double complex signal;
        double transmit_time;
};

struct receiver_ray_ribbon *init_receiver_ray_ribbon(
        struct receiver *rx,
        struct ray_ribbon *ribbon,
        struct environment *env);
void populate_receiver_ray_ribbons(struct environment *env);
bool update_receiver_ray_ribbons(struct receiver *rx,
                                struct environment *env);
void update_all_receiver_ray_ribbons(struct environment *env);
void update_receiver_ray_ribbons_signal_buffer(struct receiver *rx,
                                              struct environment *env);
void destroy_receiver_ray_ribbon(struct receiver_ray_ribbon *rrbn);
struct signal_buffer *init_signal_buffer();
void destroy_signal_buffer(struct signal_buffer *sgn);
struct signal_buffer *destroy_signal_buffer_first(struct signal_buffer *sgn);
struct ray_ribbon_array *init_ray_ribbon_array(int number);
void populate_ray_ribbon_array(struct transmitter *tx,
                               const struct perfect_reflector **patcharray,
                               struct ray_ribbon_array *rarr,
                               int num_segments, int num_reflections,
        bool single_type);
void populate_ray_ribbon_array_long(struct transmitter *tx,
                                    const struct perfect_reflector **ref_arr,
                                    struct ray_ribbon_array *rarr,
                                    int num_ref,
                                    const double phi_start,
                                    const double phi_end,
                                    const double phi_delta,
                                    const double theta_start,
                                    const double theta_end,
                                    const double theta_delta,
                                    bool single_type);
void populate_ray_ribbon_array_full(const struct transmitter *tx,
                                    const struct perfect_reflector **ref_arr,
                                    int num_ref, int num_points,
                                    const double complex *angles,
                                    struct ray_ribbon_array *rarr,
                                    bool single_type);
struct ray_ribbon *init_ray_ribbon(struct ribbon_node *rn);
struct ribbon_node *init_ribbon_node();
void destroy_ray_ribbon(struct ray_ribbon *rb);
void destroy_ray_ribbon_nodes(struct ray_ribbon *rb);
void destroy_ray_ribbon_array(struct ray_ribbon_array *array);
void destroy_ray_ribbon_array_all_but_first(struct ray_ribbon_array *array);
void destroy_chain_of_ribbon_nodes(struct ribbon_node *rn);
void destroy_ribbon_node(struct ribbon_node *rn);

bool check_same_type(const struct ray_ribbon *ray_rb1,
                     const struct ray_ribbon *ray_rb2);
bool add_ray_ribbon(struct ray_ribbon_array *array, struct ray_ribbon *rb,
                    bool single_type);
double complex compute_intersection(struct half_infinite_ray *hr,
                                    const struct perfect_reflector *pr);
void remove_ribbon_node_duplicates(struct ribbon_node *rn);
bool process_vertical_chain(struct ribbon_node *rn,
                            const struct perfect_reflector **pr,
                            int num_reflections);
void print_ray_ribbon(const struct ray_ribbon *rb);
void print_receiver_ray_ribbon(const struct receiver_ray_ribbon *rb);
void print_ray_ribbon_flattened(const struct ray_ribbon *rb);
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
struct ray_ribbon *refine_ray_ribbon_image(const struct transmitter *tx,
                                     const struct ray_ribbon *rb,
                                     const struct receiver *rx,
                                     const struct perfect_reflector **pr);
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
                                                    const struct
                                                    perfect_reflector **p,
                                                    int num_ref,
                                                    const double phi_start,
                                                    const double phi_end,
                                                    const double phi_incr,
                                                    const double thet_start,
                                                    const double thet_end,
                                                    const double thet_incr);
struct ribbon_node *get_last_ribbon_node(const struct ray_ribbon *rb);
void populate_tx_paths(struct environment *env);
void populate_env_paths(struct environment *env);
void update_env_paths_delay_dopplers(struct environment *env);
void update_ribbon_delay_dopplers(struct ray_ribbon *rb,
                                  const struct environment *env);
void update_receiver_ribbon_delay_dopplers(struct receiver_ray_ribbon *rrb,
                                  const struct environment *env);
double compute_doppler(const struct ray_ribbon *rb,
                       const struct environment *env);
double length_ribbon_node(const struct ribbon_node *rn);
void reflect(const double *pos1, const double *n1, double *vel, double *pos);
void reflection_operation(const double *v1, const double *n1, double *vref);
void readout_all_signals(struct environment *env, FILE *fpout);
void readout_all_signals_buffer(struct environment *env, FILE *fpout);
void clear_tx_paths(struct environment *env);
void clear_env_paths(struct environment *env);

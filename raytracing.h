#include "models.h"

#define NUM_TYPES 100
#define _RAYTRACING_DEBUG 0

struct half_infinite_ray {
        double point[3];
        double unit_direction[3];
        double length;
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
        const struct transmitter *start_tx;
        const struct receiver *end_rx;
};

struct ray_ribbon_array {
        struct ray_ribbon **ribbons;
        int max_len;
        int current_len;
};

struct receiver_ray_ribbon {
        struct ray_ribbon *ribbon;
        struct signal_buffer *signal;
        double delay;
        double doppler;
        // phases not normalized by 2 PI
        double integrated_doppler_phase;
        double gain;
        double reflection_phase;
        const struct transmitter *start_tx;
};

struct receiver_ray_ribbon_ll_node {
        struct receiver_ray_ribbon *rrbn;
        struct receiver_ray_ribbon_ll_node *next;
};

struct signal_buffer {
        // linked list of signal buffer
        struct signal_buffer *next;
        double complex signal;
        double transmit_time;
        double delay;
        bool receiver_read;
};

struct memory_block {
        size_t size;
        void *ptr;
        int sizehash;
        int ptrhash;
        struct memory_block *next;
};

struct memory_register {
        struct memory_block **malloc_array;
        // malloc_array contains array from which blocks will be retrieved
        struct memory_block **free_array;
        // free_array contains array of blocks which have already been served
};

void *custom_malloc(size_t size);
void *custom_calloc(size_t num, size_t size);
void custom_free(void *ptr);
unsigned long hash(unsigned char *str);
unsigned long hash_ptr(void *ptr);
unsigned long hash_int(size_t num);

struct memory_register *init_memory_register();
void destroy_memory_register(struct memory_register *mem_reg);
struct memory_block *init_memory_block(size_t size);
void destroy_memory_block(struct memory_block *mem_block);
void *get_or_free_memory(size_t size, void *ptr, int zero_for_get);
static void free_memory(void *ptr, struct memory_register *mem_reg);
static void *malloc_memory(size_t size, struct memory_register *mem_reg);
void release_all_blocks();

struct receiver_ray_ribbon *init_receiver_ray_ribbon(
        struct receiver *rx,
        const struct ray_ribbon *ribbon,
        struct environment *env);
bool populate_if_ray_ribbon_doesnt_exist
(const struct ray_ribbon *rb1, struct receiver *rx, struct environment *env);
void populate_receiver_ray_ribbons(struct environment *env);
void update_all_receiver_ray_ribbons(struct environment *env);
bool update_receiver_ray_ribbons(struct receiver *rx,
                                struct environment *env);
void update_receiver_ray_ribbons_signal_buffer(struct receiver *rx,
                                              struct environment *env);
void destroy_receiver_ray_ribbon_ll_node(struct receiver_ray_ribbon_ll_node *p);
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
void populate_ray_ribbon_array_full_malloc(const struct transmitter *tx,
                                    const struct perfect_reflector **ref_arr,
                                    int num_ref, int num_points,
                                    const double complex *angles,
                                    struct ray_ribbon_array *rarr,
                                    bool single_type);
void populate_ray_ribbon_array_full_copy(const struct transmitter *tx,
                                    const struct perfect_reflector **ref_arr,
                                    int num_ref, int num_points,
                                    const double complex *angles,
                                    struct ray_ribbon_array *rarr,
                                    bool single_type);
struct ray_ribbon *init_ray_ribbon(struct ribbon_node *rn);
struct ribbon_node *init_ribbon_node();
struct ribbon_node *init_ribbon_node_from_copy(const struct ribbon_node *rn);
struct ribbon_node *init_ribbon_node_from_points(const double *pt1,
                                                             const double *pt2);
struct ribbon_node *init_chain_of_ribbon_nodes(int length);
struct ribbon_node *copy_chain_of_ribbon_nodes(const struct ribbon_node *rn);
struct ribbon_node *copy_chain_of_ribbon_nodes_till_dest(const struct ribbon_node *rn);
struct ray_ribbon *copy_ray_ribbon(const struct ray_ribbon *rb, bool till_dest);
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
bool add_ray_ribbon_copy(struct ray_ribbon_array *array,
                         const struct ray_ribbon *rb,
                         bool single_type);
double complex compute_intersection(const struct half_infinite_ray *hr,
                                    const struct perfect_reflector *pr);
bool process_vertical_chain_nomalloc(struct ribbon_node *rn,
                                     const struct perfect_reflector **pr,
                                     int num_reflections);
void print_ray_ribbon(const struct ray_ribbon *rb);
void print_receiver_ray_ribbon(const struct receiver_ray_ribbon *rb);
void print_ray_ribbon_flattened(const struct ray_ribbon *rb);
void print_ray_ribbon_array(const struct ray_ribbon_array *rarr);
void print_vertical_strip(const struct ribbon_node *rn);
void print_ribbon_node(const struct ribbon_node *rn);
int count_segments(const struct ribbon_node *rn);
void invert_spherical_angles(const double *unit_vector, double *phi,
                             double *theta);


struct ray_ribbon_array *generate_nearby_ribbons(const struct transmitter *tx,
                                                 const struct
                                                 perfect_reflector **ref_arr,
                                                 int num_ref,
                                                 const struct ray_ribbon *rb);
int count_ribbon_nodes(const struct ribbon_node *rn);
struct ray_ribbon *refine_ray_ribbon_image(const struct transmitter *tx,
                                     const struct ray_ribbon *rb,
                                     const struct receiver *rx,
                                     const struct perfect_reflector **pr);

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

// Deprecated
void populate_ray_ribbon_array_full_malloc(const struct transmitter *tx,
                                    const struct perfect_reflector **ref_arr,
                                    int num_ref, int num_points,
                                    const double complex *angles,
                                    struct ray_ribbon_array *rarr,
                                    bool single_type);
bool process_vertical_chain(struct ribbon_node *rn,
                            const struct perfect_reflector **pr,
                            int num_reflections);
struct ray_ribbon *refine_ray_ribbon(const struct transmitter *tx,
                                     const struct ray_ribbon *rb,
                                     const struct receiver *rx,
                                     const struct perfect_reflector **pr);
bool is_close_ribbon(const struct ray_ribbon *rb, const double *point);
bool isclose(const struct ribbon_node *rn, const double *point);
void compute_average_ribbon_node(struct ribbon_node *rn,
                                const struct ray_ribbon_array *rba,
                                double *weights);
void compute_averaging_coefficients(const double *point,
                                    const struct ray_ribbon_array *rba,
                                    double *weights);
void remove_ribbon_node_duplicates(struct ribbon_node *rn);
void update_ribbon_delay_dopplers(struct ray_ribbon *rb,
                                  const struct environment *env);

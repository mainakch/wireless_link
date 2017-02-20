#include <assert.h>
#include <cblas.h>
#include <complex.h>
#include <errno.h>
#include <fftw3.h>
#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#define C 299792458
#define PI 3.141592654
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define _MODEL_DEBUG 0
#define MAX_TX 1
#define MAX_RX 1
#define MAX_RIBBON_SIZE 20
#define MAX_SURFACES 20
// if max values exceeded, may cause memory bugs

struct simulation {
        struct environment *env;
        struct file_reader *fr;
};

struct spatial_motion_model {
        double velocity[3];
        double position[3];
};

struct transmission_model {
        double power_in_dBm;
};

struct propagation_model {
        double distance;
};

struct general_node {
        struct spatial_motion_model *smm;
        struct transmission_model *tm;
        int id;
};

struct transmitter {
        struct general_node *gn;
        double complex baseband_signal;
};

struct receiver {
        struct general_node *gn;
        double recv_noise_power;
        struct receiver_ray_ribbon_ll_node *rlln;
};

struct perfect_reflector {
        double unit_normal[3];
        double unit_length_normal[3];
        double unit_width_normal[3];
        double center_point[3];
        double length, width;
};

struct environment {
        struct receiver **receivers_array;
        struct transmitter **transmitters_array;
        struct perfect_reflector **prarray;
        struct general_node **node_array;
        struct ray_ribbon_array **env_paths;
        struct ray_ribbon_array **tx_paths;

        double recv_unit_normal[3];
        double time;
        double complex *unit_power_gaussian_noise;
        double frequency;
        double wavelength;
        double delta_time;
        double end_time;
        double max_limit;
        double min_limit;
        double boundary_tolerance;

        int num_transmitters;
        int num_receivers;
        int num_reflectors;

        int sz_array_tx;
        int sz_array_rx;
        int sz_array_gn;
        int sz_array_pr;

        // flags
        bool read_in_nodes;
        bool updated_tx_paths;
        /* should be true by the time TX or RX is true */
        bool node_memory_allocated;
};

struct file_reader {
        FILE *infile;
        FILE *outfile;

        // filenames longer than 999 will be truncated
        char input_filename[1000];
        char output_filename[1000];
};

int id();
struct simulation *init_simulation();
struct spatial_motion_model *init_spatial_motion_model();
struct transmission_model *init_transmission_model();
struct propagation_model *init_propagation_model();
struct general_node *init_general_node();
struct transmitter *init_transmitter();
struct receiver *init_receiver();
struct perfect_reflector *init_perfect_reflector(const double *normal,
                                         const double *center_point,
                                         const double *length_normal,
                                         double length, double width);
struct perfect_reflector **init_perfect_reflectorarray(int number);
struct environment *init_environment();
int malloc_environment(struct environment *env);

void print_vector(const double *db);
void print_spatial_motion_model(const struct spatial_motion_model *smm);
void print_transmission_model(const struct transmission_model *tm);
void print_propagation_model(const struct propagation_model *pm);
void print_general_node(const struct general_node *gn);
void print_transmitter(const struct transmitter *tx);
void print_receiver(const struct receiver *rx);
void print_perfect_reflectors(const struct perfect_reflector *pr);
void print_environment(const struct environment *env);
void print_env_paths(const struct environment *env);
void print_tx_paths(const struct environment *env);

void destroy_spatial_motion_model(struct spatial_motion_model *smm);
void destroy_simulation(struct simulation *sim);
void destroy_transmission_model(struct transmission_model *tm);
void destroy_propagation_model(struct propagation_model *pm);
void destroy_transmitter(struct transmitter *tx);
void destroy_receiver(struct receiver *tx);
void destroy_environment(struct environment *env);
void destroy_file_reader(struct file_reader *fr);
void destroy_perfect_reflector(struct perfect_reflector *pr);
void destroy_perfect_reflectorarray(struct perfect_reflector **pr_begin);
double _gaussrand(); // http://c-faq.com/lib/gaussian.html
struct file_reader *init_file_reader(int argc, char *argv[]);
void cross_product(const double *v1, const double *v2, double *v3);
double normalize_unit_vector(double *v1);
void diff(const double *v1, const double *v2, double *v3);
int find_len(void **ptr);

void add_receiver_patch(struct environment *env, int length);
void destroy_last_reflector(struct environment *env);
double distance(const struct general_node *gn1,
                const struct general_node *gn2);

bool update_environment_from_file(struct environment *env,
                                FILE *fp);
bool update_environment_from_file_sim(struct simulation *sim);
bool handle_request(struct environment *env, FILE *fp, const char *req_type);
bool custom_fscanf(FILE *fp, const char *str, void *ptr);

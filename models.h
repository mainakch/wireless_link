#include <cblas.h>
#include <complex.h>
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

struct simulation
{
        double frequency;
        double wavelength;
        double delta_time;
        double max_limit;
        double min_limit;
        double boundary_tolerance;
        int total_time;
};

struct spatial_motion_model 
{
        double velocity[3];
        double position[3];
        double acceleration_factor;
        double distance_from_actual_source;
        bool is_static;
};

struct transmission_model
{
        double power_in_dBm;
        double start_time;
        double end_time;
        double initial_phase;
        double doppler_offset;
};

struct propagation_model
{
        double distance;
};

struct general_node
{
        struct spatial_motion_model smm;
        struct transmission_model tm;
        int id;
}; 

struct transmitter
{
        struct general_node gn;
        bool is_real_transmitter;
};

struct receiver
{
        struct general_node gn;
        double recv_noise_power;
};

struct perfect_reflector
{
        double unit_normal[3];
        double unit_length_normal[3];
        double unit_width_normal[3];
        double center_point[3];
        double length, width;
};

struct environment
{
        int _num_receivers_ctr;
        int _num_transmitters_ctr;
        int num_receivers;
        int num_transmitters;
        int num_virtual_transmitters;
        struct receiver *receivers_array;
        struct transmitter *transmitters_array;
        int time;
        struct general_node **node_array;
        
        double complex *unit_power_gaussian_noise;
        
        double frequency;
        double wavelength;
        double delta_time;
        double max_limit;
        double min_limit;
        double boundary_tolerance;
        int total_time;        
};

struct filereader
{
        FILE *infile;
        FILE *outfile;
        // filenames longer than 1000 will be rejected
        char input_filename[1000];
        char output_filename[1000];
};

int id();
void init_simulation(struct simulation *sim_not_null);
void init_spatial_motion_model(struct spatial_motion_model *smm);
void init_transmission_model(struct transmission_model *tm);
void init_propagation_model(struct propagation_model *pm);
void init_general_node(struct general_node *gn);
void init_transmitter(struct transmitter *tn);
void init_receiver(struct receiver *rc);
struct perfect_reflector *init_perfect_reflector(const double *normal,
                                         const double *center_point,
                                         const double *length_normal,
                                         double length, double width);
struct perfect_reflector **init_perfect_reflectorarray(int number);
void init_environment(struct environment *env);
void init_environment_malloc(struct environment *env);
void init_filereader(struct filereader *fr);
void destroy_environment(struct environment *env);
void destroy_filereader(struct filereader *fr);
void destroy_perfect_reflector(struct perfect_reflector *pr);
void destroy_perfect_reflectorarray(struct perfect_reflector **pr_begin);
double _gaussrand(); // http://c-faq.com/lib/gaussian.html
void interaction_scatterer(void *sc, struct transmitter *tx,
                           struct transmitter *out_array, int *number);
void parse_input(int argc, char *argv[], struct filereader *fr);
void cross_product(const double *v1, const double *v2, double *v3);

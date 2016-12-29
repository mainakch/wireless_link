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
#define MAX_TRANSMITTERS 200 // maximum number of transmitters, real or virtual
#define MAX_RECEIVERS 6400 // maximum number of receivers
#define MAX_TIME 2000 // maximum number of receivers
#define LINE_LENGTH 256


typedef struct 
{
  double frequency;
  double wavelength;
  double delta_time;
  double max_limit;
  double min_limit;
  double boundary_tolerance;
  int total_time;

} Simulation;

typedef struct 
{
  double velocity[3];
  double position[3];
  double acceleration_factor;
  double distance_from_actual_source;
  bool is_static;
} SpatialMotionModel;

typedef struct
{
  double power_in_dBm;
  double start_time;
  double end_time;
  double initial_phase;
  double doppler_offset;
} TransmissionModel;

typedef struct
{
  double distance;
} PropagationModel;

typedef struct
{
  SpatialMotionModel smm;
  TransmissionModel tm;
  int id;
} GeneralNode;

typedef struct
{
  GeneralNode gn;
  bool is_real_transmitter;
}  Transmitter;

typedef struct
{
  GeneralNode gn;
  double recv_noise_power;
}  Receiver;

typedef struct
{
  double unit_normal[3];
  double unit_length_normal[3];
  double unit_width_normal[3];
  double center_point[3];
  double length, width;
} Perfectreflector;
  
typedef struct
{
  int _num_receivers_ctr;
  int _num_transmitters_ctr;

  int num_receivers;
  int num_transmitters;
  int num_virtual_transmitters;
  Receiver * receivers_array;
  Transmitter * transmitters_array;
  int time;
  GeneralNode ** node_array;

  double complex * unit_power_gaussian_noise;

  double frequency;
  double wavelength;
  double delta_time;
  double max_limit;
  double min_limit;
  double boundary_tolerance;
  int total_time;
  
} Environment;

typedef struct
{
  FILE * infile;
  FILE * outfile;
  char input_filename[1000];
  char output_filename[1000];
} Filereader;

int id();
void init_simulation(Simulation * sim_not_null);
void init_spatial_motion_model(SpatialMotionModel * smm);
void init_transmission_model(TransmissionModel * tm);
void init_propagation_model(PropagationModel * pm);
void init_general_node(GeneralNode * gn);
void init_transmitter(Transmitter * tn);
void init_receiver(Receiver * rc);
Perfectreflector * init_perfectreflector(const double * normal,
			   const double * center_point,
			   const double * length_normal,
			   double length, double width);
Perfectreflector ** init_perfectreflectorarray(int number);
void init_environment(Environment * env);
void init_environment_malloc(Environment * env);
void init_filereader(Filereader * fr);
void destroy_environment(Environment * env);
void destroy_filereader(Filereader * fr);
void destroy_perfectreflector(Perfectreflector * pr);
void destroy_perfectreflectorarray(Perfectreflector ** pr_begin);

double _gaussrand(); // http://c-faq.com/lib/gaussian.html

//void add_node_to_environment_array(Environment * env, GeneralNode * gn);
//void print_node_path(Environment * env);
//void print_current_locations(Environment * env);
//double gaussenv(Environment * env);
void interaction_scatterer(void * sc, Transmitter * tx,
			   Transmitter * out_array, int * number);
void parse_input(int argc, char * argv[], Filereader * fr);
void cross_product(const double * v1, const double * v2, double * v3);

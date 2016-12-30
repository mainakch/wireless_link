#include "models.h"

int id()
{
        static int id = 0;
        return (id++);
}

void init_simulation(struct simulation *sim_not_null)
{
        sim_not_null->frequency = 60e9; // 60 Ghz 
        sim_not_null->wavelength = C / sim_not_null->frequency;
        sim_not_null->delta_time = 0.1;
        sim_not_null->max_limit = 1000;
        sim_not_null->min_limit = 0;
        sim_not_null->boundary_tolerance = 1;
        sim_not_null->total_time = 200;
}

void init_spatial_motion_model(struct spatial_motion_model *smm)
{
        int ctr;
        for (ctr = 0; ctr<3; ctr++) {
                smm->velocity[ctr] = 0;
                smm->position[ctr] = 5;
        }
        smm->distance_from_actual_source = 0;
        smm->acceleration_factor = 20;
        smm->is_static = false;
}

void init_transmission_model(struct transmission_model *tm)
{
        tm->power_in_dBm = -200;
        tm->start_time = 0;
        tm->end_time = 1000;
        tm->initial_phase = 0;
        tm->doppler_offset = 0;
}

void init_propagation_model(struct propagation_model *pm)
{
        pm->distance = 0;
}

void init_general_node(struct general_node *gn)
{
        init_spatial_motion_model(&(gn->smm));
        init_transmission_model(&(gn->tm));
        gn->id = id();
}

void init_transmitter(struct transmitter *tn)
{
        init_general_node(&(tn->gn));
        tn->is_real_transmitter = false;
}

void init_receiver(struct receiver *rc)
{
        init_general_node(&(rc->gn));
        rc->recv_noise_power = pow(10, -16);
}

struct perfect_reflector *init_perfect_reflector(const double *normal,
                                         const double *center_point,
                                         const double *length_normal,
                                         double length, double width)
{
        struct perfect_reflector *pr =
                calloc(1, sizeof(struct perfect_reflector));
        double norm = cblas_dnrm2(3, normal, 1);
        cblas_daxpy(3, 1 / norm, normal, 1, pr->unit_normal, 1);
        cblas_dcopy(3, center_point, 1, pr->center_point, 1);
        norm = cblas_dnrm2(3, length_normal, 1);
        cblas_daxpy(3, 1 / norm, length_normal, 1,
                    pr->unit_length_normal, 1);
        cross_product(pr->unit_length_normal,
                      pr->unit_normal, pr->unit_width_normal);
        pr->width = width;
        pr->length = length;
        return pr;
}

struct perfect_reflector **init_perfect_reflectorarray(int number)
{
        struct perfect_reflector **pr = calloc(number,
                                        sizeof(struct perfect_reflector *));
        return pr;
}

void init_environment(struct environment *env)
{
        if (!(env->num_transmitters > 0 && env->num_receivers > 0)
            || env->total_time < 1) { 
                fprintf(stderr,
                        "Please set non zero tx, rx and total time\n");
                exit(EXIT_FAILURE);
        }
        
        env->_num_transmitters_ctr = 0;
        env->_num_receivers_ctr = 0;
        env->num_virtual_transmitters = 0;
        env->time = 0;
        init_environment_malloc(env);
}

void init_environment_malloc(struct environment *env)
{
        fprintf(stderr, "Inside malloc for environment\n");
        env->receivers_array = calloc(env->num_receivers,
                                      sizeof(struct receiver));
        env->transmitters_array = calloc(env->num_transmitters,
                                         sizeof(struct transmitter));
        env->node_array = calloc((env->num_receivers +
                                  env->num_transmitters),
                                 sizeof(struct general_node *));
        env->unit_power_gaussian_noise = calloc((env->num_receivers),
                                                sizeof(complex double));
}

void init_filereader(struct filereader *fr)
{
        printf("Opening %s for writing \n", fr->output_filename);
        printf("Opening %s for reading \n", fr->input_filename); 
        
        fr->infile = fopen(fr->input_filename, "r");
        fr->outfile = fopen(fr->output_filename, "w");  
}

void destroy_environment(struct environment *env)
{
        free(env->receivers_array);
        free(env->transmitters_array);
        free(env->node_array);
        free(env->unit_power_gaussian_noise);
}

void destroy_filereader(struct filereader *fr)
{
        if (fr->infile != NULL)  fclose(fr->infile);
        if (fr->outfile != NULL)  fclose(fr->outfile);
}

void destroy_perfect_reflector(struct perfect_reflector *pr)
{
        free(pr);
}

void destroy_perfect_reflectorarray(struct perfect_reflector **pr_begin)
{
        int ctr = 0;
        while (*(pr_begin + ctr) != NULL) {
                destroy_perfect_reflector(*(pr_begin + ctr));
                ctr++;
        }
        free(pr_begin);
}

double _gaussrand() // http://c-faq.com/lib/gaussian.html
{
        static double U, V;
        static int phase = 0;
        double Z;
        
        if(phase == 0) {
                U = (rand() + 1.) / (RAND_MAX + 2.);
                V = rand() / (RAND_MAX + 1.);
                Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
        } else {
                Z = sqrt(-2 * log(U)) * cos(2 * PI * V);
        }
        
        phase = 1 - phase;
        return Z;
}

void interaction_scatterer(void *sc, struct transmitter *tx,
                           struct transmitter *out_array, int *number)
{
        // assume that sc is perfect_reflector
        struct perfect_reflector *pr = (struct perfect_reflector *) sc;
        
        // check if source is on the reflective side
        double pos[3];
        cblas_dcopy(3, tx->gn.smm.position, 1, pos, 1);
        cblas_daxpy(3, -1, pr->center_point, 1, pos, 1);
        if (cblas_ddot(3, pos, 1, pr->unit_normal, 1) > 0) {
                // facing non reflective side
                // return after setting number of virtual
                // sources added to zero
                *number = 0; 
                return;
        }
  
        // source is on the reflective side now
        // compute position transformation
        // pos contains difference of tx - center_point
        double factor0 = -2 * cblas_ddot(3, pos, 1, pr->unit_normal, 1);
        double reflected_position[3];
        cblas_dcopy(3, tx->gn.smm.position, 1, reflected_position, 1);
        cblas_daxpy(3, factor0, pr->unit_normal, 1,
                    reflected_position, 1);
        
        // compute velocity transformation
        double reflected_velocity[3];
        double factor1;
        cblas_dcopy(3, tx->gn.smm.velocity, 1, reflected_velocity, 1);
        factor1 = -2 * cblas_ddot(3, tx->gn.smm.velocity, 1,
                                pr->unit_normal, 1);
        cblas_daxpy(3, factor1, pr->unit_normal, 1,
                    reflected_velocity, 1);
  
        // compute power attenuation factor
        
        double norm = cblas_dnrm2(3, pos, 1);
        double cosangle = abs(factor0 / (2 * norm));
        double factor2 = pr->length * pr->width * cosangle/
                (4 * PI * norm * norm);
        
        // add virtual transmitter to the array
        // TODO
        *number = 1;
}

void parse_input(int argc, char *argv[], struct filereader *fr)
{
        strcpy(fr->input_filename, "/tmp/input.txt");
        strcpy(fr->output_filename, "/tmp/output.txt");
        
        struct option long_options[] = {
                {"input_filename", required_argument, 0, 'i'},
                {"output_filename", required_argument, 0, 'o'},
                {0, 0, 0, 0}
        };
        
        int c=0, option_index;
        c = getopt_long(argc, argv, "i:o:",
                        long_options, &option_index);
        while (c != -1) {
                switch (c) {
                case 0:
                        break;
                case 'i':
                        strcpy(fr->input_filename, optarg);
                        break;
                case 'o':
                        strcpy(fr->output_filename, optarg);
                        break;
                default:
                        fprintf(stderr, "Error parsing options\n");
                        exit(EXIT_FAILURE);
                }
                c = getopt_long(argc, argv, "i:o:",
                                long_options, &option_index);
        }
}

void cross_product(const double *v1, const double *v2, double *v3)
{
        v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
        v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
        v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

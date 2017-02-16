#include "raytracing.h"

static void add_reflector(struct environment *env, struct perfect_reflector *pr);
static void add_transmitter(struct environment *env, struct transmitter *tx);
static void add_receiver(struct environment *env, struct receiver *rx);
static void add_general_node(struct environment *env, struct general_node *gn);


int id()
{
        static int id = 0;
        return (id++);
}

struct simulation *init_simulation()
{
        struct simulation *sim_not_null = custom_malloc(sizeof(struct simulation));
        sim_not_null->frequency = 60e9; // 60 Ghz
        sim_not_null->wavelength = C / sim_not_null->frequency;
        sim_not_null->delta_time = 0.1;
        sim_not_null->max_limit = 1000;
        sim_not_null->min_limit = 0;
        sim_not_null->boundary_tolerance = 1;
        sim_not_null->total_time = 200;
        return sim_not_null;
}

struct spatial_motion_model *init_spatial_motion_model()
{
        struct spatial_motion_model *smm = custom_calloc(1,
                sizeof(struct spatial_motion_model));
        return smm;
}

struct transmission_model *init_transmission_model()
{
        struct transmission_model *tm = custom_malloc(
                sizeof(struct transmission_model));
        tm->power_in_dBm = -200;
        return tm;
}

struct propagation_model *init_propagation_model()
{
        struct propagation_model *pm = custom_calloc(1,
                sizeof(struct propagation_model));
        return pm;
}

struct general_node *init_general_node()
{
        struct general_node *gn = custom_malloc(sizeof(struct general_node));
        gn->smm = init_spatial_motion_model();
        gn->tm = init_transmission_model();
        gn->id = id();
        return gn;
}

struct transmitter *init_transmitter()
{
        struct transmitter *tn = custom_malloc(sizeof(struct transmitter));
        tn->gn = init_general_node();
        tn->baseband_signal = 1;
        return tn;
}

struct receiver *init_receiver()
{
        struct receiver *rc = custom_malloc(sizeof(struct receiver));
        rc->gn = init_general_node();
        rc->recv_noise_power = pow(10, -16);
        rc->rlln = 0;//custom_calloc(16, sizeof(struct receiver_ray_ribbon *));
        return rc;
}

struct perfect_reflector *init_perfect_reflector(const double *normal,
                                         const double *center_point,
                                         const double *length_normal,
                                         double length, double width)
{
        struct perfect_reflector *pr =
                custom_calloc(1, sizeof(struct perfect_reflector));
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

struct perfect_reflector *init_perfect_reflector_nine_pts_direction(
        const double *pt1, const double *pt2,
        const double *pt3, bool direction) {
        // find out normal
        double v1[3];
        double v2[3];
        double vnormal[3] = {0, 0, 0};
        double center_pt[3];
        diff(pt1, pt2, v1);
        diff(pt3, pt2, v2);
        cross_product(v1, v2, vnormal);
        normalize_unit_vector(vnormal);
        // handle user supplied input if pt1 is at origin
        double out = cblas_ddot(3, vnormal, 1, pt1, 1);
        if (out * direction > 0) {
                cblas_dscal(3, -1, vnormal, 1);
        }
        // compute center point
        cblas_dcopy(3, pt1, 1, center_pt, 1);
        cblas_daxpy(3, 1, pt3, 1, center_pt, 1);
        cblas_dscal(3, 0.5, center_pt, 1);
        double length = cblas_dnrm2(3, v1, 1);
        double width = cblas_dnrm2(3, v2, 1);
        normalize_unit_vector(v1);
        return init_perfect_reflector(vnormal, center_pt, v1, length, width);
}

struct perfect_reflector **init_perfect_reflectorarray(int number)
{
        struct perfect_reflector **pr = custom_calloc(number,
                                        sizeof(struct perfect_reflector *));
        return pr;
}

struct environment *init_environment()
{
        struct environment *env = custom_calloc(1, sizeof(struct environment));
        env->receivers_array = custom_calloc(MAX_RX + 1, sizeof(struct receiver *));
        env->transmitters_array = custom_calloc(MAX_TX + 1,
                                         sizeof(struct transmitter *));
        env->node_array = custom_calloc(MAX_TX + MAX_RX + 1,
                                 sizeof(struct general_node *));
        // env->unit_power_gaussian_noise = custom_calloc(11, sizeof(complex double));
        env->prarray = custom_calloc(MAX_SURFACES + 2, sizeof(struct perfect_reflector *));
        env->env_paths = custom_calloc(MAX_RX + 1, sizeof(struct ray_ribbon_array *));
        // size of env paths is the same as the size of receivers_array
        env->tx_paths = custom_calloc(MAX_TX + 1, sizeof(struct ray_ribbon_array *));
        // size of tx paths is the same as the size of transmitters_array

        env->sz_array_tx = MAX_TX + 1;
        env->sz_array_rx = MAX_RX + 1;
        env->sz_array_gn = MAX_TX + MAX_RX + 1;
        env->sz_array_pr = MAX_SURFACES + 2;
        env->updated_tx_paths = false;
        return env;
}

struct file_reader *init_file_reader(const char *input, const char *output)
{
        struct file_reader *fr = custom_calloc(1, sizeof(struct file_reader));

        printf("Opening %s for writing \n", output);
        printf("Opening %s for reading \n", input);

        fr->infile = fopen(input, "r");
        fr->outfile = fopen(output, "w");
        strncpy(fr->input_filename, input, 999);
        fr->input_filename[999] = 0;
        strncpy(fr->output_filename, output, 999);
        fr->output_filename[999] = 0;
        return fr;
}

void print_vector(const double *db)
{
        int ctr;
        for(ctr = 0; ctr < 3; ctr++) {
                fprintf(stderr, "%lf ", *(db + ctr));
        }
        fprintf(stderr, "\n");
}

void print_spatial_motion_model(const struct spatial_motion_model *smm) {
        fprintf(stderr, ANSI_COLOR_BLUE "Velocity: \n");
        print_vector(smm->velocity);
        fprintf(stderr, ANSI_COLOR_RESET);
        fprintf(stderr, "Position: \n");
        print_vector(smm->position);
}

void print_transmission_model(const struct transmission_model *tm) {
        fprintf(stderr, "Transmission model: \n");
        fprintf(stderr, "Power in dBm: %lf\n", tm->power_in_dBm);
}

void print_propagation_model(const struct propagation_model *pm) {
        fprintf(stderr, "%lf\n", pm->distance);
}

void print_general_node(const struct general_node *gn) {
        fprintf(stderr, "Printing general node id %d:\n", gn->id);
        print_spatial_motion_model(gn->smm);
        print_transmission_model(gn->tm);
}

void print_transmitter(const struct transmitter *tx) {
        fprintf(stderr, "Printing transmitter with id %d:\n", tx->gn->id);
        print_general_node(tx->gn);
}

void print_receiver(const struct receiver *rx) {
        fprintf(stderr, "Printing receiver with id %d:\n", rx->gn->id);
        print_general_node(rx->gn);
        fprintf(stderr, "Receiver noise power: %lf\n", rx->recv_noise_power);
}

void print_perfect_reflectors(const struct perfect_reflector *pr) {
        fprintf(stderr, "Printing perfect reflector ...\n");
        fprintf(stderr, ANSI_COLOR_RED "Printing unit normal\n");
        print_vector(pr->unit_normal);
        fprintf(stderr, ANSI_COLOR_RESET);
        fprintf(stderr, "Printing unit length normal\n");
        print_vector(pr->unit_length_normal);
        fprintf(stderr, "Printing unit width normal\n");
        print_vector(pr->unit_width_normal);
        fprintf(stderr, "Printing center point\n");
        print_vector(pr->center_point);
        fprintf(stderr, "Length: %lf, Width: %lf\n", pr->length, pr->width);
}

void print_environment(const struct environment *env) {
        fprintf(stderr, ANSI_COLOR_RED "At transmitter time %lf: \n"
                ANSI_COLOR_RESET, env->time);
        fprintf(stderr, "Printing environment\n");
        fprintf(stderr, "Printing transmitters:\n");
        int ctr = 0;
        while (*(env->transmitters_array + ctr) != 0) {
                print_transmitter(*(env->transmitters_array + ctr));
                ctr++;
        }
        ctr = 0;
        fprintf(stderr, "Printing receivers:\n");
        while (*(env->receivers_array + ctr) != 0) {
                print_receiver(*(env->receivers_array + ctr));
                ctr++;
        }
        ctr = 0;
        fprintf(stderr, "Printing perfect reflectors:\n");
        while (*(env->prarray + ctr) != 0) {
                print_perfect_reflectors(*(env->prarray + ctr));
                ctr++;
        }
        fprintf(stderr, "Wavelength: %lf, Frequency: %lf\n", env->wavelength,
                env->frequency);

        fprintf(stderr, ANSI_COLOR_GREEN "Printing tx paths \n");
        ctr = 0;
        while (*(env->tx_paths + ctr) != 0) {
                print_ray_ribbon_array(*(env->tx_paths + ctr));
                ctr++;
        }
        fprintf(stderr, ANSI_COLOR_GREEN "Printing env paths \n");
        ctr = 0;
        while (*(env->env_paths + ctr) != 0) {
                print_ray_ribbon_array(*(env->env_paths + ctr));
                ctr++;
        }
        fprintf(stderr, ANSI_COLOR_RESET);
}

void print_env_paths(const struct environment *env)
{
        int ctr = 0;
        // ctr loops over receivers
        struct ray_ribbon_array *rba = *(env->env_paths);
        fprintf(stderr, ANSI_COLOR_RED);
        while (rba != 0) {
                int ctr1 = 0;
                fprintf(stderr, "\n\n Printing for receiver %d\n\n",
                        ctr);
                // ctr1 loops over rays
                while (*(rba->ribbons + ctr1) != 0) {
                        fprintf(stderr, "Time: %lf; ", env->time);
                        fprintf(stderr, "Doppler: %lf; ",
                                (*(rba->ribbons + ctr1))->doppler);
                        print_ray_ribbon_flattened(
                                *(rba->ribbons + ctr1));
                        ctr1++;
                }
                ctr++;
                rba = *(env->env_paths + ctr);
        }
        fprintf(stderr, ANSI_COLOR_RESET);
}

void print_tx_paths(const struct environment *env)
{
        int ctr = 0;
        struct ray_ribbon_array *rba = *(env->tx_paths);
        fprintf(stderr, ANSI_COLOR_BLUE);

        while (ctr < env->num_transmitters) {
                fprintf(stderr, "\n\n Printing transmitter %d:\n\n", ctr);
                int ctr1 = 0;
                while (*(rba->ribbons + ctr1) != 0) {
                        fprintf(stderr, "Time: %lf, ", env->time);
                        fprintf(stderr, "Doppler: %lf, ",
                                (*(rba->ribbons + ctr1))->doppler);
                        print_ray_ribbon_flattened(
                                *(rba->ribbons + ctr1));
                        ctr1++;
                }
                ctr++;
                rba = *(env->tx_paths + ctr);
        }
        fprintf(stderr, ANSI_COLOR_RESET);
}

void destroy_spatial_motion_model(struct spatial_motion_model *smm)
{
        custom_free(smm);
}

void destroy_simulation(struct simulation *sim)
{
        custom_free(sim);
}

void destroy_transmission_model(struct transmission_model *tm)
{
        custom_free(tm);
}

void destroy_propagation_model(struct propagation_model *pm)
{
        custom_free(pm);
}

void destroy_general_node(struct general_node *gn)
{
        custom_free(gn->smm);
        custom_free(gn->tm);
        custom_free(gn);
}

void destroy_transmitter(struct transmitter *tn)
{
        destroy_general_node(tn->gn);
        custom_free(tn);
}

void destroy_receiver(struct receiver *rc)
{
        destroy_general_node(rc->gn);
        struct receiver_ray_ribbon_ll_node *rlln = rc->rlln;
        int ctr = 1;
        while (rlln != 0) {
                struct receiver_ray_ribbon_ll_node *rll = rlln->next;
                destroy_receiver_ray_ribbon_ll_node(rlln);
                rlln = rll;
        }
        custom_free(rc);
}

void destroy_environment(struct environment *env)
{
        int ctr = 0;
        struct transmitter *tn = *(env->transmitters_array + ctr);
        while (tn != NULL) {
                destroy_transmitter(tn);
                destroy_ray_ribbon_array(*(env->tx_paths + ctr));
                ctr++;
                tn = *(env->transmitters_array + ctr);
        }
        ctr = 0;
        struct receiver *rx = *(env->receivers_array + ctr);
        while (rx != NULL) {
                destroy_receiver(rx);
                destroy_ray_ribbon_array(*(env->env_paths + ctr));
                ctr++;
                rx = *(env->receivers_array + ctr);
        }
        ctr = 0;
        struct perfect_reflector *pr = *(env->prarray + ctr);
        while (pr != NULL) {
                destroy_perfect_reflector(pr);
                ctr++;
                pr = *(env->prarray + ctr);
        }

        custom_free(env->receivers_array);
        custom_free(env->transmitters_array);
        custom_free(env->node_array);
        custom_free(env->unit_power_gaussian_noise);
        custom_free(env->prarray);
        custom_free(env->env_paths);
        custom_free(env->tx_paths);
        custom_free(env);
}

void destroy_file_reader(struct file_reader *fr)
{
        if (fr->infile != NULL)  fclose(fr->infile);
        if (fr->outfile != NULL)  fclose(fr->outfile);
        custom_free(fr);
}

void destroy_perfect_reflector(struct perfect_reflector *pr)
{
        custom_free(pr);
}

void destroy_perfect_reflectorarray(struct perfect_reflector **pr_begin)
{
        int ctr = 0;
        while (*(pr_begin + ctr) != NULL) {
                destroy_perfect_reflector(*(pr_begin + ctr));
                ctr++;
        }
        custom_free(pr_begin);
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

struct file_reader *parse_input(int argc, char *argv[])
{
        // initialized to zero in C99
        char inputname[1000];
        char outputname[1000];
        strcpy(inputname, "/tmp/input.txt");
        strcpy(outputname, "/tmp/output.txt");

        struct option long_options[] = {
                {"input_filename", required_argument, 0, 'i'},
                {"output_filename", required_argument, 0, 'o'},
                {0, 0, 0, 0}
        };

        int c = 0, option_index;
        c = getopt_long(argc, argv, "i:o:",
                        long_options, &option_index);
        while (c != -1) {
                switch (c) {
                case 0:
                        break;
                case 'i':
                        strncpy(inputname, optarg, 999);
                        break;
                case 'o':
                        strncpy(outputname, optarg, 999);
                        break;
                default:
                        fprintf(stderr, "Error parsing options\n");
                        exit(EXIT_FAILURE);
                }
                c = getopt_long(argc, argv, "i:o:",
                                long_options, &option_index);
        }

        return init_file_reader(inputname, outputname);
}

void cross_product(const double *v1, const double *v2, double *v3)
{
        v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
        v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
        v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

double normalize_unit_vector(double *v1)
{
        double normv1 =  cblas_dnrm2(3, v1, 1);
        cblas_dscal(3, 1/normv1, v1, 1);
        return normv1;
}

void diff(const double *v1, const double *v2, double *v3)
{
        cblas_dcopy(3, v2, 1, v3, 1);
        cblas_daxpy(3, -1, v1, 1, v3, 1);
}

int find_len(void **ptr)
{
        int ctr = 0;
        while (*(ptr + ctr) != 0) {
                ctr++;
        }
        return ctr;
}

void add_receiver_patch(struct environment *env, int length)
{
        double vec[3];
        double vec_perp[3] = {0, 0, 0};
        vec[0] = rand(); vec[1] = rand(); vec[2] = rand();
        cross_product(vec, env->recv_unit_normal, vec_perp);
        normalize_unit_vector(vec_perp);
        double *position = (*(env->receivers_array))->gn->smm->position;
        //print_vector(position);
        struct perfect_reflector *new_reflector =
                init_perfect_reflector(
                        env->recv_unit_normal,
                        (*(env->receivers_array))->gn->smm->position,
                        vec_perp,
                        length, length);

        if (_MODEL_DEBUG) {
                fprintf(stderr, "Studenchdnehth\n");
                print_perfect_reflectors(new_reflector);
        }
        // make this value sufficiently large 5, 5 so that it can catch
        // all rays, setting it too high may cause problems with non linearity

        add_reflector(env, new_reflector);
}

void destroy_last_reflector(struct environment *env)
{
        struct perfect_reflector *pr = *(env->prarray + env->num_reflectors - 1);
        destroy_perfect_reflector(pr);
        env->num_reflectors--;
        *(env->prarray + env->num_reflectors) = 0;
}

static void add_reflector(struct environment *env, struct perfect_reflector *pr)
{
        /* if (env->num_reflectors >= env->sz_array_pr) { */
        /*         int new_size = 2 * env->sz_array_pr; */
        /*         void *tmp = realloc(env->prarray, */
        /*                       new_size * sizeof(struct perfect_reflector *)); */
        /*         if (tmp == 0) exit(1); */
        /*         env->prarray = tmp; */
        /*         env->sz_array_pr = new_size; */
        /* } */
        *(env->prarray + env->num_reflectors) = pr;
        env->num_reflectors++;
        *(env->prarray + env->num_reflectors) = 0;
}

static void add_transmitter(struct environment *env, struct transmitter *tx)
{
        /* if (env->num_transmitters + 1 >= env->sz_array_tx) { */
        /*         int new_size = 2 * env->sz_array_tx; */
        /*         void *tmp = realloc(env->transmitters_array, */
        /*                       new_size * sizeof(struct transmitter *)); */
        /*         if (tmp == 0) exit(1); */
        /*         env->transmitters_array = tmp; */
        /*         env->sz_array_tx = new_size; */
        /* } */
        *(env->transmitters_array + env->num_transmitters) = tx;
        env->num_transmitters++;
        *(env->transmitters_array + env->num_transmitters) = 0;
        add_general_node(env, tx->gn);
}

static void add_receiver(struct environment *env, struct receiver *rx)
{
        /* if (env->num_receivers >= env->sz_array_tx) { */
        /*         int new_size = 2 * env->sz_array_tx; */
        /*         void *tmp = realloc(env->receivers_array, */
        /*                       new_size * sizeof(struct receiver *)); */
        /*         if (tmp == 0) exit(1); */
        /*         env->receivers_array = tmp; */
        /*         *(env->receivers_array + env->num_receivers + 1) = 0; */
        /*         env->sz_array_tx = new_size; */
        /* } */
        *(env->receivers_array + env->num_receivers) = rx;
        env->num_receivers++;
        // fprintf(stderr, "Nulling %d\n", env->num_receivers);
        *(env->receivers_array + env->num_receivers) = 0;
        add_general_node(env, rx->gn);
}

static void add_general_node(struct environment *env, struct general_node *gn)
{
        int length = env->num_transmitters + env->num_receivers;
        /* if (length >= env->sz_array_gn) { */
        /*         int new_size = 2 * env->sz_array_gn; */
        /*         void *tmp = realloc(env->node_array, */
        /*                       new_size * sizeof(struct general_node *)); */
        /*         if (tmp == 0) exit(1); */
        /*         env->node_array = tmp; */
        /*         env->sz_array_gn = length; */
        /* } */
        *(env->node_array + length - 1) = gn;
        *(env->node_array + length) = 0;
}


double distance(const struct general_node *gn1, const struct general_node *gn2)
{
        double *pos1 = gn1->smm->position;
        double *pos2 = gn2->smm->position;

        double diff[3];
        cblas_dcopy(3, pos1, 1, diff, 1);
        cblas_daxpy(3, -1, pos2, 1, diff, 1);

        return cblas_dnrm2(3, diff, 1);
}

bool update_environment_from_file(struct environment *env, FILE *fp)
{
        if (fp == NULL)
        {
                fprintf(stderr, "Fatal error: input file not provided!\n");
                exit(EXIT_FAILURE); // fix potential memory errors
        }

        char mode[] = "r";
        char fmt[] = "%s"; // check for vulnerabilities
        char buff[100];
        bool eofflag = true;
        while (fscanf(fp, fmt, buff) != EOF) {
                eofflag = false;
                handle_request(env, fp, buff);
                if (strcmp(buff, "End") == 0) break;
        }
        return eofflag;
}

void handle_request(struct environment *env, FILE *fp, const char *req_type)
{
        int ctr;
        if (!strcmp(req_type, "Nodepath")) {
                int id;
                fscanf(fp, "%d", &id);
                if (_MODEL_DEBUG) fprintf(stderr, "id: %d\n", id);

                for (ctr=0; ctr<3; ctr++) {
                        fscanf(fp, "%lf",
                               &(env->node_array[id]->smm->position[ctr]));
                        fscanf(fp, "%lf",
                               &(env->node_array[id]->smm->velocity[ctr]));
                }
        } else if (!strcmp(req_type, "Gaussianrand")) {
                for (ctr=0; ctr < env->num_receivers; ctr++) {
                        double re, im;
                        fscanf(fp, "%lf", &re);
                        fscanf(fp, "%lf", &im);
                        env->unit_power_gaussian_noise[ctr] =
                                pow(0.5, -0.5)*(re + I*im);
                }
        } else if (!strcmp(req_type, "Time")) {
                fscanf(fp, "%lf", &(env->time));
        } else if (!strcmp(req_type, "End")) {
                if (!(env->read_in_nodes)) {
                        env->read_in_nodes = 1;
                        env->unit_power_gaussian_noise =
                                custom_calloc(env->num_receivers + 1,
                                       sizeof(complex double));
                }
        } else if (!strcmp(req_type, "Transmittersignal")) {
                int id;
                fscanf(fp, "%d", &id);
                double re, im;
                fscanf(fp, "%lf", &re);
                fscanf(fp, "%lf", &im);
                struct transmitter *tx;
                tx = *(env->transmitters_array + id);
                tx->baseband_signal = re + I * im;
        } else if (!strcmp(req_type, "Frequency")) {
                fscanf(fp, "%lf", &env->frequency);
                env->wavelength = C/env->frequency;
        } else if (!strcmp(req_type, "Transmitter")) {
                struct transmitter *tx = init_transmitter();
                fscanf(fp, "%d", &tx->gn->id);
                fscanf(fp, "%lf", &tx->gn->tm->power_in_dBm);
                add_transmitter(env, tx);
        } else if (!strcmp(req_type, "Receiver")) {
                struct receiver *rx = init_receiver();
                fscanf(fp, "%d", &rx->gn->id);
                double np;
                fscanf(fp, "%lf", &np);
                rx->recv_noise_power = pow(10, np/10);
                add_receiver(env, rx);
        } else if (!strcmp(req_type, "Perfectreflector")) {
                double pt1[3][3];
                double direction;
                int ctr1;
                for (ctr1 = 0; ctr1 < 3; ctr1++) {
                        for (ctr = 0; ctr < 3; ctr++) {
                                fscanf(fp, "%lf", &pt1[ctr1][ctr]);
                        }
                }
                fscanf(fp, "%lf", &direction);
                add_reflector(env, init_perfect_reflector_nine_pts_direction(
                                      pt1[0], pt1[1], pt1[2],  direction));
        } else if (!strcmp(req_type, "Frequency")) {
                fscanf(fp, "%lf", &(env->frequency));
                env->wavelength = C/env->frequency;
        } else if (!strcmp(req_type, "Receivernormal")) {
	        double v1[3];
		for(ctr = 0; ctr < 3; ctr++) {
		        fscanf(fp, "%lf", &v1[ctr]);
		}
                cblas_dcopy(3, v1, 1, env->recv_unit_normal, 1);
                normalize_unit_vector(env->recv_unit_normal);
	} else if (!strcmp(req_type, "Deltatime")) {
                fscanf(fp, "%lf", &(env->delta_time));
        } else if (!strcmp(req_type, "Endtime")) {
                fscanf(fp, "%lf", &(env->end_time));
        }
}

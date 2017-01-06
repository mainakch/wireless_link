#include "raytracing.h"

static struct ribbon_node *_refine_ribbon_node(struct ribbon_node **node_array,
                                      const double *point,
                                      struct ribbon_node *rn,
                                      bool null_node_array,
                                      const struct perfect_reflector **pr);

int counter()
{
        static int ctr = 0;
        return ctr++;
}

struct set_of_types *init_set_of_types(int number)
{
        struct set_of_types *st = malloc(sizeof(struct set_of_types));
        st->types = malloc(number * sizeof(long));

        // initialize to -1
        int ctr = 0;
        while (ctr < number) {
                *(st->types + ctr) = -1;
                ctr++;
        }
        return st;
}

void destroy_set_of_types(struct set_of_types *st)
{
        free(st->types);
        free(st);
}

void add_entry(long number, struct set_of_types *st)
{
        int ctr = 0;
        while (ctr < NUM_TYPES && *(st->types + ctr) >= 0) {
                if (*(st->types + ctr) == number) {
                        return;
                }
                ctr++;
        }

        // does not check for max length
        *(st->types + ctr) = number;
}

struct ray_ribbon_array *init_ray_ribbon_array(int number)
{
        struct ray_ribbon_array *rarr =
                calloc(1, sizeof(struct ray_ribbon_array));
        rarr->ribbons = calloc(number, sizeof(struct ray_ribbon *));
        rarr->max_len = number;
        rarr->current_len = 0;
        return rarr;
}

void populate_ray_ribbon_array(struct transmitter *tx,
                               const struct perfect_reflector **ref_arr,
                               struct ray_ribbon_array *rarr,
                               int num_divs, int num_ref, bool single_type)
{
        populate_ray_ribbon_array_long(tx, ref_arr, rarr, num_ref,
                                       -PI / 2, PI / 2, PI / num_divs,
                                       0, 2 * PI, 2 * PI / num_divs,
                                       single_type);
}

void populate_ray_ribbon_array_long(struct transmitter *tx,
                                    const struct perfect_reflector **ref_arr,
                                    struct ray_ribbon_array *rarr,
                                    int num_ref,
                                    const double phi_start,
                                    const double phi_end,
                                    const double phi_delta,
                                    const double thet_start,
                                    const double thet_end,
                                    const double thet_delt,
                                    bool single_type)
{
        int num_points = (1 + floor((phi_end - phi_start) / phi_delta)) *
                (1 + floor((thet_end - thet_start) / thet_delt));
        fprintf(stderr, "Number is %d\n", num_points);
        complex double *angles = malloc(num_points * sizeof(complex double));

        int ctr = 0;

        double phi, theta;
        for (phi = phi_start; phi < phi_end; phi += phi_delta) {
                for (theta = thet_start; theta < thet_end; theta += thet_delt) {
                        *(angles + ctr) = phi + I * theta;
                        ctr++;
                }
        }
        populate_ray_ribbon_array_full(tx, ref_arr, num_ref, ctr,
                                       angles, rarr, single_type);
        free(angles);
}

void populate_ray_ribbon_array_full(const struct transmitter *tx,
                                    const struct perfect_reflector **ref_arr,
                                    int num_ref, int num_points,
                                    const complex double *angles,
                                    struct ray_ribbon_array *rarr,
                                    bool single_type)
{
        int ctr = 0;
        double phi, theta;
        while (ctr < num_points) {
                phi = creal(*(angles + ctr));
                theta = cimag(*(angles + ctr));
                struct ribbon_node *rn = init_ribbon_node();
                cblas_dcopy(3, tx->gn->smm->position, 1,
                            rn->current->point, 1);
                double direction[3] = {cos(phi),
                                       sin(phi) * cos(theta),
                                       sin(phi) * sin(theta)};
                cblas_dcopy(3, direction, 1,
                            rn->current->unit_direction, 1);

                bool hit_des = process_vertical_chain(rn,
                                                      ref_arr,
                                                      num_ref);

                if (hit_des) {
                        struct ray_ribbon *rb = init_ray_ribbon(rn);
                        bool ribbon_added = add_ray_ribbon(rarr, rb,
                                                           single_type);
                        if (!ribbon_added) destroy_ray_ribbon(rb);

                } else {
                        destroy_chain_of_ribbon_nodes(rn);
                }
                ctr++;
        }
}

struct ray_ribbon *init_ray_ribbon(struct ribbon_node *rn)
{
        struct ray_ribbon *rb = calloc(1, sizeof(struct ray_ribbon));
        rb->head = rn;
        return rb;
}

struct ribbon_node *init_ribbon_node()
{
        struct ribbon_node *rn = malloc(sizeof(struct ribbon_node));
        rn->current = malloc(sizeof(struct half_infinite_ray));
        rn->down = 0;
        rn->ctr = counter();
        rn->num_reflections = 0;
        rn->hit_destination_patch = false;
        rn->surface_index = -1;
        return rn;
}
void destroy_ray_ribbon(struct ray_ribbon *rb)
{
        destroy_ray_ribbon_nodes(rb);
        free(rb);
}

void destroy_ray_ribbon_nodes(struct ray_ribbon *rb)
{
        if (rb == 0) return;
        struct ribbon_node *rn = rb->head;
        destroy_chain_of_ribbon_nodes(rn);
        rb->head = 0;
}

void destroy_ray_ribbon_array(struct ray_ribbon_array *array)
{
        if (array == 0) return;
        int ctr = 0;
        while (*(array->ribbons + ctr) != NULL) {
                destroy_ray_ribbon_nodes(*(array->ribbons + ctr));
                free(*(array->ribbons + ctr));
                ctr++;
        }
        free(array->ribbons);
        free(array);
}

void destroy_ray_ribbon_array_all_but_first(struct ray_ribbon_array *array)
{
        if (array == 0) return;
        int ctr = 1;
        while (*(array->ribbons + ctr) != NULL) {
                destroy_ray_ribbon_nodes(*(array->ribbons + ctr));
                free(*(array->ribbons + ctr));
                ctr++;
        }
        free(array->ribbons);
        free(array);
}

void destroy_ray_ribbon_array_ribbons(struct ray_ribbon_array *array)
{
        int ctr = 0;
        if (array == 0) return;
        while (*(array->ribbons + ctr) != NULL) {
                destroy_ray_ribbon(*(array->ribbons + ctr));
                ctr++;
        }
        free(array->ribbons);
        array->max_len = 0;
        array->current_len = 0;
        array->ribbons = 0;
}

void destroy_chain_of_ribbon_nodes(struct ribbon_node *rn)
{
        if (rn == NULL) return;
        if (rn->current != NULL) free(rn->current);
        struct ribbon_node *rndown = rn->down;
        free(rn);
        destroy_chain_of_ribbon_nodes(rndown);
}

bool check_same_type(const struct ray_ribbon *ray_rb1,
                     const struct ray_ribbon *ray_rb2)
{
        struct ribbon_node *rn1 = ray_rb1->head;
        struct ribbon_node *rn2 = ray_rb2->head;

        while(rn1 != NULL && rn2 != NULL) {
                if (rn1->surface_index != rn2->surface_index) return false;
                rn1 = rn1->down;
                rn2 = rn2->down;
        }
        return (rn1 == NULL && rn2 == NULL);
}

bool add_ray_ribbon(struct ray_ribbon_array *array, struct ray_ribbon *rb,
                    bool single_type)
{
        if (array->current_len + 1 >= array->max_len) {
                // automatic resizing
                struct ray_ribbon **tmp =
                        realloc(array->ribbons, 2 * array->max_len);
                if (tmp == NULL) {
                        return false;
                } else {
                        array->max_len *= 2;
                        array->ribbons = tmp;
                }
        }

        if (single_type) {
                int ctr = 0;
                // check for type here
                while (*(array->ribbons + ctr) != 0) {
                        if(check_same_type(*(array->ribbons + ctr), rb)) {
                                return false;
                        }
                        ctr++;
                }
        }

        *(array->ribbons + array->current_len) = rb;
        array->current_len++;
        // this enforces null termination
        *(array->ribbons + array->current_len) = 0;
        return true;

}

complex double compute_intersection(struct half_infinite_ray *hr,
                                    const struct perfect_reflector *pr)
{
        double t, sgn;
        double diff[3];
        cblas_dcopy(3, hr->point, 1, diff, 1);
        cblas_daxpy(3, -1, pr->center_point, 1, diff, 1);
        t = -cblas_ddot(3, diff, 1, pr->unit_normal, 1)
                /cblas_ddot(3, hr->unit_direction, 1, pr->unit_normal, 1);
        //fprintf(stderr, "t obtained is %lf\n", t);

        // check if t lies within the bounds of the patch

        if (t < INFINITY) {
                // Verify the signs
                cblas_daxpy(3, t, hr->unit_direction, 1, diff, 1);

                //print_vector(diff);
                //print_vector(pr->unit_length_normal);
                //print_vector(pr->unit_width_normal);

                double lengtht = cblas_ddot(3, diff, 1,
                                            pr->unit_length_normal, 1);
                double widtht = cblas_ddot(3, diff, 1,
                                           pr->unit_width_normal, 1);

                if (fabs(lengtht) > pr->length / 2 ||
                    fabs(widtht) > pr->width / 2 ||
                    t < 0) {
                        t = INFINITY;
                }
        }

        if (t < 1e4) {
                cblas_dcopy(3, hr->point, 1, hr->end_pt, 1);
                cblas_daxpy(3, t, hr->unit_direction, 1, hr->end_pt, 1);
        }

        sgn = cblas_ddot(3, hr->unit_direction, 1, pr->unit_normal, 1);
        return t + I * sgn; // if sgn positive then ray is blocked
}

bool process_vertical_chain(struct ribbon_node *rn,
                            const struct perfect_reflector **pr,
                            int num_reflections)
{
        // this function computes whether a ray can hit the
        // destination after a max num_reflections
        int ctr=0, ctrindex=-1, num_reflectors=0;
        double tmin = INFINITY, sgn = -1;

        const struct perfect_reflector *prsurf = *(pr + ctr);
        while(prsurf != NULL) {
                complex double dbl = compute_intersection(rn->current, prsurf);
                if (creal(dbl) < tmin && ctr != rn->surface_index) {
                        tmin = creal(dbl);
                        sgn = cimag(dbl);
                        ctrindex = ctr;
                }
                ctr++;
                prsurf = *(pr + ctr);
        }
        num_reflectors = ctr;

        if (sgn>0 || tmin>1e4) return false;
        if (ctrindex == num_reflectors - 1) {
                rn->hit_destination_patch = true;
                return true;
        }

        if (rn->num_reflections > num_reflections) return false;

        // only case remaining is if there is intersection with
        // reflector and number of reflections is small
        struct ribbon_node *rn_next = malloc(sizeof(struct ribbon_node));
        rn_next->down = 0;
        rn_next->current = calloc(1, sizeof(struct half_infinite_ray));
        rn_next->hit_destination_patch = false;
        rn_next->num_reflections = rn->num_reflections + 1;

        cblas_dcopy(3, rn->current->point, 1,
                    rn_next->current->point, 1);
        cblas_daxpy(3, tmin, rn->current->unit_direction, 1,
                    rn_next->current->point, 1);
        rn_next->surface_index = ctrindex;

        // next update direction
        cblas_dcopy(3, rn->current->unit_direction, 1,
                    rn_next->current->unit_direction, 1);
        const struct perfect_reflector *prsurface = pr[ctrindex];
        double factor = -2*cblas_ddot(3, rn->current->unit_direction,
                                      1, prsurface->unit_normal, 1);
        cblas_daxpy(3, factor, prsurface->unit_normal, 1,
                    rn_next->current->unit_direction, 1);

        // update pointers
        rn->down = rn_next;
        return process_vertical_chain(rn_next, pr, num_reflections);
}

void print_ray_ribbon(const struct ray_ribbon *rb)
{
        fprintf(stderr, "Printing rayribbon: \n\n");
        struct ribbon_node *rn = rb->head;
        int ctr = 0;
        while(rn != NULL) {
                fprintf(stderr, "Level %d:\n", ctr);
                print_ribbon_node(rn);
                rn = rn->down;
                ctr++;
        }
        fprintf(stderr, "Delay, doppler for rayribbon are: %lf ns, %lf\n",
                (10e9) * rb->delay, rb->doppler);

}

void print_ray_ribbon_array(const struct ray_ribbon_array *rarr)
{
        struct ray_ribbon * rb;
        rb = *(rarr->ribbons);
        int ctr = 0;
        while (rb != NULL) {
                fprintf(stderr, "Printing ribbon %d:\n", ctr);
                print_ray_ribbon(rb);
                ctr++;
                rb = *(rarr->ribbons + ctr);
        }
}

void print_vertical_strip(const struct ribbon_node *rn)
{
        int ctr = 0;
        while (rn != NULL) {
                fprintf(stderr, "Level %d\n", ctr++);
                print_ribbon_node(rn);
                rn = rn->down;
        }
}

void print_ribbon_node(const struct ribbon_node *rn)
{
        if (rn == NULL) return;
        fprintf(stderr, "Starting point: ");
        int ctr = 0;
        for(ctr = 0; ctr < 3; ctr++) {
                fprintf(stderr, "%lf ", rn->current->point[ctr]);
        }

        fprintf(stderr, "Unit direction: ");
        for(ctr = 0; ctr < 3; ctr++) {
                fprintf(stderr, "%lf ",
                        rn->current->unit_direction[ctr]);
        }

        fprintf(stderr, "Ending point: ");
        for(ctr = 0; ctr < 3; ctr++) {
                fprintf(stderr, "%lf ", rn->current->end_pt[ctr]);
        }

        fprintf(stderr,
                "Hit dest: %d, num reflec: %d, Surf index: %d",
                rn->hit_destination_patch,
                rn->num_reflections, rn->surface_index);
        fprintf(stderr, "\n");
}

int count_segments(const struct ribbon_node *rn)
{
        int ctr = 0;
        while (rn != NULL) {
                rn = rn->down;
                ctr++;
        }
        return ctr;
}

void compute_average_ribbon_node(struct ribbon_node *rn,
                                const struct ray_ribbon_array *rba,
                                double *weights)
{
        double thet, phi;
        double thetaav = 0, phiav = 0;
        int ctr = 0;
        struct ribbon_node *node;
        for (; ctr < 3; ctr++) {
                node = (*(rba->ribbons + ctr))->head;
                invert_spherical_angles(node->current->unit_direction,
                                        &phi, &thet);
                thetaav += *(weights+ctr) *thet;
                phiav += *(weights+ctr) *phi;
        }
        rn->current->unit_direction[0] = cos(phiav);
        rn->current->unit_direction[1] = sin(phiav)*cos(thetaav);
        rn->current->unit_direction[2] = sin(phiav)*sin(thetaav);
        rn->hit_destination_patch = 0;
        rn->num_reflections = 0;
        rn->ctr = 1;

        cblas_dcopy(3, node->current->point, 1,
                    rn->current->point, 1);
}

void invert_spherical_angles(const double *unit_vector, double *phi,
                             double *thet)
{
        *thet = atan(unit_vector[2]/unit_vector[1]);
        *phi = acos(unit_vector[0]);
        if (sin(*phi) * sin(*thet) / unit_vector[2] < 0 )
                *phi = 2 * PI - (*phi);
}

void compute_averaging_coefficients(const double *point,
                                    const struct ray_ribbon_array *rba,
                                    double *weights)
{
        gsl_matrix *mat = gsl_matrix_alloc(3, 2);
        gsl_permutation *perm = gsl_permutation_alloc(3);
        gsl_vector *x = gsl_vector_alloc(2);
        gsl_vector *b = gsl_vector_alloc(3);
        gsl_vector *tau = gsl_vector_alloc(2);
        gsl_vector *residual = gsl_vector_alloc(3);

        int c0, c1;
        for (c0=0; c0<3; c0++) {
                const struct ribbon_node *node =
                        get_last_ribbon_node(*(rba->ribbons));
                gsl_vector_set(b, c0, *(point + c0)
                               - node->current->end_pt[c0]);
                for (c1=0; c1<2; c1++) {
                        const struct ribbon_node *nodeset =
                                get_last_ribbon_node(*(rba->ribbons + c1 + 1));
                        gsl_matrix_set(mat, c0, c1,
                                       nodeset->current->end_pt[c0]
                                       - node->current->end_pt[c0]);
                }
        }

        gsl_linalg_QR_decomp(mat, tau);
        gsl_linalg_QR_lssolve(mat, tau, b, x, residual);

        *(weights) = 1 - gsl_vector_get(x, 0) - gsl_vector_get(x, 1);
        *(weights + 1) = gsl_vector_get(x, 0);
        *(weights + 2) = gsl_vector_get(x, 1);

        gsl_matrix_free(mat);
        gsl_permutation_free(perm);
        gsl_vector_free(x);
        gsl_vector_free(b);
        gsl_vector_free(tau);
        gsl_vector_free(residual);
}

struct ray_ribbon_array *generate_nearby_ribbons(const struct transmitter *tx,
                                                 const struct
                                                 perfect_reflector **ref_arr,
                                                 int num_ref,
                                                 const struct ray_ribbon *rb)
{
        double phi, theta;
        invert_spherical_angles(rb->head->current->unit_direction,
                                &phi, &theta);
        struct ray_ribbon_array *rarr = init_ray_ribbon_array(4);

        complex double *angles = malloc(3 * sizeof(complex double));
        *angles = (phi - 0.0000059) + I * (theta + 0.0000057);
        *(angles + 1) = phi + 0.000004 + I * (theta + 0.000009);
        *(angles + 2) = phi + 0.0000039 + I * (theta + 0.0000011);

        populate_ray_ribbon_array_full(tx, ref_arr, num_ref, 3, angles, rarr,
                                       false);
        free(angles);
        return rarr;
}

struct ray_ribbon *refine_ray_ribbon(const struct transmitter *tx,
                                     const struct ray_ribbon *rb,
                                     const struct receiver *rx,
                                     const struct perfect_reflector **pr)
{
        struct ray_ribbon_array *node_array_mod = generate_nearby_ribbons(
                tx, pr, 3, rb);
        const double *point = rx->gn->smm->position;
        fprintf(stderr, "Position\n");
        print_vector(point);
        int ctr_typ_so_far = 0, ctr = 0;
        double weights[3];
        struct ray_ribbon *rbn = *(node_array_mod->ribbons);
        struct ribbon_node *rn = 0;

        while (!is_close_ribbon(rbn, point)) {
                compute_averaging_coefficients(point, node_array_mod, weights);
                rn = init_ribbon_node();
                /* fprintf(stderr, "Weight is:\n"); print_vector(weights); */
                compute_average_ribbon_node(rn, node_array_mod, weights);
                destroy_ray_ribbon_nodes(rbn);
                bool has_hit = process_vertical_chain(rn, pr, 3);
                if (!has_hit) {
                        fprintf(stderr, "Unexpected error.\n");
                        exit(1);
                }
                (*(node_array_mod->ribbons))->head = rn;
                rbn =  *(node_array_mod->ribbons);
                // print_ray_ribbon(*(node_array_mod->ribbons));
        }

        destroy_ray_ribbon_array_all_but_first(node_array_mod);
        rbn->start_gn = tx->gn;
        rbn->end_gn = rx->gn;
        return rbn;
}


bool is_close_ribbon(const struct ray_ribbon *rb, const double *point)
{
        return isclose(rb->head, point);
}

bool isclose(const struct ribbon_node *rn, const double *point)
{
        double diff[3];
        while (rn->down != NULL) {
                rn = rn->down;
        }
        cblas_dcopy(3, rn->current->end_pt, 1, diff, 1);
        cblas_daxpy(3, -1, point, 1, diff, 1);
        if (cblas_dnrm2(3, diff, 1) < 1e-5)  return true;
        return false;
}

long type_ray_ribbon(const struct ray_ribbon *rb)
{
        int ctr=0;
        struct ribbon_node *rn = rb->head;
        while (rn != NULL) {
                ctr = MAX_SURFACES * ctr + rn->surface_index + 1;
                rn = rn->down;
        }
        return ctr;
}

struct ray_ribbon_array *throw_three_dim_ray_ribbon(struct transmitter *tn,
                                                    const struct
                                                    perfect_reflector **p,
                                                    int num_ref,
                                                    const double phi_start,
                                                    const double phi_end,
                                                    const double phi_incr,
                                                    const double thet_start,
                                                    const double thet_end,
                                                    const double thet_incr)
{
        struct ray_ribbon_array *rarr = malloc(sizeof(struct ray_ribbon_array));
        populate_ray_ribbon_array_long(tn, p, rarr,
                                       num_ref, phi_start, phi_end,
                                       phi_incr, thet_start, thet_end,
                                       thet_incr, true);
        return rarr;
}

struct set_of_types *get_unique_types(const struct ray_ribbon_array *arr)
{
        struct set_of_types *st = init_set_of_types(NUM_TYPES);
        int ctr = 0;
        const struct ray_ribbon *rb = *(arr->ribbons + ctr);
        while (*(st->types + ctr) >= 0) {
                long typenum = type_ray_ribbon(rb);
                add_entry(typenum, st);
                ctr++;
        }
        return st;
}

void print_ray_ribbon_types(const struct ray_ribbon_array *arr)
{
        int ctr = 0;
        const struct ray_ribbon *rb = *(arr->ribbons + ctr);
        while (rb != NULL) {
                long typenum = type_ray_ribbon(rb);
                fprintf(stderr, "Type of ribbon %d is %ld\n ", ctr,
                        typenum);
                ctr++;
                rb = *(arr->ribbons + ctr);
        }

}

struct ribbon_node *get_last_ribbon_node(const struct ray_ribbon *rb)
{
        struct ribbon_node *rn = rb->head;
        while (rn->down != NULL) {
                rn = rn->down;
        }
        return rn;
}

void populate_env_paths(struct environment *env)
{
        add_receiver_patch(env);

        // first sound the channel
        const struct perfect_reflector **prconst =
                (const struct perfect_reflector **) env->prarray;
        // clear existing paths
        int ctr = 0;
        while (*(env->tx_paths + ctr) != 0) {
                destroy_ray_ribbon_array(*(env->tx_paths + ctr));
                *(env->tx_paths + ctr) = 0;
                ctr++;
        }

        // now populate individual paths
        ctr = 0;
        struct transmitter *tx = *(env->transmitters_array);
        while (tx != 0) {
                struct ray_ribbon_array *rb_arr = init_ray_ribbon_array(20);
                populate_ray_ribbon_array(tx, prconst, rb_arr, 300, 3, true);
                *(env->tx_paths + ctr) = rb_arr;
                ctr++;
                tx = *(env->transmitters_array + ctr);
        }

        // now generate rayribbons for each receiver
        // clear existing ray ribbons
        ctr = 0;
        while (*(env->env_paths + ctr) != 0) {
                destroy_ray_ribbon_array(*(env->env_paths + ctr));
                *(env->env_paths + ctr) = 0;
                ctr++;
        }

        ctr = 0;
        struct receiver *rx = *(env->receivers_array);
        while (rx != 0) {
                int ctr1 = 0;
                struct ray_ribbon_array *rba = *(env->tx_paths);
                struct ray_ribbon *rb;
                struct transmitter *tx = *(env->transmitters_array);
                while (rba != 0) {
                        struct ray_ribbon_array *rb_arr =
                                init_ray_ribbon_array(5); // change this number
                        int ctr2 = 0;
                        while((rb = *(rba->ribbons + ctr2)) != 0) {
                                add_ray_ribbon(
                                        rb_arr,
                                        refine_ray_ribbon(tx, rb, rx, prconst),
                                        true);
                                ctr2++;
                                rb = *(rba->ribbons + ctr2);

                        }
                        *(env->env_paths + ctr1) = rb_arr;
                        ctr1++;
                        rba = *(env->tx_paths + ctr1);
                        tx = *(env->transmitters_array + ctr1);
                }
                ctr++;
                rx = *(env->receivers_array + ctr);
        }

        destroy_last_reflector(env);
}

void update_env_paths_delay_dopplers(struct environment *env) {
        int ctr = 0;
        struct ray_ribbon_array *rba = *(env->env_paths + ctr);
        while (rba != 0) {
                int ctr1 = 0;
                struct ray_ribbon *rb = *(rba->ribbons + ctr1);
                while (rb != NULL) {
                        update_ribbon_delay_dopplers(rb, env);
                        ctr1++;
                        rb = *(rba->ribbons + ctr1);
                }
                ctr++;
                rba = *(env->env_paths + ctr);
        }
}

void update_ribbon_delay_dopplers(struct ray_ribbon *rb,
                                  const struct environment *env) {
        if (rb == 0) return;
        // compute delay
        double dist = 0;
        struct ribbon_node *rn = rb->head;
        while (rn != NULL) {
                dist += length_ribbon_node(rn);
                rn = rn->down;
        }
        rb->delay = dist/C;
        // compute doppler
        rb->doppler = compute_doppler(rb, env);
}

double compute_doppler(const struct ray_ribbon *rb,
                       const struct environment *env) {
        struct perfect_reflector *pr;
        struct ribbon_node *rn = rb->head->down;
        double src_pos[3];
        double src_vel[3];

        cblas_dcopy(3, rb->start_gn->smm->position, 1, src_pos, 1);
        cblas_dcopy(3, rb->start_gn->smm->velocity, 1, src_vel, 1);

        while (rn != NULL) {
                pr = *(env->prarray + rn->surface_index);
                reflect(pr->center_point, pr->unit_normal, src_vel, src_pos);
                rn = rn->down;
        }
        double rel_pos[3];
        double rel_vel[3];
        diff(rb->start_gn->smm->position, rb->end_gn->smm->position, rel_pos);
        diff(rb->start_gn->smm->velocity, rb->end_gn->smm->velocity, rel_vel);
        double tmp = cblas_dnrm2(3, rel_pos, 1);
        return -cblas_ddot(3, rel_vel, 1, rel_pos, 1) / tmp / env->wavelength;
}

double length_ribbon_node(const struct ribbon_node *rn) {
        double diff_vector[3];
        diff(rn->current->point, rn->current->end_pt, diff_vector);
        return cblas_dnrm2(3, diff_vector, 1);
}

void reflect(const double *pos1, const double *n1, double *vel, double *pos) {
        double diff_vec[3];
        double r1[3];
        double v1[3];

        diff(pos1, pos, diff_vec);
        reflection_operation(diff_vec, n1, r1);
        reflection_operation(vel, n1, v1);

        // now update vel and pos
        cblas_dcopy(3, pos1, 1, pos, 1);
        cblas_daxpy(3, 1, r1, 1, pos, 1);
        cblas_dcopy(3, v1, 1, vel, 1);
}

void reflection_operation(const double *v1, const double *n1, double *vref) {
        cblas_dcopy(3, v1, 1, vref, 1);
        double tmp = cblas_ddot(3, n1, 1, vref, 1);
        cblas_daxpy(3, -2 * tmp, n1, 1, vref, 1);
}

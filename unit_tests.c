#include "raytracing.h"
#include "CUnit/CUnit.h"

int maxi(int i1, int i2)
{
        return (i1 > i2) ? i1 : i2;
}

void test_maxi(void)
{
        CU_ASSERT(maxi(0,2) == 2);
        CU_ASSERT(maxi(0,-2) == 0);
}

void test_diff(void)
{
        double v1[3] = {0, -2, 3};
        double v2[3] = {0, 2, 4};
        double v3[3];
        double v4[3] = {0, 4, 1};

        diff(v1, v2, v3);

        CU_ASSERT_DOUBLE_EQUAL(v4[0], v3[0], 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(v4[1], v3[1], 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(v4[2], v3[2], 1e-6);
}

void test_reflect(void)
{
        double pos0[3] = {0, 0, 0};
        double n0[3] = {0, 0, 1};
        double pos1[3] = {1, 0, 0};
        double n1[3] = {0, 1, 0};
        double vel[3] = {-2, -5, 9};
        double pos[3] = {2, 3, 1};
        double ref_vel[3] = {-2, -5, -9};
        double ref_pos[3] = {2, 3, -1};
	double ref_vel_[3] = {-2, 5, -9};
	double ref_pos_[3] = {2, -3, -1};

        reflect(pos0, n0, vel, pos);

        int ctr = 0;

        for (ctr = 0; ctr < 3; ctr++) {
                CU_ASSERT_DOUBLE_EQUAL(vel[ctr], ref_vel[ctr], 1e-6);
                CU_ASSERT_DOUBLE_EQUAL(pos[ctr], ref_pos[ctr], 1e-6);
        }

	reflect(pos1, n1, vel, pos);

        for (ctr = 0; ctr < 3; ctr++) {
                CU_ASSERT_DOUBLE_EQUAL(vel[ctr], ref_vel_[ctr], 1e-6);
                CU_ASSERT_DOUBLE_EQUAL(pos[ctr], ref_pos_[ctr], 1e-6);
        }
}

void test_malloc(void)
{
        double *num_bytes0 = get_or_free_memory(sizeof(double), 0, 0);
        get_or_free_memory(0, num_bytes0, 1);
	double *num_bytes1 = get_or_free_memory(sizeof(double), 0, 0);
	release_all_blocks();
        CU_ASSERT_TRUE(num_bytes0 == num_bytes1);
}

int main()
{
        CU_initialize_registry();
        struct CU_Suite *suite = CU_add_suite("uniquetest", NULL, NULL);
        CU_add_test(suite, "Test diff", test_diff);
        CU_add_test(suite, "Test reflect", test_reflect);
        CU_add_test(suite, "Test malloc", test_malloc);
        CU_basic_run_tests();
        CU_cleanup_registry();
}

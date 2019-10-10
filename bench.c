#include "quotient-filter.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

const uint32_t Q_MAX = 12;
const uint32_t R_MAX = 6;
const uint32_t ROUNDS_MAX = 1000;

static void qf_bench()
{
    quotient_filter qf;
    const uint32_t q_large = 28;
    const uint32_t q_small = 16;
    const uint32_t nlookups = 1000000;
    struct timeval tv1, tv2;
    uint64_t sec;

    /* Test random inserts + lookups */
    uint32_t ninserts = (3 * (1 << q_large) / 4);
    printf("Testing %u random inserts and %u lookups", ninserts, nlookups);
    fflush(stdout);
    qf_init(&qf, q_large, 1);
    gettimeofday(&tv1, NULL);
    while (qf.entries < ninserts) {
        assert(qf_insert(&qf, (uint64_t) rand()));
        if (qf.entries % 10000000 == 0) {
            printf(".");
            fflush(stdout);
        }
    }
    for (uint32_t i = 0; i < nlookups; ++i)
        qf_may_contain(&qf, (uint64_t) rand());
    gettimeofday(&tv2, NULL);
    sec = tv2.tv_sec - tv1.tv_sec;
    printf(" done (%lu seconds).\n", sec);
    fflush(stdout);
    qf_destroy(&qf);

    /* Create a large cluster. Test random lookups. */
    qf_init(&qf, q_small, 1);
    printf("Testing %u contiguous inserts and %u lookups", 1 << q_small,
           nlookups);
    fflush(stdout);
    gettimeofday(&tv1, NULL);
    for (uint64_t quot = 0; quot < (1 << (q_small - 1)); ++quot) {
        uint64_t hash = quot << 1;
        assert(qf_insert(&qf, hash));
        assert(qf_insert(&qf, hash | 1));
        if (quot % 2000 == 0) {
            printf(".");
            fflush(stdout);
        }
    }
    for (uint32_t i = 0; i < nlookups; ++i) {
        qf_may_contain(&qf, (uint64_t) rand());
        if (i % 50000 == 0) {
            printf(".");
            fflush(stdout);
        }
    }
    gettimeofday(&tv2, NULL);
    sec = tv2.tv_sec - tv1.tv_sec;
    printf(" done (%lu seconds).\n", sec);
    fflush(stdout);
    qf_destroy(&qf);
}

int main()
{
    srand(0);
    qf_bench();

    return 0;
}

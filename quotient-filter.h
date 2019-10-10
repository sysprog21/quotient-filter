#ifndef QUOTIENT_FILTER_H
#define QUOTIENT_FILTER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
    uint8_t qbits, rbits;
    uint8_t elem_bits;
    uint32_t entries;
    uint64_t index_mask, rmask, elem_mask;
    uint64_t max_size;
    uint64_t *table;
} quotient_filter;
typedef struct __qf_iterator qf_iterator;

/**
 * Initializes a quotient filter with capacity 2^q.
 * Increasing r improves the filter's accuracy but uses more space.
 *
 * Returns false if q == 0, r == 0, q+r > 64, or on ENOMEM.
 */
bool qf_init(quotient_filter *qf, uint32_t q, uint32_t r);

/**
 * Inserts a hash into the QF.
 * Only the lowest q+r bits are actually inserted into the QF table.
 *
 * Returns false if the QF is full.
 */
bool qf_insert(quotient_filter *qf, uint64_t hash);

/**
 * Returns true if the QF may contain the hash. Returns false otherwise.
 */
bool qf_may_contain(quotient_filter *qf, uint64_t hash);

/**
 * Removes a hash from the QF.
 *
 * Caution: If you plan on using this function, make sure that your hash
 * function emits no more than q+r bits. Consider the following scenario;
 *   insert(qf, A:X)   # X is in the lowest q+r bits.
 *   insert(qf, B:X)   # This is a no-op, since X is already in the table.
 *   remove(qf, A:X)   # X is removed from the table.
 *
 * Now, may-contain(qf, B:X) == false, which is a ruinous false negative.
 *
 * Returns false if the hash uses more than q+r bits.
 */
bool qf_remove(quotient_filter *qf, uint64_t hash);

/**
 * Resets the QF table. This function does not deallocate any memory.
 */
void qf_clear(quotient_filter *qf);

/**
 * Finds the size (in bytes) of a QF table.
 * Caution: sizeof(quotient_filter) is not included.
 */
size_t qf_table_size(uint32_t q, uint32_t r);

/**
 * Deallocates the QF table.
 */
void qf_destroy(quotient_filter *qf);

/**
 * Initialize an iterator for the QF.
 */
void qfi_start(quotient_filter *qf, qf_iterator *i);

/**
 * Returns true if there are no elements left to visit.
 */
bool qfi_done(quotient_filter *qf, qf_iterator *i);

/**
 * Returns the next (q+r)-bit fingerprint in the QF.
 * Limitation: Can not call this routine if qfi_done() == true.
 */
uint64_t qfi_next(quotient_filter *qf, qf_iterator *i);

#endif

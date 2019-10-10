#include <stdlib.h>
#include <string.h>

#include "quotient-filter.h"

#define LOW_MASK(n) ((1ULL << (n)) - 1ULL)

struct __qf_iterator {
    uint64_t index;
    uint64_t quotient;
    uint64_t visited;
};

bool qf_init(quotient_filter *qf, uint32_t q, uint32_t r)
{
    if (q == 0 || r == 0 || q + r > 64)
        return false;

    qf->qbits = q;
    qf->rbits = r;
    qf->elem_bits = qf->rbits + 3;
    qf->index_mask = LOW_MASK(q);
    qf->rmask = LOW_MASK(r);
    qf->elem_mask = LOW_MASK(qf->elem_bits);
    qf->entries = 0;
    qf->max_size = 1 << q;
    qf->table = (uint64_t *) calloc(qf_table_size(q, r), 1);
    return qf->table != NULL;
}

/* Return QF[idx] in the lower bits. */
static uint64_t get_elem(quotient_filter *qf, uint64_t idx)
{
    size_t bitpos = qf->elem_bits * idx;
    size_t tabpos = bitpos / 64;
    size_t slotpos = bitpos % 64;
    int spillbits = (slotpos + qf->elem_bits) - 64;
    uint64_t elt = (qf->table[tabpos] >> slotpos) & qf->elem_mask;
    if (spillbits > 0) {
        ++tabpos;
        uint64_t x = qf->table[tabpos] & LOW_MASK(spillbits);
        elt |= x << (qf->elem_bits - spillbits);
    }
    return elt;
}

/* Store the lower bits of elt into QF[idx]. */
static void set_elem(quotient_filter *qf, uint64_t idx, uint64_t elt)
{
    size_t bitpos = qf->elem_bits * idx;
    size_t tabpos = bitpos / 64;
    size_t slotpos = bitpos % 64;
    int spillbits = (slotpos + qf->elem_bits) - 64;
    elt &= qf->elem_mask;
    qf->table[tabpos] &= ~(qf->elem_mask << slotpos);
    qf->table[tabpos] |= elt << slotpos;
    if (spillbits > 0) {
        ++tabpos;
        qf->table[tabpos] &= ~LOW_MASK(spillbits);
        qf->table[tabpos] |= elt >> (qf->elem_bits - spillbits);
    }
}

static inline uint64_t incr(quotient_filter *qf, uint64_t idx)
{
    return (idx + 1) & qf->index_mask;
}

static inline uint64_t decr(quotient_filter *qf, uint64_t idx)
{
    return (idx - 1) & qf->index_mask;
}

/* quickly let the search algorithm determine whether a fingerprint
 * s.t Fq = x exists in the QF, where x is an index in the QF table.
 * I.e, the `is_occupied` bit is attached to a slot in the table and
 * not to any ephemeral fingerprint.
 */
static inline int is_occupied(uint64_t elt)
{
    return elt & 1;
}

static inline uint64_t set_occupied(uint64_t elt)
{
    return elt | 1;
}

static inline uint64_t clr_occupied(uint64_t elt)
{
    return elt & ~1;
}

/* Set to false for the first fingerprint in a run
 * Set to true for subsequent fingerprints in a run
 */
static inline int is_continuation(uint64_t elt)
{
    return elt & 2;
}

static inline uint64_t set_continuation(uint64_t elt)
{
    return elt | 2;
}

static inline uint64_t clr_continuation(uint64_t elt)
{
    return elt & ~2;
}

/* Set to true when a fingerprint is not in its canonical slot (i.e when
 * Fq != x, where x is the slot index). Set to false otherwise. Note that
 * this metadata bit is attached to the fingerprint, and not to a particular
 * slot.
 */
static inline int is_shifted(uint64_t elt)
{
    return elt & 4;
}

static inline uint64_t set_shifted(uint64_t elt)
{
    return elt | 4;
}

static inline uint64_t clr_shifted(uint64_t elt)
{
    return elt & ~4;
}

static inline uint64_t get_remainder(uint64_t elt)
{
    return elt >> 3;
}

static inline bool is_empty_element(uint64_t elt)
{
    return (elt & 7) == 0;
}

static inline bool is_cluster_start(uint64_t elt)
{
    return is_occupied(elt) && !is_continuation(elt) && !is_shifted(elt);
}

static inline bool is_run_start(uint64_t elt)
{
    return !is_continuation(elt) && (is_occupied(elt) || is_shifted(elt));
}

static inline uint64_t hash_to_quotient(quotient_filter *qf, uint64_t hash)
{
    return (hash >> qf->rbits) & qf->index_mask;
}

static inline uint64_t hash_to_remainder(quotient_filter *qf, uint64_t hash)
{
    return hash & qf->rmask;
}

/* Find the start index of the run for fq (given that the run exists). */
static uint64_t find_run_index(quotient_filter *qf, uint64_t fq)
{
    /* Find the start of the cluster. */
    uint64_t b = fq;
    while (is_shifted(get_elem(qf, b)))
        b = decr(qf, b);

    /* Find the start of the run for fq. */
    uint64_t s = b;
    while (b != fq) {
        do {
            s = incr(qf, s);
        } while (is_continuation(get_elem(qf, s)));

        do {
            b = incr(qf, b);
        } while (!is_occupied(get_elem(qf, b)));
    }
    return s;
}

/* Insert elt into QF[s], shifting over elements as necessary. */
static void insert_into(quotient_filter *qf, uint64_t s, uint64_t elt)
{
    uint64_t curr = elt;
    bool empty;

    do {
        uint64_t prev = get_elem(qf, s);
        empty = is_empty_element(prev);
        if (!empty) {
            /* Fix up `is_shifted' and `is_occupied'. */
            prev = set_shifted(prev);
            if (is_occupied(prev)) {
                curr = set_occupied(curr);
                prev = clr_occupied(prev);
            }
        }
        set_elem(qf, s, curr);
        curr = prev;
        s = incr(qf, s);
    } while (!empty);
}

bool qf_insert(quotient_filter *qf, uint64_t hash)
{
    if (qf->entries >= qf->max_size)
        return false;

    uint64_t fq = hash_to_quotient(qf, hash);
    uint64_t fr = hash_to_remainder(qf, hash);
    uint64_t T_fq = get_elem(qf, fq);
    uint64_t entry = (fr << 3) & ~7;

    /* Special-case filling canonical slots to simplify insert_into(). */
    if (is_empty_element(T_fq)) {
        set_elem(qf, fq, set_occupied(entry));
        ++qf->entries;
        return true;
    }

    if (!is_occupied(T_fq))
        set_elem(qf, fq, set_occupied(T_fq));

    uint64_t start = find_run_index(qf, fq);
    uint64_t s = start;

    if (is_occupied(T_fq)) {
        /* Move the cursor to the insert position in the fq run. */
        do {
            uint64_t rem = get_remainder(get_elem(qf, s));
            if (rem == fr) {
                return true;
            } else if (rem > fr) {
                break;
            }
            s = incr(qf, s);
        } while (is_continuation(get_elem(qf, s)));

        if (s == start) {
            /* The old start-of-run becomes a continuation. */
            uint64_t old_head = get_elem(qf, start);
            set_elem(qf, start, set_continuation(old_head));
        } else {
            /* The new element becomes a continuation. */
            entry = set_continuation(entry);
        }
    }

    /* Set the shifted bit if we can't use the canonical slot. */
    if (s != fq)
        entry = set_shifted(entry);

    insert_into(qf, s, entry);
    ++qf->entries;
    return true;
}

bool qf_may_contain(quotient_filter *qf, uint64_t hash)
{
    uint64_t fq = hash_to_quotient(qf, hash);
    uint64_t fr = hash_to_remainder(qf, hash);
    uint64_t T_fq = get_elem(qf, fq);

    /* If this quotient has no run, give up. */
    if (!is_occupied(T_fq))
        return false;

    /* Scan the sorted run for the target remainder. */
    uint64_t s = find_run_index(qf, fq);
    do {
        uint64_t rem = get_remainder(get_elem(qf, s));
        if (rem == fr) {
            return true;
        } else if (rem > fr) {
            return false;
        }
        s = incr(qf, s);
    } while (is_continuation(get_elem(qf, s)));
    return false;
}

/* Remove the entry in QF[s] and slide the rest of the cluster forward. */
static void delete_entry(quotient_filter *qf, uint64_t s, uint64_t quot)
{
    uint64_t curr = get_elem(qf, s);
    uint64_t sp = incr(qf, s);
    uint64_t orig = s;

    while (1) {
        uint64_t next = get_elem(qf, sp);
        bool curr_occupied = is_occupied(curr);

        if (is_empty_element(next) || is_cluster_start(next) || sp == orig) {
            set_elem(qf, s, 0);
            return;
        } else {
            /* Fix entries which slide into canonical slots. */
            uint64_t updated_next = next;
            if (is_run_start(next)) {
                do {
                    quot = incr(qf, quot);
                } while (!is_occupied(get_elem(qf, quot)));

                if (curr_occupied && quot == s) {
                    updated_next = clr_shifted(next);
                }
            }

            set_elem(qf, s,
                     curr_occupied ? set_occupied(updated_next)
                                   : clr_occupied(updated_next));
            s = sp;
            sp = incr(qf, sp);
            curr = next;
        }
    }
}

bool qf_remove(quotient_filter *qf, uint64_t hash)
{
    uint64_t highbits = hash >> (qf->qbits + qf->rbits);
    if (highbits)
        return false;

    uint64_t fq = hash_to_quotient(qf, hash);
    uint64_t fr = hash_to_remainder(qf, hash);
    uint64_t T_fq = get_elem(qf, fq);

    if (!is_occupied(T_fq) || !qf->entries)
        return true;

    uint64_t start = find_run_index(qf, fq);
    uint64_t s = start;
    uint64_t rem;

    /* Find the offending table index */
    do {
        rem = get_remainder(get_elem(qf, s));
        if (rem == fr) {
            break;
        } else if (rem > fr) {
            return true;
        }
        s = incr(qf, s);
    } while (is_continuation(get_elem(qf, s)));
    if (rem != fr) {
        return true;
    }

    uint64_t kill = (s == fq) ? T_fq : get_elem(qf, s);
    bool replace_run_start = is_run_start(kill);

    /* If we are deleting the last entry in a run, clear `is_occupied'. */
    if (is_run_start(kill)) {
        /* Write your code here */	    
    }

    delete_entry(qf, s, fq);

    if (replace_run_start) {
        /* Write your code here */	    
    }

    --qf->entries;
    return true;
}

void qf_clear(quotient_filter *qf)
{
    qf->entries = 0;
    memset(qf->table, 0, qf_table_size(qf->qbits, qf->rbits));
}

size_t qf_table_size(uint32_t q, uint32_t r)
{
    size_t bits = (1 << q) * (r + 3);
    size_t bytes = bits / 8;
    return (bits % 8) ? (bytes + 1) : bytes;
}

void qf_destroy(quotient_filter *qf)
{
    free(qf->table);
}

void qfi_start(quotient_filter *qf, qf_iterator *i)
{
    /* Mark the iterator as done. */
    i->visited = qf->entries;

    if (qf->entries == 0)
        return;

    /* Find the start of a cluster. */
    uint64_t start;
    for (start = 0; start < qf->max_size; ++start) {
        if (is_cluster_start(get_elem(qf, start)))
            break;
    }

    i->visited = 0;
    i->index = start;
}

bool qfi_done(quotient_filter *qf, qf_iterator *i)
{
    return qf->entries == i->visited;
}

uint64_t qfi_next(quotient_filter *qf, qf_iterator *i)
{
    while (!qfi_done(qf, i)) {
        uint64_t elt = get_elem(qf, i->index);

        /* Keep track of the current run. */
        if (is_cluster_start(elt)) {
            i->quotient = i->index;
        } else {
            if (is_run_start(elt)) {
                uint64_t quot = i->quotient;
                do {
                    quot = incr(qf, quot);
                } while (!is_occupied(get_elem(qf, quot)));
                i->quotient = quot;
            }
        }

        i->index = incr(qf, i->index);

        if (!is_empty_element(elt)) {
            uint64_t quot = i->quotient;
            uint64_t rem = get_remainder(elt);
            uint64_t hash = (quot << qf->rbits) | rem;
            ++i->visited;
            return hash;
        }
    }

    /* shall not reach here */
    abort();
}



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "Baseline.h"

#define MAX_ELEMENT_VALUE 1024
#define MAX_SETS 128
#define PARENT_IS_ROOT 0xFFFFFFFF


typedef struct {
    uint32_t size;
    Bitvector compressed;
} Universe;

typedef struct {
    uint32_t size;
    uint32_t parent_index;
    Bitvector compressed;
} SetNode;

typedef struct {
    uint32_t set_indices[MAX_SETS];
    uint32_t count;
} InvertedListEntry;

typedef struct {
    InvertedListEntry entries[MAX_ELEMENT_VALUE];
} InvertedList;

typedef struct {
    Universe universe;
    SetNode nodes[MAX_SETS];
    uint32_t num_sets;
} SetCollection;


// Scratch buffers used only during construction. Freed by scratch_free.
static uint8_t (*scratch_bits)[MAX_ELEMENT_VALUE] = NULL;
static uint8_t *scratch_universe = NULL;

// Allocates the temporary construction buffers.
static void scratch_alloc(void)
{
    scratch_bits = calloc(MAX_SETS, sizeof(uint8_t[MAX_ELEMENT_VALUE]));
    scratch_universe = calloc(MAX_ELEMENT_VALUE, sizeof(uint8_t));
}

// Frees the construction buffers so only the compressed form remains.
static void scratch_free(void)
{
    free(scratch_bits);
    scratch_bits = NULL;
    free(scratch_universe);
    scratch_universe = NULL;
}


// Builds the universe as the OR of all input sets and counts its size.
void build_universe(uint8_t input_sets[][MAX_ELEMENT_VALUE], uint32_t num_sets, SetCollection *col)
{
    col->universe.size = 0;

    for (uint32_t i = 0; i < MAX_ELEMENT_VALUE; i++)
    {
        scratch_universe[i] = 0;
    }

    for (uint32_t i = 0; i < num_sets; i++)
    {
        for (uint32_t j = 0; j < MAX_ELEMENT_VALUE; j++)
        {
            scratch_universe[j] |= input_sets[i][j];
        }
    }

    for (uint32_t i = 0; i < MAX_ELEMENT_VALUE; i++)
    {
        if (scratch_universe[i] == 1) col->universe.size++;
    }
}

// Copies each input set into scratch and records its size.
void load_sets(uint8_t input_sets[][MAX_ELEMENT_VALUE], uint32_t num_sets, SetCollection *col)
{
    col->num_sets = num_sets;

    for (uint32_t i = 0; i < num_sets; i++)
    {
        col->nodes[i].size = 0;

        for (uint32_t j = 0; j < MAX_ELEMENT_VALUE; j++)
        {
            scratch_bits[i][j] = input_sets[i][j];
            if (input_sets[i][j] == 1) col->nodes[i].size++;
        }

        col->nodes[i].parent_index = PARENT_IS_ROOT;
    }
}


// Sorts sets largest-first, keeping their scratch rows aligned.
void sort_sets_by_size_desc(SetCollection *col)
{
    for (uint32_t i = 1; i < col->num_sets; i++)
    {
        SetNode tmp_node = col->nodes[i];
        uint8_t tmp_bits[MAX_ELEMENT_VALUE];
        memcpy(tmp_bits, scratch_bits[i], MAX_ELEMENT_VALUE);

        int32_t j = (int32_t)i - 1;
        while (j >= 0 && col->nodes[j].size < tmp_node.size)
        {
            col->nodes[j + 1] = col->nodes[j];
            memcpy(scratch_bits[j + 1], scratch_bits[j], MAX_ELEMENT_VALUE);
            j--;
        }
        col->nodes[j + 1] = tmp_node;
        memcpy(scratch_bits[j + 1], tmp_bits, MAX_ELEMENT_VALUE);
    }
}


// For each element, lists which sets contain it.
void build_inverted_list(SetCollection *col, InvertedList *il)
{
    for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
    {
        il->entries[e].count = 0;
    }

    for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
    {
        for (uint32_t i = 0; i < col->num_sets; i++)
        {
            if (scratch_bits[i][e] == 1)
            {
                uint32_t c = il->entries[e].count;
                il->entries[e].set_indices[c] = i;
                il->entries[e].count++;
            }
        }
    }
}


// Finds each set's smallest superset and stores it as the parent.
void build_hierarchy(SetCollection *col, InvertedList *il)
{
    for (uint32_t i = 0; i < col->num_sets; i++)
    {
        uint8_t is_superset[MAX_SETS];
        for (uint32_t j = 0; j < col->num_sets; j++)
        {
            is_superset[j] = 1;
        }

        for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
        {
            if (scratch_bits[i][e] == 0) continue;

            uint8_t contains_e[MAX_SETS];
            for (uint32_t j = 0; j < col->num_sets; j++)
            {
                contains_e[j] = 0;
            }
            for (uint32_t k = 0; k < il->entries[e].count; k++)
            {
                uint32_t j = il->entries[e].set_indices[k];
                contains_e[j] = 1;
            }

            for (uint32_t j = 0; j < col->num_sets; j++)
            {
                if (contains_e[j] == 0)
                {
                    is_superset[j] = 0;
                }
            }
        }

        is_superset[i] = 0;

        uint32_t best_size = UINT32_MAX;
        uint32_t best_index = PARENT_IS_ROOT;

        for (uint32_t j = 0; j < col->num_sets; j++)
        {
            if (is_superset[j] == 1 && col->nodes[j].size < best_size)
            {
                best_size = col->nodes[j].size;
                best_index = j;
            }
        }

        col->nodes[i].parent_index = best_index;
    }
}

// Replaces each parent with the highest ancestor up to twice the set's size.
void contract_parents(SetCollection *col)
{
    for (uint32_t i = 0; i < col->num_sets; i++)
    {
        uint32_t limit = 2 * col->nodes[i].size;

        uint32_t curr = col->nodes[i].parent_index;
        if (curr == PARENT_IS_ROOT) continue;

        while (1)
        {
            uint32_t up = col->nodes[curr].parent_index;
            if (up == PARENT_IS_ROOT) break;
            if (col->nodes[up].size > limit) break;
            curr = up;
        }

        col->nodes[i].parent_index = curr;
    }
}

// Rewrites each set's bits relative to its parent's elements.
void build_relative_bitvectors(SetCollection *col)
{
    for (int32_t i = (int32_t)col->num_sets - 1; i >= 0; i--)
    {
        uint8_t *parent_bits;
        if (col->nodes[i].parent_index == PARENT_IS_ROOT)
        {
            parent_bits = scratch_universe;
        }
        else
        {
            parent_bits = scratch_bits[col->nodes[i].parent_index];
        }

        uint8_t new_bitvector[MAX_ELEMENT_VALUE] = {0};
        uint32_t k = 0;

        for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
        {
            if (parent_bits[e] == 1)
            {
                new_bitvector[k] = scratch_bits[i][e];
                k++;
            }
        }

        for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
        {
            scratch_bits[i][e] = new_bitvector[e];
        }
    }
}


// Packs a uint8 bit array into padded uint64 words for the baseline encoder.
static void pack_bits(uint8_t *src, uint32_t len, uint64_t *W, uint32_t *padded_n)
{
    uint32_t pad = BLOCK_SIZE * SUPERBLOCK_FREQ;
    uint32_t n = ((len + pad - 1) / pad) * pad;
    if (n == 0) n = pad;
    *padded_n = n;

    uint32_t words = (n + 63) / 64;
    for (uint32_t i = 0; i < words; i++) W[i] = 0;

    for (uint32_t i = 0; i < len; i++)
    {
        if (src[i] == 1) bit_set(W, i + 1);
    }
}

// Compresses the universe and every set into baseline Bitvectors.
void compress_all(SetCollection *col)
{
    uint64_t W[MAX_WORDS];
    uint32_t padded_n;

    pack_bits(scratch_universe, MAX_ELEMENT_VALUE, W, &padded_n);
    col->universe.compressed = bitvector_build(W, padded_n, BLOCK_SIZE, SUPERBLOCK_FREQ);

    for (uint32_t i = 0; i < col->num_sets; i++)
    {
        uint32_t parent_size;
        if (col->nodes[i].parent_index == PARENT_IS_ROOT)
        {
            parent_size = col->universe.size;
        }
        else
        {
            parent_size = col->nodes[col->nodes[i].parent_index].size;
        }

        pack_bits(scratch_bits[i], parent_size, W, &padded_n);
        col->nodes[i].compressed = bitvector_build(W, padded_n, BLOCK_SIZE, SUPERBLOCK_FREQ);
    }
}


// Returns the k-th element of a set by walking up the hierarchy.
uint32_t retrieve(SetCollection *col, uint32_t set_idx, uint32_t k)
{
    SetNode *node = &col->nodes[set_idx];

    if (node->parent_index == PARENT_IS_ROOT)
    {
        uint32_t p = (uint32_t)bitvector_select(node->compressed, (int)k);
        return (uint32_t)bitvector_select(col->universe.compressed, (int)p) - 1;
    }

    uint32_t p = (uint32_t)bitvector_select(node->compressed, (int)k);
    return retrieve(col, node->parent_index, p);
}

// Returns how many elements of a set are <= k by walking up the hierarchy.
uint32_t rank_set(SetCollection *col, uint32_t set_idx, uint32_t k)
{
    SetNode *node = &col->nodes[set_idx];

    uint32_t p;
    if (node->parent_index == PARENT_IS_ROOT)
    {
        p = (uint32_t)bitvector_rank(col->universe.compressed, (int)(k + 1));
    }
    else
    {
        p = rank_set(col, node->parent_index, k);
    }

    if (p == 0) return 0;
    return (uint32_t)bitvector_rank(node->compressed, (int)p);
}


// Turns a list of element values into a 0/1 bit row.
void add_set(uint8_t *row, int *elements, int n)
{
    for (int j = 0; j < MAX_ELEMENT_VALUE; j++) row[j] = 0;
    for (int i = 0; i < n; i++) row[elements[i]] = 1;
}


static SetCollection col;
static InvertedList il;
static uint8_t input[MAX_SETS][MAX_ELEMENT_VALUE];

// Value lists for the generated collection and a sorted ground-truth snapshot.
static int gen_sets[MAX_SETS][MAX_ELEMENT_VALUE];
static int gen_sizes[MAX_SETS];
static int gt_elems[MAX_SETS][MAX_ELEMENT_VALUE];
static int gt_count[MAX_SETS];
static int g_pool[MAX_ELEMENT_VALUE];


// Deterministic xorshift generator so the synthetic data is reproducible.
static uint64_t g_rng;
static void rng_seed(uint64_t s) { g_rng = s ? s : 1; }
static uint32_t rng_u32(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 7;
    g_rng ^= g_rng << 17;
    return (uint32_t)g_rng;
}
static uint32_t rng_below(uint32_t m) { return m ? rng_u32() % m : 0; }

// Picks k distinct values from src[0..src_len-1] into dst (partial shuffle).
static void sample_subset(int *dst, const int *src, int src_len, int k)
{
    for (int i = 0; i < src_len; i++) g_pool[i] = src[i];
    for (int i = 0; i < k; i++)
    {
        int r = i + (int)rng_below((uint32_t)(src_len - i));
        int tmp = g_pool[i]; g_pool[i] = g_pool[r]; g_pool[r] = tmp;
        dst[i] = g_pool[i];
    }
}

// Generates a nested collection: set 0 is a random subset of [1..U], and every
// later set is a random subset of an earlier (larger) set. This guarantees real
// containment structure, which is exactly what the hierarchy exploits.
static uint32_t generate_synthetic_collection(uint32_t num_sets,
                                              int base_size, int universe_range)
{
    static int universe_vals[MAX_ELEMENT_VALUE];
    for (int i = 0; i < universe_range; i++) universe_vals[i] = i + 1;

    sample_subset(gen_sets[0], universe_vals, universe_range, base_size);
    gen_sizes[0] = base_size;

    for (uint32_t i = 1; i < num_sets; i++)
    {
        uint32_t parent;
        do { parent = rng_below(i); } while (gen_sizes[parent] <= 1);

        int child_size = gen_sizes[parent] / 2;
        if (child_size < 1) child_size = 1;

        sample_subset(gen_sets[i], gen_sets[parent], gen_sizes[parent], child_size);
        gen_sizes[i] = child_size;
    }

    for (uint32_t i = 0; i < num_sets; i++)
        add_set(input[i], gen_sets[i], gen_sizes[i]);

    return num_sets;
}


// Number of bits used to store one class value in [0, b].
static int class_bits(int b)
{
    int bits = 0;
    for (int t = b; t > 0; t >>= 1) bits++;
    return bits ? bits : 1;
}

// Stored payload of one compressed bitvector, in bits: class array + offsets.
static uint64_t payload_bits(Bitvector bv, int cb)
{
    return bv.num_blocks * (uint64_t)cb + (bv.o_pos - 1);
}

// Snapshots each node's absolute element list (must run while scratch is still
// in universe coordinates, i.e. before build_relative_bitvectors).
static void snapshot_ground_truth(void)
{
    for (uint32_t i = 0; i < col.num_sets; i++)
    {
        gt_count[i] = 0;
        for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
            if (scratch_bits[i][e] == 1)
                gt_elems[i][gt_count[i]++] = (int)e;
    }
}

// Depth of a node, counted as the number of edges up to the root.
static uint32_t node_depth(uint32_t i)
{
    uint32_t d = 0;
    while (col.nodes[i].parent_index != PARENT_IS_ROOT)
    {
        i = col.nodes[i].parent_index;
        d++;
    }
    return d;
}

static double ms(struct timespec a, struct timespec b)
{ return (b.tv_sec - a.tv_sec) * 1000.0 + (b.tv_nsec - a.tv_nsec) / 1e6; }


// Generates one nested collection, builds the containment hierarchy and the
// standalone baseline on the same data, then reports space (bits per element),
// hierarchy depth, build time and query latency, verifying every query answer.
static int run_synthetic(uint32_t num_sets, int base_size, int universe_range,
                         int num_queries)
{
    int cb = class_bits(BLOCK_SIZE);
    struct timespec t0, t1;

    rng_seed(42);
    num_sets = generate_synthetic_collection(num_sets, base_size, universe_range);

    build_universe(input, num_sets, &col);

    // --- Standalone baseline: store each set on its own over the union
    //     universe (the same domain the hierarchy shares at its root). ---
    uint64_t baseline_bits = 0, total_n = 0;
    {
        static uint64_t W[MAX_WORDS];
        static uint8_t  ub[MAX_ELEMENT_VALUE];
        uint32_t padded_n;
        for (uint32_t i = 0; i < num_sets; i++)
        {
            uint32_t k = 0;
            for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
                if (scratch_universe[e] == 1) ub[k++] = input[i][e];

            pack_bits(ub, k, W, &padded_n);
            Bitvector bv = bitvector_build(W, padded_n, BLOCK_SIZE, SUPERBLOCK_FREQ);
            baseline_bits += payload_bits(bv, cb);
            total_n += (uint64_t)gen_sizes[i];
        }
    }

    // --- Containment hierarchy: timed construction, then compression. ---

    clock_gettime(CLOCK_MONOTONIC, &t0);
    load_sets(input, num_sets, &col);
    sort_sets_by_size_desc(&col);
    build_inverted_list(&col, &il);
    build_hierarchy(&col, &il);
    contract_parents(&col);
    snapshot_ground_truth();
    build_relative_bitvectors(&col);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ms_build = ms(t0, t1);

    compress_all(&col);

    // --- Space of the containment representation: universe + every set. ---
    uint64_t ceb_bits = payload_bits(col.universe.compressed, cb);
    for (uint32_t i = 0; i < num_sets; i++)
        ceb_bits += payload_bits(col.nodes[i].compressed, cb);

    // --- Hierarchy depth (max and mean). ---
    uint32_t max_depth = 0, sum_depth = 0;
    for (uint32_t i = 0; i < num_sets; i++)
    {
        uint32_t d = node_depth(i);
        if (d > max_depth) max_depth = d;
        sum_depth += d;
    }

    // --- Query latency and correctness for retrieve and rank. ---
    long long ret_pass = 0, ret_fail = 0, rank_pass = 0, rank_fail = 0;
    double ms_ret = 0, ms_rank = 0;

    for (int q = 0; q < num_queries; q++)
    {
        uint32_t s = rng_below(num_sets);
        if (gt_count[s] == 0) continue;

        uint32_t k = 1 + rng_below((uint32_t)gt_count[s]);
        clock_gettime(CLOCK_MONOTONIC, &t0);
        uint32_t got = retrieve(&col, s, k);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        ms_ret += ms(t0, t1);
        if ((int)got == gt_elems[s][k - 1]) ret_pass++; else ret_fail++;
    }

    for (int q = 0; q < num_queries; q++)
    {
        uint32_t s = rng_below(num_sets);
        uint32_t cut = rng_below((uint32_t)universe_range + 1);

        int expected = 0;
        for (int t = 0; t < gt_count[s]; t++)
            if (gt_elems[s][t] <= (int)cut) expected++;

        clock_gettime(CLOCK_MONOTONIC, &t0);
        uint32_t got = rank_set(&col, s, cut);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        ms_rank += ms(t0, t1);
        if ((int)got == expected) rank_pass++; else rank_fail++;
    }

    scratch_free();

    long long ret_total = ret_pass + ret_fail;
    long long rank_total = rank_pass + rank_fail;

    printf("Sets=%u  universe=%u  elements=%llu\n\n",
           num_sets, col.universe.size, (unsigned long long)total_n);
    printf("Space (baseline)    : %.3f bits/element  (%llu bits)\n",
           total_n ? (double)baseline_bits / total_n : 0.0,
           (unsigned long long)baseline_bits);
    printf("Space (containment) : %.3f bits/element  (%llu bits)\n",
           total_n ? (double)ceb_bits / total_n : 0.0,
           (unsigned long long)ceb_bits);
    printf("Space saving        : %.1f%%\n\n",
           baseline_bits ? 100.0 * (1.0 - (double)ceb_bits / baseline_bits) : 0.0);
    printf("Depth               : max %u   mean %.2f\n",
           max_depth, num_sets ? (double)sum_depth / num_sets : 0.0);
    printf("Build               : %.3f ms\n", ms_build);
    printf("Retrieve            : %lld/%lld passed   %.6f ms/query   %s\n",
           ret_pass, ret_total, ret_total ? ms_ret / ret_total : 0.0,
           ret_fail ? "FAIL" : "OK");
    printf("Rank                : %lld/%lld passed   %.6f ms/query   %s\n",
           rank_pass, rank_total, rank_total ? ms_rank / rank_total : 0.0,
           rank_fail ? "FAIL" : "OK");

    return (ret_fail || rank_fail) ? 1 : 0;
}

int main(void)
{
    baseline_init();
    scratch_alloc();
    return run_synthetic(64, 200, 512, 10000);
}

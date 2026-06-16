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

//scratch buffers used only while building
static uint8_t (*scratch_bits)[MAX_ELEMENT_VALUE] = NULL;
static uint8_t *scratch_universe = NULL;

//allocate the temporary build buffers
static void scratch_alloc(void)
{
    scratch_bits = calloc(MAX_SETS, sizeof(uint8_t[MAX_ELEMENT_VALUE]));
    scratch_universe = calloc(MAX_ELEMENT_VALUE, sizeof(uint8_t));
}

//free the temporary build buffers
static void scratch_free(void)
{
    free(scratch_bits);
    scratch_bits = NULL;
    free(scratch_universe);
    scratch_universe = NULL;
}

//build the universe as the union of all sets and count its size
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

//copy each set into scratch and record its size
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

//sort the sets from largest to smallest keeping their bits aligned
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

//for each value list which sets hold it
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

//give each set its smallest superset as a parent
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
                if (contains_e[j] == 0) is_superset[j] = 0;
            }
        }

        is_superset[i] = 0;

        uint32_t best_size = UINT32_MAX;
        uint32_t best_index = PARENT_IS_ROOT;
        for (uint32_t j = 0; j < col->num_sets; j++)
        {
            if (is_superset[j] == 1
                && col->nodes[j].size > col->nodes[i].size
                && col->nodes[j].size < best_size)
            {
                best_size = col->nodes[j].size;
                best_index = j;
            }
        }

        col->nodes[i].parent_index = best_index;
    }
}

//pull each parent up to the highest ancestor within twice the set size
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

//rewrite each set as bits inside its parent
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

//pack a bit array into padded words for the encoder
static void pack_bits(uint8_t *src, uint32_t len, uint64_t *W, uint32_t *padded_n, int b)
{
    uint32_t pad = (uint32_t)b * SUPERBLOCK_FREQ;
    uint32_t n = ((len + pad - 1) / pad) * pad;
    if (n == 0) n = pad;
    *padded_n = n;

    uint32_t words = (n + 63) / 64;
    for (uint32_t i = 0; i < words; i++)
    {
        W[i] = 0;
    }

    for (uint32_t i = 0; i < len; i++)
    {
        if (src[i] == 1) bit_set(W, i + 1);
    }
}

//compress the universe and every set at block size b
void compress_all(SetCollection *col, int b)
{
    uint64_t W[MAX_WORDS];
    uint32_t padded_n;

    pack_bits(scratch_universe, MAX_ELEMENT_VALUE, W, &padded_n, b);
    col->universe.compressed = bitvector_build(W, padded_n, b, SUPERBLOCK_FREQ);

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

        pack_bits(scratch_bits[i], parent_size, W, &padded_n, b);
        col->nodes[i].compressed = bitvector_build(W, padded_n, b, SUPERBLOCK_FREQ);
    }
}

//return the kth element of a set by walking up to the root
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

//count the elements of a set at most k by walking up to the root
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

//turn a list of values into a bit row
void add_set(uint8_t *row, int *elements, int n)
{
    for (int j = 0; j < MAX_ELEMENT_VALUE; j++)
    {
        row[j] = 0;
    }
    for (int i = 0; i < n; i++)
    {
        row[elements[i]] = 1;
    }
}

static SetCollection col;
static InvertedList il;
static uint8_t input[MAX_SETS][MAX_ELEMENT_VALUE];

//plain bitvectors one per set kept for query timing
static Bitvector base_bv[MAX_SETS];

//generated values and a sorted copy of the true sets
static int gen_sets[MAX_SETS][MAX_ELEMENT_VALUE];
static int gen_sizes[MAX_SETS];
static int gt_elems[MAX_SETS][MAX_ELEMENT_VALUE];
static int gt_count[MAX_SETS];
static int g_pool[MAX_ELEMENT_VALUE];

static uint64_t g_rng;

//seed the random generator
static void rng_seed(uint64_t s)
{
    g_rng = s ? s : 1;
}

//next random number
static uint32_t rng_u32(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 7;
    g_rng ^= g_rng << 17;
    return (uint32_t)g_rng;
}

//random number below m
static uint32_t rng_below(uint32_t m)
{
    return m ? rng_u32() % m : 0;
}

//pick k different values from src into dst
static void sample_subset(int *dst, const int *src, int src_len, int k)
{
    for (int i = 0; i < src_len; i++)
    {
        g_pool[i] = src[i];
    }
    for (int i = 0; i < k; i++)
    {
        int r = i + (int)rng_below((uint32_t)(src_len - i));
        int tmp = g_pool[i];
        g_pool[i] = g_pool[r];
        g_pool[r] = tmp;
        dst[i] = g_pool[i];
    }
}

//make a deep nested chain that is good for containment
static uint32_t generate_nested(uint32_t num_sets, int base_size, int universe_range)
{
    static int universe_vals[MAX_ELEMENT_VALUE];
    for (int i = 0; i < universe_range; i++)
    {
        universe_vals[i] = i + 1;
    }

    sample_subset(gen_sets[0], universe_vals, universe_range, base_size);
    gen_sizes[0] = base_size;

    for (uint32_t i = 1; i < num_sets; i++)
    {
        uint32_t parent = i - 1;
        int ps = gen_sizes[parent];
        int factor = 93 + (int)rng_below(6);
        int child_size = ps * factor / 100;
        if (child_size >= ps) child_size = ps - 1;
        if (child_size < 2) child_size = 2;

        sample_subset(gen_sets[i], gen_sets[parent], ps, child_size);
        gen_sizes[i] = child_size;
    }

    for (uint32_t i = 0; i < num_sets; i++)
    {
        add_set(input[i], gen_sets[i], gen_sizes[i]);
    }

    return num_sets;
}

//make random sets of similar size with no containment
static uint32_t generate_independent(uint32_t num_sets, int base_size, int universe_range)
{
    (void)base_size;
    static int universe_vals[MAX_ELEMENT_VALUE];
    for (int i = 0; i < universe_range; i++)
    {
        universe_vals[i] = i + 1;
    }

    int lo = universe_range / 3;
    int hi = universe_range / 2;
    int span = hi - lo > 0 ? hi - lo : 1;

    for (uint32_t i = 0; i < num_sets; i++)
    {
        int sz = lo + (int)rng_below((uint32_t)span);
        if (sz < 1) sz = 1;
        sample_subset(gen_sets[i], universe_vals, universe_range, sz);
        gen_sizes[i] = sz;
        add_set(input[i], gen_sets[i], gen_sizes[i]);
    }

    return num_sets;
}

//bits needed to store one class value
static int class_bits(int b)
{
    int bits = 0;
    for (int t = b; t > 0; t >>= 1)
    {
        bits++;
    }
    return bits ? bits : 1;
}

//stored bits of one compressed bitvector
static uint64_t payload_bits(Bitvector bv, int cb)
{
    return bv.num_blocks * (uint64_t)cb + (bv.o_pos - 1);
}

//logical bits needed to store a value in the range 0..maxval
static int bits_for(uint64_t maxval)
{
    int bits = 0;
    while (maxval > 0)
    {
        bits++;
        maxval >>= 1;
    }
    return bits ? bits : 1;
}

//logical bits of the query index (rank, offset-position and select arrays) of one compressed bitvector
static uint64_t index_bits(Bitvector bv)
{
    uint64_t ones = 0;
    for (uint64_t i = 0; i < bv.num_blocks; i++)
    {
        ones += bv.C[i];
    }

    uint64_t num_checkpoints = (bv.num_blocks + bv.k - 1) / bv.k;
    int rank_w = bits_for(bv.n);        
    int pos_w  = bits_for(bv.o_pos);   
    int sel_w  = bits_for(bv.n);       

    return num_checkpoints * (uint64_t)rank_w
         + num_checkpoints * (uint64_t)pos_w
         + ones * (uint64_t)sel_w;
}

//save each sets true values before they are made relative
static void record_brute_force(void)
{
    for (uint32_t i = 0; i < col.num_sets; i++)
    {
        gt_count[i] = 0;
        for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
        {
            if (scratch_bits[i][e] == 1) gt_elems[i][gt_count[i]++] = (int)e;
        }
    }
}

//number of steps from a set up to the root
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

//milliseconds between two times
static double ms(struct timespec a, struct timespec b)
{
    return (b.tv_sec - a.tv_sec) * 1000.0 + (b.tv_nsec - a.tv_nsec) / 1e6;
}

//bits needed to choose k things out of n
static double log2_choose(uint32_t n, uint32_t k)
{
    if (k > n) return 0.0;
    return (lgamma((double)n + 1.0) - lgamma((double)k + 1.0)
            - lgamma((double)(n - k) + 1.0)) / log(2.0);
}

//baseline kth element found through the universe
static uint32_t baseline_retrieve(Bitvector set_bv, Bitvector uni_bv, uint32_t k)
{
    uint32_t p = (uint32_t)bitvector_select(set_bv, (int)k);
    return (uint32_t)bitvector_select(uni_bv, (int)p) - 1;
}

//baseline count at most k found through the universe
static uint32_t baseline_rank(Bitvector set_bv, Bitvector uni_bv, uint32_t k)
{
    uint32_t p = (uint32_t)bitvector_rank(uni_bv, (int)(k + 1));
    if (p == 0) return 0;
    return (uint32_t)bitvector_rank(set_bv, (int)p);
}

#define MAX_QUERIES 5000

//query streams shared by both solutions so the test is fair
static int q_ret_s[MAX_QUERIES], q_ret_k[MAX_QUERIES], q_ret_exp[MAX_QUERIES];
static int q_rank_s[MAX_QUERIES], q_rank_cut[MAX_QUERIES], q_rank_exp[MAX_QUERIES];
static uint8_t amask[MAX_ELEMENT_VALUE];

//build one plain bitvector per set against the universe
static void build_baseline(uint32_t num_sets, int b)
{
    static uint64_t W[MAX_WORDS];
    static uint8_t ub[MAX_ELEMENT_VALUE];
    uint32_t padded_n;

    for (uint32_t i = 0; i < num_sets; i++)
    {
        for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
        {
            amask[e] = 0;
        }
        for (int t = 0; t < gt_count[i]; t++)
        {
            amask[gt_elems[i][t]] = 1;
        }

        uint32_t k = 0;
        for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
        {
            if (scratch_universe[e] == 1) ub[k++] = amask[e];
        }

        pack_bits(ub, k, W, &padded_n, b);
        base_bv[i] = bitvector_build(W, padded_n, b, SUPERBLOCK_FREQ);
    }
}

//make random queries and their correct answers
static void make_queries(uint32_t num_sets, int universe_range, int num_queries, int *nq_ret, int *nq_rank)
{
    int r = 0;
    for (int q = 0; q < num_queries; q++)
    {
        uint32_t s = rng_below(num_sets);
        if (gt_count[s] == 0) continue;
        q_ret_s[r] = (int)s;
        q_ret_k[r] = 1 + (int)rng_below((uint32_t)gt_count[s]);
        q_ret_exp[r] = gt_elems[s][q_ret_k[r] - 1];
        r++;
    }

    int n = 0;
    for (int q = 0; q < num_queries; q++)
    {
        uint32_t s = rng_below(num_sets);
        int cut = (int)rng_below((uint32_t)universe_range + 1);
        int expected = 0;
        for (int t = 0; t < gt_count[s]; t++)
        {
            if (gt_elems[s][t] <= cut) expected++;
        }
        q_rank_s[n] = (int)s;
        q_rank_cut[n] = cut;
        q_rank_exp[n] = expected;
        n++;
    }

    *nq_ret = r;
    *nq_rank = n;
}

//run one setup and write its row of results
static int run_config(int dataset, uint64_t seed, int b, int use_contraction,
                      uint32_t num_sets, int base_size, int universe_range,
                      int num_queries, FILE *csv)
{
    int cb = class_bits(b);
    struct timespec t0, t1;
    if (num_queries > MAX_QUERIES) num_queries = MAX_QUERIES;

    rng_seed(seed);
    if (dataset == 0)
    {
        num_sets = generate_nested(num_sets, base_size, universe_range);
    }
    else
    {
        num_sets = generate_independent(num_sets, base_size, universe_range);
    }

    build_universe(input, num_sets, &col);

    clock_gettime(CLOCK_MONOTONIC, &t0);
    load_sets(input, num_sets, &col);
    sort_sets_by_size_desc(&col);
    build_inverted_list(&col, &il);
    build_hierarchy(&col, &il);
    if (use_contraction) contract_parents(&col);
    record_brute_force();
    build_relative_bitvectors(&col);
    compress_all(&col, b);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ceb_build_ms = ms(t0, t1);

    clock_gettime(CLOCK_MONOTONIC, &t0);
    build_baseline(num_sets, b);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double base_build_ms = ms(t0, t1);

    uint64_t ceb_bits = payload_bits(col.universe.compressed, cb);
    uint64_t ceb_index = index_bits(col.universe.compressed);
    double L_bits = 0.0;
    for (uint32_t i = 0; i < num_sets; i++)
    {
        ceb_bits += payload_bits(col.nodes[i].compressed, cb);
        ceb_index += index_bits(col.nodes[i].compressed);
        uint32_t psz = col.nodes[i].parent_index == PARENT_IS_ROOT
                     ? col.universe.size
                     : col.nodes[col.nodes[i].parent_index].size;
        L_bits += log2_choose(psz, col.nodes[i].size);
    }

    uint64_t baseline_bits = 0, base_index = 0, total_n = 0;
    double H_bits = 0.0;
    for (uint32_t i = 0; i < num_sets; i++)
    {
        baseline_bits += payload_bits(base_bv[i], cb);
        base_index += index_bits(base_bv[i]);
        total_n += (uint64_t)gt_count[i];
        H_bits += log2_choose(col.universe.size, (uint32_t)gt_count[i]);
    }

    uint32_t max_depth = 0, sum_depth = 0;
    for (uint32_t i = 0; i < num_sets; i++)
    {
        uint32_t d = node_depth(i);
        if (d > max_depth) max_depth = d;
        sum_depth += d;
    }

    int nq_ret, nq_rank;
    make_queries(num_sets, universe_range, num_queries, &nq_ret, &nq_rank);

    long long ret_pass = 0, rank_pass = 0;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int q = 0; q < nq_ret; q++)
    {
        if ((int)retrieve(&col, q_ret_s[q], q_ret_k[q]) == q_ret_exp[q]) ret_pass++;
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ms_ret = ms(t0, t1);

    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int q = 0; q < nq_rank; q++)
    {
        if ((int)rank_set(&col, q_rank_s[q], q_rank_cut[q]) == q_rank_exp[q]) rank_pass++;
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ms_rank = ms(t0, t1);

    long long bret_pass = 0, brank_pass = 0;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int q = 0; q < nq_ret; q++)
    {
        if ((int)baseline_retrieve(base_bv[q_ret_s[q]], col.universe.compressed,
                                   (uint32_t)q_ret_k[q]) == q_ret_exp[q]) bret_pass++;
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ms_bret = ms(t0, t1);

    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int q = 0; q < nq_rank; q++)
    {
        if ((int)baseline_rank(base_bv[q_rank_s[q]], col.universe.compressed,
                               (uint32_t)q_rank_cut[q]) == q_rank_exp[q]) brank_pass++;
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ms_brank = ms(t0, t1);

    double base_bpe = (double)baseline_bits / total_n;
    double ceb_bpe = (double)ceb_bits / total_n;
    double H_bpe = H_bits / total_n;
    double L_bpe = L_bits / total_n;
    double saving = 100.0 * (1.0 - (double)ceb_bits / baseline_bits);
    double base_total_bpe = (double)(baseline_bits + base_index) / total_n;
    double ceb_total_bpe = (double)(ceb_bits + ceb_index) / total_n;
    double total_saving = 100.0 * (1.0 - (double)(ceb_bits + ceb_index)
                                         / (baseline_bits + base_index));
    double mean_depth = (double)sum_depth / num_sets;
    double ret_us = 1000.0 * ms_ret / nq_ret;
    double rank_us = 1000.0 * ms_rank / nq_rank;
    double bret_us = 1000.0 * ms_bret / nq_ret;
    double brank_us = 1000.0 * ms_brank / nq_rank;

    int all_ok = (ret_pass == nq_ret) && (rank_pass == nq_rank)
              && (bret_pass == nq_ret) && (brank_pass == nq_rank);

    const char *name = dataset == 0 ? "nested" : "independent";

    fprintf(csv,
        "%s,%llu,%d,%d,%u,%u,%llu,%llu,%llu,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d\n",
        name, (unsigned long long)seed, b, use_contraction, num_sets, col.universe.size,
        (unsigned long long)total_n,
        (unsigned long long)baseline_bits, (unsigned long long)ceb_bits,
        base_bpe, ceb_bpe, H_bpe, L_bpe, saving,
        base_total_bpe, ceb_total_bpe, total_saving,
        ceb_build_ms, base_build_ms, max_depth, mean_depth,
        ret_us, rank_us, bret_us, brank_us,
        all_ok, nq_ret, nq_rank);

    printf("%-11s seed=%-6llu b=%2d con=%d  base=%5.2f ceb=%5.2f bpe  save=%5.1f%%  depth=%2u/%.2f  build ceb=%.3f base=%.3f ms  %s\n",
        name, (unsigned long long)seed, b, use_contraction,
        base_bpe, ceb_bpe, saving, max_depth, mean_depth,
        ceb_build_ms, base_build_ms, all_ok ? "OK" : "FAIL");

    return all_ok ? 0 : 1;
}

//pass through every setup and write the results file
int main(void)
{
    baseline_init();
    scratch_alloc();

    const char *path = "results.csv";
    FILE *csv = fopen(path, "w");
    if (!csv)
    {
        perror("fopen");
        return 1;
    }

    fprintf(csv,
        "dataset,seed,b,contraction,num_sets,universe,elements,"
        "baseline_bits,ceb_bits,baseline_bpe,ceb_bpe,H_bpe,L_bpe,saving_pct,"
        "baseline_total_bpe,ceb_total_bpe,total_saving_pct,"
        "ceb_build_ms,base_build_ms,max_depth,mean_depth,"
        "ceb_retrieve_us,ceb_rank_us,base_retrieve_us,base_rank_us,"
        "all_ok,num_ret,num_rank\n");

    uint64_t seeds[] = {1, 7, 42, 101, 1337, 2024, 90210, 555};
    int bs[] = {4, 8, 16, 32, 64};
    int num_seeds = (int)(sizeof(seeds) / sizeof(seeds[0]));
    int num_bs = (int)(sizeof(bs) / sizeof(bs[0]));
    int fails = 0;

    for (int d = 0; d < 2; d++)
    {
        for (int si = 0; si < num_seeds; si++)
        {
            for (int bi = 0; bi < num_bs; bi++)
            {
                for (int con = 0; con < 2; con++)
                {
                    fails += run_config(d, seeds[si], bs[bi], con, 64, 480, 512, 2000, csv);
                }
            }
        }
    }

    fclose(csv);
    scratch_free();

    printf("\nWrote %s   %s\n", path, fails ? "SOME FAILS" : "ALL OK");
    return fails ? 1 : 0;
}

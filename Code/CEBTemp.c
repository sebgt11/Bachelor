#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "Baseline.h"

#define MAX_ELEMENT_VALUE 32
#define MAX_SETS 16
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

//print the stored arrays of one compressed bitvector
static void print_bitvector(const char *name, Bitvector bv)
{
    uint64_t ones = 0;
    for (uint64_t i = 0; i < bv.num_blocks; i++) ones += bv.C[i];
    int num_sb = (int)((bv.num_blocks + bv.k - 1) / bv.k);

    printf("%s  (n=%llu, b=%d, k=%d, ones=%llu)\n",
           name, (unsigned long long)bv.n, bv.b, bv.k, (unsigned long long)ones);

    printf("  C = [ ");
    for (uint64_t i = 0; i < bv.num_blocks; i++) printf("%d ", bv.C[i]);
    printf("]\n");

    printf("  O = ");
    for (uint64_t i = 1; i < bv.o_pos; i++) printf("%d", bit_read(bv.O, i));
    printf("\n");

    printf("  L = [ ");
    for (int i = 0; i <= bv.b; i++) printf("%llu ", (unsigned long long)bv.L[i]);
    printf("]\n");

    printf("  R = [ ");
    for (int i = 0; i < num_sb; i++) printf("%llu ", (unsigned long long)bv.R[i]);
    printf("]\n");

    printf("  P = [ ");
    for (int i = 0; i < num_sb; i++) printf("%llu ", (unsigned long long)bv.P[i]);
    printf("]\n");

    printf("  S = [ ");
    for (uint64_t i = 0; i < ones; i++) printf("%llu ", (unsigned long long)bv.S[i]);
    printf("]\n\n");
}

//build the example hierarchy and print its stored structure
int main(void)
{
    scratch_alloc();
    baseline_init();

    //the three nested example sets from the report
    int s1[] = {1, 3, 7, 9, 10, 11, 12, 17, 18, 19};
    int s2[] = {3, 10, 12, 17, 19};
    int s3[] = {3, 12, 19};
    add_set(input[0], s1, 10);
    add_set(input[1], s2, 5);
    add_set(input[2], s3, 3);
    uint32_t num_sets = 3;

    //run the construction pipeline with no timing or data collection
    build_universe(input, num_sets, &col);
    load_sets(input, num_sets, &col);
    sort_sets_by_size_desc(&col);
    build_inverted_list(&col, &il);
    build_hierarchy(&col, &il);
    contract_parents(&col);
    build_relative_bitvectors(&col);
    compress_all(&col, BLOCK_SIZE);

    //print the hierarchy that was built
    printf("hierarchy (%u sets, universe size %u)\n", num_sets, col.universe.size);
    for (uint32_t i = 0; i < num_sets; i++)
    {
        if (col.nodes[i].parent_index == PARENT_IS_ROOT)
            printf("  set %u  size %u  parent = universe\n", i, col.nodes[i].size);
        else
            printf("  set %u  size %u  parent = set %u\n",
                   i, col.nodes[i].size, col.nodes[i].parent_index);
    }
    printf("\n");

    //print the stored arrays of the universe and every set
    print_bitvector("universe", col.universe.compressed);
    for (uint32_t i = 0; i < num_sets; i++)
    {
        char label[16];
        snprintf(label, sizeof(label), "set %u", i);
        print_bitvector(label, col.nodes[i].compressed);
    }

    //demonstrate retrieve and rank on the smallest set
    uint32_t s = num_sets - 1;
    printf("queries on set %u (size %u)\n", s+1, col.nodes[s].size);
    printf("  retrieve k = 1, 2 %u : ", col.nodes[s].size);
    for (uint32_t k = 1; k <= col.nodes[s].size; k++)
        printf("%u ", retrieve(&col, s, k));
    printf("\n");
    printf("  rank at cutoff 6  = %u\n", rank_set(&col, s, 6));
    printf("  rank at cutoff 12 = %u\n", rank_set(&col, s, 12));

    scratch_free();
    return 0;
}

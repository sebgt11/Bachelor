// Containment entropy with baseline-compressed bitvectors.
// Construction uses a temporary uint8_t scratch buffer which is freed after
// compress_all runs. After that, the only stored form of each set is the
// compressed Bitvector from Baseline.h.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

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

static void scratch_alloc(void)
{
    scratch_bits = calloc(MAX_SETS, sizeof(uint8_t[MAX_ELEMENT_VALUE]));
    scratch_universe = calloc(MAX_ELEMENT_VALUE, sizeof(uint8_t));
}

static void scratch_free(void)
{
    free(scratch_bits);
    scratch_bits = NULL;
    free(scratch_universe);
    scratch_universe = NULL;
}


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


// Pack the first len bits of src into a uint64_t word array W. Baseline
// requires the length to be a multiple of BLOCK_SIZE * SUPERBLOCK_FREQ, so we
// pad up to the next multiple with zeros.
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


void add_set(uint8_t *row, int *elements, int n)
{
    for (int j = 0; j < MAX_ELEMENT_VALUE; j++) row[j] = 0;
    for (int i = 0; i < n; i++) row[elements[i]] = 1;
}


static SetCollection col;
static InvertedList il;
static uint8_t input[MAX_SETS][MAX_ELEMENT_VALUE];

int main(void)
{
    baseline_init();
    scratch_alloc();

    int big[] = {1, 3, 5, 7, 9};
    int medium[] = {3, 5, 9};
    int small[] = {5, 9};

    add_set(input[0], big, 5);
    add_set(input[1], medium, 3);
    add_set(input[2], small, 2);

    build_universe(input, 3, &col);
    load_sets(input, 3, &col);
    sort_sets_by_size_desc(&col);
    build_inverted_list(&col, &il);
    build_hierarchy(&col, &il);
    contract_parents(&col);
    build_relative_bitvectors(&col);
    compress_all(&col);
    scratch_free();

    uint32_t big_idx = 0;
    uint32_t medium_idx = 1;
    uint32_t small_idx = 2;

    printf("rank(Big, 6) = %u (expected 3)\n", rank_set(&col, big_idx, 6));
    printf("rank(Medium, 6) = %u (expected 2)\n", rank_set(&col, medium_idx, 6));
    printf("rank(Small, 6) = %u (expected 1)\n", rank_set(&col, small_idx, 6));

    printf("retrieve(Big, 4) = %u (expected 7)\n", retrieve(&col, big_idx, 4));
    printf("retrieve(Medium, 3) = %u (expected 9)\n", retrieve(&col, medium_idx, 3));
    printf("retrieve(Small, 1) = %u (expected 5)\n", retrieve(&col, small_idx, 1));

    return 0;
}

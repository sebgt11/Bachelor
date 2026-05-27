#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

//General structure
#define MAX_ELEMENT_VALUE 1024
#define MAX_SETS 128

#define PARENT_IS_ROOT 0xFFFFFFFF

typedef struct {
    uint8_t bits[MAX_ELEMENT_VALUE];
    uint32_t size;
} Universe;

typedef struct {
    uint32_t size;
    uint8_t bitvector[MAX_ELEMENT_VALUE]; //
    uint32_t parent_index;
} SetNode;

typedef struct {
    Universe universe;
    SetNode nodes[MAX_SETS];
    uint32_t num_sets;
} SetCollection;

typedef struct {
    uint32_t set_indices[MAX_SETS];  // the sets that contain that element
    uint32_t count; //how many sets are in list
} InvertedListEntry;

typedef struct {
    InvertedListEntry entries[MAX_ELEMENT_VALUE]; //list of InvertedListEntries
} InvertedList;


//checks if element is in set
int element_in_set(SetNode *node, uint32_t elem) {
    return node->bitvector[elem] == 1;
}

int is_strict_subset(SetNode *a, SetNode *b) {

    // cant be strict subset if same size or larger
    if (a->size >= b->size) return 0;

    for (uint32_t i = 0; i < MAX_ELEMENT_VALUE; i++) {
        if (a->bitvector[i] == 1 && b->bitvector[i] == 0) 
        {
            return 0;  // a has an element b doesnt
        }
    }

    return 1;
}

void build_universe(uint8_t input_sets[][MAX_ELEMENT_VALUE],
                    uint32_t num_sets, SetCollection *col)
{

    // start with all zeros
    col->universe.size = 0;

    for (uint32_t i = 0; i < MAX_ELEMENT_VALUE; i++) 
    {
        col->universe.bits[i] = 0;
    }

    // OR every set into the universe
    for (uint32_t i = 0; i < num_sets; i++) 
    {
        for (uint32_t j = 0; j < MAX_ELEMENT_VALUE; j++) 
        {
            col->universe.bits[j] |= input_sets[i][j];
        }
    }

    // count how many 1s ended up in U
    for (uint32_t i = 0; i < MAX_ELEMENT_VALUE; i++) 
    {
        if (col->universe.bits[i] == 1) 
        {
            col->universe.size++;
        }
    }
}

//Copies input bitvectors into the collection and computes each sets size
void load_sets(uint8_t input_sets[][MAX_ELEMENT_VALUE],
               uint32_t num_sets,
               SetCollection *col) {

    col->num_sets = num_sets;

    for (uint32_t i = 0; i < num_sets; i++) 
    {

        // reset size
        col->nodes[i].size = 0;

        //copy bits and count
        for (uint32_t j = 0; j < MAX_ELEMENT_VALUE; j++) 
        {
            col->nodes[i].bitvector[j] = input_sets[i][j];

            if (input_sets[i][j] == 1) 
            {
                col->nodes[i].size++;
            }
        }

        //parent is unknown so we do this for now
        col->nodes[i].parent_index = PARENT_IS_ROOT;
    }
}

//sorts set by descending order it sues insertion sort
void sort_sets_by_size_desc(SetCollection *col)
{
    for (uint32_t i = 1; i < col->num_sets; i++)
    {
        SetNode tmp = col->nodes[i];
        int32_t j = (int32_t)i - 1;
        while (j >= 0 && col->nodes[j].size < tmp.size)
        {
            col->nodes[j + 1] = col->nodes[j];
            j--;
        }
        col->nodes[j + 1] = tmp;
    }
}

//	For each element it records the indices of all sets containing it
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
            if (col->nodes[i].bitvector[e] == 1)
            {
                uint32_t c = il->entries[e].count;
                il->entries[e].set_indices[c] = i;
                il->entries[e].count++;
            }
        }
    }
}

//Finds each sets smallest superset by intersecting inverted list entries
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
            if (col->nodes[i].bitvector[e] == 0) continue;

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

        uint32_t best_size  = UINT32_MAX;
        uint32_t best_index = PARENT_IS_ROOT;

        for (uint32_t j = 0; j < col->num_sets; j++)
        {
            if (is_superset[j] == 1 && col->nodes[j].size < best_size)
            {
                best_size  = col->nodes[j].size;
                best_index = j;
            }
        }

        col->nodes[i].parent_index = best_index;
    }
}

//Replaces each parent with the highest ancestor whose size is at most 2x
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

//makes bitvectors relative to parent
void build_relative_bitvectors(SetCollection *col)
{
    for (int32_t i = (int32_t)col->num_sets - 1; i >= 0; i--)
    {
        uint8_t *parent_bits;
        if (col->nodes[i].parent_index == PARENT_IS_ROOT)
        {
            parent_bits = col->universe.bits;
        }
        else
        {
            parent_bits = col->nodes[col->nodes[i].parent_index].bitvector;
        }

        uint8_t new_bitvector[MAX_ELEMENT_VALUE];
        uint32_t k = 0;

        for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
        {
            if (parent_bits[e] == 1)
            {
                new_bitvector[k] = col->nodes[i].bitvector[e];
                k++;
            }
        }

        for (uint32_t e = 0; e < MAX_ELEMENT_VALUE; e++)
        {
            col->nodes[i].bitvector[e] = new_bitvector[e];
        }
    }
}

//Naive scan returning the position of the k-th 1-bit
uint32_t select_in_bitvector(uint8_t *bv, uint32_t len, uint32_t k)
{
    uint32_t count = 0;
    for (uint32_t i = 0; i < len; i++)
    {
        if (bv[i] == 1)
        {
            count++;
            if (count == k) return i;
        }
    }
    return UINT32_MAX;
}

//Recursively walks up the hierarchy to translate a position into a universe element
uint32_t retrieve(SetCollection *col, uint32_t set_idx, uint32_t k)
{
    SetNode *node = &col->nodes[set_idx];

    if (node->parent_index == PARENT_IS_ROOT)
    {
        uint32_t p = select_in_bitvector(node->bitvector, col->universe.size, k);
        return select_in_bitvector(col->universe.bits, MAX_ELEMENT_VALUE, p + 1);
    }

    uint32_t parent_size = col->nodes[node->parent_index].size;
    uint32_t p = select_in_bitvector(node->bitvector, parent_size, k);
    return retrieve(col, node->parent_index, p + 1);
}

//Counts the number of 1-bits in bitvector
uint32_t rank_in_bitvector(uint8_t *bv, uint32_t i)
{
    uint32_t count = 0;
    for (uint32_t j = 0; j <= i; j++)
    {
        if (bv[j] == 1) count++;
    }
    return count;
}
//Recursively translates the cutoff k through the hierarchy and counts elements <= k
uint32_t rank_set(SetCollection *col, uint32_t set_idx, uint32_t k)
{
    SetNode *node = &col->nodes[set_idx];

    uint32_t p;
    if (node->parent_index == PARENT_IS_ROOT)
    {
        p = rank_in_bitvector(col->universe.bits, k);
    }
    else
    {
        p = rank_set(col, node->parent_index, k);
    }

    if (p == 0) return 0;
    return rank_in_bitvector(node->bitvector, p - 1);
}

//Helper that converts a list of element values into a 0/1 bitvector row
void add_set(uint8_t *row, int *elements, int n)
{
    for (int j = 0; j < MAX_ELEMENT_VALUE; j++) row[j] = 0;
    for (int i = 0; i < n; i++) row[elements[i]] = 1;
}

int main(void)
{
    int big[] = {1, 3, 5, 7, 9};
    int medium[] = {3, 5, 9};
    int small[] = {5, 9};

    uint8_t input[3][MAX_ELEMENT_VALUE];
    add_set(input[0], big, 5);
    add_set(input[1], medium, 3);
    add_set(input[2], small, 2);

    SetCollection col;
    InvertedList  il;

    build_universe(input, 3, &col);
    load_sets(input, 3, &col);
    sort_sets_by_size_desc(&col);
    build_inverted_list(&col, &il);
    build_hierarchy(&col, &il);
    contract_parents(&col);
    build_relative_bitvectors(&col);

    uint32_t big_idx = 0, medium_idx = 1, small_idx = 2;

    printf("rank(Big, 6) = %u  (expected 3)\n", rank_set(&col, big_idx, 6));
    printf("rank(Medium, 6) = %u  (expected 2)\n", rank_set(&col, medium_idx, 6));
    printf("rank(Small, 6) = %u  (expected 1)\n", rank_set(&col, small_idx, 6));

    printf("retrieve(Big 4) = %u  (expected 7)\n", retrieve(&col, big_idx, 4));
    printf("retrieve(Medium 3) = %u  (expected 9)\n", retrieve(&col, medium_idx, 3));
    printf("retrieve(Small 1) = %u  (expected 5)\n", retrieve(&col, small_idx, 1));

    return 0;
}
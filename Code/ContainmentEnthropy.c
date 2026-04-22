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
    uint8_t bits[MAX_ELEMENT_VALUE]; //element number i is in this set.
    uint32_t size;
    uint8_t bitvector[MAX_ELEMENT_VALUE]; // element number j of the parent is also in the current set
    uint32_t parent_index;
} SetNode;

typedef struct {
    Universe universe;
    SetNode nodes[MAX_SETS];
    uint32_t num_sets;
} SetCollection;


//checks if element is in set
int element_in_set(SetNode *node, uint32_t elem) {
    return node->bits[elem] == 1;
}

int is_strict_subset(SetNode *a, SetNode *b) {

    // cant be strict subset if same size or larger
    if (a->size >= b->size) return 0;

    for (uint32_t i = 0; i < MAX_ELEMENT_VALUE; i++) {
        if (a->bits[i] == 1 && b->bits[i] == 0) {
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

    for (uint32_t i = 0; i < MAX_ELEMENT_VALUE; i++) {
        col->universe.bits[i] = 0;
    }

    // OR every set into the universe
    for (uint32_t i = 0; i < num_sets; i++) {
        for (uint32_t j = 0; j < MAX_ELEMENT_VALUE; j++) {
            col->universe.bits[j] |= input_sets[i][j];
        }
    }

    // count how many 1s ended up in U
    for (uint32_t i = 0; i < MAX_ELEMENT_VALUE; i++) {
        if (col->universe.bits[i] == 1) {
            col->universe.size++;
        }
    }
}

void load_sets(uint8_t input_sets[][MAX_ELEMENT_VALUE],
               uint32_t num_sets,
               SetCollection *col) {

    col->num_sets = num_sets;

    for (uint32_t i = 0; i < num_sets; i++) {

        // reset size
        col->nodes[i].size = 0;

        //copy bits and count
        for (uint32_t j = 0; j < MAX_ELEMENT_VALUE; j++) {
            col->nodes[i].bits[j] = input_sets[i][j];

            if (input_sets[i][j] == 1) {
                col->nodes[i].size++;
            }
        }

        //parent is unknown so we do this for now
        col->nodes[i].parent_index = PARENT_IS_ROOT;
    }
}





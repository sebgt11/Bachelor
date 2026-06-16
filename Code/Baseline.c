#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define BLOCK_SIZE      4
#define SUPERBLOCK_FREQ 4
#define MAX_N           4096
#define WORD_SIZE       64

#define MAX_BLOCKS      (MAX_N / BLOCK_SIZE)
#define MAX_SUPERBLOCKS (MAX_N / (BLOCK_SIZE * SUPERBLOCK_FREQ))
#define MAX_WORDS       ((MAX_N + 63) / 64)


typedef struct {
    uint8_t   C[MAX_BLOCKS];
    uint64_t  O[MAX_WORDS];
    uint64_t  P[MAX_SUPERBLOCKS];
    uint64_t  R[MAX_SUPERBLOCKS];
    uint64_t  S[MAX_N];
    uint64_t  K[65][65];
    uint64_t  L[65];
    uint64_t  n;
    int       b;
    int       k;
    uint64_t  num_blocks;
    uint64_t  o_pos;
} Bitvector;


//read the bit at position i (1-indexed) of the packed words
static int bit_read(uint64_t *W, uint64_t i)
{
    i--;
    return (W[i / WORD_SIZE] >> (WORD_SIZE - 1 - (i % WORD_SIZE))) & 1;
}

static void bit_set(uint64_t *W, uint64_t i)
{
    i--;
    W[i / WORD_SIZE] |= (1ULL << (WORD_SIZE - 1 - (i % WORD_SIZE)));
}

//write the lowest len bits of val into W starting at position start
static void bits_write(uint64_t *W, uint64_t start, uint64_t len, uint64_t val)
{
    for (uint64_t i = start; i <= start + len - 1; i++)
        if ((val >> (start + len - 1 - i)) & 1) bit_set(W, i);
}

//Pascal's triangle of binomial coefficients
static void populateK(uint64_t K[65][65])
{
    for (int i = 0; i <= 64; i++)
    {
        K[i][0] = 1;
        K[i][i] = 1;
        for (int j = 1; j < i; j++)
            K[i][j] = K[i-1][j-1] + K[i-1][j];
    }
}

//number of offset bits needed for each class
static void createL(int b, uint64_t K[65][65], uint64_t L[65])
{
    for (int i = 0; i <= b; i++)
    {
        uint64_t combos = K[b][i];
        int bits = 0;
        uint64_t tmp = combos - 1;
        while (tmp > 0) { bits++; tmp >>= 1; }
        L[i] = combos <= 1 ? 0 : bits;
    }
}

//turn one block into its class and its offset within that class
static void encode(uint64_t *W, uint64_t start, int b, uint64_t K[65][65],
                   int *c_out, uint64_t *o_out)
{
    *c_out = 0;
    for (int i = 1; i <= b; i++)
        if (bit_read(W, start + i - 1) == 1) (*c_out)++;

    *o_out = 0;
    int ctemp = *c_out;
    for (int j = 1; j <= b; j++)
        if (bit_read(W, start + j - 1) == 1)
        {
            *o_out += K[b - j][ctemp];
            ctemp--;
        }
}

//store the position of every 1 bit so select is a direct lookup
static void CreateS(Bitvector *bv, uint64_t *W)
{
    uint64_t count = 0;
    for (uint64_t i = 1; i <= bv->n; i++)
        if (bit_read(W, i) == 1) bv->S[count++] = i;
}

//build all the arrays for one bitvector
static Bitvector bitvector_build(uint64_t *W, uint64_t n, int b, int k)
{
    Bitvector bv = {0};
    bv.n = n;
    bv.b = b;
    bv.k = k;
    bv.num_blocks = n / b;

    populateK(bv.K);
    createL(bv.b, bv.K, bv.L);

    uint64_t o_pos = 1;
    uint64_t rank_count = 0;

    for (uint64_t i = 0; i < bv.num_blocks; i++)
    {
        uint64_t start = i * b + 1;

        if (i % k == 0)                 //checkpoint every k blocks
        {
            bv.P[i / k] = o_pos;
            bv.R[i / k] = rank_count;
        }

        int c;
        uint64_t o;
        encode(W, start, b, bv.K, &c, &o);

        bv.C[i] = c;

        int len = bv.L[c];
        if (len > 0)
        {
            bits_write(bv.O, o_pos, len, o);
            o_pos += len;
        }

        rank_count += c;
    }

    bv.o_pos = o_pos;
    CreateS(&bv, W);
    return bv;
}

int main(void)
{
    //the example set S1 from the report, as a bitvector
    const char *bits = "10100010111100001110";

    int n = (int)strlen(bits);
    int b = BLOCK_SIZE;
    int k = SUPERBLOCK_FREQ;

    //pack the 0/1 text into words
    uint64_t W[MAX_WORDS] = {0};
    int ones = 0;
    for (int i = 0; i < n; i++)
        if (bits[i] == '1') { bit_set(W, i + 1); ones++; }

    Bitvector bv = bitvector_build(W, n, b, k);

    //print the set, the bitvector, and every array
    printf("set S1     = {1,3,7,9,10,11,12,17,18,19}\n");
    printf("bitvector  = %s   (n = %d, b = %d, k = %d)\n\n", bits, n, b, k);

    printf("C = [ ");
    for (uint64_t i = 0; i < bv.num_blocks; i++) printf("%d ", bv.C[i]);
    printf("]\n");

    printf("O = ");
    for (uint64_t i = 1; i < bv.o_pos; i++) printf("%d", bit_read(bv.O, i));
    printf("\n");

    printf("L = [ ");
    for (int i = 0; i <= b; i++) printf("%llu ", (unsigned long long)bv.L[i]);
    printf("]\n");

    int num_sb = (int)((bv.num_blocks + k - 1) / k);
    printf("R = [ ");
    for (int i = 0; i < num_sb; i++) printf("%llu ", (unsigned long long)bv.R[i]);
    printf("]\n");

    printf("P = [ ");
    for (int i = 0; i < num_sb; i++) printf("%llu ", (unsigned long long)bv.P[i]);
    printf("]\n");

    printf("S = [ ");
    for (int i = 0; i < ones; i++) printf("%llu ", (unsigned long long)bv.S[i]);
    printf("]\n");

    return 0;
}

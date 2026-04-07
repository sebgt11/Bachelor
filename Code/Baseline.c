#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define WORD_SIZE 64

#define MAX_BITS 1024
#define MAX_BLOCKS 256
#define MAX_WORDS 16

typedef struct {
    uint64_t *W;  //bits
    uint8_t   C[MAX_BLOCKS];  //class of each block
    uint64_t  O[MAX_WORDS];  // offset
    uint64_t  P[MAX_BLOCKS];  // starting position of bit
    uint64_t  R[MAX_BLOCKS];  //rank for each block
    uint64_t  K[65][65];  //Pascal precomputed
    uint64_t  L[65];  //bits per class
    uint64_t  n;   // total bits
    int       b;    // block size
    int       k;  // how often P and K
    uint64_t  num_blocks;   // total # of blocks
    uint64_t  o_pos;      // current pos in O
} Bitvector;


//bit functionality:

int static inline bit_read(uint64_t *W, uint64_t i)
{
    i--;
    return (W[i / WORD_SIZE] >> (WORD_SIZE - 1 - (i % WORD_SIZE))) & 1;
}

static inline void bit_set(uint64_t *W, uint64_t i)
{
    i--;
    W[i / WORD_SIZE] |= (1ULL << (WORD_SIZE - 1 - (i % WORD_SIZE)));
}

static inline void bit_clear(uint64_t *W, uint64_t i)
{
    i--;
    W[i / WORD_SIZE] &= ~(1ULL << (WORD_SIZE - 1 - (i % WORD_SIZE)));
}

static inline void bits_write(uint64_t *W, uint64_t start, uint64_t len, uint64_t val) 
{ 
    for (uint64_t i = start; i <= start + len - 1; i++) 
    { 
        if ((val >> (start + len - 1 - i)) & 1) 
            bit_set(W, i); 
        else 
            bit_clear(W, i); 
    } 
}

static inline uint64_t bits_read_simple(uint64_t *W, uint64_t start, uint64_t stop)
{
    uint64_t temp = 0;                   

    for(int i = start; i <= (stop); i++)   
    {   
        int bit = bit_read(W, i);         

        temp = (temp << 1) | bit;         
    }

    return temp;
}



//There is a more complicated and quicker way of does bits read, implement this later
//this is more readable and simple


static inline int popcount(uint64_t x)
{
    return __builtin_popcountll(x);
}

static void populateK(uint64_t K[65][65])
{
    for(int i = 0; i <= 64; i++)
    {
        K[i][0] = 1;      
        K[i][i] = 1;      

        for(int j = 1; j < i; j++)  
        {
            K[i][j] = K[i-1][j-1] + K[i-1][j];
        }
    }
}

static void createL(int b, uint64_t K[65][65], uint64_t L[65])
{
    for(int i =0; i<=b; i++)
    {
        uint64_t combos = K[b][i];
        if(combos <= 1)
        {
            L[i] = 0;
        }
        else
        {
            int bits = 0;
            uint64_t tmp = combos - 1;
            while(tmp > 0)
            {
                bits++;
                tmp >>= 1;
            }
            L[i] = bits;
        }
    }
}

static void encode(
    uint64_t *W,        
    uint64_t start,     
    int b,             
    uint64_t K[65][65], 
    int *c_out,         
    uint64_t *o_out     
)
{
    *c_out = 0;                         
    for(int i = 1; i <= b; i++)
    {
        if(bit_read(W, start + i - 1) == 1)  
        {
            (*c_out)++;                  
        }
    }

    *o_out = 0;                         
    int ctemp = *c_out;                 

    for(int j = 1; j <= b; j++)         
    {
        if(bit_read(W, start + j - 1) == 1)  
        {
            *o_out += K[b - j][ctemp];  
            ctemp--;
        }
    }
}


static uint64_t decode(int c, uint64_t o, int b, uint64_t K[65][65])
{
    uint64_t B = 0;   

    for(int j = 1; j <= b; j++)
    {
        if(o >= K[b-j][c])        
        {
            bit_set(&B, j);       
            o -= K[b-j][c];       
            c--;                   
        }
    }

    return B;
}

static Bitvector bitvector_build(uint64_t *W, uint64_t n, int b, int k)
{
    Bitvector bv;

    bv.n = n;
    bv.b= b;
    bv.W = W;
    bv.k = k;
    bv.num_blocks = n/b;

    populateK(bv.K);
    createL(bv.b, bv.K, bv.L);



    int c = 0;
    uint64_t o = 0;

    uint64_t o_pos = 1; 
    uint64_t rank_count = 0;   

    for(uint64_t i = 0; i < bv.num_blocks; i++)
    {
        uint64_t start = i * b + 1;   

        if(i % k == 0)
        {
            bv.P[i / k] = o_pos;      
            bv.R[i / k] = rank_count;
        }

        int c;
        uint64_t o;
        encode(W, start, b, bv.K, &c, &o);  

    
        bv.C[i] = c;

        
        int len = bv.L[c];
        if(len > 0)
        {
            bits_write(bv.O, o_pos, o_pos + len, o);
            o_pos += len;
        }

        
        rank_count += c;  
    }

    bv.o_pos = o_pos;   

    return bv;   
}

static int bitvector_access(Bitvector bv, int i)                             
{
    uint64_t is = (uint64_t)ceil((double)i/(bv.k*bv.b));                 
    uint64_t p = bv.P[is-1];                                                 

    for(uint64_t t = (is-1)*bv.k+1; t<= (uint64_t)ceil((double)i/bv.b)-1; t++) 
    {
        p+= bv.L[bv.C[t-1]];                                                 
    }
    uint64_t littleC = bv.C[(int)ceil((double)i/bv.b)-1];                    
    uint64_t o = bits_read_simple(bv.O, p, p+bv.L[littleC]-1);
    uint64_t blk = decode(littleC, o, bv.b, bv.K);
    return bit_read(&blk, ((i-1)%bv.b)+1);                                   
}

static int bitvector_rank(Bitvector bv, int i)
{
    uint64_t is = (i - 1) / (bv.k * bv.b);
    uint64_t p = bv.P[is];
    uint64_t r = bv.R[is];

    for(uint64_t t = is * bv.k; t < (uint64_t)(i - 1) / bv.b; t++)
    {
        int c = bv.C[t];
        p += bv.L[c];
        r += c;
    }


    int      c = bv.C[(i - 1) / bv.b];
    uint64_t o = bits_read_simple(bv.O, p, p + bv.L[c] - 1);
    uint64_t blk = decode(c, o, bv.b, bv.K);
    uint64_t partial = bits_read_simple(&blk, 1, ((i - 1) % bv.b) + 1);
    return r + popcount(partial);
}

int main()
{
    uint64_t W[1] = {0};
    const char *bits = "1000001000000110000010100000101110000001";
    int n = 40;

    for(int i = 0; i < n; i++)
    {
        if(bits[i] == '1')
        {
            bit_set(W, i + 1);
        }
            
    }
        

    Bitvector bv = bitvector_build(W, n, 4, 4);

    printf("C: ");
    for(int i = 0; i < (int)bv.num_blocks; i++)
        printf("%d ", bv.C[i]);
    printf("\n");

    printf("P: ");
    for(int i = 0; i < (int)(bv.num_blocks + bv.k - 1) / bv.k; i++)
        printf("%llu ", bv.P[i]);
    printf("\n");

    printf("R: ");
    for(int i = 0; i < (int)(bv.num_blocks + bv.k - 1) / bv.k; i++) 
        printf("%llu ", bv.R[i]);
    printf("\n");

    printf("\n rank tests\n");
    printf("%d (expected 4)\n",  bitvector_rank(bv, 16));
    printf("%d (expected 8)\n",  bitvector_rank(bv, 31));
    printf("%d (expected 9)\n",  bitvector_rank(bv, 32));
    printf("%d (expected 11)\n", bitvector_rank(bv, 40));

    return 0;
}
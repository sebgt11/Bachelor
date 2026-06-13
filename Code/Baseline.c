#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#define BLOCK_SIZE      4
#define SUPERBLOCK_FREQ 4
#define CSV_PATH        "bitvectors.csv"
#define MAX_N           4096

#define WORD_SIZE 64

#define MAX_BLOCKS      (MAX_N / BLOCK_SIZE)
#define MAX_SUPERBLOCKS (MAX_N / (BLOCK_SIZE * SUPERBLOCK_FREQ))
#define MAX_BITS        MAX_N
#define MAX_OFFSET_BITS MAX_N
#define MAX_WORDS       ((MAX_OFFSET_BITS + 63) / 64)



typedef struct {
    uint8_t   C[MAX_BLOCKS];
    uint64_t  O[MAX_WORDS];
    uint64_t  P[MAX_SUPERBLOCKS];
    uint64_t  R[MAX_SUPERBLOCKS];
    uint64_t  S[MAX_BITS];
    uint64_t  K[65][65];
    uint64_t  L[65];
    uint64_t  n;
    int       b;
    int       k;
    uint64_t  num_blocks;
    uint64_t  o_pos;
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

static void CreateS(Bitvector *bv, uint64_t *W)
{
    uint64_t count = 0;
    for (uint64_t i = 1; i <= bv->n; i++)
    {
        if (bit_read(W, i) == 1)
        {
            bv->S[count] = i;
            count++;
        }
    }
}

static Bitvector bitvector_build(uint64_t *W, uint64_t n, int b, int k)
{
    Bitvector bv = {0};

    bv.n = n;
    bv.b= b;
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
            bits_write(bv.O, o_pos, len, o);
            o_pos += len;
        }

        
        rank_count += c;  
    }

    bv.o_pos = o_pos;   

    CreateS(&bv, W);

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

static int bitvector_select(Bitvector bv, int j) { return (int)bv.S[j-1]; }


// De her laver sådan set bare hel normal Brute force gennemgang af bits og laver et count
//det blir brugt til at validate resultater
static int gt_rank(const char *s, int i)
    { int r=0; for(int j=0;j<i;j++) r+=(s[j]=='1'); return r; }
static int gt_select(const char *s, int n, int j)
    { int cnt=0; for(int i=1;i<=n;i++) if(s[i-1]=='1'&&++cnt==j) return i; return -1; }
 
static double ms(struct timespec a, struct timespec b)
    { return (b.tv_sec-a.tv_sec)*1000.0 + (b.tv_nsec-a.tv_nsec)/1e6; }
 


// Number of bits used to store one class value in [0, b].
static int class_bits(int b)
{
    int bits = 0;
    for (int t = b; t > 0; t >>= 1) bits++;
    return bits ? bits : 1;
}

// Stored payload of one compressed bitvector, in bits: the class array plus
// the packed offsets. R, P and S are query-acceleration structures and are
// reported separately, so this figure is the part that approaches the entropy.
static uint64_t bitvector_payload_bits(Bitvector bv, int cb)
{
    return bv.num_blocks * (uint64_t)cb + (bv.o_pos - 1);
}

// Runs every set in the CSV through the encoder at one block size b, measuring
// space (bits per element), build time, and rank/select speed, and verifying
// every answer against a brute-force ground truth. Each bitvector is zero-padded
// up to a multiple of b*k so the same rows are used for every block size, which
// keeps the sweep a fair comparison.
static int run_block_size(const char *path, int b, int k)
{
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return 1; }

    static char     bvs[MAX_N+2];
    static uint64_t W[MAX_N/WORD_SIZE];
    char line[MAX_N+64];
    fgets(line, sizeof(line), fp);

    int cb  = class_bits(b);
    int pad = b * k;

    long long rows=0, rank_pass=0, rank_fail=0, sel_pass=0, sel_fail=0;
    double ms_build=0, ms_rank=0, ms_sel=0;
    uint64_t total_payload_bits=0, total_ones=0;
    struct timespec t0, t1;

    while (fgets(line, sizeof(line), fp))
    {
        strtok(line, ","); strtok(NULL, ",");
        char *tok = strtok(NULL, ",\r\n");
        if (!tok) continue;
        strncpy(bvs, tok, MAX_N); bvs[MAX_N] = '\0';
        int n = (int)strlen(bvs);
        if (n == 0) continue;

        int np = ((n + pad - 1) / pad) * pad;   // pad up to a whole superblock
        if (np > MAX_N) continue;

        memset(W, 0, sizeof(W));
        for (int i = 0; i < n; i++) if (bvs[i]=='1') bit_set(W, i+1);

        //Test 1: build time
        clock_gettime(CLOCK_MONOTONIC, &t0);
        Bitvector bv = bitvector_build(W, np, b, k);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        ms_build += ms(t0, t1);
        rows++;

        int ones = gt_rank(bvs, n);

        //Space: stored payload bits for this set
        total_payload_bits += bitvector_payload_bits(bv, cb);
        total_ones += ones;

        //Test 2: rank at every position
        clock_gettime(CLOCK_MONOTONIC, &t0);
        for (int i = 1; i <= n; i++) {
            int fast = bitvector_rank(bv, i);
            clock_gettime(CLOCK_MONOTONIC, &t1);
            ms_rank += ms(t0, t1);

            int slow = gt_rank(bvs, i);

            if (fast == slow) rank_pass++;
            else rank_fail++;

            clock_gettime(CLOCK_MONOTONIC, &t0);
        }

        //Test 3: select for every 1-bit
        clock_gettime(CLOCK_MONOTONIC, &t0);
        for (int j = 1; j <= ones; j++) {
            int fast = bitvector_select(bv, j);
            clock_gettime(CLOCK_MONOTONIC, &t1);
            ms_sel += ms(t0, t1);

            int slow = gt_select(bvs, n, j);

            if (fast == slow) sel_pass++;
            else sel_fail++;

            clock_gettime(CLOCK_MONOTONIC, &t0);
        }
    }
    fclose(fp);

    if (rows == 0) { fprintf(stderr, "No usable rows in %s\n", path); return 1; }

    long long total_ranks = rank_pass + rank_fail;
    long long total_sels  = sel_pass + sel_fail;
    int ok = !(rank_fail || sel_fail);

    printf("b=%-3d k=%d  rows=%-6lld  %7.3f bits/elem   build %.4f ms/row   "
           "rank %.6f ms/q   select %.6f ms/q   %s\n",
           b, k, rows,
           total_ones ? (double)total_payload_bits/total_ones : 0.0,
           ms_build/rows,
           total_ranks ? ms_rank/total_ranks : 0.0,
           total_sels ? ms_sel/total_sels : 0.0,
           ok ? "OK" : "FAIL");

    return ok ? 0 : 1;
}

// Sweeps a range of block sizes so the effect of b on space and query time can
// be read off directly. The select-helper array S caps usable bit-lengths at
// MAX_BITS, and the class array at MAX_BLOCKS = MAX_N / 4, so the smallest block
// size in the sweep must be at least 4 to stay within the static bounds.
static int run_baseline(const char *path)
{
    int sizes[] = {4, 8, 16, 32};
    int k = SUPERBLOCK_FREQ;
    int status = 0;

    printf("Block-size sweep over %s  (k = %d)\n\n", path, k);
    for (int s = 0; s < (int)(sizeof(sizes) / sizeof(sizes[0])); s++)
        status |= run_block_size(path, sizes[s], k);

    return status;
}

int main(void)
{
    return run_baseline(CSV_PATH);
}


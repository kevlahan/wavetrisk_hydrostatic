#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdint.h>
#include <unistd.h>
#include <zlib.h>

/* One stream in each MPI process. No background worker or subprocess exists.
 * Stage325 writes exact metadata and finite real-value digests.
 * The Fortran caller enforces the uncompressed128MiB bound and lifecycle. */
static gzFile stream322 = NULL;
static int durable_fd322 = -1;

int audit_compress322_open(const char *path)
{
    int fd;
    if (stream322 != NULL || durable_fd322 != -1 || path == NULL) return 1;
    fd = open(path, O_WRONLY | O_CREAT | O_EXCL, 0600);
    if (fd < 0) return 2;
    durable_fd322 = dup(fd);
    if (durable_fd322 < 0) { close(fd); return 3; }
    stream322 = gzdopen(fd, "wb1");
    if (stream322 == NULL) {
        close(fd); close(durable_fd322); durable_fd322 = -1; return 4;
    }
    if (gzbuffer(stream322, 1024U * 1024U) != 0) {
        gzclose(stream322); stream322 = NULL;
        close(durable_fd322); durable_fd322 = -1; return 5;
    }
    return 0;
}

int audit_compress322_write(const void *data, int64_t length)
{
    const unsigned char *next = (const unsigned char *)data;
    if (stream322 == NULL || length < 0 || (length > 0 && data == NULL)) return 1;
    while (length > 0) {
        unsigned int count = length > INT_MAX ? (unsigned int)INT_MAX : (unsigned int)length;
        int written = gzwrite(stream322, next, count);
        if (written <= 0 || (unsigned int)written != count) return 2;
        next += count;
        length -= count;
    }
    return 0;
}

int audit_compress322_close(void)
{
    int failed = 0;
    int result;
    if (stream322 == NULL || durable_fd322 < 0) return 1;
    if (gzclose(stream322) != Z_OK) failed = 2;
    stream322 = NULL;
    do { result = fsync(durable_fd322); } while (result < 0 && errno == EINTR);
    if (result < 0) failed = 3;
    if (close(durable_fd322) < 0) failed = 4;
    durable_fd322 = -1;
    return failed;
}

/* Stage325 qualification: SHA-256 of exact IEEE payload bytes, with finiteness
 * checked before hashing. This code never changes the numerical arrays. */
#include <string.h>
#include <math.h>

typedef struct {
    uint32_t h[8];
    uint64_t length;
    size_t used;
    unsigned char block[64];
} sha325_state;

static uint32_t rotr325(uint32_t x, unsigned n) { return (x >> n) | (x << (32U-n)); }

static void sha325_block(sha325_state *s, const unsigned char *p)
{
    static const uint32_t k[64] = {
        0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
        0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
        0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
        0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
        0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
        0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
        0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
        0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2 };
    uint32_t w[64],a,b,c,d,e,f,g,h,t1,t2;
    unsigned i;
    for(i=0;i<16;i++) w[i]=((uint32_t)p[4*i]<<24)|((uint32_t)p[4*i+1]<<16)|((uint32_t)p[4*i+2]<<8)|p[4*i+3];
    for(i=16;i<64;i++) {
        uint32_t x=w[i-15],y=w[i-2];
        w[i]=w[i-16]+(rotr325(x,7)^rotr325(x,18)^(x>>3))+w[i-7]+(rotr325(y,17)^rotr325(y,19)^(y>>10));
    }
    a=s->h[0];b=s->h[1];c=s->h[2];d=s->h[3];e=s->h[4];f=s->h[5];g=s->h[6];h=s->h[7];
    for(i=0;i<64;i++) {
        t1=h+(rotr325(e,6)^rotr325(e,11)^rotr325(e,25))+((e&f)^((~e)&g))+k[i]+w[i];
        t2=(rotr325(a,2)^rotr325(a,13)^rotr325(a,22))+((a&b)^(a&c)^(b&c));
        h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
    }
    s->h[0]+=a;s->h[1]+=b;s->h[2]+=c;s->h[3]+=d;s->h[4]+=e;s->h[5]+=f;s->h[6]+=g;s->h[7]+=h;
}

static void sha325_update(sha325_state *s,const unsigned char *p,size_t n)
{
    while(n) {
        size_t part=64-s->used;
        if(part>n) part=n;
        memcpy(s->block+s->used,p,part);s->used+=part;p+=part;n-=part;s->length+=part;
        if(s->used==64) { sha325_block(s,s->block);s->used=0; }
    }
}

/* Exported separately so the gate can check known vectors without a stream. */
int audit_sha325(const void *data,int64_t length,unsigned char out[32])
{
    sha325_state s={{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19},0,0,{0}};
    unsigned char pad[128]={0x80},bits[8];
    uint64_t nbits;
    unsigned i;
    size_t padding;
    if(length<0 || (length && data==NULL) || out==NULL || (uint64_t)length>UINT64_MAX/8) return 1;
    sha325_update(&s,(const unsigned char *)data,(size_t)length);
    nbits=s.length*8;padding=s.used<56 ? 56-s.used : 120-s.used;
    for(i=0;i<8;i++) bits[7-i]=(unsigned char)(nbits>>(8*i));
    sha325_update(&s,pad,padding);sha325_update(&s,bits,8);
    for(i=0;i<8;i++) { out[4*i]=(unsigned char)(s.h[i]>>24);out[4*i+1]=(unsigned char)(s.h[i]>>16);out[4*i+2]=(unsigned char)(s.h[i]>>8);out[4*i+3]=(unsigned char)s.h[i]; }
    return 0;
}

int audit_digest325_write(const void *data,int64_t length)
{
    unsigned char digest[32];
    const unsigned char *p=(const unsigned char *)data;
    int64_t i;
    if(sizeof(double)!=8 || length<0 || length%8 || (length && data==NULL)) return 1;
    for(i=0;i<length;i+=8) {
        double x;
        memcpy(&x,p+i,8);
        if(!isfinite(x)) return 2;
    }
    if(audit_sha325(data,length,digest)) return 3;
    return audit_compress322_write(digest,32);
}

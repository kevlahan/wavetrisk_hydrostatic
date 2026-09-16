#define _GNU_SOURCE
/* Identical observer clock on macOS and Linux; no adjustable wall-clock epoch. */
#define _POSIX_C_SOURCE 200809L
#include <time.h>

double workflow_clock199(void)
{
    static struct timespec origin;
    static int initialized;
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0)
        return -1.0;
    if (!initialized) {
        origin = now;
        initialized = 1;
    }
    return (double)(now.tv_sec - origin.tv_sec)
         + (double)(now.tv_nsec - origin.tv_nsec) * 1e-9;
}

#ifdef __linux__
#include <sched.h>
int workflow_affinity200(void) {
    cpu_set_t cpus;
    if (sched_getaffinity(0, sizeof(cpus), &cpus) || CPU_COUNT(&cpus) != 1) return -1;
    for (int i=0; i<CPU_SETSIZE; ++i) if (CPU_ISSET(i,&cpus)) return i;
    return -1;
}
#else
int workflow_affinity200(void) { return -2; }
#endif

/* Stage201: thread-local PMU counts and absolute CLOCK_MONOTONIC intervals.
 * Four simultaneous user-space events, no multiplex scaling, no hot-loop hooks.
 * No allocation or output within measured timesteps. See event provenance doc.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <inttypes.h>
#include <unistd.h>
#ifdef __linux__
#include <linux/perf_event.h>
#include <sys/syscall.h>
#include <sys/ioctl.h>
#include <cpuid.h>
#endif
#define MAX201 16384
#define EVENTS201 4
struct snapshot201 { uint64_t nr, enabled, running, value[EVENTS201]; };
struct interval201 { int step, context; uint64_t start, end, cpu, enabled, running, value[EVENTS201]; };
static struct interval201 rows201[MAX201];
static struct snapshot201 prev201;
static uint64_t start201, cpu201;
static int mode201, fd201[EVENTS201]={-1,-1,-1,-1}, n201;
static const char *names201[3][EVENTS201]={
    {"none","none","none","none"},
    {"cycles","instructions","ex_no_retire.not_complete","ex_no_retire.load_not_complete"},
    {"cycles","instructions","l2_cache_req_stat.ls_rd_blk_c","ls_dmnd_fills_from_sys.dram_io_near_or_far"}};
static void fail201(const char *where) { fprintf(stderr,"Stage201 %s: %s\n",where,strerror(errno)); exit(201); }
static uint64_t ns201(clockid_t clock) {
    struct timespec t; if(clock_gettime(clock,&t)) fail201("clock_gettime");
    return (uint64_t)t.tv_sec*1000000000ULL+(uint64_t)t.tv_nsec;
}
static struct snapshot201 read201(void) {
    struct snapshot201 s={0};
#ifdef __linux__
    if(mode201 && (read(fd201[0],&s,sizeof(s))!=(ssize_t)sizeof(s) || s.nr!=EVENTS201)) fail201("PMU group read");
#endif
    return s;
}
void workflow_hw_init201(void) {
    const char *s=getenv("WAVETRISK_PMU201");
    if(s && strcmp(s,"0") && strcmp(s,"1") && strcmp(s,"2")) {errno=EINVAL;fail201("mode");}
    mode201=s?atoi(s):0;
    if(!mode201) return;
#ifdef __linux__
    unsigned int a,b,c,d, family,model; char vendor[13];
    __cpuid(0,a,b,c,d);memcpy(vendor,&b,4);memcpy(vendor+4,&d,4);memcpy(vendor+8,&c,4);vendor[12]=0;
    __cpuid(1,a,b,c,d);family=((a>>8)&15)+((a>>20)&255);model=((a>>4)&15)|((a>>12)&240);
    if(strcmp(vendor,"AuthenticAMD") || family!=25 || model!=160) {errno=ENOTSUP;fail201("requires validated EPYC 9754 family25/model160");}
    uint64_t raw[2]={mode201==1?0x02d6:0x0864,mode201==1?0xa2d6:0x4843};
    for(int i=0;i<EVENTS201;i++) {
        struct perf_event_attr attr;memset(&attr,0,sizeof(attr));attr.size=sizeof(attr);
        attr.type=i<2?PERF_TYPE_HARDWARE:PERF_TYPE_RAW;
        attr.config=i==0?PERF_COUNT_HW_CPU_CYCLES:i==1?PERF_COUNT_HW_INSTRUCTIONS:raw[i-2];
        attr.disabled=i==0;attr.exclude_kernel=1;attr.exclude_hv=1;
        attr.read_format=PERF_FORMAT_GROUP|PERF_FORMAT_TOTAL_TIME_ENABLED|PERF_FORMAT_TOTAL_TIME_RUNNING;
        fd201[i]=(int)syscall(__NR_perf_event_open,&attr,0,-1,i?fd201[0]:-1,PERF_FLAG_FD_CLOEXEC);
        if(fd201[i]<0) fail201(names201[mode201][i]);
    }
    if(ioctl(fd201[0],PERF_EVENT_IOC_RESET,PERF_IOC_FLAG_GROUP) || ioctl(fd201[0],PERF_EVENT_IOC_ENABLE,PERF_IOC_FLAG_GROUP)) fail201("enable PMU group");
#else
    errno=ENOTSUP;fail201("PMU requires Linux AMD cluster; mode0 is portable");
#endif
}
void workflow_hw_begin201(void) {
    prev201=read201();cpu201=ns201(CLOCK_THREAD_CPUTIME_ID);start201=ns201(CLOCK_MONOTONIC);
}
void workflow_hw_tick201(int step,int context) {
    if(n201>=MAX201 || step<1 || context<0 || context>8) {errno=EOVERFLOW;fail201("interval overflow/invalid context");}
    struct snapshot201 now=read201();uint64_t cpu=ns201(CLOCK_THREAD_CPUTIME_ID),end=ns201(CLOCK_MONOTONIC);
    struct interval201 *r=&rows201[n201++];r->step=step;r->context=context;r->start=start201;r->end=end;r->cpu=cpu-cpu201;
    r->enabled=now.enabled-prev201.enabled;r->running=now.running-prev201.running;
    for(int i=0;i<EVENTS201;i++) {if(now.value[i]<prev201.value[i]) {errno=ERANGE;fail201("counter reversal");} r->value[i]=now.value[i]-prev201.value[i];}
    start201=end;cpu201=cpu;prev201=now;
}
void workflow_hw_finish201(int rank) {
    char name[80];snprintf(name,sizeof(name),"hardware201-rank-%d.txt",rank);
    FILE *f=fopen(name,"wx");if(!f)fail201("open hardware report");
    fprintf(f,"schema 1 %d %d %d\n",mode201,rank,(int)getpid());
    fprintf(f,"events %s %s %s %s\n",names201[mode201][0],names201[mode201][1],names201[mode201][2],names201[mode201][3]);
    for(int j=0;j<n201;j++) {
        struct interval201 *r=&rows201[j];
        fprintf(f,"interval %d %d %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64,r->step,r->context,r->start,r->end,r->cpu,r->enabled,r->running);
        for(int i=0;i<EVENTS201;i++)fprintf(f," %"PRIu64,r->value[i]);
        fputc('\n',f);
    }
    if(fclose(f))fail201("close report");
    for(int i=0;i<EVENTS201;i++)if(fd201[i]>=0)close(fd201[i]);
}
#ifdef STANDALONE201
int main(void) {
    workflow_hw_init201();workflow_hw_begin201();
    volatile uint64_t sum=1;
    for(int i=0;i<100000000;i++)sum=sum*1664525+1013904223;
    workflow_hw_tick201(1,0);workflow_hw_finish201(0);
    if(mode201 && (!rows201[0].value[0] || !rows201[0].value[1] || !rows201[0].enabled || (double)rows201[0].running/rows201[0].enabled<0.95)) {
        errno=EAGAIN;fail201("PMU unscheduled/under 95 percent running");
    }
    (void)sum; return 0;
}
#endif

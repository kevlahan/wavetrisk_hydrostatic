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

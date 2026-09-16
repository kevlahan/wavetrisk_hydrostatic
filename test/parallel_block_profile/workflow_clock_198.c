/* Identical observer clock on macOS and Linux; no adjustable wall-clock epoch. */
#define _POSIX_C_SOURCE 200809L
#include <time.h>

double workflow_clock198(void)
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

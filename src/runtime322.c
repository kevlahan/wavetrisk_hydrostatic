#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#ifdef __linux__
#include <sched.h>
#include <dirent.h>
#include <sys/syscall.h>
#endif

/* Qualification-only snapshots; bare runs take only pre/post-window snapshots.
 * All output is deferred until the existing final checkpoint has completed. */
#define OBS_MAX 7
#define THREAD_MAX 64
#define CPU_MAX 1024
#ifdef __linux__
_Static_assert(CPU_SETSIZE == CPU_MAX, "Stage322 CPU affinity bound mismatch");
#endif
struct thread322 { long tid; int count; uint64_t mask[CPU_MAX/64]; };
struct observation322 {
  int stage, cpu, policy, count;
  long main_tid;
  struct thread322 threads[THREAD_MAX];
};
static struct observation322 observations322[OBS_MAX];
static int count322;
static char host322[256], namespace322[256], children322[256], boot322[256];

int64_t monotonic322(void) {
  struct timespec ts;
  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) return -1;
  return (int64_t)ts.tv_sec * INT64_C(1000000000) + ts.tv_nsec;
}

int observe322(int stage) {
#ifdef __linux__
  if (count322 >= OBS_MAX || stage < 0 || stage > 6) return 1;
  struct observation322 *o = &observations322[count322];
  memset(o, 0, sizeof *o);
  o->stage = stage; o->main_tid = syscall(SYS_gettid); o->cpu = sched_getcpu();
  if (o->main_tid != (long)getpid() || o->cpu < 0) return 2;
  if (syscall(SYS_get_mempolicy, &o->policy, NULL, 0, NULL, 0) != 0) return 3;
  DIR *directory = opendir("/proc/self/task");
  if (!directory) return 4;
  struct dirent *entry;
  while ((entry = readdir(directory)) != NULL) {
    char *end; long tid = strtol(entry->d_name, &end, 10);
    if (tid <= 0 || *end != '\0') continue;
    if (o->count == THREAD_MAX) { closedir(directory); return 5; }
    struct thread322 *t = &o->threads[o->count];
    t->tid = tid;
    cpu_set_t cpus;
    if (sched_getaffinity((pid_t)tid, sizeof cpus, &cpus) != 0) {
      closedir(directory); return 6;
    }
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
      if (CPU_ISSET(cpu, &cpus)) {
        if (t->count == CPU_MAX) { closedir(directory); return 7; }
        t->mask[cpu/64] |= UINT64_C(1) << (cpu%64);
        ++t->count;
      }
    }
    if (t->count == 0) { closedir(directory); return 8; }
    ++o->count;
  }
  closedir(directory);
  const char *wanted = getenv("WAVETRISK_EXPECT_CPU200");
  if (!wanted) return 9;
  char *end; long expected = strtol(wanted, &end, 10);
  if (end == wanted || *end != '\0' || expected < 0 || expected >= CPU_MAX) return 10;
  int found = 0;
  for (int i = 0; i < o->count; ++i) {
    if (o->threads[i].tid == o->main_tid) {
      if (o->threads[i].count != 1 || (o->threads[i].mask[expected/64] & (UINT64_C(1) << (expected%64))) == 0) return 11;
      found = 1;
    }
  }
  if (!found || o->cpu != expected || o->policy != 4) return 12;
  if (count322 == 0) {
    if (gethostname(host322, sizeof host322) != 0) return 13;
    ssize_t n = readlink("/proc/self/ns/time", namespace322, sizeof namespace322 - 1);
    if (n < 0) return 14;
    namespace322[n] = '\0';
    n = readlink("/proc/self/ns/time_for_children", children322, sizeof children322 - 1);
    if (n < 0) return 15;
    children322[n] = '\0';
    FILE *f = fopen("/proc/sys/kernel/random/boot_id", "r");
    if (!f) return 16;
    if (!fgets(boot322, sizeof boot322, f)) { fclose(f); return 17; }
    fclose(f); boot322[strcspn(boot322, "\r\n")] = '\0';
  }
  ++count322;
  return 0;
#else
  (void)stage;
  return 18;
#endif
}

int flush_runtime322(int rank, int64_t begin, int64_t end) {
  char filename[128];
  if (count322 == 0) return 1;
  if (snprintf(filename, sizeof filename, "runtime322-rank-%d.json", rank) < 0) return 2;
  FILE *f = fopen(filename, "wx");
  if (!f) return 3;
  /* Hostnames and Linux namespace/boot identifiers cannot contain JSON escapes. */
  fprintf(f, "{\"rank\":%d,\"pid\":%ld,\"host\":\"%s\",\"clock_namespace\":\"%s\","
    "\"time_namespace_for_children\":\"%s\",\"boot_id\":\"%s\",\"begin_ns\":%lld,"
    "\"end_ns\":%lld,\"observations\":[", rank, (long)getpid(), host322,
    namespace322, children322, boot322, (long long)begin, (long long)end);
  for (int i = 0; i < count322; ++i) {
    struct observation322 *o = &observations322[i];
    fprintf(f, "%s{\"stage\":%d,\"main_tid\":%ld,\"cpu\":%d,\"memory_policy\":%d,\"threads\":[",
      i ? "," : "", o->stage, o->main_tid, o->cpu, o->policy);
    for (int j = 0; j < o->count; ++j) {
      struct thread322 *t = &o->threads[j];
      fprintf(f, "%s{\"tid\":%ld,\"mask\":[", j ? "," : "", t->tid);
      int printed = 0;
      for (int k = 0; k < CPU_MAX; ++k)
        if (t->mask[k/64] & (UINT64_C(1) << (k%64)))
          fprintf(f, "%s%d", printed++ ? "," : "", k);
      fprintf(f, "]}");
    }
    fprintf(f, "]}");
  }
  fprintf(f, "]}\n");
  int failed = ferror(f);
  if (fclose(f) != 0) failed = 1;
  return failed ? 4 : 0;
}

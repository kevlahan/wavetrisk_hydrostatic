#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdint.h>
#include <unistd.h>
#include <zlib.h>

/* One stream in each MPI process. No background worker or subprocess exists.
 * The compressed stream contains the original Stage316 bytes without changes.
 * The caller enforces its uncompressed512MiB bound and lifecycle. */
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

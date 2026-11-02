#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

static uint64_t parse_u64(char **cursor) {
    while (**cursor == ' ' || **cursor == '\t') (*cursor)++;
    uint64_t value = 0;
    while (**cursor >= '0' && **cursor <= '9') {
        value = value * 10 + (uint64_t)(**cursor - '0');
        (*cursor)++;
    }
    return value;
}

static void die(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
    if (argc != 5) {
        fprintf(stderr, "usage: %s matrix.mtx.gz target_rows.tsv library_u64.bin target_counts_u32.bin\n", argv[0]);
        return EXIT_FAILURE;
    }
    const char *matrix_path = argv[1];
    const char *targets_path = argv[2];
    const char *library_path = argv[3];
    const char *counts_path = argv[4];

    gzFile input = gzopen(matrix_path, "rb");
    if (!input) die("cannot open gzip matrix");
    gzbuffer(input, 8 * 1024 * 1024);
    char line[256];
    uint64_t n_rows = 0, n_cols = 0, n_entries = 0;
    while (gzgets(input, line, sizeof(line))) {
        if (line[0] == '%') continue;
        char *p = line;
        n_rows = parse_u64(&p);
        n_cols = parse_u64(&p);
        n_entries = parse_u64(&p);
        break;
    }
    if (!n_rows || !n_cols || !n_entries) die("invalid Matrix Market header");

    int32_t *row_to_target = malloc((n_rows + 1) * sizeof(int32_t));
    if (!row_to_target) die("cannot allocate row map");
    for (uint64_t i = 0; i <= n_rows; i++) row_to_target[i] = -1;

    FILE *targets = fopen(targets_path, "r");
    if (!targets) die("cannot open target row file");
    uint64_t row_id = 0;
    int32_t target_id = 0;
    while (fgets(line, sizeof(line), targets)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char *p = line;
        row_id = parse_u64(&p);
        if (row_id == 0 || row_id > n_rows) die("target row outside matrix dimensions");
        row_to_target[row_id] = target_id++;
    }
    fclose(targets);
    if (target_id == 0) die("no target rows provided");

    uint64_t *library = calloc(n_cols, sizeof(uint64_t));
    uint32_t *selected = calloc((uint64_t)target_id * n_cols, sizeof(uint32_t));
    if (!library || !selected) die("cannot allocate output arrays");

    uint64_t processed = 0;
    while (gzgets(input, line, sizeof(line))) {
        char *p = line;
        uint64_t row = parse_u64(&p);
        uint64_t col = parse_u64(&p);
        uint64_t value = parse_u64(&p);
        if (row == 0 || row > n_rows || col == 0 || col > n_cols) die("matrix coordinate outside dimensions");
        library[col - 1] += value;
        int32_t index = row_to_target[row];
        if (index >= 0) {
            uint64_t offset = (uint64_t)index * n_cols + (col - 1);
            uint64_t updated = (uint64_t)selected[offset] + value;
            if (updated > UINT32_MAX) die("selected count overflow");
            selected[offset] = (uint32_t)updated;
        }
        processed++;
        if (processed % 50000000ULL == 0) {
            fprintf(stderr, "processed=%" PRIu64 "/%" PRIu64 "\n", processed, n_entries);
        }
    }
    gzclose(input);
    if (processed != n_entries) {
        fprintf(stderr, "entry count mismatch: processed=%" PRIu64 " expected=%" PRIu64 "\n", processed, n_entries);
        return EXIT_FAILURE;
    }

    FILE *library_out = fopen(library_path, "wb");
    FILE *counts_out = fopen(counts_path, "wb");
    if (!library_out || !counts_out) die("cannot open output binaries");
    if (fwrite(library, sizeof(uint64_t), n_cols, library_out) != n_cols) die("failed writing library output");
    uint64_t selected_n = (uint64_t)target_id * n_cols;
    if (fwrite(selected, sizeof(uint32_t), selected_n, counts_out) != selected_n) die("failed writing selected counts output");
    fclose(library_out);
    fclose(counts_out);
    fprintf(stderr, "complete rows=%" PRIu64 " cols=%" PRIu64 " entries=%" PRIu64 " targets=%d\n", n_rows, n_cols, n_entries, target_id);

    free(row_to_target);
    free(library);
    free(selected);
    return EXIT_SUCCESS;
}

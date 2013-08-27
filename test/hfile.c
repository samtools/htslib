#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <sys/stat.h>

#include "hfile.h"

void fail(const char *format, ...)
{
    int err = errno;
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, ": %s\n", strerror(err));
    exit(EXIT_FAILURE);
}

char *slurp(const char *filename)
{
    char *text;
    struct stat sbuf;
    size_t filesize;
    FILE *f = fopen(filename, "r");
    if (f == NULL) fail("fopen(\"%s\", \"r\")", filename);
    if (fstat(fileno(f), &sbuf) != 0) fail("fstat(\"%s\")", filename);
    filesize = sbuf.st_size;

    text = (char *) malloc(filesize + 1);
    if (text == NULL) fail("malloc(text)");

    if (fread(text, 1, filesize, f) != filesize) fail("fread");
    fclose(f);

    text[filesize] = '\0';
    return text;
}

hFILE *fin = NULL;
hFILE *fout = NULL;

void reopen(const char *infname, const char *outfname)
{
    if (fin) { if (hclose(fin) != 0) fail("hclose(input)"); }
    if (fout) { if (hclose(fout) != 0) fail("hclose(output)"); }

    fin = hopen(infname, "r");
    if (fin == NULL) fail("hopen(\"%s\")", infname);

    fout = hopen(outfname, "w");
    if (fout == NULL) fail("hopen(\"%s\")", outfname);
}

int main() {
    static const int size[] = { 1, 13, 403, 999, 30000 };

    char buffer[40000];
    char *original;
    int c, i;
    ssize_t n;

    reopen("vcf.c", "test/hfile1.tmp");
    while ((c = hgetc(fin)) != EOF) {
        if (hputc(c, fout) == EOF) fail("hputc");
    }
    if (herrno(fin)) { errno = herrno(fin); fail("hgetc"); }

    reopen("test/hfile1.tmp", "test/hfile2.tmp");
    if (hpeek(fin, buffer, 50) < 0) fail("hpeek");
    while ((n = hread(fin, buffer, 17)) > 0) {
        if (hwrite(fout, buffer, n) != n) fail("hwrite");
    }
    if (n < 0) fail("hread");

    reopen("test/hfile2.tmp", "test/hfile3.tmp");
    while ((n = hread(fin, buffer, sizeof buffer)) > 0) {
        if (hwrite(fout, buffer, n) != n) fail("hwrite");
        if (hpeek(fin, buffer, 700) < 0) fail("hpeek");
    }
    if (n < 0) fail("hread");

    reopen("test/hfile3.tmp", "test/hfile4.tmp");
    i = 0;
    while ((n = hread(fin, buffer, size[i++ % 5])) > 0) {
        buffer[n] = '\0';
        if (hputs(buffer, fout) == EOF) fail("hputs");
        if ((n = hpeek(fin, buffer, size[(i+3) % 5])) < 0) fail("hpeek");
    }
    if (n < 0) fail("hread");

    if (hflush(fout) == EOF) fail("hflush");

    original = slurp("vcf.c");
    for (i = 1; i <= 4; i++) {
        char *text;
        sprintf(buffer, "test/hfile%d.tmp", i);
        text = slurp(buffer);
        if (strcmp(original, text) != 0) {
            fprintf(stderr, "%s differs from vcf.c\n", buffer);
            return EXIT_FAILURE;
        }
        free(text);
    }
    free(original);

    if (hclose(fin) != 0) fail("hclose(input)");
    if (hclose(fout) != 0) fail("hclose(output)");

    return EXIT_SUCCESS;
}

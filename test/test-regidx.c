/*  test/test-regidx.c -- Regions index test harness.

    gcc -g -Wall -O0 -I. -I../htslib/ -L../htslib regidx.c -o test-regidx test-regidx.c -lhts

    Copyright (C) 2014,2016,2018, 2020, 2023 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

#include <config.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#include "../htslib/kstring.h"
#include "../htslib/regidx.h"
#include "../htslib/hts_defs.h"
#include "../textutils_internal.h"

static int verbose = 0;

HTS_FORMAT(HTS_PRINTF_FMT, 1, 2)
static void debug(const char *format, ...)
{
    if ( verbose<2 ) return;
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
}

HTS_FORMAT(HTS_PRINTF_FMT, 1, 2)
static void info(const char *format, ...)
{
    if ( verbose<1 ) return;
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
}

HTS_NORETURN HTS_FORMAT(HTS_PRINTF_FMT, 1, 2)
static void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

int custom_parse(const char *line, char **chr_beg, char **chr_end, hts_pos_t *beg, hts_pos_t *end, void *payload, void *usr)
{
    // Use the standard parser for CHROM,FROM,TO
    int i, ret = regidx_parse_tab(line,chr_beg,chr_end,beg,end,NULL,NULL);
    if ( ret!=0 ) return ret;

    // Skip the fields that were parsed above
    char *ss = (char*) line;
    while ( *ss && isspace_c(*ss) ) ss++;
    for (i=0; i<3; i++)
    {
        while ( *ss && !isspace_c(*ss) ) ss++;
        if ( !*ss ) return -2;  // wrong number of fields
        while ( *ss && isspace_c(*ss) ) ss++;
    }
    if ( !*ss ) return -2;

    // Parse the payload
    char *se = ss;
    while ( *se && !isspace_c(*se) ) se++;
    char **dat = (char**) payload;
    *dat = (char*) malloc(se-ss+1);
    memcpy(*dat,ss,se-ss+1);
    (*dat)[se-ss] = 0;
    return 0;
}
void custom_free(void *payload)
{
    char **dat = (char**)payload;
    free(*dat);
}

void test_sequential_access(void)
{
    // Init index with no file name, we will insert the regions manually
    regidx_t *idx = regidx_init(NULL,custom_parse,custom_free,sizeof(char*),NULL);
    if ( !idx ) error("init failed\n");

    // Insert regions
    kstring_t str = {0,0,0};
    int i, n = 10;
    for (i=0; i<n; i++)
    {
        int beg = 10*(i+1);
        str.l = 0;
        ksprintf(&str,"1\t%d\t%d\t%d",beg,beg,beg);
        if ( regidx_insert(idx,str.s)!=0 ) error("insert failed: %s\n",str.s);
    }

    // Test
    regitr_t *itr = regitr_init(idx);
    i = 0;
    while ( regitr_loop(itr) )
    {
        if ( itr->beg!=itr->end || itr->beg+1!=10*(i+1) ) error("listing failed, expected %d, found %"PRIhts_pos"\n",10*(i+1),itr->beg+1);
        str.l = 0;
        ksprintf(&str,"%"PRIhts_pos, itr->beg+1);
        if ( strcmp(regitr_payload(itr,char*),str.s) ) error("listing failed, expected payload \"%s\", found \"%s\"\n",str.s,regitr_payload(itr,char*));
        i++;
    }
    if ( i!=n ) error("Expected %d regions, listed %d\n", n,i);
    debug("ok: listed %d regions\n", n);

    // Clean up
    regitr_destroy(itr);
    regidx_destroy(idx);
    free(str.s);
}

void test_custom_payload(void)
{
    // Init index with no file name, we will insert the regions manually
    regidx_t *idx = regidx_init(NULL,custom_parse,custom_free,sizeof(char*),NULL);
    if ( !idx ) error("init failed\n");

    // Insert regions
    char *line;
    line = "1 10000000 10000000 1:10000000-10000000"; if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);
    line = "1 20000000 20000001 1:20000000-20000001"; if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);
    line = "1 20000002 20000002 1:20000002-20000002"; if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);
    line = "1 30000000 30000000 1:30000000-30000000"; if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);
    line = "1 8000000000 8000000000 1:8000000000-8000000000"; if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);

    // Test
    regitr_t *itr = regitr_init(idx);
    hts_pos_t from, to;

    from = to = 10000000;
    if ( !regidx_overlap(idx,"1",from-1,to-1,itr) ) error("query failed: 1:%"PRIhts_pos"-%"PRIhts_pos"\n",from,to);
    if ( strcmp("1:10000000-10000000",regitr_payload(itr,char*)) ) error("query failed: 1:%"PRIhts_pos"-%"PRIhts_pos" vs %s\n", from,to,regitr_payload(itr,char*));
    if ( !regidx_overlap(idx,"1",from-2,to-1,itr) ) error("query failed: 1:%"PRIhts_pos"-%"PRIhts_pos"\n",from-1,to);
    if ( !regidx_overlap(idx,"1",from-2,to+3,itr) ) error("query failed: 1:%"PRIhts_pos"-%"PRIhts_pos"\n",from-1,to+2);
    if ( regidx_overlap(idx,"1",from-2,to-2,itr) ) error("query failed: 1:%"PRIhts_pos"-%"PRIhts_pos"\n",from-1,to-1);

    from = to = 20000000;
    if ( !regidx_overlap(idx,"1",from-1,to-1,itr) ) error("query failed: 1:%"PRIhts_pos"-%"PRIhts_pos"\n",from,to);

    from = to = 20000002;
    if ( !regidx_overlap(idx,"1",from-1,to-1,itr) ) error("query failed: 1:%"PRIhts_pos"-%"PRIhts_pos"\n",from,to);

    from = to = 30000000;
    if ( !regidx_overlap(idx,"1",from-1,to-1,itr) ) error("query failed: 1:%"PRIhts_pos"-%"PRIhts_pos"\n",from,to);

    from = to = 8000000000;
    if ( !regidx_overlap(idx,"1",from-1,to-1,itr) ) error("query failed: 1:%"PRIhts_pos"-%"PRIhts_pos"\n",from,to);

    // This shouldn't bring anything back
    from &= 0xffffffffU;
    to   &= 0xffffffffU;
    if ( regidx_overlap(idx,"1",from-1,to-1,itr) ) error("query should not succeed: 1:%"PRIhts_pos"-%"PRIhts_pos"\n",from,to);

    // Clean up
    regitr_destroy(itr);
    regidx_destroy(idx);
}

void get_random_region(uint32_t min, uint32_t max, uint32_t *beg, uint32_t *end)
{
    uint64_t b = rand(), e = rand();
    *beg = min + (b * (max-min)) / RAND_MAX;
    *end = *beg + (e * (max-*beg)) / RAND_MAX;
}

void test_random(int nregs, uint32_t min, uint32_t max)
{
    min--;
    max--;

    // Init index with no file name, we will insert the regions manually
    regidx_t *idx = regidx_init(NULL,custom_parse,custom_free,sizeof(char*),NULL);
    if ( !idx ) error("init failed\n");

    // Test region
    uint32_t beg,end;
    get_random_region(min,max,&beg,&end);

    // Insert regions
    int i, nexp = 0;
    kstring_t str = {0,0,0};
    for (i=0; i<nregs; i++)
    {
        uint32_t b,e;
        get_random_region(min,max,&b,&e);
        str.l = 0;
        ksprintf(&str,"1\t%"PRIu32"\t%"PRIu32"\t1:%"PRIu32"-%"PRIu32"",b+1,e+1,b+1,e+1);
        if ( regidx_insert(idx,str.s)!=0 ) error("insert failed: %s\n", str.s);
        if ( e>=beg && b<=end ) nexp++;
    }

    // Test
    regitr_t *itr = regitr_init(idx);
    int nhit = 0, ret = regidx_overlap(idx,"1",beg,end,itr);
    if ( nexp && !ret ) error("query failed, expected %d overlap(s), found none: %d-%d\n", nexp,beg+1,end+1);
    if ( !nexp && ret ) error("query failed, expected no overlaps, found some: %d-%d\n", beg+1,end+1);
    while ( ret && regitr_overlap(itr) )
    {
        str.l = 0;
        ksprintf(&str,"1:%"PRIhts_pos"-%"PRIhts_pos"",itr->beg+1,itr->end+1);
        if ( strcmp(str.s,regitr_payload(itr,char*)) )
            error("query failed, incorrect payload: %s vs %s (%d-%d)\n",str.s,regitr_payload(itr,char*),beg+1,end+1);
        if ( itr->beg > end || itr->end < beg )
            error("query failed, incorrect hit: %d-%d vs %"PRIhts_pos"-%"PRIhts_pos", payload %s\n", beg+1,end+1,itr->beg+1,itr->end+1,regitr_payload(itr,char*));
        nhit++;
    }
    if ( nexp!=nhit ) error("query failed, expected %d overlap(s), found %d: %d-%d\n",nexp,nhit,beg+1,end+1);
    debug("ok: found %d overlaps\n", nexp);

    // Clean up
    regitr_destroy(itr);
    regidx_destroy(idx);
    free(str.s);
}
void test_explicit(char *tgt, char *qry, char *exp)
{
    regidx_t *idx = regidx_init(NULL,regidx_parse_reg,NULL,0,NULL);

    char *beg = tgt, *end, *exp_ori = exp;
    kstring_t str = {0,0,0};
    while ( *beg )
    {
        end = tgt;
        while ( *end && *end!=';' ) end++;
        str.l = 0;
        kputsn(beg, end-beg, &str);
        debug("insert: %s\n", str.s);
        if ( regidx_insert(idx,str.s)!=0 ) error("insert failed: %s\n", str.s);
        beg = *end ? end + 1 : end;
    }

    beg = qry;
    while ( *beg )
    {
        end = qry;
        while ( *end && *end!=';' ) end++;
        str.l = 0;
        kputsn(beg, end-beg, &str);
        beg = *end ? end + 1 : end;

        char *chr_beg, *chr_end;
        hts_pos_t reg_beg, reg_end;
        if ( regidx_parse_reg(str.s, &chr_beg, &chr_end, &reg_beg, &reg_end, NULL, NULL)!=0 ) error("could not parse: %s in %s\n", str.s, qry);
        chr_end[1] = 0;
        int hit = regidx_overlap(idx,chr_beg,reg_beg,reg_end,NULL);
        if ( *exp=='1' )
        {
            if ( !hit )
            {
                error("query failed, there should be a hit .. %s:%"PRIhts_pos"-%"PRIhts_pos"\n",chr_beg, reg_beg+1, reg_end+1);
            }
            else
            {
                debug("ok: overlap found for %s:%"PRIhts_pos"-%"PRIhts_pos"\n",chr_beg,reg_beg+1,reg_end+1);
            }
        }
        else if ( *exp=='0' )
        {
            if ( hit )
            {
                error("query failed, there should be no hit .. %s:%"PRIhts_pos"-%"PRIhts_pos"\n",chr_beg,reg_beg+1,reg_end+1);
            }
            else
            {
                debug("ok: no overlap found for %s:%"PRIhts_pos"-%"PRIhts_pos"\n",chr_beg,reg_beg+1,reg_end+1);
            }
        }
        else error("could not parse: %s\n", exp_ori);
        exp++;
    }

    free(str.s);
    regidx_destroy(idx);
}

void create_line_bed(char *line, size_t size, char *chr, int start, int end)
{
    snprintf(line,size,"%s\t%d\t%d\n",chr,start-1,end);
}
void create_line_tab(char *line, size_t size, char *chr, int start, int end)
{
    snprintf(line,size,"%s\t%d\t%d\n",chr,start,end);
}
void create_line_reg(char *line, size_t size, char *chr, int start, int end)
{
    snprintf(line,size,"%s:%d-%d\n",chr,start,end);
}

typedef void (*set_line_f)(char *line, size_t size, char *chr, int start, int end);

void test(set_line_f set_line, regidx_parse_f parse)
{
    regidx_t *idx = regidx_init(NULL,parse,NULL,0,NULL);
    if ( !idx ) error("init failed\n");

    char line[250], *chr = "1";
    int i, n = 10, start, end, nhit;
    for (i=1; i<n; i++)
    {
        start = end = 10*i;
        set_line(line,sizeof(line),chr,start,end);
        debug("insert: %s", line);
        if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);

        start = end = 10*i + 1;
        set_line(line,sizeof(line),chr,start,end);
        debug("insert: %s", line);
        if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);

        start = 20000*i; end = start + 2000;
        set_line(line,sizeof(line),chr,start,end);
        debug("insert: %s", line);
        if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);
    }

    regitr_t *itr = regitr_init(idx);
    for (i=1; i<n; i++)
    {
        // no hit
        start = end = 10*i - 1;
        if ( regidx_overlap(idx,chr,start-1,end-1,itr) ) error("query failed, there should be no hit: %s:%d-%d\n",chr,start,end);
        debug("ok: no overlap found for %s:%d-%d\n",chr,start,end);


        // one hit
        start = end = 10*i;
        if ( !regidx_overlap(idx,chr,start-1,end-1,itr) ) error("query failed, there should be a hit: %s:%d-%d\n",chr,start,end);
        debug("ok: overlap(s) found for %s:%d-%d\n",chr,start,end);
        nhit = 0;
        while ( regitr_overlap(itr) )
        {
            if ( itr->beg > end-1 || itr->end < start-1 ) error("query failed, incorrect region: %"PRIhts_pos"-%"PRIhts_pos" for %d-%d\n",itr->beg+1,itr->end+1,start,end);
            debug("\t %"PRIhts_pos"-%"PRIhts_pos"\n",itr->beg+1,itr->end+1);
            nhit++;
        }
        if ( nhit!=1 ) error("query failed, expected one hit, found %d: %s:%d-%d\n",nhit,chr,start,end);


        // one hit
        start = end = 10*i+1;
        if ( !regidx_overlap(idx,chr,start-1,end-1,itr) ) error("query failed, there should be a hit: %s:%d-%d\n",chr,start,end);
        debug("ok: overlap(s) found for %s:%d-%d\n",chr,start,end);
        nhit = 0;
        while ( regitr_overlap(itr) )
        {
            if ( itr->beg > end-1 || itr->end < start-1 ) error("query failed, incorrect region: %"PRIhts_pos"-%"PRIhts_pos" for %d-%d\n",itr->beg+1,itr->end+1,start,end);
            debug("\t %"PRIhts_pos"-%"PRIhts_pos"\n",itr->beg+1,itr->end+1);
            nhit++;
        }
        if ( nhit!=1 ) error("query failed, expected one hit, found %d: %s:%d-%d\n",nhit,chr,start,end);


        // two hits
        start = 10*i; end = start+1;
        if ( !regidx_overlap(idx,chr,start-1,end-1,itr) ) error("query failed, there should be a hit: %s:%d-%d\n",chr,start,end);
        debug("ok: overlap(s) found for %s:%d-%d\n",chr,start,end);
        nhit = 0;
        while ( regitr_overlap(itr) )
        {
            if ( itr->beg > end-1 || itr->end < start-1 ) error("query failed, incorrect region: %"PRIhts_pos"-%"PRIhts_pos" for %d-%d\n",itr->beg+1,itr->end+1,start,end);
            debug("\t %"PRIhts_pos"-%"PRIhts_pos"\n",itr->beg+1,itr->end+1);
            nhit++;
        }
        if ( nhit!=2 ) error("query failed, expected two hits, found %d: %s:%d-%d\n",nhit,chr,start,end);

        // fully contained interval, one hit
        start = 20000*i - 5000; end = 20000*i + 3000;
        set_line(line,sizeof(line),chr,start,end);
        if ( !regidx_overlap(idx,chr,start-1,end-1,itr) ) error("query failed, there should be a hit: %s:%d-%d\n",chr,start,end);
        debug("ok: overlap(s) found for %s:%d-%d\n",chr,start,end);
        nhit = 0;
        while ( regitr_overlap(itr) )
        {
            if ( itr->beg > end-1 || itr->end < start-1 ) error("query failed, incorrect region: %"PRIhts_pos"-%"PRIhts_pos" for %d-%d\n",itr->beg+1,itr->end+1,start,end);
            debug("\t %"PRIhts_pos"-%"PRIhts_pos"\n",itr->beg+1,itr->end+1);
            nhit++;
        }
        if ( nhit!=1 ) error("query failed, expected one hit, found %d: %s:%d-%d\n",nhit,chr,start,end);
    }
    regitr_destroy(itr);
    regidx_destroy(idx);
}

static void usage(void)
{
    fprintf(stderr, "Usage: test-regidx [OPTIONS]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -h, --help          this help message\n");
    fprintf(stderr, "   -s, --seed <int>    random seed\n");
    fprintf(stderr, "   -v, --verbose       increase verbosity by giving multiple times\n");

    exit(1);
}

int main(int argc, char **argv)
{
    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"verbose",0,0,'v'},
        {"seed",1,0,'s'},
        {0,0,0,0}
    };
    int c;
    int seed = (int)time(NULL);
    while ((c = getopt_long(argc, argv, "hvs:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 's': seed = atoi(optarg); break;
            case 'v': verbose++; break;
            default: usage(); break;
        }
    }

    info("Testing sequential access\n");
    test_sequential_access();

    info("Testing TAB\n");
    test(create_line_tab,regidx_parse_tab);

    info("Testing REG\n");
    test(create_line_reg,regidx_parse_reg);

    info("Testing BED\n");
    test(create_line_bed,regidx_parse_bed);

    info("Testing custom payload\n");
    test_custom_payload();

    info("Testing cases encountered in past\n");
    test_explicit("12:2064519-2064763","12:2064488-2067434","1");

    int i, ntest = 1000, nreg = 50;
    srand(seed);
    info("%d randomized tests, %d regions per test. Random seed is %d\n", ntest,nreg,seed);
    for (i=0; i<ntest; i++) test_random(nreg,1,1000);

    return 0;
}



#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "vcf.h"

int main_samview(int argc, char *argv[]);
int main_vcfview(int argc, char *argv[]);
int main_bamidx(int argc, char *argv[]);
int main_bcfidx(int argc, char *argv[]);
int main_bamshuf(int argc, char *argv[]);
int main_bam2fq(int argc, char *argv[]);
int main_tabix(int argc, char *argv[]);
int main_abreak(int argc, char *argv[]);
int main_bam2bed(int argc, char *argv[]);
int main_vcfcheck(int argc, char *argv[]);
int main_vcfisec(int argc, char *argv[]);
int main_vcfmerge(int argc, char *argv[]);
int main_vcfquery(int argc, char *argv[]);
int main_vcffilter(int argc, char *argv[]);
int main_vcfgtcheck(int argc, char *argv[]);

typedef struct
{
    int (*func)(int, char*[]);
    const char *alias, *help, *sep;
}
cmd_t;

// By symlinking, the htscmd program can be invoked in multiple
//  ways. For example, by creating a symlink 'htsvcf', the command
//  'htscmd vcfcheck' can be run as 'htsvcf check'. By creating
//  a symlink 'vcf', the same command can be run as 'vcf check'.
static cmd_t cmds[] =
{
    { .func  = main_samview,  
      .alias = "htscmd samview", 
      .help  = "SAM<->BAM conversion",
      .sep   = "General and indexing tools:"
    },
    { .func  = main_vcfview,  
      .alias = "htscmd vcfview, htsvcf view, vcf view", 
      .help  = "VCF<->BCF conversion",
      .sep   = NULL
    },
    { .func  = main_tabix,    
      .alias = "htscmd tabix, htsvcf tabix, vcf tabix",
      .help  = "tabix for BGZF'd BED, GFF, SAM, VCF and more",
      .sep   = NULL
    },
    { .func  = main_bamidx,   
      .alias = "htscmd bamidx",
      .help  = "index BAM", 
      .sep   = NULL
    },
    { .func = main_bcfidx,   
      .alias = "htscmd bcfidx, htsvcf idx, vcf idx",
      .help = "index BCF",
      .sep   = NULL
    },
    { .func  = main_bamshuf,  
      .alias = "htscmd bamshuf", 
      .help  = "shuffle BAM and group alignments by query name", 
      .sep   = "Alignment tools:"
    },
    { .func  = main_bam2fq,   
      .alias = "htscmd bam2fq",  
      .help  = "convert name grouped BAM to interleaved fastq",
      .sep   = NULL
    },
    { .func  = main_abreak,   
      .alias = "htscmd abreak",  
      .help  = "summarize assembly break points",
      .sep   = NULL
    },
    { .func  = main_bam2bed,  
      .alias = "htscmd bam2bed", 
      .help  = "BAM->BED conversion",
      .sep   = NULL
    },
    { .func  = main_vcfcheck, 
      .alias = "htscmd vcfcheck, htsvcf check, vcf check",
      .help  = "produce VCF stats",
      .sep   = "VCF/BCF tools:"
    },
    { .func  = main_vcffilter, 
      .alias = "htscmd vcffilter, htsvcf filter, vcf filter",
      .help  = "filter VCF files",
      .sep   = NULL
    },
    { .func  = main_vcfgtcheck, 
      .alias = "htscmd vcfgtcheck, htsvcf gtcheck, vcf gtcheck",
      .help  = "check sample identity",
    },
    { .func  = main_vcfisec,  
      .alias = "htscmd vcfisec, htsvcf isec, vcf isec", 
      .help  = "intersections of VCF files",
      .sep   = NULL
    },
    { .func  = main_vcfmerge, 
      .alias = "htscmd vcfmerge, htsvcf merge, vcf merge",
      .help  = "merge VCF files",
      .sep   = NULL
    },
    { .func  = main_vcfquery, 
      .alias = "htscmd vcfquery, htsvcf query, vcf query",
      .help  = "transform VCF into user-defined formats",
      .sep   = NULL
    },
    { .func  = NULL,
      .alias = NULL,
      .help  = NULL,
      .sep   = NULL
    }
};

static int set_alias(cmd_t *cmd, const char *argv0, char **buf, int *nbuf)
{
    int n = strlen(argv0);
    const char *b = cmd->alias;
    while ( *b )
    {
        while ( *b && isspace(*b) ) b++;    // skip leading spaces
        if ( !strncmp(argv0,b,n) && isspace(b[n]) )
        {
            b += n+1;
            break;  // found
        }

        while ( *b && !isspace(*b) ) b++;   // skip the argv0 string
        while ( *b && *b!=',' ) b++;        // skip the cmd string
        if ( *b==',' ) b++;
    }

    if ( !*b ) return 0;
    const char *e = b;
    while ( *e && *e!=',' ) e++;
    if ( !*buf || e-b+1 > *nbuf ) { *nbuf = e-b+1; *buf = (char*) realloc(*buf, *nbuf); }
    memcpy(*buf,b,e-b);
    (*buf)[e-b] = 0;
    return 1;
}


static int known_alias(char *argv0)
{
    int i = 0, nbuf = 0, known = 0;
    char *buf = NULL;
    while (cmds[i].func)
    {
        if ( (known = set_alias(&cmds[i], argv0, &buf, &nbuf)) ) break;
        i++;
    }
    if ( buf ) free(buf);
    return known;
}

static int usage(char *argv0)
{
	fprintf(stderr, "\nUsage:   %s <command> <argument>\n", argv0);
	fprintf(stderr, "Commands:\n");

    if ( !known_alias(argv0) ) argv0 = "htscmd";

    int i = 0, nbuf = 0;
    char *buf = NULL;
    const char *sep = NULL;
    while (cmds[i].func)
    {
        if ( cmds[i].sep ) sep = cmds[i].sep;
        if ( set_alias(&cmds[i], argv0, &buf, &nbuf) )
        {
            if ( sep )
            {
                printf("\n -- %s\n", sep);
                sep = NULL;
            }
            printf("\t%-15s %s\n", buf, cmds[i].help);
        }
        i++;
    }
    if ( buf ) free(buf);

    fprintf(stderr,"\n");
	return 1;
}

int main(int argc, char *argv[])
{
    char *a0 = argv[0] + strlen(argv[0]) - 1;
    while ( a0>argv[0] && a0[-1]!='/' && a0[-1]!='\\' ) a0--;
    
	if (argc < 2) return usage(a0);
    if ( !known_alias(a0) ) a0 = "htscmd";

    int i = 0, nbuf = 0;
    char *buf = NULL;
    while (cmds[i].func)
    {
        if ( set_alias(&cmds[i], a0, &buf, &nbuf) )
        {
            if ( !strcmp(argv[1],buf) ) 
            {
                free(buf);
                return cmds[i].func(argc-1,argv+1);
            }
        }
        i++;
    }
    if ( buf ) free(buf);
    fprintf(stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
    return 1;
}


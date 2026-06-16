#    Copyright (C) 2017, 2024, 2026 Genome Research Ltd.
#
#    Author: Robert Davies <rmd@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# First field:
#   INIT = initialisation, not counted in testing
#   P = expected to pass (zero return; expected output matches, if present)
#   N = expected to return non-zero
#   F = expected to fail
#
# Second field (P/N/F only):
#   Filename of expected output.  If '.', output is not checked
#
# Rest:
#   Command to execute.  $bgzip and $tabix are replaced with the path to
# bgzip and tabix. 

# TBI index on VCF
INIT $bgzip -c vcf_file.vcf > vcf_file.tbi.tmp.vcf.gz
P . $tabix -f -p vcf vcf_file.tbi.tmp.vcf.gz
P vcf_file.1.3000151.out $tabix vcf_file.tbi.tmp.vcf.gz 1:3000151-3000151
P vcf_file.1.3000151.out $tabix -u vcf_file.tbi.tmp.vcf.gz 1:3000151-3000151
P vcf_file.2.3199812.out $tabix vcf_file.tbi.tmp.vcf.gz 2:3199812-3199812
P vcf_file.2.3199812.out $tabix -u vcf_file.tbi.tmp.vcf.gz 2:3199812-3199812

# CSI index on VCF
INIT $bgzip -c vcf_file.vcf > vcf_file.csi.tmp.vcf.gz
P . $tabix -f -C -p vcf vcf_file.csi.tmp.vcf.gz
P vcf_file.1.3000151.out $tabix vcf_file.csi.tmp.vcf.gz 1:3000151-3000151
P vcf_file.1.3000151.out $tabix -u vcf_file.csi.tmp.vcf.gz 1:3000151-3000151
P vcf_file.2.3199812.out $tabix vcf_file.csi.tmp.vcf.gz 2:3199812-3199812
P vcf_file.2.3199812.out $tabix -u vcf_file.csi.tmp.vcf.gz 2:3199812-3199812

# CSI index on BCF
INIT cp -f vcf_file.bcf vcf_file.tmp.bcf
P . $tabix -f vcf_file.tmp.bcf
P vcf_file.1.3000151.out $tabix vcf_file.tmp.bcf 1:3000151-3000151
P vcf_file.1.3000151.out $tabix -u vcf_file.tmp.bcf 1:3000151-3000151
P vcf_file.2.3199812.out $tabix vcf_file.tmp.bcf 2:3199812-3199812
P vcf_file.2.3199812.out $tabix -u vcf_file.tmp.bcf 2:3199812-3199812

# VCF file with chromosome > 2^29-1 bases long
# TBI cannot index this file, so building the index should fail
INIT $bgzip -c large_chr.vcf > large_chr.tmp.vcf.gz
N . $tabix -f -p vcf large_chr.tmp.vcf.gz

# CSI can handle positions > 2^29-1, so building should work
P . $tabix -f -C -p vcf large_chr.tmp.vcf.gz
P large_chr.20.1.2147483647.out $tabix large_chr.tmp.vcf.gz chr20:1-2147483647
P large_chr.20.1.2147483647.out $tabix -u large_chr.tmp.vcf.gz chr20:1-2147483647

# VCF file with overlaps due to deletions
INIT $bgzip -c vcf_overlaps.vcf > vcf_overlaps.tmp.vcf.gz
INIT printf 'chr1 1000002 1000002\nchr1 1000007 1000007\nchr2 1000001 1000001\nchr2 1000009 1000009\nchr3 1000002 1000002\nchr3 1000006 1000006\n' > vcf_overlaps.tmp.regions
P . $tabix -f -p vcf vcf_overlaps.tmp.vcf.gz
P vcf_overlaps.no_u.out $tabix vcf_overlaps.tmp.vcf.gz chr1:1000002-1000002 chr1:1000007-1000007 chr2:1000001-1000001 chr2:1000009-1000009 chr3:1000002-1000002 chr3:1000006-1000006
P vcf_overlaps.unique.out $tabix -u vcf_overlaps.tmp.vcf.gz chr1:1000002-1000002 chr1:1000007-1000007 chr2:1000001-1000001 chr2:1000009-1000009 chr3:1000002-1000002 chr3:1000006-1000006
P vcf_overlaps.no_u.out $tabix -R vcf_overlaps.tmp.regions vcf_overlaps.tmp.vcf.gz
P vcf_overlaps.unique.out $tabix -u -R vcf_overlaps.tmp.regions vcf_overlaps.tmp.vcf.gz
P vcf_overlaps.unique.out $tabix -T vcf_overlaps.tmp.regions vcf_overlaps.tmp.vcf.gz

# BCF file with overlaps due to deletions
# ( BCF file made using test_view -b -l 0 vcf_overlaps.vcf )
INIT cp -f vcf_overlaps.bcf vcf_overlaps.tmp.bcf
P . $tabix -f vcf_overlaps.tmp.bcf
P vcf_overlaps.no_u.out $tabix vcf_overlaps.tmp.bcf chr1:1000002-1000002 chr1:1000007-1000007 chr2:1000001-1000001 chr2:1000009-1000009 chr3:1000002-1000002 chr3:1000006-1000006
P vcf_overlaps.unique.out $tabix -u vcf_overlaps.tmp.bcf chr1:1000002-1000002 chr1:1000007-1000007 chr2:1000001-1000001 chr2:1000009-1000009 chr3:1000002-1000002 chr3:1000006-1000006
P vcf_overlaps.no_u.out $tabix -R vcf_overlaps.tmp.regions vcf_overlaps.tmp.bcf
P vcf_overlaps.unique.out $tabix -u -R vcf_overlaps.tmp.regions vcf_overlaps.tmp.bcf
P vcf_overlaps.unique.out $tabix -T vcf_overlaps.tmp.regions vcf_overlaps.tmp.bcf

# TBI index on BED
INIT $bgzip -c bed_file.bed > bed_file.tbi.tmp.bed.gz
P . $tabix -f -p bed bed_file.tbi.tmp.bed.gz
P bed_file.Y.100200.out $tabix bed_file.tbi.tmp.bed.gz Y:100200-100200
P bed_file.Y.100200.out $tabix -u bed_file.tbi.tmp.bed.gz Y:100200-100200
P bed_file.unique.out $tabix -u bed_file.tbi.tmp.bed.gz Y:100201-100201 Y:100801-100801

# TBI index on GFF3
INIT $bgzip -c gff_file.gff > gff_file.tbi.tmp.gff.gz
P . $tabix -f -p gff gff_file.tbi.tmp.gff.gz
P gff_file.X.2934832.2935190.out $tabix gff_file.tbi.tmp.gff.gz X:2934832-2935190
P gff_file.X.2934832.2935190.out $tabix -u gff_file.tbi.tmp.gff.gz X:2934832-2935190
P gff_file.unique.out $tabix -u gff_file.tbi.tmp.gff.gz X:2935191-2936741 X:2936864-2938094 X:2943200-2945997

# tabix with --separate-regions
P bed_file.separate.out $tabix --separate-regions bed_file.tbi.tmp.bed.gz X:1100-1400 Y:100000-100550 Z:100000-100005

# Using threads with tabix
P . $tabix -f -p bed bed_file.tbi.tmp.bed.gz -@ 2
P vcf_file.1.3000151.out $tabix vcf_file.tbi.tmp.vcf.gz 1:3000151-3000151 --threads 2
P vcf_file.1.3000151.out $tabix -u vcf_file.tbi.tmp.vcf.gz 1:3000151-3000151 --threads 2

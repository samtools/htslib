# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.")
endif()

# choose the architecture 
set(DEPLOY_PLATFORM "unknown" CACHE STRING "choose the instruction set architecture")
set_property(CACHE DEPLOY_PLATFORM PROPERTY STRINGS unknown x86 arm64-v8a armeabi-v7a)
IF (${DEPLOY_PLATFORM} STREQUAL "unknown")
	message(FATAL_ERROR "choose the DEPLOY_PLATFORM") 
	return() # This is to stop proceeding further and to stop opencv getting set to the default ANDROID_ABI
ENDIF()


# Set minimum cmake version + cmake settings
cmake_minimum_required (VERSION 3.9)
cmake_policy(SET CMP0048 OLD)

# Specify the project, and set the use of c++ and c
project(htslib LANGUAGES CXX C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fPIC -rdynamic -Wall -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -rdynamic -Wall -O2 -std=c++11 ")


include_directories(${PROJECT_SOURCE_DIR})


set (bcf_sr_sort bcf_sr_sort.c)
set (bgzf bgzf.c)
# set (bgzip bgzip.c) # has a main()
# set (hfile_s3_write hfile_s3_write.c) # only if libcurl (and maybe s3) is enabled 
set (hts hts.c)
# set (htsfile htsfile.c) # has a main()
set (hts_os hts_os.c)
set (kfunc kfunc.c)
set (knetfile knetfile.c)
set (kstring kstring.c)
set (md5 md5.c)
# set (os_rand os/rand.c) # refer hts_os.c 
set (pluging plugin.c)
set (probaln probaln.c)
set (realn realn.c)
set (regidx regidx.c)
set (region region.c)
set (sam sam.c)
set (synced_bcf_reader synced_bcf_reader.c)
# set (tabix tabix.c) # has a main()
set (tbx tbx.c)
set (textutils textutils.c)
set (thread_pool thread_pool.c)
set (vcf vcf.c)
set (vcfutils vcfutils.c)
set (vcf_sweep vcf_sweep.c)
set (errmod errmod.c)
set (faidx faidx.c)
set (header header.c)
set (hfile hfile.c)
set (hfile_net hfile_net.c)
set (cram_stats cram/cram_stats.c)
set (cram_codecs cram/cram_codecs.c)
set (cram_decode cram/cram_decode.c)
set (cram_encode cram/cram_encode.c)
set (cram_external cram/cram_external.c)
set (cram_index cram/cram_index.c)
set (cram_io cram/cram_io.c)
set (cram_samtools cram/cram_samtools.c)
set (mFILE cram/mFILE.c)
set (open_trace_file cram/open_trace_file.c)
set (pooled_alloc cram/pooled_alloc.c)
set (rANS_static cram/rANS_static.c)
set (string_alloc cram/string_alloc.c)
# set (hfile_s3 hfile_s3.c) # only if s3 is enabled
set (multipart multipart.c)


# Build a shared lib
add_library(hts SHARED ${bcf_sr_sort} ${bgzf} ${hts} ${hts_os} ${kfunc} 
	${knetfile} ${kstring} ${md5} ${pluging} ${probaln} ${realn} ${regidx} ${region} ${sam} 
	${synced_bcf_reader} ${tbx} ${textutils} ${thread_pool} ${vcf} ${vcfutils} ${vcf_sweep} 
	${errmod} ${faidx} ${header} ${hfile} ${hfile_net} ${cram_stats} ${cram_codecs} ${cram_decode} 
	${cram_encode} ${cram_external} ${cram_index} ${cram_io} ${cram_samtools} ${mFILE} ${open_trace_file} 
	${pooled_alloc} ${rANS_static} ${string_alloc} ${multipart} )

# Build a static lib
# add_library(hts STATIC ${bcf_sr_sort} ${bgzf} ${hts} ${hts_os} ${kfunc} 
# 	${knetfile} ${kstring} ${md5} ${pluging} ${probaln} ${realn} ${regidx} ${region} ${sam} 
# 	${synced_bcf_reader} ${tbx} ${textutils} ${thread_pool} ${vcf} ${vcfutils} ${vcf_sweep} 
# 	${errmod} ${faidx} ${header} ${hfile} ${hfile_net} ${cram_stats} ${cram_codecs} ${cram_decode} 
# 	${cram_encode} ${cram_external} ${cram_index} ${cram_io} ${cram_samtools} ${mFILE} ${open_trace_file} 
# 	${pooled_alloc} ${rANS_static} ${string_alloc} ${multipart} )


IF (${DEPLOY_PLATFORM} STREQUAL "x86")
	target_link_libraries(hts -lz -ldl -lm -lpthread)
ELSE()
	target_link_libraries(hts -lz -ldl -lm)
ENDIF()
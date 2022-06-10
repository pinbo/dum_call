# dum_call
An example of using HTS lib to read bam files and call SNPs.

This tool can check SNPs, indels, big deletions and inversions from bam files mapped by `bwa mem` for CRISPR editing NGS check, same as [editcall](https://github.com/pinbo/editcall).


## Compile

```sh
# Step 1: get the htslib; I downloaded version 1.15.1 from the release
wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2
tar -xf htslib-1.15.1.tar.bz2
cd htslib-1.15.1

# Step 2: compile it based on the instruction
./configure --disable-libcurl    # you can disable more. Refer the file INSTALL
make

# Step 3: get dum_call
git clone https://github.com/pinbo/dum_call.git
cd dub_call
gcc -Wall -g -O2 dum_call.c -o dum_call ../libhts.a -lz -lm -llzma -lbz2 -lpthread

# to see the help
./dum_call
```

## Update

- 2022-06-10: added dum_callm.c to call multiple bam/sam files at one time.
- 2022-06-09: added dum_callv2.c that can read just a region of an indexed bam.
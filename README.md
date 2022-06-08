# dum_call
An example of using HTS lib to read bam files and call SNPs


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
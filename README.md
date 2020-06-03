# RabbitMash (experimental)

RabbitMash is an efficient highly optimized implementation of [Mash](https://github.com/marbl/Mash) which can take full advantage of modern hardware including multi-threading, vectorization, and fast I/O.



## Installation

The easiest way to use RabbitMash is to grab a binary release from [here](https://github.com/ZekunYin/RabbitMash/releases). We provide four versions of RabbitMash for different x86_64 architectures including:

| Version     | CPU flags                 |
| ----------- | ------------------------- |
| mash_nosimd | no requirement            |
| mash_sse4   | sse4_1                    |
| mash_avx2   | avx2                      |
| mash_avx512 | avx512f avx512bw avx512vl |

You can check the CPU Flags by `lscpu` to select corresponding binary.

All binaries have been tested on our server with **Intel Xeon Gold CPUs** running **CentOS 7.5/8.1**.

The binaries should also work on other linuxes. But if you have problem running these binaries, feel free to iussue us.



## Notes

**Process gzipped files**

RabbitMash supports plain FASTQ/FASTA and gzipped FASTQ/FASTA file formats.  When processing gziped files, the performance of `sketch -i` and `screen` operations are limited by decompression speed. Instead of processing gziped files directly, a more efficient strategy is to process these files by two steps: (i) decompress gziped files to FASTQ format by [libdeflate](https://github.com/ebiggers/libdeflate) or [pugz](https://github.com/Piezoid/pugz), and (ii) process FASTQ/FASTA files by RabbitMash. 

But when you need to sketch large dataset by files, there won't be much performance penalty to process gzipped files.



## Build

**Dependencies:**

   - Cap'n Proto ( https://capnproto.org/ )
   - GNU Scientific Library ( http://www.gnu.org/software/gsl/ )
   - Zlib ( included with most Linuxes, http://www.zlib.net ) 
   - GCC >= 5 (C++14 required)

**Build:**

```bash
git clone https://github.com/ZekunYin/RabbitMash.git
cd RabbitMash
./bootstrap.sh
./configure [--prefix=...] [--with-capnp=...] [--with-gsl=...] \
            [--with-simd=yes/no]
make -j4
#optional
make install
#optional
make test
```

**Build dependency-free binary:**

```bash
git clone https://github.com/ZekunYin/RabbitMash.git
cd RabbitMash
./bootstrap.sh
./configure [--prefix=...] [--with-capnp=...] [--with-gsl=...] \
            [--with-simd=yes/no] [--enable-static-gsl=yes]     \
            [--enable-static-cpp=yes]
make -j4
#optional
make install
#optional
make test
```

You can also type `./configure -h` for configure help information.

**Install dependency on CentOS 8.1 (root user):**

```bash
sudo dnf install capnproto capnproto-devel gsl gsl-devel
```

If you are not a root user, you need to build the dependecies from source code.



## Simple Usage

**sketch:**

```bash
./mash sketch test/genome1.fna -o test/genome1.fna.msh
./mash sketch test/genome2.fna -o test/genome2.fna.msh
```

**dist:**

```bash
 ./mash dist test/genome1.fna.msh test/genome2.fna.msh -o dist.bin
 #optional
 ./mash dumpdist test/genome1.fna.msh test/genome2.fna.msh dist.bin -o dist.txt
```

**triangle:**

```bash
./mash triangle test/genome1.fna.msh -o tri.bin
#optional
./mash dumptri test/genome1.fna.msh tri.bin -o tri.txt
```

**screen:**

```bash
./mash screen test/genome1.fna.msh test/reads1.fastq > scr.out
```



## Document

RabbitMash is based on [Mash](https://github.com/marbl/Mash) . All functions and most parameters of RabbitMash is the same with Mash.  Just type `mash` for command information and type `mash <command_name>` for help information.

See Mash's document  ([http://mash.readthedocs.org](http://mash.readthedocs.org)) for more information.



## Different Commands or Parameters to Mash

#### New parameter 

**sketch:**

```bash
-fw #Create mutiple msh files to keep low memory footprint for sketching massive sequences.
```

**dist:**

```bash
-o <text> #Create binary format result file for better performance. If -o is not specified, text results will be written to stdout.
```

**triangle:**

```bash
-o <text> #Create binary format result file for better performance. If -o is not specified, text results will be written to stdout.
```

#### New Command

```bash
mash dumpdist #Convert binary dist results to human-readable texts.
```

```bash
mash dumptri  #Convert binary triangle results to human-readable texts.
```



## Bug Report

All bug reports, comments and suggestions are welcome.

Feel free to open a new issue, normally I can make a response in one day if I'm not on vacation. 

## Cite

RabbitMash paper is under review now.

## Limitations

- OSX version has not been tested.
- x86 version has not been tested.

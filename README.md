# RabbitMash (experimental)

RabbitMash is an efficient highly optimized implementation of [Mash](https://github.com/marbl/Mash) which can take full advantage of modern hardware including multi-threading, vectorization, and fast I/O.



## Installation

The easiest way to use RabbitMash is to grab a binary release from [here](https://github.com/ZekunYin/RabbitMash/releases). We provide four versions of RabbitMash for different architectures including:

| Version      | CPU flags                 |
| ------------ | ------------------------- |
| mash_nosimd  | no requirement            |
| mash_sse4    | sse4_2                    |
| mash_avx2    | avx2                      |
| mash_avx512  | avx512f avx512bw avx512vl |

You can check the CPU Flags by `lscpu` to select corresponding binary.

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
#optimal
make install
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
#optimal
make install
```

You can also type `./configure -h` for configure help information.

**Install dependency on CentOS 8.1 (root user):**

```bash
sudo dnf install capnproto capnproto-devel gsl gsl-devel
```

If you are not a root user, you can build the dependecies from source code.



## Simple Usage

**sketch:**

```bash
./mash sketch test/genome1.fna -o test/genome1.fna.msh
./mash sketch test/genome2.fna -o test/genome2.fna.msh
```

**dist:**

```bash
 ./mash dist test/genome1.fna.msh test/genome2.fna.msh
```

**triangle:**

```bash
./mash triangle test/genome1.fna.msh
```

**screen:**

```bash
./mash screen test/genome1.fna.msh test/reads1.fastq
```



## Document

RabbitMash is based on [Mash](https://github.com/marbl/Mash) . All functions and most parameters of RabbitMash is the same with Mash. 

See Mash's document  ([http://mash.readthedocs.org](http://mash.readthedocs.org)) for more information.



## Different Commands or Parameters to Mash

#### New parameter 

**sketch:**

```bash
-fw #create mutiple msh files to keep low memory footprint for mass sequences.
```

**dist:**

```bash
-o #output file name in binary format for better performance
```

**triangle:**

```bash
-o #output to the binary phylip format for better performance.
```

#### New Command

```bash
outputBin: dump binary phylip format to phylip format
```



## Bug Report

All bug reports, comments and suggestions are welcome.

Feel free to open a new issue, normally I can make a response in one day if I'm not on vacation. 

## Cite

RabbitMash paper is under review now.

## Limitations

- OSX version is not tested.
- x86 version is not tested.

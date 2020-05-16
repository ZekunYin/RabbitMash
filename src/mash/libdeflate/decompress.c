/*
 * decompress.c - a file decompression apis from file to memory buffer
 *
 * Copyright 2016 Eric Biggers
 *
 * Last modified 5/17/2020 by Zekun Yin
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "prog_util.h"

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
#  include <sys/utime.h>
#else
#  include <sys/time.h>
#  include <unistd.h>
#  include <utime.h>
#endif

#include <stdint.h>

struct options {
	bool to_stdout;
	bool decompress;
	bool force;
	bool keep;
	int compression_level;
	const tchar *suffix;
};

static u32
load_u32_gzip(const u8 *p)
{
	return ((u32)p[0] << 0) | ((u32)p[1] << 8) |
		((u32)p[2] << 16) | ((u32)p[3] << 24);
}

static const tchar *
get_suffix(const tchar *path, const tchar *suffix)
{
	size_t path_len = tstrlen(path);
	size_t suffix_len = tstrlen(suffix);
	const tchar *p;

	if (path_len <= suffix_len)
		return NULL;
	p = &path[path_len - suffix_len];
	if (tstrxcmp(p, suffix) == 0)
		return p;
	return NULL;
}

static tchar *
append_suffix(const tchar *path, const tchar *suffix)
{
	size_t path_len = tstrlen(path);
	size_t suffix_len = tstrlen(suffix);
	tchar *suffixed_path;

	suffixed_path = xmalloc((path_len + suffix_len + 1) * sizeof(tchar));
	if (suffixed_path == NULL)
		return NULL;
	tmemcpy(suffixed_path, path, path_len);
	tmemcpy(&suffixed_path[path_len], suffix, suffix_len + 1);
	return suffixed_path;
}

static int
stat_file(struct file_stream *in, stat_t *stbuf, bool allow_hard_links)
{
	if (tfstat(in->fd, stbuf) != 0) {
		msg("%"TS": unable to stat file", in->name);
		return -1;
	}

	if (!S_ISREG(stbuf->st_mode) && !in->is_standard_stream) {
		msg("%"TS" is %s -- skipping",
		    in->name, S_ISDIR(stbuf->st_mode) ? "a directory" :
							"not a regular file");
		return -2;
	}

	if (stbuf->st_nlink > 1 && !allow_hard_links) {
		msg("%"TS" has multiple hard links -- skipping "
		    "(use -f to process anyway)", in->name);
		return -2;
	}

	return 0;
}

static char *
do_decompress(struct libdeflate_decompressor *decompressor,
	      struct file_stream *in, int *ret, int64_t * out_size)
{
	const u8 *compressed_data = in->mmap_mem;
	size_t compressed_size = in->mmap_size;
	void *uncompressed_data = NULL;
	size_t uncompressed_size;
	size_t actual_in_nbytes;
	size_t actual_out_nbytes;
	enum libdeflate_result result;
	//int ret = 0;
	char *out = NULL;
	int nblocks = 0;
	if (compressed_size < sizeof(u32)) {
	       msg("%"TS": not in gzip format", in->name);
	       *ret = -1;
	       goto output;
	}

	/*
	 * Use the ISIZE field as a hint for the decompressed data size.  It may
	 * need to be increased later, however, because the file may contain
	 * multiple gzip members and the particular ISIZE we happen to use may
	 * not be the largest; or the real size may be >= 4 GiB, causing ISIZE
	 * to overflow.  In any case, make sure to allocate at least one byte.
	 */
	uncompressed_size = load_u32_gzip(&compressed_data[compressed_size - 4]);
	if (uncompressed_size == 0)
		uncompressed_size = 1;

	do {
		if (uncompressed_data == NULL) {
			uncompressed_data = xmalloc(uncompressed_size);
			if (uncompressed_data == NULL) {
				msg("%"TS": file is probably too large to be "
				    "processed by this program", in->name);
				*ret = -1;
				goto output;
			}
		}

		result = libdeflate_gzip_decompress_ex(decompressor,
						       compressed_data,
						       compressed_size,
						       uncompressed_data,
						       uncompressed_size,
						       &actual_in_nbytes,
						       &actual_out_nbytes);

		if (result == LIBDEFLATE_INSUFFICIENT_SPACE) {
			if (uncompressed_size * 2 <= uncompressed_size) {
				msg("%"TS": file corrupt or too large to be "
				    "processed by this program", in->name);
				*ret = -1;
				goto output;
			}
			uncompressed_size *= 2;
			free(uncompressed_data);
			uncompressed_data = NULL;
			continue;
		}

		if (result != LIBDEFLATE_SUCCESS) {
			msg("%"TS": file corrupt or not in gzip format",
			    in->name);
			*ret = -1;
			goto output;
		}

		if (actual_in_nbytes == 0 ||
		    actual_in_nbytes > compressed_size ||
		    actual_out_nbytes > uncompressed_size) {
			msg("Bug in libdeflate_gzip_decompress_ex()!");
			*ret = -1;
			goto output;
		}

		//yzk
		//ret = full_write(out, uncompressed_data, actual_out_nbytes);
		if (*ret != 0)
			goto output;

		compressed_data += actual_in_nbytes;
		compressed_size -= actual_in_nbytes;
		nblocks ++;

	} while (compressed_size != 0);
	if(nblocks > 1)
		msg("%"TS": contains multiple blocks. WARNING!!!",
		    in->name);
output:
	//free(uncompressed_data);
	out = (char *)uncompressed_data;
	*out_size = uncompressed_size;
	return out;
}

static char *
decompress_file(struct libdeflate_decompressor *decompressor, const tchar *path,
		const struct options *options, int *ret, int64_t *out_size)
{
	tchar *oldpath = (tchar *)path;
	tchar *newpath = NULL;
	struct file_stream in;
	//struct file_stream out;
	stat_t stbuf;
	//int ret;
	char *out = NULL;
	//int ret2;

	if (path != NULL) {
		const tchar *suffix = get_suffix(path, options->suffix);
		if (suffix == NULL) {
			/*
			 * Input file is unsuffixed.  If the file doesn't exist,
			 * then try it suffixed.  Otherwise, if we're not
			 * writing to stdout, skip the file with warning status.
			 * Otherwise, go ahead and try to open the file anyway
			 * (which will very likely fail).
			 */
			if (tstat(path, &stbuf) != 0 && errno == ENOENT) {
				oldpath = append_suffix(path, options->suffix);
				if (oldpath == NULL){
					*ret = -1;
					return NULL;
				}
				if (!options->to_stdout)
					newpath = (tchar *)path;
			} else if (!options->to_stdout) {
				msg("\"%"TS"\" does not end with the %"TS" "
				    "suffix -- skipping",
				    path, options->suffix);
				//return -2;
				*ret = -2;
				return NULL;
			}
		} else if (!options->to_stdout) {
			/*
			 * Input file is suffixed, and we're not writing to
			 * stdout.  Strip the suffix to get the path to the
			 * output file.
			 */
			newpath = xmalloc((suffix - oldpath + 1) *
					  sizeof(tchar));
			if (newpath == NULL){
				*ret = -1;
				return NULL;
				//return -1;
			}
			tmemcpy(newpath, oldpath, suffix - oldpath);
			newpath[suffix - oldpath] = '\0';
		}
	}

	*ret = xopen_for_read(oldpath, options->force || options->to_stdout,
			     &in);
	if (*ret != 0)
		goto out_free_paths;

	if (!options->force && isatty(in.fd)) {
		msg("Refusing to read compressed data from terminal.  "
		    "Use -f to override.\nFor help, use -h.");
		*ret = -1;
		goto out_close_in;
	}

	*ret = stat_file(&in, &stbuf, options->force || options->keep ||
			oldpath == NULL || newpath == NULL);
	if (*ret != 0)
		goto out_close_in;

	//yzk
	//ret = xopen_for_write(newpath, options->force, &out); 
	//if (ret != 0) 
	//	goto out_close_in;

	/* TODO: need a streaming-friendly solution */
	*ret = map_file_contents(&in, stbuf.st_size);
	if (*ret != 0)
		goto out_close_in;//yzk
	//yzk
	out = do_decompress(decompressor, &in, ret, out_size);
	if (*ret != 0)
		goto out_close_in;//yzk

	//yzk
	//if (oldpath != NULL && newpath != NULL)
	//	restore_metadata(&out, newpath, &stbuf);
	*ret = 0;
//out_close_out:
	//yzk
	//ret2 = xclose(&out);
	//if (ret == 0)
	//	ret = ret2;
	//if (ret != 0 && newpath != NULL)
	//	tunlink(newpath);
out_close_in:
	xclose(&in);
//yzk we do not need to delete any files
//	if (ret == 0 && oldpath != NULL && newpath != NULL && !options->keep)
//		tunlink(oldpath);
out_free_paths:
	if (newpath != path)
		free(newpath);
	if (oldpath != path)
		free(oldpath);
	return out;
}
char * decompress_libdeflate(const char * filename, int *ret, int64_t * out_size){

	char *out = NULL;
	struct options options;
	options.to_stdout = false;
	options.decompress = true;
	options.force = false;
	options.keep = true;//keep original files
	options.compression_level = 6;//useless in this function
	options.suffix = T(".gz");

	struct libdeflate_decompressor *d;

	d = alloc_decompressor();
	if (d == NULL){
		*ret = 1;
		return NULL;
	}
	out = decompress_file(d, filename, &options, ret, out_size);

	libdeflate_free_decompressor(d);

	return out;
	
}

//int main()
//{	
//	const char *file = "test.fna.gz";
//	char *file_data = NULL;
//	int64_t file_size;
//	int ret;
//	
//	file_data = decompress_libdeflate(file, &ret, &file_size);
//
//	if(ret != 0)
//		fprintf(stderr, "decompress failed");
//	else
//		fprintf(stderr, "size of %s is %ld\n", file, file_size);
//
//	if(file_data != NULL){
//		FILE *out = fopen("test.fna", "w");
//		fwrite(file_data, 1, file_size, out);
//		fclose(out);
//		free(file_data);
//	}
//	else
//		fprintf(stderr, "no file data\n");
//
//	return 0;
//}

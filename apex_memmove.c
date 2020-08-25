/*
 * Apex memmove of T. Herselman.
 * The fastest memcpy/memmove on x86/x64 ... EVER,
 * written in C (according to the author).
 *
 * The original version has some fallback functions to use SSE/SSE2/SSE3,
 * I only keep the SSE4.2 version since most computers now support SSE4.2
 */


#include <emmintrin.h>   /* Intel/AMD SSE intrinsics */

#include "apex_memmove.h"


void *apex_kryptonite(void *dst, const void *src, size_t size)
{
  /*
   * based on memmove09 for "size <= 112" and memmove41 for "size > 112";
   * memmove09's "size <= 112" proved fastest overall (weighted),
   * even on Core i5!
   */
  if (size <= 112) {

    if (size >= 16) {
      const __m128i xmm0 = _mm_loadu_si128( (__m128i*)src );

      if (size > 16) {

        if (size >= 32) {
          const __m128i xmm1 = _mm_loadu_si128( (__m128i*)((char*)src + 16) );

          if (size > 32) {
            long long rax = *(long long*)((char*)src + size - 16);
            long long rcx = *(long long*)((char*)src + size - 8);

            if (size > 48) {
              const __m128i xmm2 = _mm_loadu_si128( (__m128i*)((char*)src + 32) );

              if (size > 64) {
                const __m128i xmm3 = _mm_loadu_si128( (__m128i*)((char*)src + 48) );

                if (size > 80) {
                  const __m128i xmm4 = _mm_loadu_si128( (__m128i*)((char*)src + 64) );

                  if (size > 96) {
                    const __m128i xmm5 = _mm_loadu_si128( (__m128i*)((char*)src + 80) );
                    *(long long*)((char*)dst + size - 16) = rax;
                    *(long long*)((char*)dst + size - 8) = rcx;
                    _mm_storeu_si128( (__m128i*)dst, xmm0 );
                    _mm_storeu_si128( (__m128i*)((char*)dst + 16), xmm1 );
                    _mm_storeu_si128( (__m128i*)((char*)dst + 32), xmm2 );
                    _mm_storeu_si128( (__m128i*)((char*)dst + 48), xmm3 );
                    _mm_storeu_si128( (__m128i*)((char*)dst + 64), xmm4 );
                    _mm_storeu_si128( (__m128i*)((char*)dst + 80), xmm5 );
                    return dst;
                  }

                  *(long long*)((char*)dst + size - 16) = rax;
                  *(long long*)((char*)dst + size - 8) = rcx;
                  _mm_storeu_si128( (__m128i*)dst, xmm0 );
                  _mm_storeu_si128( (__m128i*)((char*)dst + 16), xmm1 );
                  _mm_storeu_si128( (__m128i*)((char*)dst + 32), xmm2 );
                  _mm_storeu_si128( (__m128i*)((char*)dst + 48), xmm3 );
                  _mm_storeu_si128( (__m128i*)((char*)dst + 64), xmm4 );
                  return dst;
                }

                *(long long*)((char*)dst + size - 16) = rax;
                *(long long*)((char*)dst + size - 8) = rcx;
                _mm_storeu_si128( (__m128i*)dst, xmm0 );
                _mm_storeu_si128( (__m128i*)((char*)dst + 16), xmm1 );
                _mm_storeu_si128( (__m128i*)((char*)dst + 32), xmm2 );
                _mm_storeu_si128( (__m128i*)((char*)dst + 48), xmm3 );
                return dst;
              }

              *(long long*)((char*)dst + size - 16) = rax;
              *(long long*)((char*)dst + size - 8) = rcx;
              _mm_storeu_si128( (__m128i*)dst, xmm0 );
              _mm_storeu_si128( (__m128i*)((char*)dst + 16), xmm1 );
              _mm_storeu_si128( (__m128i*)((char*)dst + 32), xmm2 );
              return dst;
            }

            *(long long*)((char*)dst + size - 16) = rax;
            *(long long*)((char*)dst + size - 8) = rcx;
          }

          _mm_storeu_si128( (__m128i*)dst, xmm0 );
          _mm_storeu_si128( (__m128i*)((char*)dst + 16), xmm1 );
          return dst;
        }

        long long rax = *(long long*)((char*)src + size - 16);
        long long rcx = *(long long*)((char*)src + size - 8);
        *(long long*)((char*)dst + size - 16) = rax;
        *(long long*)((char*)dst + size - 8) = rcx;
      }

      _mm_storeu_si128( (__m128i*)dst, xmm0 );
      return dst;
    }

    if (size >= 8) {
      long long rax = *(long long*)src;

      if (size > 8) {
        long long rcx = *(long long*)((char*)src + size - 8);
        *(long long*)dst = rax;
        *(long long*)((char*)dst + size - 8) = rcx;
      }
      else *(long long*)dst = rax;
    }
    else if (size >= 4) {
      int eax = *(int*)src;

      if (size > 4) {
        int ecx = *(int*)((char*)src + size - 4);
        *(int*)dst = eax;
        *(int*)((char*)dst + size - 4) = ecx;
      }
      else *(int*)dst = eax;
    }
    else if (size >= 1) {
      char al = *(char*)src;

      if (size > 1) {
        short cx = *(short*)((char*)src + size - 2);
        *(char*)dst = al;
        *(short*)((char*)dst + size - 2) = cx;
      }
      else *(char*)dst = al;
    }

    return dst;
  }
  else {
    void * const ret = dst;

    if (((size_t)dst - (size_t)src) >= size) {

      if (size < (1024 * 256)) {
        /* "Round down to nearest multiple of 64" */
        long long offset = (long long)(size & -0x20);
        dst = (char*)dst + offset;  /* "Point to the end" */
        src = (char*)src + offset;  /* "Point to the end" */
        size -= offset;             /* "Remaining data after loop" */
        offset = -offset;           /* "Negative index from the end" */

        do {
          const __m128i xmm0 = _mm_loadu_si128( (__m128i*)((char*)src + offset) );
          const __m128i xmm1 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 16) );
          _mm_storeu_si128( (__m128i*)((char*)dst + offset), xmm0 );
          _mm_storeu_si128( (__m128i*)((char*)dst + offset + 16), xmm1 );
        } while (offset += 32);

        if (size >= 16) {

          if (size > 16) {
            const __m128i xmm7 = _mm_loadu_si128( (__m128i*)((char*)src + size - 16) );
            const __m128i xmm0 = _mm_loadu_si128( (__m128i*)src );
            _mm_storeu_si128( (__m128i*)((char*)dst + size - 16), xmm7 );
            _mm_storeu_si128( (__m128i*)dst, xmm0 );
            return ret;
          }

          _mm_storeu_si128( (__m128i*)dst, _mm_loadu_si128( (__m128i*)src ) );
          return ret;
        }
      }
      /* do forward streaming copy/move */
      else {
        /* We MUST do prealignment on streaming copies! */
        const size_t prealign = -(size_t)dst & 0xf;

        if (prealign) {

          if (prealign >= 8) {
            long long rax = *(long long*)src;

            if (prealign > 8) {
              long long rcx = *(long long*)((char*)src + prealign - 8);
              *(long long*)dst = rax;
              *(long long*)((char*)dst + prealign - 8) = rcx;
            }
            else *(long long*)dst = rax;
          }
          else if (prealign >= 4) {
            int eax = *(int*)src;

            if (prealign > 4) {
              int ecx = *(int*)((char*)src + prealign - 4);
              *(int*)dst = eax;
              *(int*)((char*)dst + prealign - 4) = ecx;
            }
            else *(int*)dst = eax;
          }
          else {
            char al = *(char*)src;

            if (prealign > 1) {
              short cx = *(short*)((char*)src + prealign - 2);
              *(char*)dst = al;
              *(short*)((char*)dst + prealign - 2) = cx;
            }
            else *(char*)dst = al;
          }

          src = (char*)src + prealign;
          dst = (char*)dst + prealign;
          size -= prealign;
        }

        /* Begin prefetching upto 4KB */
        for (long long offset = 0; offset < 4096; offset += 256) {
          _mm_prefetch( ((char*)src + offset), _MM_HINT_NTA );
          _mm_prefetch( ((char*)src + offset + 64), _MM_HINT_NTA );
          _mm_prefetch( ((char*)src + offset + 128), _MM_HINT_NTA );
          _mm_prefetch( ((char*)src + offset + 192), _MM_HINT_NTA );
        }

        /* "Round down to nearest multiple of 64" */
        long long offset = (long long)(size & -0x40);
        size -= offset;             /* "Remaining data after loop" */
        offset -= 4096;             /* stage 1 INCLUDES prefetches */
        dst = (char*)dst + offset;  /* "Point to the end" */
        src = (char*)src + offset;  /* "Point to the end" */
        offset = -offset;           /* "Negative index from the end" */

        /* stage 1 ~~ WITH prefetching */
        do {
          _mm_prefetch( (char*)src + offset + 4096, _MM_HINT_NTA );
          const __m128i xmm0 = _mm_loadu_si128( (__m128i*)((char*)src + offset) );
          const __m128i xmm1 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 16) );
          const __m128i xmm2 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 32) );
          const __m128i xmm3 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 48) );
          _mm_stream_si128( (__m128i*)((char*)dst + offset), xmm0 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 16), xmm1 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 32), xmm2 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 48), xmm3 );
        } while (offset += 64);

        offset = -4096;
        dst = (char*)dst + 4096;
        src = (char*)src + 4096;

        /* prefetch the final tail section */
        _mm_prefetch( ((char*)src + size - 64), _MM_HINT_NTA );

        /* stage 2 ~~ WITHOUT further prefetching */
        do {
          const __m128i xmm0 = _mm_loadu_si128( (__m128i*)((char*)src + offset) );
          const __m128i xmm1 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 16) );
          const __m128i xmm2 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 32) );
          const __m128i xmm3 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 48) );
          _mm_stream_si128( (__m128i*)((char*)dst + offset), xmm0 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 16), xmm1 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 32), xmm2 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 48), xmm3 );
        } while (offset += 64);

        if (size >= 16) {
          const __m128i xmm0 = _mm_loadu_si128( (__m128i*)src );

          if (size > 16) {

            if (size > 32) {
              const __m128i xmm1 = _mm_loadu_si128( (__m128i*)((char*)src + 16) );
              const __m128i xmm6 = _mm_loadu_si128( (__m128i*)((char*)src + size - 32) );
              const __m128i xmm7 = _mm_loadu_si128( (__m128i*)((char*)src + size - 16) );
              _mm_stream_si128( (__m128i*)dst, xmm0 );
              _mm_stream_si128( (__m128i*)((char*)dst + 16), xmm1 );
              _mm_storeu_si128( (__m128i*)((char*)dst + size - 32), xmm6 );
              _mm_storeu_si128( (__m128i*)((char*)dst + size - 16), xmm7 );
              return ret;
            }

            const __m128i xmm7 = _mm_loadu_si128( (__m128i*)((char*)src + size - 16) );
            _mm_stream_si128( (__m128i*)dst, xmm0 );
            _mm_storeu_si128( (__m128i*)((char*)dst + size - 16), xmm7 );
            return ret;
          }

          _mm_stream_si128( (__m128i*)dst, xmm0 );
          return ret;
        }
      }

      if (size >= 8) {
        long long rax = *(long long*)src;

        if (size > 8) {
          long long rcx = *(long long*)((char*)src + size - 8);
          *(long long*)dst = rax;
          *(long long*)((char*)dst + size - 8) = rcx;
        }
        else *(long long*)dst = rax;
      }
      else if (size >= 4) {
        int eax = *(int*)src;

        if (size > 4) {
          int ecx = *(int*)((char*)src + size - 4);
          *(int*)dst = eax;
          *(int*)((char*)dst + size - 4) = ecx;
        }
        else *(int*)dst = eax;
      }
      else if (size >= 1) {
        char al = *(char*)src;

        if (size > 1) {
          short cx = *(short*)((char*)src + size - 2);
          *(char*)dst = al;
          *(short*)((char*)dst + size - 2) = cx;
        }
        else *(char*)dst = al;
      }

      return ret;
    }
    /* src < dst ... do reverse copy */
    else {
      src = (char*)src + size;
      dst = (char*)dst + size;

      if (size < (1024 * 256)) {
        /* "Round down to nearest multiple of 64" */
        long long offset = (long long)(size & -0x20);
        /* "Point to the end" ... actually, we point to the start! */
        dst = (char*)dst - offset;
        /* "Point to the end" ... actually, we point to the start! */
        src = (char*)src - offset;
        /* "Remaining data after loop" */
        size -= offset;
        /* "Negative index from the end" ... not when doing reverse copy/move! */
        /* offset = -offset; */

        offset -= 32;

        do {
          const __m128i xmm2 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 16) );
          const __m128i xmm3 = _mm_loadu_si128( (__m128i*)((char*)src + offset) );
          _mm_storeu_si128( (__m128i*)((char*)dst + offset + 16), xmm2 );
          _mm_storeu_si128( (__m128i*)((char*)dst + offset), xmm3 );
        } while ((offset -= 32) >= 0);

        if (size >= 16) {

          if (size > 16) {
            size = -size;
            /*
             * The order has been mixed so the compiler
             * will not re-order the statements!
             */
            const __m128i xmm7 = _mm_loadu_si128( (__m128i*)((char*)src + size) );
            const __m128i xmm0 = _mm_loadu_si128( (__m128i*)((char*)src - 16) );
            _mm_storeu_si128( (__m128i*)((char*)dst + size), xmm7 );
            _mm_storeu_si128( (__m128i*)((char*)dst - 16), xmm0 );
            return ret;
          }

          _mm_storeu_si128( (__m128i*)((char*)dst - 16), _mm_loadu_si128( (__m128i*)((char*)src - 16) ) );
          return ret;
        }
      }
      /* do reversed streaming copy/move */
      else {
        /* We MUST do prealignment on streaming copies! */
        const size_t prealign = (size_t)dst & 0xf;

        if (prealign) {
          src = (char*)src - prealign;
          dst = (char*)dst - prealign;
          size -= prealign;

          if (prealign >= 8) {
            long long rax = *(long long*)((char*)src + prealign - 8);

            if (prealign > 8) {
              long long rcx = *(long long*)src;
              *(long long*)((char*)dst + prealign - 8) = rax;
              *(long long*)dst = rcx;
            }
            /*
             * different on purpose, because we know the exact size now,
             * which is 8, and "dst" has already been aligned!
             */
            else *(long long*)dst = rax;
          }
          else if (prealign >= 4) {
            int eax = *(int*)((char*)src + prealign - 4);

            if (prealign > 4) {
              int ecx = *(int*)src;
              *(int*)((char*)dst + prealign - 4) = eax;
              *(int*)dst = ecx;
            }
            else *(int*)dst = eax;  /* different on purpose! */
          }
          else {
            char al = *(char*)((char*)src + prealign - 1);

            if (prealign > 1) {
              short cx = *(short*)src;
              *(char*)((char*)dst + prealign - 1) = al;
              *(short*)dst = cx;
            }
            else *(char*)dst = al;  /* different on purpose! */
          }
        }

        /* Begin prefetching upto 4KB */
        for (long long offset = 0; offset > -4096; offset -= 256) {
          _mm_prefetch( ((char*)src + offset - 64), _MM_HINT_NTA );
          _mm_prefetch( ((char*)src + offset - 128), _MM_HINT_NTA );
          _mm_prefetch( ((char*)src + offset - 192), _MM_HINT_NTA );
          _mm_prefetch( ((char*)src + offset - 256), _MM_HINT_NTA );
        }

        /* "Round down to nearest multiple of 64" */
        long long offset = (long long)(size & -0x40);
        /* "Remaining data after loop" */
        size -= offset;
        /* stage 1 INCLUDES prefetches */
        offset -= 4096;
        /* "Point to the end" ... actually, we point to the start! */
        dst = (char*)dst - offset;
        /* "Point to the end" ... actually, we point to the start! */
        src = (char*)src - offset;
        /* "Negative index from the end" ... not when doing reverse copy/move! */
        /* offset = -offset; */

        offset -= 64;

        /* stage 1 ~~ WITH prefetching */
        do {
          _mm_prefetch( (char*)src + offset - 4096, _MM_HINT_NTA );
          const __m128i xmm0 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 48) );
          const __m128i xmm1 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 32) );
          const __m128i xmm2 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 16) );
          const __m128i xmm3 = _mm_loadu_si128( (__m128i*)((char*)src + offset) );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 48), xmm0 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 32), xmm1 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 16), xmm2 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset), xmm3 );
        } while ((offset -= 64) >= 0);

        offset = 4096;
        dst = (char*)dst - 4096;
        src = (char*)src - 4096;

        /* prefetch the final tail section */
        _mm_prefetch( ((char*)src - 64), _MM_HINT_NTA );

        offset -= 64;

        /* stage 2 ~~ WITHOUT further prefetching */
        do {
          const __m128i xmm0 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 48) );
          const __m128i xmm1 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 32) );
          const __m128i xmm2 = _mm_loadu_si128( (__m128i*)((char*)src + offset + 16) );
          const __m128i xmm3 = _mm_loadu_si128( (__m128i*)((char*)src + offset) );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 48), xmm0 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 32), xmm1 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset + 16), xmm2 );
          _mm_stream_si128( (__m128i*)((char*)dst + offset), xmm3 );
        } while ((offset -= 64) >= 0);

        if (size >= 16) {
          const __m128i xmm0 = _mm_loadu_si128( (__m128i*)((char*)src - 16) );

          if (size > 16) {
            if (size > 32) {
              size = -size;
              const __m128i xmm1 = _mm_loadu_si128( (__m128i*)((char*)src - 32) );
              const __m128i xmm6 = _mm_loadu_si128( (__m128i*)((char*)src + size + 16) );
              const __m128i xmm7 = _mm_loadu_si128( (__m128i*)((char*)src + size) );
              _mm_stream_si128( (__m128i*)((char*)dst - 16), xmm0 );
              _mm_stream_si128( (__m128i*)((char*)dst - 32), xmm1 );
              _mm_storeu_si128( (__m128i*)((char*)dst + size + 16), xmm6 );
              _mm_storeu_si128( (__m128i*)((char*)dst + size), xmm7 );
              return ret;
            }

            size = -size;
            const __m128i xmm7 = _mm_loadu_si128( (__m128i*)((char*)src + size) );
            _mm_stream_si128( (__m128i*)((char*)dst - 16), xmm0 );
            _mm_storeu_si128( (__m128i*)((char*)dst + size), xmm7 );
            return ret;
          }

          _mm_stream_si128( (__m128i*)((char*)dst - 16), xmm0 );
          return ret;
        }
      }

      if (size >= 8) {
        long long rax = *(long long*)((char*)src - 8);

        if (size > 8) {
          /*
           * that's right, we're converting an unsigned value to a negative,
           * saves 2 clock cycles!
           */
          size = -size;
          long long rcx = *(long long*)((char*)src + size);
          *(long long*)((char*)dst - 8) = rax;
          *(long long*)((char*)dst + size) = rcx;
        }
        else *(long long*)((char*)dst - 8) = rax;
      }
      else if (size >= 4) {
        int eax = *(int*)((char*)src - 4);

        if (size > 4) {
          size = -size;
          int ecx = *(int*)((char*)src + size);
          *(int*)((char*)dst - 4) = eax;
          *(int*)((char*)dst + size) = ecx;
        }
        else *(int*)((char*)dst - 4) = eax;
      }
      else if (size >= 1) {
        char al = *((char*)src - 1);

        if (size > 1) {
          size = -size;
          short cx = *(short*)((char*)src + size);
          *((char*)dst - 1) = al;
          *(short*)((char*)dst + size) = cx;
        }
        else *((char*)dst - 1) = al;
      }

      return ret;
    }
  }
}

void *apex_memmove_dispatcher(void *dst, const void *src, size_t size)
{
  apex_memcpy = apex_memmove = &apex_kryptonite;

  return apex_memmove(dst, src, size);
}


void *(*apex_memcpy)(void *dst, const void *src, size_t size) = apex_memmove_dispatcher;
void *(*apex_memmove)(void *dst, const void *src, size_t size) = apex_memmove_dispatcher;

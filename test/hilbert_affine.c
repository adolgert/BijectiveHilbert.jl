/*
 * hilbert_affine.c
 *
 * Anisotropic (activation) Hilbert encode/decode using the affine-map
 * formulation from lean/refine_affine_aniso.md.
 *
 * - Axes are labeled 0..n-1.
 * - At level s (MSB-first), active axes are those with m_j >= s, ordered by
 *   priority (m_j, j).
 * - The per-level state is the affine map S_{e,delta}(x) = rotL(delta) x XOR e
 *   on the active list (delta = d+1).
 * - Digits are variable-width (k_s bits) and are packed MSB-first into the
 *   Hilbert index.
 *
 * Coordinate type is uint32_t, so each m_j must be in [0, 32].
 * Hilbert indices use __uint128_t; sum(m_j) must be <= 128.
 */

#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#define MAX_DIMS 32
#define MAX_LEVELS 32
#define MAX_INDEX_BITS 128

typedef __uint128_t hindex_t;
typedef uint32_t coord_t;

typedef struct {
  uint32_t e; /* entry mask in the current active list */
  uint32_t d; /* direction index in the current active list */
} hilbert_state_t;

static inline uint32_t mask_bits(uint32_t bits) {
  return (bits >= 32u) ? 0xFFFFFFFFu : ((1u << bits) - 1u);
}

static inline uint32_t rotl_bits(uint32_t x, uint32_t r, uint32_t bits) {
  if (bits == 0u) return x;
  if (bits == 32u) {
    r &= 31u;
    return (r == 0u) ? x : (uint32_t)((x << r) | (x >> (32u - r)));
  }
  const uint32_t mask = mask_bits(bits);
  x &= mask;
  r %= bits;
  if (r == 0u) return x;
  return (uint32_t)(((x << r) | (x >> (bits - r))) & mask);
}

static inline uint32_t rotr_bits(uint32_t x, uint32_t r, uint32_t bits) {
  if (bits == 0u) return x;
  if (bits == 32u) {
    r &= 31u;
    return (r == 0u) ? x : (uint32_t)((x >> r) | (x << (32u - r)));
  }
  const uint32_t mask = mask_bits(bits);
  x &= mask;
  r %= bits;
  if (r == 0u) return x;
  return (uint32_t)(((x >> r) | (x << (bits - r))) & mask);
}

static inline uint32_t gray_code(uint32_t x) {
  return x ^ (x >> 1);
}

static inline uint32_t gray_decode(uint32_t g) {
  uint32_t x = g;
  x ^= x >> 1;
  x ^= x >> 2;
  x ^= x >> 4;
  x ^= x >> 8;
  x ^= x >> 16;
  return x;
}

static inline uint32_t trailing_ones(uint32_t x) {
  uint32_t c = 0;
  while ((x & 1u) != 0u) {
    c++;
    x >>= 1;
  }
  return c;
}

/* Hamilton entry sequence e(w) for a k-dimensional cube. */
static inline uint32_t child_entry(uint32_t w) {
  if (w == 0u) return 0u;
  return gray_code((w - 1u) & ~1u);
}

/* Hamilton direction sequence d(w) for a k-dimensional cube. */
static inline uint32_t child_dir(uint32_t w, uint32_t k) {
  if (w == 0u) return 0u;
  if ((w & 1u) != 0u) return trailing_ones(w) % k;
  return trailing_ones(w - 1u) % k;
}

/* S_{e,delta}(x) = rotL(delta) x XOR e, with delta = d+1. */
static inline uint32_t affine_apply(uint32_t x, uint32_t e, uint32_t d, uint32_t k) {
  return (rotl_bits(x, d + 1u, k) ^ e) & mask_bits(k);
}

/* S^{-1}(y) = rotR(delta) (y XOR e), with delta = d+1. */
static inline uint32_t affine_apply_inv(uint32_t y, uint32_t e, uint32_t d, uint32_t k) {
  return rotr_bits((y ^ e), d + 1u, k) & mask_bits(k);
}

typedef struct {
  int axis;
  int exp;
} axis_exp_t;

/* Sort axes by priority (m_j, j). */
static void sort_axes_by_priority(const int *m, int n, int *order) {
  axis_exp_t items[MAX_DIMS];
  for (int i = 0; i < n; i++) {
    items[i].axis = i;
    items[i].exp = m[i];
  }
  for (int i = 1; i < n; i++) {
    axis_exp_t key = items[i];
    int j = i - 1;
    while (j >= 0 && (items[j].exp > key.exp ||
           (items[j].exp == key.exp && items[j].axis > key.axis))) {
      items[j + 1] = items[j];
      j--;
    }
    items[j + 1] = key;
  }
  for (int i = 0; i < n; i++) {
    order[i] = items[i].axis;
  }
}

static bool build_active_axes(const int *m, int n,
                              int *mmax_out, int *total_bits_out,
                              int axes_level[MAX_LEVELS + 1][MAX_DIMS],
                              int k_level[MAX_LEVELS + 1],
                              int pos_level[MAX_LEVELS + 1][MAX_DIMS]) {
  if (!m || n <= 0 || n > MAX_DIMS) return false;

  int mmax = 0;
  int total_bits = 0;
  for (int i = 0; i < n; i++) {
    if (m[i] < 0 || m[i] > MAX_LEVELS) return false;
    if (m[i] > mmax) mmax = m[i];
    total_bits += m[i];
  }
  if (total_bits > MAX_INDEX_BITS) return false;

  int order[MAX_DIMS];
  sort_axes_by_priority(m, n, order);

  k_level[0] = 0;
  for (int ax = 0; ax < n; ax++) pos_level[0][ax] = -1;

  for (int s = 1; s <= mmax; s++) {
    int k = 0;
    for (int i = 0; i < n; i++) {
      int ax = order[i];
      if (m[ax] >= s) {
        axes_level[s][k++] = ax;
      }
    }
    k_level[s] = k;
    for (int ax = 0; ax < n; ax++) pos_level[s][ax] = -1;
    for (int j = 0; j < k; j++) {
      pos_level[s][axes_level[s][j]] = j;
    }
  }

  if (mmax_out) *mmax_out = mmax;
  if (total_bits_out) *total_bits_out = total_bits;
  return true;
}

static void embed_state(const int *A_old, int k_old, const int *pos_new,
                        uint32_t e_old, uint32_t d_old,
                        uint32_t *e_new, uint32_t *d_new) {
  uint32_t e = 0u;
  for (int j = 0; j < k_old; j++) {
    if ((e_old >> j) & 1u) {
      int new_pos = pos_new[A_old[j]];
      e |= 1u << new_pos;
    }
  }
  int dir_axis = A_old[d_old];
  *e_new = e;
  *d_new = (uint32_t)pos_new[dir_axis];
}

int hilbert_affine_index_bits(const int *m, int n) {
  if (!m || n <= 0 || n > MAX_DIMS) return 0;
  int total = 0;
  for (int i = 0; i < n; i++) {
    if (m[i] < 0) return 0;
    total += m[i];
  }
  return total;
}

hindex_t hilbert_affine_encode(const coord_t *point, const int *m, int n) {
  int mmax = 0;
  int total_bits = 0;
  int axes_level[MAX_LEVELS + 1][MAX_DIMS];
  int pos_level[MAX_LEVELS + 1][MAX_DIMS];
  int k_level[MAX_LEVELS + 1];

  if (!point || !m) return (hindex_t)0;
  if (!build_active_axes(m, n, &mmax, &total_bits, axes_level, k_level, pos_level)) {
    return (hindex_t)0;
  }
  if (mmax == 0) return (hindex_t)0;

  hilbert_state_t st = {0u, 0u};
  hindex_t h = 0;

  for (int s = mmax; s >= 1; s--) {
    const int *A = axes_level[s];
    int k = k_level[s];
    if (k == 0) return (hindex_t)0;

    const uint32_t mask = mask_bits((uint32_t)k);
    st.e &= mask;
    st.d %= (uint32_t)k;

    uint32_t plane = 0u;
    for (int j = 0; j < k; j++) {
      int ax = A[j];
      plane |= ((point[ax] >> (s - 1)) & 1u) << j;
    }
    plane &= mask;

    uint32_t pre = affine_apply_inv(plane, st.e, st.d, (uint32_t)k);
    uint32_t w = gray_decode(pre) & mask;

    h = (h << k) | (hindex_t)w;

    uint32_t entry = child_entry(w) & mask;
    st.e = (st.e ^ rotl_bits(entry, st.d + 1u, (uint32_t)k)) & mask;
    st.d = (st.d + child_dir(w, (uint32_t)k) + 1u) % (uint32_t)k;

    if (s > 1 && k_level[s - 1] > k) {
      uint32_t e_new = 0u;
      uint32_t d_new = 0u;
      embed_state(A, k, pos_level[s - 1], st.e, st.d, &e_new, &d_new);
      st.e = e_new;
      st.d = d_new;
    }
  }

  return h;
}

void hilbert_affine_decode(hindex_t h, const int *m, int n, coord_t *point) {
  if (!point || !m || n <= 0 || n > MAX_DIMS) return;
  memset(point, 0, (size_t)n * sizeof(coord_t));

  int mmax = 0;
  int total_bits = 0;
  int axes_level[MAX_LEVELS + 1][MAX_DIMS];
  int pos_level[MAX_LEVELS + 1][MAX_DIMS];
  int k_level[MAX_LEVELS + 1];

  if (!build_active_axes(m, n, &mmax, &total_bits, axes_level, k_level, pos_level)) {
    return;
  }
  if (mmax == 0) return;

  int bit_pos = total_bits;
  hilbert_state_t st = {0u, 0u};

  for (int s = mmax; s >= 1; s--) {
    const int *A = axes_level[s];
    int k = k_level[s];
    if (k == 0) return;

    const uint32_t mask = mask_bits((uint32_t)k);
    st.e &= mask;
    st.d %= (uint32_t)k;

    bit_pos -= k;
    uint32_t w = (uint32_t)((h >> bit_pos) & (hindex_t)mask);

    uint32_t g = gray_code(w) & mask;
    uint32_t plane = affine_apply(g, st.e, st.d, (uint32_t)k);

    for (int j = 0; j < k; j++) {
      int ax = A[j];
      point[ax] |= ((plane >> j) & 1u) << (s - 1);
    }

    uint32_t entry = child_entry(w) & mask;
    st.e = (st.e ^ rotl_bits(entry, st.d + 1u, (uint32_t)k)) & mask;
    st.d = (st.d + child_dir(w, (uint32_t)k) + 1u) % (uint32_t)k;

    if (s > 1 && k_level[s - 1] > k) {
      uint32_t e_new = 0u;
      uint32_t d_new = 0u;
      embed_state(A, k, pos_level[s - 1], st.e, st.d, &e_new, &d_new);
      st.e = e_new;
      st.d = d_new;
    }
  }
}

uint64_t hilbert_affine_encode_64(const coord_t *point, const int *m, int n) {
  return (uint64_t)hilbert_affine_encode(point, m, n);
}

void hilbert_affine_decode_64(uint64_t h, const int *m, int n, coord_t *point) {
  hilbert_affine_decode((hindex_t)h, m, n, point);
}

void hilbert_affine_encode_128(const coord_t *point, const int *m, int n,
                               uint64_t *h_lo, uint64_t *h_hi) {
  hindex_t h = hilbert_affine_encode(point, m, n);
  *h_lo = (uint64_t)h;
  *h_hi = (uint64_t)(h >> 64);
}

void hilbert_affine_decode_128(uint64_t h_lo, uint64_t h_hi, const int *m, int n,
                               coord_t *point) {
  hindex_t h = ((hindex_t)h_hi << 64) | h_lo;
  hilbert_affine_decode(h, m, n, point);
}

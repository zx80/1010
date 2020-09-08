/*
 * $Id: 1010.c 144 2016-08-19 07:14:26Z coelho $
 *
 * COPYRIGHT
 *   (c) 2016 Fabien Coelho <1010 dot bang at coelho dot net>
 *
 * LICENSE
 *   This code is free software.
 *   It is licensed under the GNU General Public License v3 or better.
 *   For details about your rights: https://www.gnu.org/licenses/gpl-3.0.html
 *   Summary: use at your own risks.
 *
 * CAVEATS
 *   All software have bugs, this is a software, it must have some bugs.
 *   Beware, you may loose your friends or your hairs because of this software.
 *
 * COMPILATION
 *   sh> cc -DPARALLEL=6 -Ofast -march=native -o 1010 1010.c -lpthread -lm
 *    - If you have plenty of cores and you want to reduce latency further,
 *      consider -DPARALLEL=12
 *    - For better, althought not reproducible random, enable -DMY_RANDOM
 *    - For using one random call per round instead of per piece -DONE_RANDOM
 *    - For also trying to optimize for higher score, try -DMAX_SCORE
 *
 * DOCUMENTATION
 *   Read the source code below:-)
 *
 */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

// #define PARALLEL 12
// #undef PARALLEL

#ifdef PARALLEL
#include <pthread.h>
#ifndef PARALLEL_DEGREE
#define PARALLEL_DEGREE ((PARALLEL + 5) / 6)
#endif // PARALLEL_DEGREE
#else // ! PARALLEL
#undef PARALLEL_DEGREE
#endif // PARALLEL

#define VERSION_ID "$Id: 1010.c 144 2016-08-19 07:14:26Z coelho $"

typedef enum { false, true } bool;

/******************************************************************** RANDOM */

#ifdef MY_RANDOM
#define WITH_MY_RANDOM "yes"

/* local random implementation: probably much better than random, but
 * runs are not reproducible, which is a drawback.
 */

#include <unistd.h>

// could be something else, but 1-3 random numbers are drawn per round.
#ifndef RANDOM_SOURCE
#define RANDOM_SOURCE "/dev/urandom"
#endif // RANDOM_SOURCE

static FILE * random_source = NULL;

// random initialization (must be called)
static void my_srandom(unsigned int seed)
{
  // forward seeding
  srandom(seed);
  // use external source
  fprintf(stdout, "using \"%s\" as random source\n", RANDOM_SOURCE);
  random_source = fopen(RANDOM_SOURCE, "r");
  assert(random_source != NULL);
}

// return a positive integer
static long int my_random(void)
{
  assert(random_source != NULL);
  int32_t r;
  ssize_t nr = read(fileno(random_source), &r, sizeof(r));
  assert(nr == sizeof(r));
  // combine external source & random sequential output
  return (r ^ random()) & 0x7fffffff;
}

// on the fly local random & srandom replacement.
#define random my_random
#define srandom my_srandom

#else // ! MY_RANDOM
#define WITH_MY_RANDOM "no"
#endif // MY_RANDOM

/************************************************************************ T7 */

// gcc seems to have 128-bit integers (unsigned __int128),
// but the doc is scarce, eg how to write constants for instance.

typedef struct {
  uint64_t f1, f2;   // two-part 128 bit bitfield
  uint32_t score;    // current score
  uint32_t points;   // awarded points on last round
} t7_t;

#define t7_rshifteq(i, s)                                           \
  i.f2 = s == 0 ? i.f2 :                                            \
    s < 64 ? i.f2 >> (s) | (i.f1 << (64-(s))) : i.f1 >> ((s) - 64), \
    i.f1 = s == 0 ? i.f1 : s < 64 ? i.f1 >> (s) : 0
#define t7_oreq(i, j) i.f1 |= j.f1, i.f2 |= j.f2
#define t7_xoreq(i, j) i.f1 ^= j.f1, i.f2 ^= j.f2
#define t7_andeq(i, j) i.f1 &= j.f1, i.f2 &= j.f2
#define t7_noteq(i) i.f1 = ~i.f1, i.f2 = ~i.f2
#define t7_free(i, j) ((i.f1 & j.f1) == 0 && (i.f2 & j.f2) == 0)
#define t7_in(m, i) ((i.f1 & m.f1) == m.f1 && (i.f2 & m.f2) == m.f2)

#define T7(f1, f2) (t7_t) { f1, f2, 0, 0 }

#define t7_ZER T7(0x0000000000000000L, 0x0000000000000000L)
#define t7_DOT T7(0x8000000000000000L, 0x0000000000000000L)
#define t7_TOP T7(0xffc0000000000000L, 0x0000000000000000L)
#define t7_CBT T7(0x003fffffffffffffL, 0xfffffffff0000000L) // center-bottom
#define t7_RO1 T7(0x003ff00000000000L, 0x0000000000000000L)
#define t7_RO2 T7(0x00000ffc00000000L, 0x0000000000000000L)
#define t7_RO3 T7(0x00000003ff000000L, 0x0000000000000000L)
#define t7_RO4 T7(0x0000000000ffc000L, 0x0000000000000000L)
#define t7_RO5 T7(0x0000000000003ff0L, 0x0000000000000000L)
#define t7_RO6 T7(0x000000000000000fL, 0xfc00000000000000L)
#define t7_RO7 T7(0x0000000000000000L, 0x03ff000000000000L)
#define t7_RO8 T7(0x0000000000000000L, 0x0000ffc000000000L)
#define t7_BOT T7(0x0000000000000000L, 0x0000003ff0000000L)
#define t7_LEF T7(0x8020080200802008L, 0x0200802000000000L)
#define t7_CRI T7(0x7fdff7fdff7fdff7L, 0xfdff7fdff0000000L) // center-right
#define t7_CO1 T7(0x4010040100401004L, 0x0100401000000000L)
#define t7_CO2 T7(0x2008020080200802L, 0x0080200800000000L)
#define t7_CO3 T7(0x1004010040100401L, 0x0040100400000000L)
#define t7_CO4 T7(0x0802008020080200L, 0x8020080200000000L)
#define t7_CO5 T7(0x0401004010040100L, 0x4010040100000000L)
#define t7_CO6 T7(0x0200802008020080L, 0x2008020080000000L)
#define t7_CO7 T7(0x0100401004010040L, 0x1004010040000000L)
#define t7_CO8 T7(0x0080200802008020L, 0x0802008020000000L)
#define t7_RIT T7(0x0040100401004010L, 0x0401004010000000L)
#define t7_100 T7(0xffffffffffffffffL, 0xfffffffff0000000L)

typedef struct {
  char name;      // one-char name of patter
  int8_t w, x, y; // weight, x and y sizes (could be computed with recent C?)
  t7_t pat;       // pattern definition (score & points unused)
} pattern_t;

#define Pattern(n, f1, w, x, y) { n, w, x, y, T7(f1, 0x0000000000000000L) }

#define npats 20
static const pattern_t PAT[npats] = {
  // standard pieces
  Pattern('.', 0x8000000000000000L, 1, 1, 1),
  Pattern('i', 0x8020000000000000L, 2, 1, 2),
  Pattern('-', 0xc000000000000000L, 2, 2, 1),
  Pattern('I', 0x8020080000000000L, 3, 1, 3),
  Pattern('_', 0xe000000000000000L, 3, 3, 1),
  Pattern('r', 0xc020000000000000L, 3, 2, 2),
  Pattern('l', 0x8030000000000000L, 3, 2, 2),
  Pattern('j', 0x4030000000000000L, 3, 2, 2),
  Pattern('t', 0xc010000000000000L, 3, 2, 2),
  Pattern('o', 0xc030000000000000L, 4, 2, 2),
  Pattern('v', 0x8020080200000000L, 4, 1, 4),
  Pattern('h', 0xf000000000000000L, 4, 4, 1),
  Pattern('V', 0x8020080200800000L, 5, 1, 5),
  Pattern('H', 0xf800000000000000L, 5, 5, 1),
  Pattern('R', 0xe020080000000000L, 5, 3, 3),
  Pattern('L', 0x80200e0000000000L, 5, 3, 3),
  Pattern('J', 0x20080e0000000000L, 5, 3, 3),
  Pattern('T', 0xe008020000000000L, 5, 3, 3),
  Pattern('O', 0xe0380e0000000000L, 9, 3, 3),
  // non existing large square
  Pattern('X', 0xf83e0f83e0f80000L, 25, 5, 5)
};

/* return pattern identified by its one character name
 */
static const pattern_t * get_pattern(char c)
{
  static bool cache_initialized = false;
  static const pattern_t * cache[128];
  if (!cache_initialized)
  {
    for (int i = 0; i < 128; i++)
      cache[i] = NULL;
    for (int i = 0; i < npats; i++)
      cache[(int) PAT[i].name] = &PAT[i];
    cache_initialized = true;
  }
  // assert(0 <= c && (int) c < 128);
  assert(0 <= c);
  return cache[(int) c];
}

/* count bits in f
 * complexity: 3 ops
 * plenty algorithms are available, or hardware instructions
 */
static int t7_count(t7_t f)
{
  return __builtin_popcountll(f.f1) + __builtin_popcountll(f.f2);
}

#define naligns 20
static t7_t ALIGNS[naligns] = {
  t7_TOP, t7_RO1, t7_RO2, t7_RO3, t7_RO4,
  t7_RO5, t7_RO6, t7_RO7, t7_RO8, t7_BOT,
  t7_LEF, t7_CO1, t7_CO2, t7_CO3, t7_CO4,
  t7_CO5, t7_CO6, t7_CO7, t7_CO8, t7_RIT
};

/* remove full rows & columns + update the score
 * complexity: 3 2n bitops
 */
static t7_t simplify(t7_t f)
{
  t7_t aligns = t7_ZER;
  int n = 0;

  // first check for completed alignements
  for (int i = 0; i < naligns; i++)
    if (t7_in(ALIGNS[i], f))
      n++, t7_oreq(aligns, ALIGNS[i]);

  // then remove if needed
  if (n != 0)
  {
    f.points += 5 * n * (n+1);
    t7_noteq(aligns);
    t7_andeq(f, aligns);
  }

  return f;
}

/* show field, short
 */
static void t7_dump(FILE * out, const char * name, t7_t i)
{
  if (name)
    fprintf(out, "%s: ", name);
  fprintf(out, "%d 0x%016lx%016lx %d (%d)\n",
          t7_count(i), i.f1, i.f2, i.score, i.points);
}

/* show field, nicer
 */
static void t7_show(FILE * out, const char * name, t7_t i)
{
  t7_dump(out, name, i);
  fputs("/  0 1 2 3 4 5 6 7 8 9\n", out);
  for (int r = 0; r < 10 ; r++)
  {
    putc('0' + r, out);
    putc(' ', out);
    for (int c = 0; c < 10 ; c++)
    {
      t7_t p = t7_DOT;
      t7_rshifteq(p, 10*r+c);
      putc(' ', out);
      putc(t7_in(p, i)? 'x': '.', out);
    }
    putc('\n', out);
  }
}

/* return whether pattern p may be put anywhere in field f
 * complexity: n² (worst case)
 */
static bool can_put(t7_t f, const pattern_t * p)
{
  for (int r = 0; r <= 10 - p->y; r++)
    for (int c = 0; c <= 10 - p->x; c++)
      if (t7_free(p->pat, f))
        return true;
  return false;
}

/* compute filling surface compaction
 * complexity: initial 4 n², bitops ~ 14
 */
static double compact(t7_t f)
{
  t7_t t1 = f, t2 = f;
  t7_rshifteq(t1, 1);
  t7_andeq(t2, t7_LEF);
  t7_andeq(t1, t7_CRI);
  t7_oreq(t1, t2);
  t7_xoreq(t1, f);
  int surface = t7_count(t1);

  t1 = f, t2 = f;
  t7_rshifteq(t1, 10);
  t7_andeq(t2, t7_TOP);
  t7_andeq(t1, t7_CBT);
  t7_oreq(t1, t2);
  t7_xoreq(t1, f);
  surface += t7_count(t1);

  return (180.0 - surface) / 180.0; // normalize
}

/* compute alignment score
 * complexity: initial n², bitops 5 2n
 */
static double alignments(t7_t f)
{
  double a = 0.0;

  for (int i = 0; i < naligns; i++)
  {
    t7_t t = f;
    t7_andeq(t, ALIGNS[i]);
    int n = t7_count(t);
    a += n * n;
  }

  return a / (9.0 * 9.0 * naligns); // normalize
}

/* large pieces, including special 5x5
 */
static const pattern_t * O, * H, * V, * X;

static void init_patterns(void)
{
  O = get_pattern('O');
  H = get_pattern('H');
  V = get_pattern('V');
  X = get_pattern('X');
}

/* default decision weights */
static double
  w_free = 1.0,
  w_X = 0.0,
  w_OVH  = 0.0,
  w_a  = 2.0,
  w_c  = 5.0;

#define npieces 3
typedef struct {
  double eval;       // evaluation of move
  char c[npieces];   // pieces to be placed, fixed order
  int move[npieces]; // chosen position for each piece
  t7_t field;        // initial then final field
} move_t;

#ifdef PARALLEL
static int partial = PARALLEL_DEGREE;
#else
#define partial 1
#endif // PARALLEL

// stats: count the number of evaluations on each round
static uint64_t eval_count = 0;

// very small to avoid 0.0
#define epsilon (1.0E-308)

/* try to place recursively 3 pieces at every possible places
 * worst case complexity: (n²)³ n² = n^8
 */
static move_t rec_play(t7_t f, int i, move_t move, move_t best,
                       int every, uint32_t * nevals)
{
  assert(i <= 3);
  const pattern_t * pc = get_pattern(move.c[i]);
  assert(pc != NULL);

  for (int r = 0; r <= 10 - pc->y; r++)
  {
    for (int c = 0; c <= 10 - pc->x; c++)
    {
      int m = 10 * r + c;

#ifdef PARALLEL

      // share work between threads on first level
      if (i == 0 && m % partial != every)
        continue;

#endif // PARALLEL

      move.move[i] = m;
      t7_t fc = pc->pat;
      t7_rshifteq(fc, m);
      if (t7_free(fc, f))
      {
        t7_t nf = f;
        nf.points += pc->w;
        t7_oreq(nf, fc);
        nf = simplify(nf);
        if (i < 2)
          // recursion
          best = rec_play(nf, i+1, move, best, 0, nevals);
        else
        {
          // final evaluation
          (*nevals) ++;

          move.eval =
#ifdef MAX_SCORE
#define WITH_MAX_SCORE "yes"
            // choose better score
            // round points max <= (9+5*6*(6+1))*3 == 657, in practice < 200
            0.0000001 * nf.points +
#else // ! MAX_SCORE
#define WITH_MAX_SCORE "no"
            // for 00000: evaluated is better than nothing
            epsilon +
#endif // MAX_SCORE
            // free space, integer number of cells
            w_free * (100.0 - t7_count(nf)) +
            // check for an artificial 5x5 piece, X * 0/1
            (w_X > 0.0 ? w_X * can_put(nf, X): 0.0) +
            // check space for 3 equaly probable large pieces, OVH * 0/1/2/3
            (w_OVH > 0.0 ?
             w_OVH * (can_put(nf, O) + can_put(nf, H) + can_put(nf, V)) : 0.0) +
            // the more alignments the better, a * [0,1] (1/1620)
            w_a * alignments(nf) +
            // compaction, to avoid holes c * [0,1] (1/180)
            w_c * compact(nf);

          if (move.eval > best.eval)
          {
            best = move;
            best.field = nf;
          }
        }
      }
    }
  }

  return best;
}

typedef struct
{
  pthread_t tid;     // thread identifier
  int every;         // work selector for thread
  move_t move;       // best move
  uint32_t nevals;   // local counter
} search_move_t;

static void * search_move(void * p)
{
  search_move_t * pm = (search_move_t *) p;
  pm->move = rec_play(pm->move.field, 0, pm->move, pm->move,
                      pm->every, &pm->nevals);
  return NULL;
}

#define MOVE(p, i0, i1, i2, f)                            \
  { 0.0, { p[i0], p[i1], p[i2] }, { -1, -1, -1 }, f }

/* return best move from field & provided pieces
 */
static move_t play(FILE * out, t7_t f, char pieces[3])
{
  move_t move[6] = {
    // try every order
    MOVE(pieces, 0, 1, 2, f), MOVE(pieces, 0, 2, 1, f),
    MOVE(pieces, 1, 0, 2, f), MOVE(pieces, 1, 2, 0, f),
    MOVE(pieces, 2, 0, 1, f), MOVE(pieces, 2, 1, 0, f)
  };

  // do not run duplicates
  bool run[6] = { true, true, true, true, true, true };
  for (int i = 1; i < 6; i++)
  {
    // check whether it is equal to one before
    for (int j = 0; j < i && run[i]; j++)
    {
      if (move[i].c[0] == move[j].c[0] &&
          move[i].c[1] == move[j].c[1] &&
          move[i].c[2] == move[j].c[2])
        run[i] = false;
    }
  }

  search_move_t search[6 * partial];

  fputs("play: ", out);

#ifdef PARALLEL

  fputs("threaded version", out);

  for (int t = 0; t < 6; t++)
  {
    if (run[t])
    {
      for (int p = 0; p < partial; p++)
      {
        int i = partial * t + p;
        search[i] = (search_move_t) { .every = p, .move = move[t], .nevals = 0 };
        int e = pthread_create(&search[i].tid, NULL, search_move, &search[i]);
        assert(e == 0);
      }
    }
  }

  for (int t = 0; t < 6; t++)
  {
    if (run[t])
    {
      for (int p = 0; p < partial; p++)
      {
        int i = partial * t + p ;
        int e = pthread_join(search[i].tid, NULL);
        assert(e == 0);
        eval_count += search[i].nevals;
        fputc('.', out);
      }
    }
  }

  fputs(" done.", out);

#else // ! PARALLEL

  fputs("single thread version", out);
  assert(partial == 1);

  for (int i = 0; i < 6; i++)
  {
    if (run[i])
    {
      search[i] = (search_move_t) { .every = 0, .move = move[i] };
      search_move(&search[i]);
      eval_count += search[i].nevals;
      fputc('.', out);
    }
  }

#endif // PARALLEL

  fputc('\n', out);

  // find best ply
  move_t best;
  for (int t = 0; t < 6; t++)
  {
    if (run[t])
    {
      for (int p = 0; p < partial; p++)
      {
        int i = t * partial + p;
        if (i == 0 || search[i].move.eval > best.eval)
          best = search[i].move;
      }
    }
  }

  return best;
}

/* interactive play (the computer plays:-)
 */
static void do_play(FILE * out, __attribute__((unused)) FILE * err, t7_t init)
{
  t7_t f1010 = init, previous = init;
  char pieces[3];
  while (true)
  {
    t7_show(out, "current", f1010);
    fprintf(out, "3 pieces (.i-I_rtjlhvoHVRTJLO, bbb=back): ");
    int i = 0;
    while (i < 3)
      if (fscanf(stdin, "%c", &pieces[i]) == 1 &&
          (get_pattern(pieces[i]) != NULL || pieces[i] == 'b'))
        i++;
    if (pieces[0] == 'b')
      f1010 = previous;
    else
    {
      f1010.points = 0;
      move_t best = play(out, f1010, pieces);
      // move_t best = play_once(f1010, pieces);
      fprintf(out, "best %c%d %c%d %c%d: %.2f\n",
              best.c[0], best.move[0],
              best.c[1], best.move[1],
              best.c[2], best.move[2],
              best.eval
        );
      previous = f1010;
      f1010 = best.field;
      f1010.score += f1010.points;
    }
  }
}

/* return time as a double, for easy arithmetics
 */
static double dtime(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 0.0000001;
}

// probable intended distribution
#define RAND_1010 "..iii---III___rrttjjllhhvvooooooHHVVRTJLOO"

// biased distribution (one j missing)
// #define RAND_1010 "..iii---III___rrttjllhhvvooooooHHVVRTJLOO"

// flat
// #define RAND_1010 ".i-I_rtjlhvoHVRTJLO"

// bad, should not last much...
// #define RAND_1010 "OVH"

/* play alone
 */
static void auto_play(FILE * out, FILE * err, t7_t init)
{
  t7_t f1010 = init;
  const int len = strlen(RAND_1010);

  // init stats
  const double start = dtime();
  int rounds = 0;
  int64_t cel_sum = 0;
  int64_t cel_sum2 = 0;
  int64_t inputs_sum = 0;
  int maxcels = 0;
  int nzeros = 0;
  int distrib[91];  // worst case is full minus diagonal
  for (int i = 0; i < 91; i++)
    distrib[i] = 0;

  while (true)
  {
    t7_show(out, "current", f1010);
    char pc[3];

#ifdef ONE_RANDOM
#define WITH_ONE_RANDOM "yes"
    int r = random();
    pc[0] = RAND_1010[r % len];
    r /= len;
    pc[1] = RAND_1010[r % len];
    r /= len;
    pc[2] = RAND_1010[r % len];
#else // ! ONE_RANDOM
#define WITH_ONE_RANDOM "no"
    pc[0] = RAND_1010[random() % len];
    pc[1] = RAND_1010[random() % len];
    pc[2] = RAND_1010[random() % len];
#endif // ONE_RANDOM

    const int w = get_pattern(pc[0])->w + get_pattern(pc[1])->w +
      get_pattern(pc[2])->w;

    fprintf(err, "%d: trying %c%c%c %d\n", rounds, pc[0], pc[1], pc[2], w);

    f1010.points = 0;
    const move_t best = play(out, f1010, pc);
    if (best.eval <= 0.0)
      break;

    fprintf(out, "best %c%d %c%d %c%d: %.3f\n",
            best.c[0], best.move[0],
            best.c[1], best.move[1],
            best.c[2], best.move[2],
            best.eval
      );

    rounds ++;
    f1010 = best.field;
    f1010.score += f1010.points;

    // collect & show stats
#define A 18
    const int cels = t7_count(f1010) - A;
    cel_sum += cels;
    cel_sum2 += cels * cels;
    inputs_sum += w;
    const int ncels = t7_count(f1010);
    if (maxcels < ncels)
      maxcels = ncels;
    if (ncels == 0)
      nzeros ++;
    distrib[ncels] ++;

    fprintf(err,
            "per round: cels=%.2f+-%.2f inputs=%.2f score=%.2f "
            "max=%d zeros=%d evals=%.0f %.3fs",
            A + (1.0 * cel_sum / rounds),
            sqrt((cel_sum2 - (cel_sum * cel_sum) / rounds) / rounds),
            1.0 * inputs_sum / rounds,
            1.0 * f1010.score / rounds,
            maxcels, nzeros, 1.0 * eval_count / rounds,
            (dtime() - start) / rounds);

    fputc('\n', err);
  }

  // final report
  fprintf(out, "stopped after %d rounds\n", rounds);
  fprintf(err, "distribution: ");
  int last_zero = 91;
  while (distrib[last_zero-1] == 0)
    last_zero--;
  // five numbers
  int cum = 0;
  bool x_min = true, x_q1 = true, x_med = true, x_q3 = true;
  for (int i = 0; i < last_zero; i++)
  {
    cum += distrib[i];
    if (x_min && cum != 0)
      fprintf(out, "[ %d ", i), x_min = false;
    else if (x_q1 && cum >= rounds/4)
      fprintf(out, "%d ", i), x_q1 = false;
    else if (x_med && cum >= rounds/2)
      fprintf(out, "%d ", i), x_med = false;
    else if (x_q3 && cum >= rounds*3/4)
      fprintf(out, "%d ", i), x_q3 = false;
    else if (cum == rounds)
      fprintf(out, "%d ]\n", i);
  }
  // detailed distribution
  for (int i = 0; i < last_zero; i++)
    fprintf(err, "%2d: %d\n", i, distrib[i]);
}

/* put pieces on any position */
static t7_t free_play(FILE * out, FILE * err)
{
  t7_t f1010 = t7_ZER;

  while (true)
  {
    t7_show(out, "current", f1010);
    fputs(".i-I_rljtovhVHRLJTO + pos or q\n", out);
    char m;
    int s;
    int n = fscanf(stdin, "%c%u", &m, &s);
    if (n >= 1 && m == 'q')
      return f1010;
    if (n == 2)
    {
      const pattern_t * p = get_pattern(m);
      if (p != NULL)
      {
        const int r = s / 10, c = s % 10;
        if (0 <= r && r <= 10-p->y &&
            0 <= c && c <= 10-p->x)
        {
          t7_t pat = p->pat;
          t7_rshifteq(pat, s);
          if (t7_free(pat, f1010))
          {
            t7_oreq(f1010, pat);
            f1010 = simplify(f1010);
          }
          else
            fprintf(err, "cannot put %c at %d: occupied\n", m, s);
        }
        else
            fprintf(err, "cannot put %c at %d: illegal\n", m, s);
      }
      else
        fprintf(err, "pattern %c not found\n", m);
    }
    else
      fprintf(err, "unexpected input, please retry\n");
  }

  return f1010;
}

// hmmm...
static char * version_id(void)
{
  static char version[20];
  char * v = version, * c = strchr(strchr(VERSION_ID, ' ') + 1, ' ') + 1;
  while (*c != ' ')
    *v++ = *c++;
  *v++ = '\0';
  return version;
}

int main(int argc, char * argv[])
{
  fprintf(stdout,
          "%s - 1010! simulator & assistant v%s\n"
          "Copyright: 2016 Fabien Coelho <1010 dot bang at coelho dot net>\n"
          "License: GPLv3+\n", argv[0], version_id());

  if (argc <= 1)
  {
    fprintf(stderr,
            "usage: %s (auto|play|free|fraut|fredo) [seed [weight]]\n",
            argv[0]);
    return 0;
  }

  init_patterns();

  unsigned int seed = argc >= 3 ? atoi(argv[2]): time(NULL);
  // seed==0 <=> seed==1 on my Linux implementation...
  if (seed == 0) seed = time(NULL);
  fprintf(stdout, "seeding random with %u one=" WITH_ONE_RANDOM "\n", seed);
  srandom(seed);

  // rough weight parameter...
  if (argc >= 4)
  {
#define BASE 36
    int w = strtol(argv[3], NULL, BASE);
    w_free = w / (BASE*BASE*BASE*BASE) % BASE;
    w_X = w / (BASE*BASE*BASE) % BASE;
    w_OVH = w / (BASE*BASE) % BASE;
    w_a = w / BASE % BASE;
    w_c = w % BASE;
  }

  fprintf(stdout,
          "weights: free=%.0f X=%.0f OVH=%.0f a=%.0f c=%.0f "
          "score=" WITH_MAX_SCORE "\n",
          w_free, w_X, w_OVH, w_a, w_c);

  if (strcmp(argv[1], "free")==0)
    free_play(stdout, stderr);
  else if (strcmp(argv[1], "auto")==0)
    auto_play(stdout, stdout, t7_ZER);
  else if (strcmp(argv[1], "play")==0)
    do_play(stdout, stderr, t7_ZER);
  else if (strcmp(argv[1], "fraut")==0)
  {
    t7_t f = free_play(stdout, stderr);
    auto_play(stdout, stdout, f);
  }
  else if (strcmp(argv[1], "fredo")==0)
  {
    t7_t f = free_play(stdout, stderr);
    do_play(stdout, stdout, f);
  }

  return 0;
}

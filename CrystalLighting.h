// C++11
#include <sys/time.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>

const bool USE_HALFSOLUTIONS = false;
const float FUTURE_FACTOR = 0.5f;

// ---- SHORTCUTS ----
using namespace std;

using VS = vector<string>;
using small = int16_t;  // signed

#define BE(X) std::begin(X), std::end(X)
#define FOR(VAR, FROM, OP_TO) for (auto VAR = (FROM); VAR OP_TO; ++(VAR))
#define FORBACK(VAR, FROM, OP_TO) for (auto VAR = (FROM); VAR OP_TO; --(VAR))
#define FORBE(VAR, X) for (auto(VAR) = begin(X); (VAR) != end(X); ++(VAR))
#define CHECK(x) assert(x)
// ----

using time_point = chrono::high_resolution_clock::time_point;
using chrono::duration;
time_point getTime(){
  return chrono::high_resolution_clock::now();
}

struct uninitialized_t
{};
constexpr const uninitialized_t uninitialized;

template <class T, int Capacity>
class InlineVector
{
public:
  using array_t = array<T, Capacity>;
  using iterator = typename array_t::iterator;
  using const_iterator = typename array_t::const_iterator;

  // Default constructor aggregate-initializes underlying array
  InlineVector() : a{} {}

  explicit InlineVector(uninitialized_t) {}

  InlineVector(int n, uninitialized_t) : s(n) { CHECK(n <= Capacity); }

  InlineVector(int n, const T& x) : s(n)
  {
    CHECK(n <= Capacity);
    for (int i = 0; i < n; ++i)
      a[i] = x;
  }

  explicit InlineVector(std::initializer_list<T> x) : s(x.size())
  {
    CHECK(x.size() <= Capacity);
    std::copy(BE(x), a.begin());
  }

  template <class C>
  void operator=(const C& x)
  {
    CHECK(x.size() <= Capacity);
    std::copy(BE(x), a.begin());
    s = x.size();
  }

  template <class C>
  void operator=(std::initializer_list<C> x)
  {
    CHECK(x.size() <= Capacity);
    std::copy(BE(x), a.begin());
    s = x.size();
  }

  template <class U, size_t N>
  void operator=(const std::array<U, N>& x)
  {
    static_assert(N <= Capacity, "");
    std::copy(BE(x), a.begin());
    s = x.size();
  }

  int size() const { return s; }
  constexpr int capacity() const { return Capacity; }
  bool empty() const { return s == 0; }
  T& operator[](int x)
  {
    assert(0 <= x && x < s);
    return a[x];
  }
  const T& operator[](int x) const
  {
    assert(0 <= x && x < s);
    return a[x];
  }
  T& front()
  {
    assert(s > 0);
    return a[0];
  }
  const T& front() const
  {
    assert(s > 0);
    return a[0];
  }
  T& back()
  {
    assert(s > 0);
    return a[s - 1];
  }
  const T& back() const
  {
    assert(s > 0);
    return a[s - 1];
  }

  iterator begin() { return a.begin(); }
  iterator end() { return a.begin() + s; }
  const_iterator begin() const { return a.begin(); }
  const_iterator end() const { return a.begin() + s; }
  void push_back(const T& x)
  {
    assert(s < Capacity);
    a[s++] = x;
  }
  void pop_back()
  {
    assert(s > 0);
    --s;
  }
  void erase(const_iterator it)
  {
    int idx = it - a.begin();
    CHECK(0 <= idx && idx < s);
    for (int i = idx; i + 1 < s; ++i)
      a[i] = std::move(a[i + 1]);
    --s;
  }

  void clear() { s = 0; }

  void resize(int i, uninitialized_t)
  {
    CHECK(0 <= i && i <= Capacity);
    s = i;
  }

  void resize(int i, const T& value = T())
  {
    CHECK(0 <= i && i <= Capacity);
    for (int j = s; j < i; ++j)
      a[j] = value;
    s = i;
  }

  auto data() { return a.data(); }
  auto data() const { return a.data(); }

private:
  int s = 0;
  array_t a;
};

void printBoard(const VS& b)
{
  for (auto& r : b) {
    FOR (i, 0, < b.size()) {
      printf("%c", r[i]);
    }
    printf("\n");
  }
  printf("\n");
}

struct Gain
{
  // 'realized' gain is, for a single move, includes
  // - cost of lanterns, incorretly lit crystals (negative)
  // - points from correctly lit crystals
  int realized = 0;
  // 'future' gain is an additional gain that can be realized later by completing the lighting of
  // crystals that has been partly lit in this move. It includes:
  // - reversing the penalty for incorrect lighting
  //- adding the points for the now fully lit crystals
  //-minimum cost of additional lanterns
  int future = 0;

  float weighted_combined() const { return realized + FUTURE_FACTOR * future; }
  int best_combined() const { return realized + future; }

  void set_int_min()
  {
    realized = INT_MIN;
    future = 0;
  }
};

int lantern_color_or_0(char c)
{
  switch (c) {
    case 'B':
      return 1;
    case 'Y':
      return 2;
    case 'R':
      return 4;
  }
  return 0;
}

int crystal_color_or_0(char c)
{
  if (!isdigit(c))
    return 0;
  int color = c - '0';
  assert(1 <= color && color <= 6);
  return color;
}

char make_lantern_char(int color)
{
  switch (color) {
    case 1:
      return 'B';
    case 2:
      return 'Y';
    case 4:
      return 'R';
  }
  assert(false);
  return 0;
}

struct RC
{
  small r, c;
};
bool operator==(RC x, RC y)
{
  return x.r == y.r && x.c == y.c;
}

bool within_table(RC rc, int H, int W)
{
  return 0 <= rc.r && rc.r < H && 0 <= rc.c && rc.c < W;
}

template <class T>
struct Matrix
{
  const int H, W;
  vector<T> d;  // data, row-major

  Matrix(int H, int W) : H(H), W(W), d(H * W) {}
  void valid_coords(RC rc) const { assert(within_table(rc, H, W)); }
  T& operator()(RC rc)
  {
    valid_coords(rc);
    return d[rc.r * W + rc.c];
  }
  const T& operator()(RC rc) const
  {
    valid_coords(rc);
    return d[rc.r * W + rc.c];
  }
};

enum class ColorsKind
{
  primary,
  secondary
};

ColorsKind colors_kind(int colors)
{
  assert(1 <= colors && colors <= 6);
  static const ColorsKind cks[6] = {ColorsKind::primary,   ColorsKind::primary,
                                    ColorsKind::secondary, ColorsKind::primary,
                                    ColorsKind::secondary, ColorsKind::secondary};
  return cks[colors - 1];
}

int num_colors(int colors)
{
  assert(0 <= colors && colors <= 7);
  static int ws[8] = {0, 1, 1, 2, 1, 2, 2, 3};
  return ws[colors];
}

struct Crystal
{
  Crystal(RC rc, int colors) : rc(rc), colors(colors), rcvd_colors(0) {}
  int missing_colors() const { return colors & (~rcvd_colors); }
  bool lit_incorrectly_unfixable() const { return (rcvd_colors & (~colors)) != 0; }
  bool lit_incorrectly_fixable() const
  {
    return !lit_incorrectly_unfixable() && missing_colors() != 0 && rcvd_colors != 0;
  }
  bool are_any_colors_not_needed_here(int colors_arg) const { return (colors_arg & ~colors) != 0; }
  bool do_colors_make_difference(int colors_arg)
  {
    if (lit_incorrectly_unfixable())
      return false;
    if (are_any_colors_not_needed_here(colors_arg))
      return true;
    return (missing_colors() & colors) != 0;
  }
  RC rc;
  small colors;
  small rcvd_colors;
};

struct Lantern
{
  Lantern() : color(0) {}
  Lantern(RC rc, int color) : rc(rc), color(color) {}
  RC rc;
  small color;
};

// 00 ROW+
// 01 ROW-
// 10 COL+
// 11 COL-
enum Dir
{
  DIR_DOWN = 0,
  DIR_UP,
  DIR_RIGHT,
  DIR_LEFT
};

// Bumping into / mirror: left, right, up, down = 3, 2, 1, 0 = 3 - dir
// Bumping into \ mirror: right, left, down, up = 2, 3, 0, 1 = dir ^ 2
int turn_at_slash_mirror(int dir)
{
  return 3 - dir;
}
int turn_at_backslash_mirror(int dir)
{
  return dir ^ 2;
}

const int NOTHING = -1;

struct CellX
{
  struct Towards
  {
    small target_crystal_ix = NOTHING;  // cyrstal idx | NOTHING
    small crystal_dist = NOTHING;       // Number of empty cells plus one between this
                                        // cell and the target crystal
    small lit_by_lantern_ix = NOTHING;  // lantern idx | NOTHING
    small lantern_dist = NOTHING;
  };
  array<Towards, 4> towards;
};

RC move_towards(RC rc, int dir)
{
  switch (dir) {
    case 0:
      ++rc.r;
      break;
    case 1:
      --rc.r;
      break;
    case 2:
      ++rc.c;
      break;
    case 3:
      --rc.c;
      break;
    default:
      assert(false);
  }
  return rc;
}

int opposite_dir(int dir)
{
  return dir ^ 1;
}

class RayEnumerator
{
  RC rc0;
  int initial_dir_counter;
  int curdir;
  RC currc;
  const Matrix<char>& board;
  int curdist;
  const int first_dir;

public:
  bool accumulate_lanterns = false;
  int accumulated_lantern_colors = 0;

  RayEnumerator(RC rc, const Matrix<char>& board, int first_dir = 0)
      : rc0(rc),
        initial_dir_counter(0),
        curdir(first_dir),
        currc{rc0},
        board(board),
        curdist(0),
        first_dir(first_dir)
  {
  }
  // Returns cell position, direction this cell has been entered along and distance (mirrors not
  // count in distance)
  tuple<RC, int, int> next()
  {
    assert(curdir < 4);
    for (;;) {
      currc = move_towards(currc, curdir);
      if (within_table(currc, board.H, board.W)) {
        char c = board(currc);
        if (c == '.') {
          ++curdist;
          return make_tuple(currc, curdir, curdist);
        }
        if (c == '/') {
          curdir = turn_at_slash_mirror(curdir);
          continue;
        }
        if (c == '\\') {
          curdir = turn_at_backslash_mirror(curdir);
          continue;
        }
        if (accumulate_lanterns) {
          int lantern_color = lantern_color_or_0(c);
          if (lantern_color != 0)
            accumulated_lantern_colors |= lantern_color;
        }
      }
      ++initial_dir_counter;
      if (initial_dir_counter >= 4)
        return make_tuple(RC{-1, -1}, -1, -1);
      currc = rc0;
      curdir = (initial_dir_counter + first_dir) % 4;
    }
  }
};

struct BoardX
{
  const int H, W;
  Matrix<char> board;
  Matrix<CellX> boardx;
  vector<Crystal> crystals;
  vector<Lantern> lanterns;

  BoardX(int H, int W) : H(H), W(W), board(H, W), boardx(H, W)
  {
    crystals.reserve(H * W / 4);  // max 25% crystal probability.
    lanterns.reserve(crystals.size());
  }

  void rebuild_boardx()
  {  // Prepare info about cells.
    FOR (crix, 0, < crystals.size()) {
      auto& cry = crystals[crix];
      // Enumerate all affected cells
      RayEnumerator re(cry.rc, board);
      re.accumulate_lanterns = true;
      for (;;) {
        RC rc;
        int dir, dist;
        tie(rc, dir, dist) = re.next();
        if (dir < 0)
          break;
        auto& cx = boardx(rc);
        int odir = opposite_dir(dir);
        auto& tw = cx.towards[odir];
        tw.crystal_dist = dist;
        tw.target_crystal_ix = crix;
      }
      cry.rcvd_colors = re.accumulated_lantern_colors;
    }
    FOR (lix, 0, < lanterns.size()) {
      auto& lantern = lanterns[lix];
      RayEnumerator re(lantern.rc, board);
      for (;;) {
        RC rc;
        int dir, dist;
        tie(rc, dir, dist) = re.next();
        if (dir < 0)
          break;
        auto& cx = boardx(rc);
        auto& tw = cx.towards[dir];
        tw.lit_by_lantern_ix = lix;
        tw.lantern_dist = dist;
      }
    }
  }
};

bool can_see_each_other(RC x, RC y, const Matrix<char>& board)
{
  if (x == y)
    return true;
  int first_dir = 0;
  if (x.r == y.r) {
    first_dir = x.c < y.c ? DIR_RIGHT : DIR_LEFT;
  } else if (x.c == y.c) {
    first_dir = x.r < y.r ? DIR_DOWN : DIR_UP;
  }
  RayEnumerator re(x, board, first_dir);
  for (;;) {
    RC rc;
    int dir;
    tie(rc, dir, ignore) = re.next();
    if (dir < 0)
      return false;
    if (rc == y)
      return true;
  }
}
struct Option
{
  RC rc;
  int colors;
};

// Which crystals will be affected by putting a lantern of a specific color
struct PutLanternEffects
{
  struct CrixColors
  {
    small crix;  // crystal index
    small colors;
  };
  InlineVector<CrixColors, 8> crystals_affected;

  // 'first time' means the caller is sure that this crystal_ix has not been added yet
  void light_crystal_first_time(int crystal_ix, int color)
  {
    crystals_affected.push_back(CrixColors{(small)crystal_ix, (small)color});
    assert(crystals_affected.size() <= 4);
  }
  bool is_empty() const { return crystals_affected.empty(); }
  // Return gain and which crystals got only one light out of
  // two.
  tuple<Gain, int> evaluate(int CL,
                            const vector<Crystal>& crystals,
                            vector<int>* crystals_needing_one_more_light) const
  {
    if (crystals_needing_one_more_light)
      crystals_needing_one_more_light->clear();
    Gain gain;
    int colors_used = 0;
    int colors_missing = 0;
    int num_crystals_going_incorrect = 0;
    for (auto cc : crystals_affected) {
      colors_used |= cc.colors;
      auto& cry = crystals[cc.crix];
      assert(!cry.lit_incorrectly_unfixable());  // Crystal should not be here if it's already
                                                 // incorrect.
      if (cry.are_any_colors_not_needed_here(cc.colors)) {
        if (cry.lit_incorrectly_fixable()) {
          assert(cry.missing_colors() != 0);
        } else {
          gain.realized -= 10;
          ++num_crystals_going_incorrect;
          if (cry.missing_colors() == 0)
            gain.realized -= colors_kind(cry.colors) == ColorsKind::secondary ? 30 : 20;
          else {
            assert(cry.rcvd_colors == 0);
          }
        }
      } else if ((cry.missing_colors() & cc.colors) != 0) {
        if ((cry.missing_colors() & ~cc.colors) == 0) {
          // This crystal would be done.
          bool needs_secondary = colors_kind(cry.colors) == ColorsKind::secondary;
          gain.realized += needs_secondary ? 30 : 20;
          if (needs_secondary && num_colors(cry.missing_colors()) == 1)
            gain.realized += 10;  // previously it's been incorrect, remove the penalty
        } else {
          // 'colors' contain missing colors but not all. Missing colors must be 2 colors
          // while colors a single one.
          assert(colors_kind(cc.colors) == ColorsKind::primary &&
                 colors_kind(cry.missing_colors()) == ColorsKind::secondary);
          colors_missing |= (cry.missing_colors() & ~cc.colors);
          if (crystals_needing_one_more_light)
            crystals_needing_one_more_light->push_back(cc.crix);
          gain.realized -= 10;
          gain.future += 10 + 30;
        }
      }
    }
    gain.realized -= num_colors(colors_used) * CL;
    gain.future -= num_colors(colors_missing) * CL;
    return make_tuple(gain, num_crystals_going_incorrect);
  }  // evaluate

  static PutLanternEffects unite(const PutLanternEffects& x, const PutLanternEffects& y)
  {
    PutLanternEffects result(x);
    for (auto y_cc : y.crystals_affected) {
      bool found = false;
      // find y_cc.crix in result
      for (auto& r_cc : result.crystals_affected) {
        if (r_cc.crix == y_cc.crix) {
          r_cc.colors |= y_cc.colors;
          found = true;
          break;
        }
      }
      if (!found) {
        result.crystals_affected.push_back(y_cc);
        assert(result.crystals_affected.size() <= 8);
      }
    }
    return result;
  }
};
bool thisisit(RC u, int uc, RC v, int vc, RC x, int xc, RC y, int yc)
{
  return (u == x && uc == xc && v == y && vc == yc) || (v == x && vc == xc && u == y && uc == yc);
}
struct LanternJury
{
  LanternJury(vector<int>* tmpvector) : crystals_needing_one_more_light(tmpvector) {}
  struct SingleLanternSolution
  {
    SingleLanternSolution() { gain.set_int_min(); }
    RC rc{-1, -1};
    int color = -1;
    Gain gain;
  };
  struct HalfSolution
  {
    RC rc;
    int color;
    PutLanternEffects ple;
    int num_crystals_going_incorrect;
  };
  SingleLanternSolution best_single_lantern_solution;
  unordered_map<int, vector<HalfSolution>> halfsolutions_by_crystal;
  vector<int>* crystals_needing_one_more_light;
  void add_option(RC rc,
                  int color,
                  const PutLanternEffects& ple,
                  int CL,
                  const vector<Crystal>& crystals)
  {
    // At this point the ple contains a single color
    assert(!ple.is_empty());
    Gain gain;
    int num_crystals_going_incorrect;
    tie(gain, num_crystals_going_incorrect) =
        ple.evaluate(CL, crystals, crystals_needing_one_more_light);
    if (gain.best_combined() > 0 &&
        gain.weighted_combined() > best_single_lantern_solution.gain.weighted_combined()) {
      best_single_lantern_solution.rc = rc;
      best_single_lantern_solution.color = color;
      best_single_lantern_solution.gain = gain;
    }
    if (USE_HALFSOLUTIONS) {
      for (auto& c : *crystals_needing_one_more_light) {
        halfsolutions_by_crystal[c].emplace_back(
            HalfSolution{rc, color, ple, num_crystals_going_incorrect});
      }
    }
  }
  struct NextMove
  {
    NextMove() { gain.set_int_min(); }

    array<Lantern, 2> lanterns;
    Gain gain;
  };
  NextMove get_best_next_move(const vector<Crystal>& crystals, int CL, const Matrix<char>& board)
  {
    NextMove nm;
    if (best_single_lantern_solution.gain.best_combined() > 0) {
      nm.lanterns[0].rc = best_single_lantern_solution.rc;
      nm.lanterns[0].color = best_single_lantern_solution.color;
      nm.gain = best_single_lantern_solution.gain;
    }
    if (USE_HALFSOLUTIONS) {
      for (auto& kv : halfsolutions_by_crystal) {
        auto& hss = kv.second;
        int N = hss.size();
        FOR (i, 0, < N - 1) {
          auto& u = hss[i];
          FOR (j, i + 1, < N) {
            auto& v = hss[j];
            if (u.color == v.color || u.rc == v.rc)
              continue;
            if (can_see_each_other(u.rc, v.rc, board))
              continue;
            // Unite the solutions u and v
            Gain gain;
            tie(gain, ignore) =
                PutLanternEffects::unite(u.ple, v.ple).evaluate(CL, crystals, nullptr);
            if (gain.best_combined() > 0 &&
                gain.weighted_combined() > nm.gain.weighted_combined()) {
              nm.lanterns[0].rc = u.rc;
              nm.lanterns[0].color = u.color;
              nm.lanterns[1].rc = v.rc;
              nm.lanterns[1].color = v.color;
              nm.gain = gain;
            }
          }
        }
      }
    }
    return nm;
  }
};

class CrystalLighting
{
  int CL, CM, CO, MM, MO, H, W;  // Game constants

public:
  vector<string> placeItems(const vector<string>& targetBoard,
                            int costLantern,
                            int costMirror,
                            int costObstacle,
                            int maxMirrors,
                            int maxObstacles)
  {
    auto t0 = getTime();

    CL = costLantern;
    CM = costMirror;
    CO = costObstacle;
    MM = maxMirrors;
    MO = maxObstacles;
    H = (int)targetBoard.size();
    W = (int)targetBoard[0].size();

    // printBoard(targetBoard);

    BoardX boardx(H, W);
    FOR (r, 0, < H) {
      FOR (c, 0, < W) {
        boardx.board(RC{(small)r, (small)c}) = targetBoard[r][c];
      }
    }

    // Collect crystals.
    FOR (r, 0, < H) {
      FOR (c, 0, < W) {
        RC rc{(small)r, (small)c};
        int crystal_color = crystal_color_or_0(boardx.board(rc));
        if (crystal_color != 0) {
          boardx.crystals.emplace_back(rc, (small)crystal_color);
        }
      }
    }

    int score = 0;
    vector<int> lantern_jury_tmpvector;
    FOR (counter, 0, < INT_MAX) {
      boardx.rebuild_boardx();

      LanternJury lantern_jury(&lantern_jury_tmpvector);
      FOR (r, 0, < H) {
        FOR (c, 0, < W) {
          RC rc{(small)r, (small)c};

          if (boardx.board(rc) != '.')
            continue;

          auto& cell = boardx.boardx(rc);
          bool skip_this = false;
          FOR (i, 0, < 4) {
            if (cell.towards[i].lit_by_lantern_ix >= 0) {
              skip_this = true;
              break;
            }
          }
          if (skip_this)
            continue;

          array<PutLanternEffects, 3> put_lantern_effects_by_colorix;
          FOR (i, 0, < 4) {
            auto& twi = cell.towards[i];
            if (twi.target_crystal_ix < 0)
              continue;
            auto& cry = boardx.crystals[twi.target_crystal_ix];
            if (cry.lit_incorrectly_unfixable())
              continue;
            FOR (colorix, 0, < 3) {
              int color = 1 << colorix;
              if (cry.do_colors_make_difference(color)) {
                put_lantern_effects_by_colorix[colorix].light_crystal_first_time(
                    twi.target_crystal_ix, color);
              }
            }
          }
          FOR (colorix, 0, < 3) {
            if (put_lantern_effects_by_colorix[colorix].is_empty())
              continue;
            lantern_jury.add_option(rc, 1 << colorix, put_lantern_effects_by_colorix[colorix], CL,
                                    boardx.crystals);
          }
        }
      }
      auto nm = lantern_jury.get_best_next_move(boardx.crystals, CL, boardx.board);
      if (nm.gain.best_combined() <= 0)
        break;
      assert(nm.lanterns[0].color != 0);
      //      fprintf(stderr, "Counter: %d\n", counter);
      for (auto& l : nm.lanterns) {
        if (l.color == 0)
          break;
        assert(boardx.board(l.rc) == '.');
        boardx.board(l.rc) = make_lantern_char(l.color);
        boardx.lanterns.emplace_back(Lantern{l.rc, l.color});
        //      fprintf(stderr, "\t %d %d %d\n", l.rc.r, l.rc.c, l.color);
      }
      score += nm.gain.realized;
      //      fprintf(stderr, "-> Score = %d\n", score);
    }
    vector<string> result;
    result.reserve(boardx.lanterns.size());
    char s[100];
    for (auto& l : boardx.lanterns) {
      sprintf(s, "%d %d %c", l.rc.r, l.rc.c, l.color + '0');
      result.emplace_back(s);
    }
    fprintf(stderr, "SCORE should be %d\n", score);
    auto t1 = getTime();
    fprintf(stderr, "Elapsed: %f\n", duration<double>(t1 - t0).count());
    return result;
  }
};

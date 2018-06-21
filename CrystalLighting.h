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

const float FUTURE_FACTOR = 0.5f;
#define SOMEVECTOR_IS_INLINEVECTOR 0

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

using hrc = chrono::high_resolution_clock;
using time_point = hrc::time_point;
using chrono::duration;
time_point getTime()
{
  return hrc::now();
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

#if SOMEVECTOR_IS_INLINEVECTOR
template <class T, int N>
using SomeVector = InlineVector<T, N>;
#else
template <class T, int N>
using SomeVector = vector<T>;
#endif

template <class T, int MaxH, int MaxW>
struct Matrix
{
  const int H, W;
  SomeVector<T, MaxH * MaxW> d;  // data, row-major

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

template <class T>
using BoardMatrix = Matrix<T, 100, 100>;

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

enum CellContent : int8_t
{
  o_empty,
  o_obstacle,
  o_mirror_slash,
  o_mirror_backslash,
  o_lantern_blue,
  o_lantern_yellow,
  o_lantern_red,
  o_crystal_blue,
  o_crystal_yellow,
  o_crystal_green,
  o_crystal_red,
  o_crystal_violet,
  o_crystal_orange,
  o_invalid
};

bool is_lantern(CellContent cc)
{
  return o_lantern_blue <= cc && cc <= o_lantern_red;
}
int lantern_color(CellContent cc)
{
  assert(is_lantern(cc));
  return 1 << (cc - o_lantern_blue);
}
bool is_crystal(CellContent cc)
{
  return o_crystal_blue <= cc && cc <= o_crystal_orange;
}
int crystal_colors(CellContent cc)
{
  assert(o_crystal_blue <= cc && cc <= o_crystal_orange);
  return cc - o_crystal_blue + 1;
}

// This could be stored on 2.5 bytes, so 5 bytes per 2 pixels but complicates storage.
struct Cell
{
  CellContent cc;
  uint8_t lights;
  uint16_t freecells;

  void set_freecells(int dir, int count)
  {
    assert(0 <= count && count <= 15);
    int shift = 4 * dir;
    freecells = (freecells & ~(0xf << shift)) | (count << shift);
  }
  int freecells_towards(int dir) const
  {
    int shift = 4 * dir;
    return (freecells >> shift);
  }
  int color_towards(int dir) const
  {
    int a = (lights >> (2 * dir)) & 0x3;
    return a == 0 ? 0 : 1 << (a - 1);
  }
  void set_color_towards(int dir, int color)
  {
    assert(colors_kind(color) == ColorsKind::primary);
    int color_code = (color_code <= 2) ? color : 4;
    int shift = 2 * dir;
    lights = (lights & ~(0x03 << shift)) | (color_code << shift);
  }
};

Cell make_lantern(int color)
{
  return Cell{(CellContent)(o_lantern_blue - 1 + color), 0, 0};
}

static_assert(sizeof(Cell) == 4, "");

struct Crystal
{
  explicit Crystal(Cell cell) : cc(cell.cc), colors(crystal_colors(cell.cc))
  {
    rcvd_colors = 0;
    FOR (dir, 0, < 4) {
      rcvd_colors |= cell.color_towards(dir);
    }
  }
  int missing_colors() const { return colors & ~rcvd_colors; }
  bool lit_incorrectly_unfixable() const { return (rcvd_colors & ~colors) != 0; }
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
  CellContent cc;
  int colors;
  int rcvd_colors;
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

int turn_ccw(int dir)
{
  static const int d[4] = {2, 3, 1, 0};
  return d[dir];
}

RC move_towards(RC rc, int dir, int count = 1)
{
  switch (dir) {
    case 0:
      rc.r += count;
      break;
    case 1:
      rc.r -= count;
      break;
    case 2:
      rc.c += count;
      break;
    case 3:
      rc.c -= count;
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
  int initial_dir_counter = 0;
  int curdir;
  RC currc;
  const BoardMatrix<Cell>& board;
  int curdist = 0;
  const int first_dir;
  const bool single_dir;
  Cell cell0, curcell;
  bool done_with_this_dir = false;

public:
  // With the default settings nothing will be reported.
  bool report_mirrors = false;
  bool report_empty = false;
  bool report_lanterns = false;
  bool report_crystals = false;
  RayEnumerator(RC rc, const BoardMatrix<Cell>& board, int first_dir = 0, bool single_dir = false)
      : rc0(rc),
        curdir(first_dir),
        currc{rc0},
        board(board),
        first_dir(first_dir),
        single_dir(single_dir),
        cell0(board(rc)),
        curcell(cell0)
  {
  }
  int current_dir() const { return curdir; }
  RC current_rc() const { return currc; }
  int current_dist() const { return curdist; }
  int initial_dir() const { return (initial_dir_counter + first_dir) % 4; }

  // Skips mirrors if report_mirrors is false.
  // Reports turned direction at mirrors.
  Cell next()
  {
    assert(curdir < 4);
    for (;;) {
      if (done_with_this_dir) {
        done_with_this_dir = false;
      } else {
        int dist;
        if (!report_empty && curcell.cc == o_empty) {
          dist = curcell.freecells_towards(curdir) + 1;
        } else {
          dist = 1;
        }
        currc = move_towards(currc, curdir, dist);
        curdist += dist;
        curcell = board(currc);
        if (within_table(currc, board.H, board.W)) {
          if (curcell.cc == o_empty) {
            if (report_empty)
              return curcell;
            continue;
          }
          if (curcell.cc == o_mirror_slash) {
            curdir = turn_at_slash_mirror(curdir);
            if (report_mirrors)
              return curcell;
            continue;
          }
          if (curcell.cc == o_mirror_backslash) {
            curdir = turn_at_slash_mirror(curdir);
            if (report_mirrors)
              return curcell;
            continue;
          }
          // curcell is crystal or lantern or obstacle
          if ((report_lanterns && is_lantern(curcell.cc)) ||
              (report_crystals && is_crystal(curcell.cc))) {
            done_with_this_dir = true;
            return curcell;
          }
        }
      }
      ++initial_dir_counter;
      if (single_dir || initial_dir_counter >= 4) {
        curdir = -1;
        currc = RC{-1, -1};
        curdist = -1;
        curcell = Cell{o_invalid, 0, 0};
        return curcell;
      }
      currc = rc0;
      curdir = (initial_dir_counter + first_dir) % 4;
      curcell = cell0;
    }
  }
};

int tolix(RC rc, int W)
{
  return rc.r * W + rc.c;
}

RC torc(int lix, int W)
{
  return RC{(small)(lix / W), (small)(lix % W)};
}

// Stores list of integers in max N1 buckets with N2 integers per bucket
template <int N1, int N2>
struct Int16Set
{
  SomeVector<InlineVector<int16_t, N2>, N1> vv;
  SomeVector<tuple<int16_t, int16_t>, N2>
      idcs;  // Stores the <first ordinal, first number> of buckets
  int total;
  void init(int16_t from, int16_t to)
  {
    assert(from <= to && to - from < N1 * N2);
    vv.clear();
    idcs.clear();
    total = from - to;
    InlineVector<int16_t, N2>* vb = nullptr;
    int n = 0;
    FOR (c, from, < to) {
      if (!vb || vb->size() == N1) {
        vv.push_back(InlineVector<int16_t, N2>{});
        vb = &vv.back();
        idcs.emplace_back(n, c);
      }
      vb->push_back(c);
    }
  }
  void remove_existing(int16_t x)
  {
    auto it = lower_bound(
        BE(idcs), x, [](const tuple<int16_t, int16_t>& a, int16_t b) { return get<1>(a) < b; });
    assert(it != idcs.end());
    auto i = it - idcs.begin();
    auto& vvi = vv[i];
    auto it2 = lower_bound(BE(vvi), x);
    assert(it2 != vvi.end());
    vvi.erase(it2);
    for (++it; it != idcs.end(); ++it) {
      --get<0>(*it);
    }
    --total;
  }
  int16_t get_nth(int16_t n)
  {
    auto it = lower_bound(
        BE(idcs), n, [](const tuple<int16_t, int16_t>& a, int16_t b) { return get<0>(a) < b; });
    assert(it != idcs.end());
    int d = n - get<0>(*it);
    auto& vvi = vv[it - idcs.begin()];
    assert(0 <= d && d < vvi.size());
    return vvi[d];
  }
};

// CellList represents an ordered set of cell coordinates.
// You can remove one cell or ask for the nth cell.
struct CellList
{
  Int16Set<100, 100> il;
  const int H, W;
  CellList(int H, int W) : H(H), W(W) { il.init(0, H * W); }
  void remove_existing(RC rc)
  {
    int lix = tolix(rc, W);
    il.remove_existing(lix);
  }
  RC get_nth(int n)
  {
    auto lix = il.get_nth(n);
    return torc(lix, W);
  }
  int count() const { return il.total; }
};

struct BoardX
{
  const int H, W;
  BoardMatrix<Cell> ccm;  // For each cell one of the 13 possible cell content.
  CellList freecells;
  BoardX(int H, int W) : H(H), W(W), ccm(H, W), freecells(H, W) {}

  void init(const vector<string>& board)
  {
    assert(board.size() == H);
    FOR (r, 0, < H) {
      auto& boardr = board[r];
      assert(boardr.size() == W);
      for (RC rc{(small)r, (small)0}; rc.c < W; ++rc.c) {
        CellContent cc;
        auto ch = boardr[rc.c];
        switch (ch) {
          case '.':
            cc = o_empty;
            break;
          case 'X':
            cc = o_obstacle;
            break;
          case '/':
            cc = o_mirror_slash;
            break;
          case '\\':
            cc = o_mirror_backslash;
            break;
          default:
            assert('1' <= ch && ch <= '6');
            {
              int color = ch - '0';
              cc = (CellContent)(o_crystal_blue - 1 + color);
            }
        }
        ccm(rc) = Cell{cc, 0, 0};
      }
    }
    // Free cell information.
    const RC startrc[4] = {RC{(small)0, (small)0}, RC{(small)(H - 1), (small)(W - 1)},
                           RC{(small)(H - 1), (small)(0)}, RC{(small)(0), (small)(W - 1)}};
    FOR (dir, 0, > 4) {
      auto ccwdir = turn_ccw(dir);
      auto oppdir = opposite_dir(dir);
      for (RC rc0 = startrc[dir]; within_table(rc0, H, W); rc0 = move_towards(rc0, ccwdir)) {
        int c = 0;
        for (RC rc = rc0; within_table(rc0, H, W); rc = move_towards(rc, dir)) {
          auto& cell = ccm(rc);
          if (cell.cc == o_empty) {
            cell.set_freecells(oppdir, std::min(c, 15));
            ++c;
          } else {
            c = 0;
          }
        }
      }
    }
  }
  // Returns the coords of lanterns, crystals, obstacles in all four directions.
  // Returns RC{-1,-1} if ray goes off board.
  tuple<array<Cell, 4>, array<RC, 4>> raytrace_in_four_dirs(RC rc0) const
  {
    tuple<array<Cell, 4>, array<RC, 4>> result;
    auto& cells = get<0>(result);
    auto& rcs = get<1>(result);
    FOR (initial_dir, 0, < 4) {
      auto rc = move_towards(rc0, initial_dir);
      int dir = initial_dir;
      rcs[dir] = RC{(small)-1, (small)-1};
      cells[dir].cc = o_invalid;
      for (; within_table(rc, H, W);) {
        auto cell = ccm(rc);
        if (cell.cc == o_empty) {
          rc = move_towards(rc, dir, cell.freecells_towards(dir) + 1);
        } else if (cell.cc == o_mirror_slash) {
          dir = turn_at_slash_mirror(dir);
          rc = move_towards(rc, dir);
        } else if (cell.cc == o_mirror_backslash) {
          dir = turn_at_backslash_mirror(dir);
          rc = move_towards(rc, dir);
        } else {
          assert(is_lantern(cell.cc) || is_crystal(cell.cc) || cell.cc == o_obstacle);
          rcs[dir] = rc;
          cells[dir] = cell;
          break;
        }
      }
    }
    return result;
  }
  void raytrace_in_four_dirs_get_all(RC rc0, array<vector<tuple<Cell, RC>>, 4>& result) const
  {
    FOR (initial_dir, 0, < 4) {
      auto& rid = result[initial_dir];
      rid.clear();
    }
  }
  array<int, 4> get_lights_towards_at_empty(RC rc) const
  {
    assert(ccm(rc).cc == o_empty);
    RayEnumerator re(rc, ccm, 0, false);
    re.report_lanterns = true;
    array<int, 4> lights_towards;
    lights_towards.fill(0);
    for (;;) {
      auto cell = re.next();
      if (cell.cc == o_invalid)
        break;
      assert(is_lantern(cell.cc));
      lights_towards[opposite_dir(re.initial_dir())] = lantern_color(cell.cc);
    }
    return lights_towards;
  }
  void set_over_empty(RC rc, CellContent cell_content)
  {
    assert(ccm(rc).cc == o_empty && cell_content != o_empty);
    Cell cell{cell_content, 0, 0};
    //- remove lights, if any
    auto lights_towards = get_lights_towards_at_empty(rc);
    FOR (dir, 0, < 4) {
      if (lights_towards[dir] == 0)
        continue;
      RayEnumerator re(rc, ccm, dir, true);
      re.report_mirrors = true;
      re.report_crystals = true;
      for (;;) {
        auto cell = re.next();
        if (cell.cc == o_invalid)
          break;
        ccm(re.current_rc()).set_color_towards(re.current_dir(), 0);
      }
    }
    //- update freecells
    FOR (dir, 0, < 4) {
      RayEnumerator re(rc, ccm, dir, true);
      re.report_empty = true;
      int odir = opposite_dir(dir);
      for(int counter = 0;;++counter) {
        auto cell = re.next();
        if (cell.cc == o_invalid||re.current_dir()!=dir) break;
        ccm(re.current_rc()).set_freecells(odir, std::min(counter, 15));
      }
    }
    // - set lights
    array<int, 4> new_lights_towards;
    if (cell.
      FOR(dir, 0, <4) {
        int new_dir = cell.cc==o_mirror_slash?turn_at_slash_mirror(dir):turn_at_backslash_mirror(dir);
        new_lights_towards[new_dir] = lights_towards[dir];
      }
    } else if (is_lantern(cell.cc)) {
      assert(all_of(BE(lights_towards), [](auto x){return x==0}));
      new_lights_towards.fill(lantern_color(cell.cc));
    } else if (is_crystal(cell.cc)) {
      new_lights_towards = lights_towards;
    } else assert(false);}
               //- set cell
               ccm(rc) = cell;
};

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
    RC rc;
    Crystal crystal;
    int colors;
  };
  InlineVector<CrixColors, 8> crystals_affected;

  // 'first time' means the caller is sure that this crystal_ix has not been added yet
  void light_crystal_first_time(RC rc, Crystal crystal, int color)
  {
    crystals_affected.push_back(CrixColors{rc, crystal, color});
    assert(crystals_affected.size() <= 4);
  }
  bool is_empty() const { return crystals_affected.empty(); }
  // Return gain and which crystals got only one light out of
  // two.
  tuple<Gain, int> evaluate(int CL) const
  {
    Gain gain;
    int colors_used = 0;
    int colors_missing = 0;
    int num_crystals_going_incorrect = 0;
    for (auto cc : crystals_affected) {
      colors_used |= cc.colors;
      assert(
          !cc.crystal.lit_incorrectly_unfixable());  // Crystal should not be here if it's already
                                                     // incorrect.
      if (cc.crystal.are_any_colors_not_needed_here(cc.colors)) {
        if (cc.crystal.lit_incorrectly_fixable()) {
          assert(cc.crystal.missing_colors() != 0);
        } else {
          gain.realized -= 10;
          ++num_crystals_going_incorrect;
          if (cc.crystal.missing_colors() == 0)
            gain.realized -= colors_kind(cc.crystal.colors) == ColorsKind::secondary ? 30 : 20;
          else {
            assert(cc.crystal.rcvd_colors == 0);
          }
        }
      } else if ((cc.crystal.missing_colors() & cc.colors) != 0) {
        if ((cc.crystal.missing_colors() & ~cc.colors) == 0) {
          // This crystal would be done.
          bool needs_secondary = colors_kind(cc.crystal.colors) == ColorsKind::secondary;
          gain.realized += needs_secondary ? 30 : 20;
          if (needs_secondary && num_colors(cc.crystal.missing_colors()) == 1)
            gain.realized += 10;  // previously it's been incorrect, remove the penalty
        } else {
          // 'colors' contain missing colors but not all. Missing colors must be 2 colors
          // while colors a single one.
          assert(colors_kind(cc.colors) == ColorsKind::primary &&
                 colors_kind(cc.crystal.missing_colors()) == ColorsKind::secondary);
          colors_missing |= (cc.crystal.missing_colors() & ~cc.colors);
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
        if (r_cc.rc == y_cc.rc) {
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
  struct SingleLanternSolution
  {
    SingleLanternSolution() { gain.set_int_min(); }
    RC rc{-1, -1};
    int color = -1;
    Gain gain;
  };
  SingleLanternSolution best_single_lantern_solution;
  void add_option(RC rc, int color, const PutLanternEffects& ple, int CL)
  {
    // At this point the ple contains a single color
    assert(!ple.is_empty());
    Gain gain;
    int num_crystals_going_incorrect;
    tie(gain, num_crystals_going_incorrect) = ple.evaluate(CL);
    if (gain.best_combined() > 0 &&
        gain.weighted_combined() > best_single_lantern_solution.gain.weighted_combined()) {
      best_single_lantern_solution.rc = rc;
      best_single_lantern_solution.color = color;
      best_single_lantern_solution.gain = gain;
    }
  }
  struct NextMove
  {
    NextMove() { gain.set_int_min(); }

    array<Lantern, 2> lanterns;
    Gain gain;
  };
  NextMove get_best_next_move(int CL)
  {
    NextMove nm;
    if (best_single_lantern_solution.gain.best_combined() > 0) {
      nm.lanterns[0].rc = best_single_lantern_solution.rc;
      nm.lanterns[0].color = best_single_lantern_solution.color;
      nm.gain = best_single_lantern_solution.gain;
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
    boardx.init(targetBoard);

    int score = 0;
    hrc::duration d_mainloop = hrc::duration::zero();
    FOR (counter, 0, < INT_MAX) {
      LanternJury lantern_jury;
      auto d0 = getTime();
      FOR (freecellcounter, 0, < boardx.freecells.count()) {
        auto rc = boardx.freecells.get_nth(freecellcounter);

        assert(boardx.ccm(rc).cc == o_empty);

        auto neighbors = boardx.raytrace_in_four_dirs(rc);
        auto& nbcells = get<0>(neighbors);
        bool skip = false;
        bool has_crystal = false;
        for (auto& nbc : nbcells) {
          if (is_lantern(nbc.cc)) {
            skip = true;
            break;
          }
          if (is_crystal(nbc.cc))
            has_crystal = true;
        }
        if (skip || !has_crystal)
          continue;

        array<PutLanternEffects, 3> put_lantern_effects_by_colorix;
        auto& nbrcs = get<1>(neighbors);
        FOR (dir, 0, < 4) {
          auto cell = nbcells[dir];
          if (!is_crystal(cell.cc))
            continue;
          auto cry = Crystal(cell);
          if (cry.lit_incorrectly_unfixable())
            continue;
          FOR (colorix, 0, < 3) {
            int color = 1 << colorix;
            if (cry.do_colors_make_difference(color)) {
              put_lantern_effects_by_colorix[colorix].light_crystal_first_time(nbrcs[dir], cry,
                                                                               color);
            }
          }
        }
        FOR (colorix, 0, < 3) {
          if (put_lantern_effects_by_colorix[colorix].is_empty())
            continue;
          lantern_jury.add_option(rc, 1 << colorix, put_lantern_effects_by_colorix[colorix], CL);
        }
      }
      d_mainloop += getTime() - d0;
      auto nm = lantern_jury.get_best_next_move(CL);
      if (nm.gain.best_combined() <= 0)
        break;
      assert(nm.lanterns[0].color != 0);
      //      fprintf(stderr, "Counter: %d\n", counter);
      for (auto& l : nm.lanterns) {
        if (l.color == 0)
          break;
        boardx.set_over_empty(l.rc, make_lantern(l.color));
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
    fprintf(stderr, "Rebuild: %f\n", duration<double>(d_rebuild).count());
    fprintf(stderr, "Mainloop: %f\n", duration<double>(d_mainloop).count());
    return result;
  }
};

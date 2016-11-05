#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <sstream>

#define CHECK(expr)                               \
  if (expr) {                                     \
  } else {                                        \
    fprintf(stderr, "CHECK Failed (%s:%d): %s\n", \
            __FILE__, __LINE__, #expr);           \
    exit(EXIT_FAILURE);                           \
  }

#define REP2(i, m, n) for (long long i = (long long)(m); i < (long long)(n); i++)
#define REP(i, n) REP2(i, 0, n)
#define ALL(S) (S).begin(), (S).end()
#define DEBUG(x) do { cerr << #x << ": " << (x) << endl; } while (0)
#define fst first
#define snd second
typedef std::pair<int, int> PI;

template <typename S, typename T> std::ostream &operator<<(std::ostream &out, const std::pair<S, T> &p) {
  out << "(" << p.first << ", " << p.second << ")";
  return out;
}

template <typename T> std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  out << "[";
  REP(i, v.size()){
    if (i > 0) out << ", ";
    out << v[i];
  }
  out << "]";
  return out;
}

#endif /* COMMON_H */

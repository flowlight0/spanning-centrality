#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <sstream>

#define fst first
#define snd second
typedef std::pair<int, int> PI;

template <typename S, typename T> std::ostream &operator<<(std::ostream &out, const std::pair<S, T> &p) {
  out << "(" << p.first << ", " << p.second << ")";
  return out;
}

template <typename T> std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  out << "[";
  for (size_t i = 0; i < v.size(); i++) {
    if (i > 0) out << ", ";
    out << v[i];
  }
  out << "]";
  return out;
}

#endif /* COMMON_H */

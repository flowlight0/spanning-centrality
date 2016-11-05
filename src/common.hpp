#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <fstream>

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


bool ReadGraph(const std::string &graph_file, std::vector<PI> &es){
  es.clear();
  
  std::ifstream ifs(graph_file);
  if (!ifs.good()){
    std::cerr << "Error: open graph_file." << std::endl;
    return false;
  }
    
  for (int u, v; ifs >> u >> v;) {
    if (u != v) es.emplace_back(u, v);
  }
  ifs.close();
  return true;
}


void ConvertToUndirectedGraph(std::vector<PI> &es){
  // remove self loop.
  size_t m = 0;
  for (size_t i = 0; i < es.size(); i++){
    if (es[i].fst != es[i].snd) es[m++] = es[i];
  }
  es.resize(m);
    
  // remove redundant edges.
  for (auto &e : es){
    if (e.fst > e.snd) std::swap(e.fst, e.snd);
  }
  std::sort(es.begin(), es.end());
  es.erase(unique(es.begin(), es.end()), es.end());
}


#endif /* COMMON_H */

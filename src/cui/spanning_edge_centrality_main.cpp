#include "algorithm/spanning_centrality.hpp"
#include "common.hpp"
#include "gflags/gflags.h"
#include <cstdio>
#include <map>
using namespace std;
DEFINE_string(graph_file, "-", "Input graph file.");
DEFINE_int32(num_samples, 1000, "");

int main(int argc, char *argv[])
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  vector<PI> es;
  ReadGraph(FLAGS_graph_file, es);
  vector<double> sc = spanning_centrality::EstimateEdgeCentrality(es, FLAGS_num_samples);
  
  map<pair<int, int> , double> centrality_map;
  for (size_t i = 0; i < es.size(); i++) {
    const auto &e = es[i];
    centrality_map[make_pair(e.fst, e.snd)] = sc[i];
  }
  
  for (const auto &c : centrality_map){
    printf("%d\t%d\t%lf\n", c.fst.fst, c.fst.snd, c.snd);
  }
  return 0;
}

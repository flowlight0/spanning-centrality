#include "algorithm/spanning_centrality.hpp"
#include "gflags/gflags.h"
#include <cstdio>
#include <map>
using namespace std;
DEFINE_string(graph_file, "-", "Input graph file.");
DEFINE_int32(num_samples, 1000, "");

int main(int argc, char *argv[])
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  vector<PI> es(ReadGraph(FLAGS_graph_file));
  vector<double> sc = spanning_centrality::EstimateVertexCentrality(es, FLAGS_num_samples);
  
  map<int, double> centrality_map;
  for (const auto &e : es){
    centrality_map[e.fst] = sc[e.fst];
    centrality_map[e.snd] = sc[e.snd];
  }

  for (const auto &c : centrality_map){
    printf("%d\t%lf\n", c.fst, c.snd);
  }
  return 0;
}

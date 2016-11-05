#include "algorithm/spanning_centrality.hpp"
#include "gflags/gflags.h"
#include <set>
#include <cstdio>
using namespace std;
DEFINE_string(graph_file, "-", "Input graph file.");
DEFINE_int32(num_samples, 10000, "");
DEFINE_bool(aggregated, false, "Output aggregated centrality values or not");

int main(int argc, char *argv[])
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  spanning_centrality::SpanningCentrality spanning_centrality;
  if (spanning_centrality.Construct(FLAGS_graph_file, FLAGS_num_samples)) {
    set<int> vertex_set;
    for (size_t i = 0; i < spanning_centrality.GetNumEdges(); i++) {
      int u = spanning_centrality.GetEdge(i).first;
      int v = spanning_centrality.GetEdge(i).second;
      vertex_set.insert(u);
      vertex_set.insert(v);
    }

    for (int v : vertex_set) {
      double centrality = FLAGS_aggregated ?
        spanning_centrality.GetAggregatedCentrality(v) : 
        spanning_centrality.GetVertexCentrality(v);
      printf("%d\t%.3lf\n", v, centrality);
    }
  } else {
    cerr << "Failed to estimate centrality values from some reason." << endl;
  }
  return 0;
}

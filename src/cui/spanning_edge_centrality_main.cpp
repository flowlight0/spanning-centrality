#include "algorithm/spanning_centrality.hpp"
#include "common.hpp"
#include "gflags/gflags.h"
#include <cstdio>
#include <map>
using namespace std;
DEFINE_string(graph_file, "-", "Input graph file.");
DEFINE_int32(num_samples, 10000, "");

int main(int argc, char *argv[])
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  spanning_centrality::SpanningCentrality spanning_centrality;
  if (spanning_centrality.Construct(FLAGS_graph_file, FLAGS_num_samples)) {
    for (size_t i = 0; i < spanning_centrality.GetNumEdges(); i++) {
      double c = spanning_centrality.GetEdgeCentrality(i);
      int u = spanning_centrality.GetEdge(i).first;
      int v = spanning_centrality.GetEdge(i).second;
      printf("%d\t%d\t%.3lf\n", u, v, c);
    }
  } else {
    cerr << "Failed to estimate centrality values from some reason." << endl;
  }
  return 0;
}

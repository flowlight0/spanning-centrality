#include "common.hpp"
#include "algorithm/articulation_point_decomposition.hpp"
#include "algorithm/spanning_centrality.hpp"
#include "gflags/gflags.h"
#include "jlog/jlog.h"
#include <vector>
using namespace std;

DEFINE_string(graph_file, "-", "Input graph file");
DEFINE_int32(num_trials, 5, "");
DEFINE_int32(num_trees, 1000, "");

namespace {
  size_t CountVertices(const vector<PI> &es){
    std::unordered_set<int> vs;
    for (const auto &e : es){
      vs.insert(e.fst);
      vs.insert(e.snd);
    }
    return vs.size();
  }
  
  void OutputLCCSize(const vector<vector<int>  > &edge_groups, const vector<PI> &es){
    size_t lcc_es = 0;
    size_t lcc_vs = 0;
    for (const auto &edge_group : edge_groups){
      vector<PI> tmp_es;
      for (int id : edge_group){
        tmp_es.push_back(es[id]);
      }

      size_t num_vs = CountVertices(tmp_es);
      if (num_vs > lcc_vs){
        lcc_vs = num_vs;
        lcc_es = tmp_es.size();
      }
    }
    JLOG_PUT("lcc_es", lcc_es);
    JLOG_PUT("lcc_vs", lcc_vs);
  }

  double Average(const std::vector<double> &vec){
    return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
  }
}

int main(int argc, char *argv[])
{
  JLOG_INIT(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  int V = 0;
  vector<PI> es;
  ReadGraph(FLAGS_graph_file, es);
  ConvertToUndirectedGraph(es);
  for (const auto &e : es){
    V = max(V, max(e.fst, e.snd) + 1);
  }
  
  JLOG_PUT("setting.dataset", FLAGS_graph_file);
  JLOG_PUT("setting.num_vs", V);
  JLOG_PUT("setting.num_es", es.size());
  JLOG_PUT("setting.num_ts", FLAGS_num_trees);
  
  JLOG_OPEN("results") JLOG_PUT_BENCHMARK("total_time"){
    ArticulationPointDecomposition decomp(es);
    OutputLCCSize(decomp.edge_groups, es);

    vector<double> elapsed_times;
    for (int q = 0; q < FLAGS_num_trials; q++) {
      double start = jlog_internal::get_current_time_sec();
      spanning_centrality::EstimateEdgeCentrality(es, FLAGS_num_trees);
      elapsed_times.push_back(jlog_internal::get_current_time_sec() - start);
      JLOG_ADD("sampling.elapsed_time", elapsed_times.back());
    }
    JLOG_PUT("sampling.elapsed_time_average", Average(elapsed_times));
  }

  return 0;
}


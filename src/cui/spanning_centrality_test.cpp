#include "algorithm/spanning_centrality.hpp"
#include "gtest/gtest.h"
using namespace spanning_centrality;
using namespace std;

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

const int NUM_SAMPLES = 100000;
const double TOLERANCE = 0.01;

void CheckEdgeCentrality(const vector<pair_int> &es, const vector<double> &answer) {
  vector<double> centrality = EstimateEdgeCentrality(es, NUM_SAMPLES);
  // In this test, es must not contrain duplicated edges. 
  for (size_t i = 0; i < es.size(); i++) {
    ASSERT_NEAR(centrality[i], answer[i], TOLERANCE)
      << es << " " << centrality << " " << answer << " " << i;
  }
}

void CheckVertexCentrality(const vector<pair_int> &es, const vector<double> &answer) {
  vector<double> centrality = EstimateVertexCentrality(es, NUM_SAMPLES);
  // In this test, es must not contrain duplicated edges.
  for (size_t i = 0; i < centrality.size(); i++) {
    ASSERT_NEAR(centrality[i], answer[i], TOLERANCE)
      << es << " " << centrality << " " << answer << " " << i;
  }
}


void CheckAggregatedCentrality(const vector<pair_int> &es, const vector<double> &answer) {
  vector<double> centrality = EstimateAggregatedCentrality(es, NUM_SAMPLES);
  // In this test, es must not contrain duplicated edges.
  for (size_t i = 0; i < centrality.size(); i++) {
    ASSERT_NEAR(centrality[i], answer[i], TOLERANCE)
      << es << " " << centrality << " " << answer << " " << i;
  }
}


TEST(SPANNING_EDGE_CENTRALITY, SMALL0){
  vector<pair_int> es = {{0, 1}, {1, 2}, {1, 3}, {2, 3}, {3, 4}};
  vector<double> ans = {1.00, 0.66, 0.66, 0.66, 1.00};
  CheckEdgeCentrality(es, ans);
}

TEST(SPANNING_EDGE_CENTRALITY, SMALL1){
  vector<pair_int> es = {{0, 1}, {1, 2}, {1, 3}, {2, 3}, {4, 5}};
  vector<double> ans = {1.00, 0.66, 0.66, 0.66, 1.00};
  CheckEdgeCentrality(es, ans);
}

TEST(SPANNING_EDGE_CENTRALITY, SMALL2){
  vector<pair_int> es = {{0, 1}, {0, 2}, {2, 3}};
  vector<double> ans = {1.00, 1.00, 1.00};
  CheckEdgeCentrality(es, ans);
}

TEST(SPANNING_VERTEX_CENTRALITY, SMALL0){
  vector<pair_int> es = {{0, 1}, {1, 2}, {1, 3}, {2, 3}, {3, 4}};
  vector<double> ans = {0.00, 1.00, 0.33, 1.00, 0.00};
  CheckVertexCentrality(es, ans);
}

TEST(SPANNING_VERTEX_CENTRALITY, SMALL1){
  vector<pair_int> es = {{0, 1}, {1, 2}, {1, 3}, {2, 3}, {4, 5}};
  vector<double> ans = {0.00, 1.00, 0.33, 0.33, 0.00, 0.00};
  CheckVertexCentrality(es, ans);
}


TEST(SPANNING_AGGREGATED_CENTRALITY, TREE){
  vector<pair_int> es = {{0, 1}, {0, 2}, {0, 3}, {1, 4}};
  vector<double> ans = {3.0, 2.0, 1.0, 1.0, 1.0};
  CheckAggregatedCentrality(es, ans);
}

TEST(SPANNING_AGGREGATED_CENTRALITY, CLIQUE){
  const int n = 20;
  vector<pair_int> es;
  for (int v = 0; v < n; v++) {
    for (int u = 0; u < v; u++) {
      es.emplace_back(u, v);
    }
  }
  
  vector<double> ans(n, (n - 1) * 2.0 / n);
  CheckAggregatedCentrality(es, ans);
}

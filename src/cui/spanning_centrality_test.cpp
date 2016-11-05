#include "algorithm/spanning_centrality.hpp"
#include "gtest/gtest.h"
using namespace spanning_centrality;
using namespace std;


TEST(SPANNING_EDGE_CENTRALITY, SMALL0){
  const double tol = 0.01;
  vector<PI> es = {{0, 1}, {1, 2}, {1, 3}, {2, 3}, {3, 4}};
  vector<double> ans = {1.00, 0.66, 0.66, 0.66, 1.00};
  vector<double> sc = EstimateEdgeCentrality(es, 10000);

  for (size_t i = 0; i < sc.size(); i++) {
    ASSERT_NEAR(sc[i], ans[i], tol * 2);
  }
}

TEST(SPANNING_EDGE_CENTRALITY, SMALL1){
  const double tol = 0.01;  
  vector<PI> es = {{0, 1}, {1, 2}, {1, 3}, {2, 3}, {4, 5}};
  vector<double> ans = {1.00, 0.66, 0.66, 0.66, 1.00};
  vector<double> sc = EstimateEdgeCentrality(es, 10000);
  for (size_t i = 0; i < sc.size(); i++) {
    ASSERT_NEAR(sc[i], ans[i], tol * 2);
  }
}

TEST(SPANNING_EDGE_CENTRALITY, SMALL2){
  const double tol = 0.01;
  vector<PI> es = {{0, 1}, {0, 2}, {2, 3}};
  ConvertToUndirectedGraph(es);
  vector<double> ans = {1.00, 1.00, 1.00};
  vector<double> sc = EstimateEdgeCentrality(es, 10000);
  for (size_t i = 0; i < sc.size(); i++) {
    ASSERT_NEAR(sc[i], ans[i], tol * 2);
  }
}

TEST(SPANNING_VERTEX_CENTRALITY, SMALL0){
  const double tol = 0.01;
  vector<PI> es = {{0, 1}, {1, 2}, {1, 3}, {2, 3}, {3, 4}};
  vector<double> ans = {0.00, 1.00, 0.33, 1.00, 0.00};
  vector<double> sc = EstimateVertexCentrality(es, 10000);
  
  for (size_t i = 0; i < sc.size(); i++) {
    ASSERT_NEAR(sc[i], ans[i], tol * 2) << sc << " " << ans << endl;
  }
}

TEST(SPANNING_VERTEX_CENTRALITY, SMALL1){
  const double tol = 0.01;  
  vector<PI> es = {{0, 1}, {1, 2}, {1, 3}, {2, 3}, {4, 5}};
  vector<double> ans = {0.00, 1.00, 0.33, 0.33, 0.00, 0.00};
  vector<double> sc = EstimateVertexCentrality(es, 10000);
  for (size_t i = 0; i < sc.size(); i++) {
    ASSERT_NEAR(sc[i], ans[i], tol * 2) << sc << " " << ans << endl;
  }
}

TEST(SPANNING_EDGE_CENTRALITY, KARATE_CLUB){
  vector<PI> es = {
    {1 , 2 }, 
    {1 , 3 },
    {2 , 3 },  
    {1 , 4 },
    {2 , 4 },  
    {3 , 4 },  
    {1 , 5 },  
    {1 , 6 },  
    {1 , 7 },  
    {5 , 7 },  
    {6 , 7 },  
    {1 , 8 },  
    {2 , 8 },  
    {3 , 8 },  
    {4 , 8 },  
    {1 , 9 },  
    {3 , 9 },  
    {3 , 10},  
    {1 , 11},  
    {5 , 11},  
    {6 , 11},  
    {1 , 12},  
    {1 , 13},  
    {4 , 13},  
    {1 , 14},  
    {2 , 14},  
    {3 , 14},  
    {4 , 14},  
    {6 , 17},  
    {7 , 17},  
    {1 , 18},  
    {2 , 18},  
    {1 , 20}, 
    {2 , 20},  
    {1 , 22},  
    {2 , 22},  
    {24, 26},  
    {25, 26},  
    {3 , 28},  
    {24, 28},  
    {25, 28},  
    {3 , 29},  
    {24, 30},
    {27, 30},  
    {2 , 31},  
    {9 , 31},  
    {1 , 32},  
    {25, 32},  
    {26, 32},  
    {29, 32},  
    {3 , 33},  
    {9 , 33},  
    {15, 33},  
    {16, 33},  
    {19, 33},  
    {21, 33},  
    {23, 33},  
    {24, 33},  
    {30, 33},  
    {31, 33},  
    {32, 33},  
    {9 , 34},  
    {10, 34},  
    {14, 34},  
    {15, 34},  
    {16, 34},  
    {19, 34},  
    {20, 34},  
    {21, 34},  
    {23, 34},  
    {24, 34},  
    {27, 34},  
    {28, 34},  
    {29, 34},  
    {30, 34},  
    {31, 34},  
    {32, 34},  
    {33, 34}
  };

  ConvertToUndirectedGraph(es);
  vector<double> sc = EstimateEdgeCentrality(es, 10000);
  
  for (size_t i = 0; i < es.size(); i++){
    ASSERT_LT(sc[i], 1.001) << es[i].fst << '\t' << es[i].snd << '\t' << sc[i] << endl;
    cout << es[i].fst << '\t' << es[i].snd << '\t' << sc[i] << endl;
  }
}


TEST(SPANNING_VERTEX_CENTRALITY_AGGREGATED, TREE){
  const double tol = 0.01;
  vector<PI> es = {{0, 1}, {0, 2}, {0, 3}, {1, 4}};
  vector<double> ans = {3.0, 2.0, 1.0, 1.0, 1.0};
  vector<double> sc = EstimateAggregatedCentrality(es, 100000);

  for (size_t i = 0; i < sc.size(); i++){
    ASSERT_NEAR(sc[i], ans[i], tol * 2) << sc << " " << ans << endl;
  }
}

TEST(SPANNING_VERTEX_CENTRALITY_AGGREGATED, CLIQUE){
  const double tol = 0.01;
  const int n = 20;
  vector<PI> es;
  for (int v = 0; v < n; v++) {
    for (int u = 0; u < v; u++) {
      es.emplace_back(u, v);
    }
  }
  
  vector<double> ans(n, (n - 1) * 2.0 / n);
  vector<double> sc = EstimateAggregatedCentrality(es, 100000);

  for (size_t i = 0; i < sc.size(); i++){
    ASSERT_NEAR(sc[i], ans[i], tol * 2) << sc << " " << ans << endl;
  }
}

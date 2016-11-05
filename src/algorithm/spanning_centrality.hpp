#ifndef SPANNING_CENTRALITY_H
#define SPANNING_CENTRALITY_H

#include "common.hpp"
#include "articulation_point_decomposition.hpp"
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <cassert>
#include <queue>

struct Edge {
  int  from;
  int    to;
  size_t id;
  double centrality;
  Edge(int f, int t, size_t i) : from(f), to(t), id(i), centrality(0){}
};


std::vector<PI> ReadGraph(const std::string &graph_file){
  std::vector<PI> es;
  
  std::ifstream ifs(graph_file);
  if (!ifs.good()){
    std::cerr << "Error: open graph_file." << std::endl;
    exit(EXIT_FAILURE);
  }
    
  for (int u, v; ifs >> u >> v;) {
    if (u != v) es.emplace_back(u, v);
  }
  ifs.close();
  return es;
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
  


namespace {
  std::vector<int> DistanceOrdering(const std::vector<std::vector<PI> > &g){
    std::vector<int> degree(g.size());
    for (size_t v = 0; v < g.size(); v++) {
      degree[v] = g[v].size();
    }
    
    int root = max_element(degree.begin(), degree.end()) - degree.begin();

    std::queue<int>  que;
    std::vector<int> dist(g.size(), std::numeric_limits<int>::max());
    que.push(root);
    dist[root] = 0;
    while (!que.empty()){
      int v = que.front(); que.pop();
      for (const auto &p : g[v]){
        int w = p.fst;
        if (dist[w] > dist[v] + 1){
          dist[w] = dist[v] + 1;
          que.push(w);
        }
      }
    }
    
    std::vector<PI> order;
    for (size_t v = 0; v < g.size(); v++) {
      order.push_back(std::make_pair(dist[v], v));
    }
    std::sort(order.begin(), order.end());
    std::vector<int> res;
    for (size_t v = 0; v < g.size(); v++) {
      res.push_back(order[v].snd);
    }
    assert(res[0] == root);
    return res;
  }
  
  std::vector<std::vector<PI> > BuildCompressedGraph(std::vector<Edge> &es){
    size_t V = 0;
    std::unordered_map<int, int> vertex2id;
    for (auto &e : es){
      if (vertex2id.count(e.from) == 0)  vertex2id[e.from] = V++;
      if (vertex2id.count(e.to  ) == 0)  vertex2id[e.to  ] = V++;
      e.from  = vertex2id[e.from];
      e.to    = vertex2id[e.to  ];
    }

    std::vector<std::vector<PI> >  G(V);
    for (size_t i = 0; i < es.size(); i++){
      G[es[i].from].emplace_back(es[i].to  , i);
      G[es[i].to  ].emplace_back(es[i].from, i);
    }
    return G;
  }
  
  void SampleSpanningTree(std::vector<int > &next_edges,
                          std::vector<bool> &visit,
                          const std::vector<std::vector<PI> > &G,
                          const std::vector<int> &order)
  {
    size_t V = G.size();

    visit[order[0]] = true;

    for (size_t s = 1; s < V; s++) {
      if (visit[order[s]]) continue;
      
      int u = order[s];
      while (!visit[u]){
        int next = rand() % G[u].size();
        next_edges[u] = next;
        u = G[u][next].fst;
      }
        
      u = order[s];
      while (!visit[u]){
        visit[u] = true;
        u = G[u][next_edges[u]].fst;
      }
    }
  }

  void EstimateEdgeCentrality_(std::vector<Edge> &es, int num_samples){
    if (!es.empty()) {
      std::vector<std::vector<PI> > g(BuildCompressedGraph(es));
      
      size_t       V     = g.size();
      std::vector<int>  order = DistanceOrdering(g);
      std::vector<bool> visit(V, false);
      std::vector<int>  next_edges(V, -1);

      for (int trial = 0; trial < num_samples; trial++) {
        fill(visit.begin(), visit.end(), false);
        SampleSpanningTree(next_edges, visit, g, order);
        for (size_t s = 1; s < V; s++) {
          int v = order[s];
          es[g[v][next_edges[v]].snd].centrality += 1.0 / num_samples;
        }
      }
    }
  }
}

std::vector<double> EstimateEdgeCentrality(const std::vector<PI> &es, int num_samples){

  std::vector<double> centrality(es.size(), 0);
  const std::vector<std::vector<int> > edge_groups = ArticulationPointDecomposition(es).edge_groups;

  for (const auto &edge_group : edge_groups) {
    if (edge_group.size() > 1){
      std::vector<Edge> ccomp_es;
      for (int eid : edge_group){
        ccomp_es.emplace_back(es[eid].fst, es[eid].snd, eid);
      }

      EstimateEdgeCentrality_(ccomp_es, num_samples);
      for (const Edge &e : ccomp_es){
        centrality[e.id] = e.centrality;
      }
    } else if (!edge_group.empty()){
      centrality[edge_group[0]] = 1.0;
    }
  }
  return centrality;
}

std::vector<double> EstimateVertexCentrality(const std::vector<PI> &es, int num_samples){
  int V = 0;
  for (const auto &e : es){
    assert(e.fst != e.snd);
    V = std::max({V, e.fst + 1, e.snd + 1});
  }

  ArticulationPointDecomposition decomp(es);
  std::vector<double> centrality(V);
  std::vector<std::vector<int> > edge_groups = decomp.edge_groups;
  std::vector<int> articulation_points = decomp.articulation_points;

  // std::cout << "ESTIMATE: " << std::endl;

  for (const auto &edge_group: edge_groups) {
    if (edge_group.empty()) continue;
    
    std::vector<Edge> ccomp_es;
    for (int eid : edge_group){
      ccomp_es.emplace_back(es[eid].fst, es[eid].snd, eid);
    }
    
    std::vector<std::vector<PI> > g(BuildCompressedGraph(ccomp_es));
    std::vector<int>  order = DistanceOrdering(g);
    // 簡単のため全体の頂点数分の領域を確保している.二重連結成分が多いと効率が悪くなる.
    std::vector<bool> visit(V, false);
    std::vector<int>  next_edges(V, -1);
    std::vector<int>  degree(V);

    for (int trial = 0; trial < num_samples; trial++) {
      fill(visit.begin(), visit.end(), false);
      fill(degree.begin(), degree.end(), 0);
      SampleSpanningTree(next_edges, visit, g, order);
      
      for (size_t s = 1; s < g.size(); s++) {
        int v = order[s];
        int e = ccomp_es[g[v][next_edges[v]].snd].id;
        
        if (++degree[es[e].fst] == 2){
          centrality[es[e].fst] += 1.0 / num_samples;
        }

        if (++degree[es[e].snd] == 2){
          centrality[es[e].snd] += 1.0 / num_samples;
        }
      }
    }
  }

  for (int v : articulation_points) {
    centrality[v] = 1.0;
  }

  for (int v = 0; v < V; v++) {
    assert(centrality[v] < 1 + 1e-9);
  }
  return centrality;
}

std::vector<double> EstimateVertexAggregatedCentrality(const std::vector<PI> &es, int num_samples){
  int V = 0;
  for (const auto &e : es){
    assert(e.fst != e.snd);
    V = std::max({V, e.fst + 1, e.snd + 1});
  }

  std::vector<double> centrality(V);
  std::vector<std::vector<int> > edge_groups = ArticulationPointDecomposition(es).edge_groups;

  // std::cout << "ESTIMATE: " << std::endl;

  std::unordered_set<int> usd;

  for (const auto &edge_group : edge_groups) {
    if (edge_group.empty()) continue;
    
    std::vector<Edge> ccomp_es;
    for (int eid : edge_group){
      ccomp_es.emplace_back(es[eid].fst, es[eid].snd, eid);
    }
    std::vector<std::vector<PI> > g(BuildCompressedGraph(ccomp_es));
    std::vector<int>  order = DistanceOrdering(g);
    // 簡単のため全体の頂点数分の領域を確保している.二重連結成分が多いと効率が悪くなる.
    std::vector<bool> visit(V, false);
    std::vector<int>  next_edges(V, -1);
    
    // std::cout << "START ITERATION: " << num_samples << std::endl;
    for (int trial = 0; trial < num_samples; trial++) {
      fill(visit.begin(), visit.end(), false);
      SampleSpanningTree(next_edges, visit, g, order);

      for (size_t s = 1; s < g.size(); s++) {
        int v = order[s];
        int e = ccomp_es[g[v][next_edges[v]].snd].id;
        centrality[es[e].fst] += 1.0 / num_samples;
        centrality[es[e].snd] += 1.0 / num_samples;
      }
    }
    // std::cout << "FINISH_ITERATION" << std::endl;
  }
  // std::cout << "FINISH_ESTIMATION" << std::endl;
  // std::cout << centrality << std::endl;
  return centrality;
}


#endif /* SPANNING_CENTRALITY_H */

#ifndef SPANNING_CENTRALITY_H
#define SPANNING_CENTRALITY_H

#include "common.hpp"
#include "articulation_point_decomposition.hpp"
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <cassert>
#include <random>
#include <queue>

namespace spanning_centrality {
  
  namespace {
    
    struct Edge {
      int src;
      int dst;
      size_t edge_id;
      double centrality;
      Edge(int src_, int dst_, size_t edge_id_) : src(src_), dst(dst_), edge_id(edge_id_), centrality(0) {}
    };
    
    
    std::vector<std::vector<PI> > BuildCompressedGraph(std::vector<Edge> &es){
      size_t V = 0;
      std::unordered_map<int, int> vertex2id;
      for (auto &e : es){
        if (vertex2id.count(e.src) == 0) vertex2id[e.src] = V++;
        if (vertex2id.count(e.dst) == 0) vertex2id[e.dst] = V++;
        e.src = vertex2id[e.src];
        e.dst = vertex2id[e.dst];
      }

      std::vector<std::vector<PI> > g(V);
      for (size_t i = 0; i < es.size(); i++){
        g[es[i].src].emplace_back(es[i].dst, i);
        g[es[i].dst].emplace_back(es[i].src, i);
      }
      return g;
    }


    std::vector<int> DistanceOrdering(const std::vector<std::vector<PI> > &g){
      int root = 0;
      for (size_t v = 0; v < g.size(); v++) {
        if (g[v].size() > g[root].size()) root = v;
      }
      
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
    
    
    void SampleSpanningTree(std::vector<int > &next_edges,
                            std::vector<bool> &visit,
                            const std::vector<std::vector<PI> > &g,
                            const std::vector<int> &order,
                            std::mt19937 &mt)
    {
      visit[order[0]] = true;
      
      for (size_t s = 1; s < g.size(); s++) {
        if (visit[order[s]]) continue;
      
        int u = order[s];
        while (!visit[u]){
          int next = mt() % g[u].size();
          next_edges[u] = next;
          u = g[u][next].fst;
        }
        
        u = order[s];
        while (!visit[u]){
          visit[u] = true;
          u = g[u][next_edges[u]].fst;
        }
      }
    }
    
    
    void EstimateEdgeCentrality_(std::vector<Edge> &es, int num_samples, std::mt19937 &mt){
      if (es.empty()) return;
      
      std::vector<std::vector<PI> > g(BuildCompressedGraph(es));
        
      size_t V = g.size();
      std::vector<int>  order = DistanceOrdering(g);
      std::vector<bool> visit(V, false);
      std::vector<int>  next_edges(V, -1);

      for (int trial = 0; trial < num_samples; trial++) {
        fill(visit.begin(), visit.end(), false);
        SampleSpanningTree(next_edges, visit, g, order, mt);
        for (size_t s = 1; s < V; s++) {
          int v = order[s];
          es[g[v][next_edges[v]].snd].centrality += 1.0 / num_samples;
        }
      }
    }
  }
  
  
  std::vector<double> EstimateEdgeCentrality(const std::vector<PI> &es, int num_samples){
    std::random_device rd;
    std::mt19937 mt(rd());

    
    std::vector<double> centrality(es.size(), 0);
    const std::vector<std::vector<int> > edge_groups = ArticulationPointDecomposition(es).edge_groups;
    
    for (const auto &edge_group : edge_groups) {
      if (edge_group.size() > 1){
        std::vector<Edge> ccomp_es;
        for (int edge_id : edge_group){
          ccomp_es.emplace_back(es[edge_id].fst, es[edge_id].snd, edge_id);
        }

        EstimateEdgeCentrality_(ccomp_es, num_samples, mt);
        for (const Edge &e : ccomp_es){
          centrality[e.edge_id] = e.centrality;
        }
      } else if (!edge_group.empty()){
        centrality[edge_group[0]] = 1.0;
      }
    }
    return centrality;
  }
  
  
  std::vector<double> EstimateVertexCentrality(const std::vector<PI> &es, int num_samples){
    std::random_device rd;
    std::mt19937 mt(rd());
    
    int V = 0;
    for (const auto &e : es){
      assert(e.fst != e.snd);
      V = std::max({V, e.fst + 1, e.snd + 1});
    }

    ArticulationPointDecomposition decomp(es);
    std::vector<double> centrality(V);
    std::vector<std::vector<int> > edge_groups = decomp.edge_groups;
    std::vector<int> articulation_points = decomp.articulation_points;

    for (const auto &edge_group: edge_groups) {
      if (edge_group.empty()) continue;
    
      std::vector<Edge> ccomp_es;
      for (int edge_id : edge_group){
        ccomp_es.emplace_back(es[edge_id].fst, es[edge_id].snd, edge_id);
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
        SampleSpanningTree(next_edges, visit, g, order, mt);
      
        for (size_t s = 1; s < g.size(); s++) {
          int v = order[s];
          int e = ccomp_es[g[v][next_edges[v]].snd].edge_id;
          
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
    
    std::random_device rd;
    std::mt19937 mt(rd());
    
    int V = 0;
    for (const auto &e : es){
      assert(e.fst != e.snd);
      V = std::max({V, e.fst + 1, e.snd + 1});
    }

    std::vector<double> centrality(V);
    std::vector<std::vector<int> > edge_groups = ArticulationPointDecomposition(es).edge_groups;
    
    for (const auto &edge_group : edge_groups) {
      if (edge_group.empty()) continue;
      
      std::vector<Edge> ccomp_es;
      for (int edge_id : edge_group){
        ccomp_es.emplace_back(es[edge_id].fst, es[edge_id].snd, edge_id);
      }
      std::vector<std::vector<PI> > g(BuildCompressedGraph(ccomp_es));
      std::vector<int>  order = DistanceOrdering(g);
      // 簡単のため全体の頂点数分の領域を確保している.二重連結成分が多いと効率が悪くなる.
      std::vector<bool> visit(V, false);
      std::vector<int>  next_edges(V, -1);
    
      // std::cout << "START ITERATION: " << num_samples << std::endl;
      for (int trial = 0; trial < num_samples; trial++) {
        fill(visit.begin(), visit.end(), false);
        SampleSpanningTree(next_edges, visit, g, order, mt);

        for (size_t s = 1; s < g.size(); s++) {
          int v = order[s];
          int e = ccomp_es[g[v][next_edges[v]].snd].edge_id;
          centrality[es[e].fst] += 1.0 / num_samples;
          centrality[es[e].snd] += 1.0 / num_samples;
        }
      }
    }
    return centrality;
  }
}


#endif /* SPANNING_CENTRALITY_H */

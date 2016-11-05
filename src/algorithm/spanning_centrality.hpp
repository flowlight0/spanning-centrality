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
  }

  
  class SpanningCentrality {
  private:
    std::vector<double> edge_centrality;
    std::vector<double> vertex_centrality;
    std::vector<double> aggregated_centrality;
    std::vector<std::pair<int, int> > original_edges;
    std::mt19937  mt;

  public:

    SpanningCentrality() {
      std::random_device rd;
      this->mt = std::mt19937(rd());
    }
      
    inline double GetEdgeCentrality(size_t edge_id) const {
      return edge_centrality.at(edge_id);
    }
      
    inline double GetVertexCentrality(size_t vertex_id) const {
      return vertex_centrality.at(vertex_id);
    }

    inline double GetAggregatedCentrality(size_t vertex_id) const {
      return aggregated_centrality.at(vertex_id);
    }

    inline size_t GetNumVertices() const {
      return vertex_centrality.size();
    }

    inline size_t GetNumEdges() const {
      return edge_centrality.size();
    }

    inline std::pair<int, int> GetEdge(size_t edge_id) const {
      return original_edges.at(edge_id);
    }

    bool Construct(const std::string &graph_file, int num_samples) {
      std::vector<std::pair<int, int> > es;
      if (ReadGraph(graph_file, es)) {
        ConvertToUndirectedGraph(es);
        return Construct(es, num_samples);
      } else {
        return false;
      }
    }

    bool Construct(const std::vector<std::pair<int, int> > &es, int num_samples) {
      this->original_edges = es;
      
      int V = 0;
      for (const auto &e : es){
        assert(e.fst != e.snd);
        V = std::max({V, e.fst + 1, e.snd + 1});
      }
        
      this->edge_centrality = std::vector<double>(es.size());
      this->vertex_centrality = std::vector<double>(V);
      this->aggregated_centrality = std::vector<double>(V);

      const ArticulationPointDecomposition decomposition = ArticulationPointDecomposition(es);
      const std::vector<std::vector<int> > edge_groups = decomposition.edge_groups;
      const std::vector<int> articulation_points = decomposition.articulation_points;
      
      // process each bi-connected component one by one
      for (const auto &edge_group : edge_groups) {
        if (edge_group.empty()) continue;

        // build a subgraph from an edge group and shuffle vertices on distances.
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
          
          SampleSpanningTree(next_edges, visit, g, order);
            
          for (size_t s = 1; s < g.size(); s++) {
            int v = order[s];
            int e = ccomp_es[g[v][next_edges[v]].snd].edge_id;

            this->edge_centrality[e]+= 1.0 / num_samples;

            if (++degree[es[e].fst] == 2){
              this->vertex_centrality[es[e].fst] += 1.0 / num_samples;
            }
            if (++degree[es[e].snd] == 2){
              this->vertex_centrality[es[e].snd] += 1.0 / num_samples;
            }
            
            this->aggregated_centrality[es[e].fst] += 1.0 / num_samples;
            this->aggregated_centrality[es[e].snd] += 1.0 / num_samples;
          }
        }
      }

      for (int v : articulation_points) {
        vertex_centrality[v] = 1.0;
      }

      for (int v = 0; v < V; v++) {
        assert(vertex_centrality[v] < 1 + 1e-9);
      }
      return true;
    }

  private:
    void SampleSpanningTree(std::vector<int > &next_edges,
                            std::vector<bool> &visit,
                            const std::vector<std::vector<PI> > &g,
                            const std::vector<int> &order)
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
  };
  
  
  
  std::vector<double> EstimateEdgeCentrality(const std::vector<PI> &es, int num_samples){
    SpanningCentrality spanning_centrality;
    spanning_centrality.Construct(es, num_samples);

    size_t E = spanning_centrality.GetNumEdges();
    std::vector<double> centrality(E);
    for (size_t e = 0; e < E; e++) {
      centrality[e] = spanning_centrality.GetEdgeCentrality(e);
    }
    return centrality;
  }
  
  
  std::vector<double> EstimateVertexCentrality(const std::vector<PI> &es, int num_samples){
    SpanningCentrality spanning_centrality;
    spanning_centrality.Construct(es, num_samples);

    size_t V = spanning_centrality.GetNumVertices();
    std::vector<double> centrality(V);
    for (size_t v = 0; v < V; v++) {
      centrality[v] = spanning_centrality.GetVertexCentrality(v);
    }
    return centrality;
  }
  

  std::vector<double> EstimateAggregatedCentrality(const std::vector<PI> &es, int num_samples){
    SpanningCentrality spanning_centrality;
    spanning_centrality.Construct(es, num_samples);

    size_t V = spanning_centrality.GetNumVertices();
    std::vector<double> centrality(V);
    for (size_t v = 0; v < V; v++) {
      centrality[v] = spanning_centrality.GetAggregatedCentrality(v);
    }
    return centrality;
  }
}

#endif /* SPANNING_CENTRALITY_H */

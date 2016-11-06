# Spanning Tree Centrality
Spanning tree centrality is one of centrality measures defined from randomly sampled spanning trees. This library provides a fast approximate implementation of spanning tree centrality of vertices and edges. It can also estimate the aggregated version of spanning tree centrality of vertices. 

## Usage
Given a graph and the number of spanning trees sampled for estimating centrality values, we compute approximate values of spanning tree centrality values. The more spanning trees we sample, the more accurate estimation we can make.

### From CLI
* Execute `./waf configure` to check if you have libraries required to use our library. 
* Execute `./waf` to build and test our library. 
* Execute `bin/spanning_edge_centrality` to estimate spanning tree centrality of edges. You need to specify the following options.  
* Execute `bin/spanning_vertex_centrality` to estimate spanning tree centrality of vertices.
* Execute `bin/spanning_vertex_centrality` with `--aggregated` option to estimate the aggregated version of spanning tree centrality of vertices. 

#### Example
    $ ./waf configure 
    $ ./waf 
    $ bin/spanning_edge_centrality --graph_file=sample-graph.txt --num_samples=10000 # estimate edge centrality
    0       1       1.000
    1       2       0.733
    1       4       0.736
    2       3       0.729
    3       4       0.581
    3       6       0.716
    4       5       0.715
    5       6       0.579
    5       7       0.703
    6       8       0.698
    7       8       0.541
    7       9       0.641
    8       9       0.627
    9       10      1.000
    $ bin/spanning_edge_centrality --graph_file=sample-graph.txt --num_samples=10000 # estimate vertex centrality
    0       0.000
    1       1.000
    2       0.460
    3       0.779
    4       0.788
    5       0.761
    6       0.764
    7       0.709
    8       0.704
    9       1.000
    10      0.000
    $ bin/spanning_vertex_centrality --aggregated --graph_file=sample-graph.txt --num_samples=10000 # estimate aggregated vertex centrality
    0       1.000
    1       2.458
    2       1.460
    3       2.047
    4       2.043
    5       1.985
    6       1.984
    7       1.881
    8       1.876
    9       2.267
    10      1.000
    
### From Your Program
You can use this library from your C++ program by including two header files `src/algorithm/articulation_point_decomposition.hpp` and `src/algorithm/spanning_centrality.hpp`. This library can be used as follows.
    
    spanning_centrality::SpanningCentrality spanning_centrality;
    spanning_centrality.Construct(FLAGS_graph_file, FLAGS_num_samples); // you can also pass vector<pair<int, int> > as a set of edges. 
    for (size_t i = 0; i < spanning_centrality.GetNumEdges(); i++) {
      double c = spanning_centrality.GetEdgeCentrality(i);
      int u = spanning_centrality.GetEdge(i).first;
      int v = spanning_centrality.GetEdge(i).second;
      printf("%d\t%d\t%.3lf\n", u, v, c);
    }
    
## Reference
Takanori Hayashi, Takuya Akiba, and Yuichi Yoshida, **[Efficient Algorithms for Spanning Tree Centrality](http://www.ijcai.org/Proceedings/16/Papers/525.pdf)**.
In IJCAI 2016.

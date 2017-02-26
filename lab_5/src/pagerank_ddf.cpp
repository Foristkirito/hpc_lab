//
// Created by 夏鑫 on 07/01/2017.
//

/*
Copyright (c) 2014-2015 Xiaowei Zhu, Tsinghua University

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include <string>
#include "core/graph.hpp"

#define DAMPING_FACTOR 0.85

int main(int argc, char ** argv) {
    if (argc<3) {
        fprintf(stderr, "usage: pagerank_ddf [path] [iterations] [memory budget in GB]\n");
        exit(-1);
    }
    std::string path = argv[1];
    int iterations = atoi(argv[2]);
    long memory_bytes = (argc>=4)?atol(argv[3])*1024l*1024l*1024l:8l*1024l*1024l*1024l;

    Graph graph(path);
    graph.set_memory_bytes(memory_bytes);
    BigVector<VertexId> degree(graph.path+"/degree_ddf", graph.vertices);
    BigVector<float> pagerank(graph.path+"/pagerank_ddf", graph.vertices);
    BigVector<float> delta(graph.path+"/delta_ddf", graph.vertices);
    BigVector<float> new_delta(graph.path+"/new_delta_ddf", graph.vertices);
    BigVector<float> tmp;

    long vertex_data_bytes = (long)graph.vertices * ( sizeof(VertexId) + sizeof(float) + sizeof(float) + sizeof(float));
    graph.set_vertex_data_bytes(vertex_data_bytes);

    double begin_time = get_time();

    degree.fill(0);
    graph.stream_edges<VertexId>(
        [&](Edge & e){
          write_add(&degree[e.source], 1);
          return 0;
        }, nullptr, 0, 0
    );
    printf("degree calculation used %.2f seconds\n", get_time() - begin_time);
    fflush(stdout);
    new_delta.fill(0);
    graph.hint(pagerank, delta);
    float init_rank = 1.f - DAMPING_FACTOR;
    graph.stream_vertices<VertexId>(
        [&](VertexId i){
          pagerank[i] = init_rank;
          delta[i] = 1.f / degree[i];
          return 0;
        }, nullptr, 0 ,
        [&](std::pair<VertexId,VertexId> vid_range){
          pagerank.load(vid_range.first, vid_range.second);
          delta.load(vid_range.first, vid_range.second);
        },
        [&](std::pair<VertexId,VertexId> vid_range){
          pagerank.save();
          delta.save();
        }
    );

    for (int iter=0;iter<iterations;iter++) {
        printf("iter : %d\n", iter);
        graph.hint(new_delta);
        graph.stream_edges<VertexId>(
            [&](Edge & e){
              write_add(&new_delta[e.target], delta[e.source]);
              return 0;
            }, nullptr, 0, 1,
            [&](std::pair<VertexId,VertexId> source_vid_range){
              new_delta.lock(source_vid_range.first, source_vid_range.second);
            },
            [&](std::pair<VertexId,VertexId> source_vid_range){
              new_delta.unlock(source_vid_range.first, source_vid_range.second);
            }
        );
        graph.hint(pagerank, delta);
        if (iter==iterations-1) {
            graph.stream_vertices<VertexId>(
                [&](VertexId i){
                  pagerank[i] = 0.15f + 0.85f * new_delta[i];
                  return 0;
                }, nullptr, 0,
                [&](std::pair<VertexId,VertexId> vid_range){
                  pagerank.load(vid_range.first, vid_range.second);
                },
                [&](std::pair<VertexId,VertexId> vid_range){
                  pagerank.save();
                }
            );
        } else {
            graph.stream_vertices<VertexId>(
                [&](VertexId i){
                  pagerank[i] = 0.15f + 0.85f * new_delta[i];
                  delta[i] = 0;
                  return 0;
                }, nullptr, 0,
                [&](std::pair<VertexId,VertexId> vid_range){
                  pagerank.load(vid_range.first, vid_range.second);
                  delta.load(vid_range.first, vid_range.second);
                },
                [&](std::pair<VertexId,VertexId> vid_range){
                  pagerank.save();
                  delta.save();
                }
            );
        }
        tmp = new_delta;
        new_delta = delta;
        delta = tmp;
    }

    double end_time = get_time();
    printf("%d iterations of pagerank_ddf took %.2f seconds\n", iterations, end_time - begin_time);

}

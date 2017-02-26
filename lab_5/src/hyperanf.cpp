//
// Created by 夏鑫 on 10/01/2017.
//

#include <string>
#include <core/graph.hpp>
#include "hyperloglog.h"

#define RSD 0.1
#define M 128

struct __attribute__((__packed__)) counter {
  unsigned char hll_ctr[M] HLL_ALIGN; //cal according to RSD
};

struct __attribute__((__packed__)) update_counter {
  bool change;
  unsigned char hll_ctr[M] HLL_ALIGN; //cal according to RSD
};



int main(int argc, char ** argv) {
    if (argc<3) {
        fprintf(stderr, "usage: hyperanf [path] [iterations] [memory budget in GB]\n");
        exit(-1);
    }
    int parallelism = std::thread::hardware_concurrency();
    std::string path = argv[1];
    int iterations = atoi(argv[2]);
    long memory_bytes = (argc>=4)?atol(argv[3])*1024l*1024l*1024l:8l*1024l*1024l*1024l;

    Graph graph(path);
    graph.set_memory_bytes(memory_bytes);
    hyper_log_log_params loglog_para;
    setup_hll_params(&loglog_para, RSD);
    int counter_size = sizeof_hll_counter(&loglog_para);
    BigVector<counter> llh(graph.path+"/llh", graph.vertices);
    BigVector<update_counter> update_llh(graph.path+"/update_llh", graph.vertices);
    BigVector<double> count(graph.path+"/count_llh", graph.vertices);
    printf("vertices: %d, edges: %d \n", graph.vertices, graph.edges);

    long vertex_data_bytes = (long)graph.vertices * (sizeof(counter) + sizeof(update_counter) + sizeof(double));
    graph.set_vertex_data_bytes(vertex_data_bytes);

    double begin_time = get_time();
    counter ct;
    memset(ct.hll_ctr, 0, M * sizeof(unsigned char));
    update_counter uct;
    memset(uct.hll_ctr, 0, M * sizeof(unsigned char));
    uct.change = true;
    llh.fill(ct);
    update_llh.fill(uct);
    count.fill(0.0);
    graph.hint(llh);
    graph.stream_vertices<VertexId>(
        [&](VertexId i){
          unsigned long hash = jenkins(i, 0xdeadbeef);
          unsigned char *llh_ct = llh[i].hll_ctr;
          unsigned char *update_ct = update_llh[i].hll_ctr;
          add_hll_counter(&loglog_para, llh_ct, hash);
          add_hll_counter(&loglog_para, update_ct, hash);
          return 0;
        }, nullptr, 0
    );
    double init_time = get_time();
    printf("init hll time: %.10f \n", init_time - begin_time);

    for (int iter=0;iter<iterations;iter++) {
        printf("itr : %d ---------------------\n", iter);
        graph.hint(update_llh);
        graph.stream_edges<VertexId>(
            [&](Edge & e){
              unsigned char *counter_src = llh[e.source].hll_ctr;
              unsigned char *counter_des = update_llh[e.target].hll_ctr;
              bool change = false;
              /*
              #pragma omp parallel for schedule(dynamic) num_threads(parallelism)
              for(unsigned int i=0;i<= loglog_para.m_minus_1;i++) {
                  if(counter_des[i] < counter_src[i]) {
                      change = true;
                      write_max(&counter_des[i], counter_src[i]);
                  }
              }
              */
              change = hll_union(counter_des, counter_src, &loglog_para);
              //printf("one edge done\n");
              update_llh[e.target].change = change;
              return 0;
            }, nullptr, 0, 1,
            [&](std::pair<VertexId,VertexId> source_vid_range){
              update_llh.lock(source_vid_range.first, source_vid_range.second);
            },
            [&](std::pair<VertexId,VertexId> source_vid_range){
              update_llh.unlock(source_vid_range.first, source_vid_range.second);
            }
        );
        graph.hint(llh, count);

        printf("stream edges done\n", iter);
        graph.stream_vertices<VertexId>(
            [&](VertexId i){
              if (update_llh[i].change){
                  memcpy(llh[i].hll_ctr, update_llh[i].hll_ctr, M * sizeof(unsigned char));
              }
              count[i] = count_hll_counter(&loglog_para, llh[i].hll_ctr);

              return 0;
            }, nullptr, 0,
            [&](std::pair<VertexId,VertexId> vid_range){
              llh.load(vid_range.first, vid_range.second);
              count.load(vid_range.first, vid_range.second);
            },
            [&](std::pair<VertexId,VertexId> vid_range){
              llh.save();
              count.save();
            }
        );
        printf("stream vertices done\n", iter);
    }

    double end_time = get_time();
    printf("%d iterations of hyperanf took %.2f seconds\n", iterations, end_time - begin_time);

}

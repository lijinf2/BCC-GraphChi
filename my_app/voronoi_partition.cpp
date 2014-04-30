#include <string>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

#define GRAPHCHI_DISABLE_COMPRESSION


#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"

using namespace graphchi;

bool        scheduler = false;

typedef int VertexDataType;
typedef int EdgeDataType;
#define SEED 0
#define NOT_SEED -1
#define NO_COLOR -2
inline bool HAS_COLOR(int x) { return x>0; }
static int numColors = 0;

struct VoronoiSampler : public GraphChiProgram<VertexDataType, EdgeDataType> {
    double sampleRate;
    VoronoiSampler(double _sampleRate) :sampleRate(_sampleRate){ }
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
        // sample itself
        if(rand()*1.0/RAND_MAX < sampleRate) {
            v.set_data(SEED);
            numColors++;
        } else {
            v.set_data(NOT_SEED);
        }
        for(int i=0; i < v.num_edges(); i++) {
            v.edge(i)->set_data(NO_COLOR);
        }
    }
};

struct VoronoiPartitioner : public GraphChiProgram<VertexDataType, EdgeDataType> {

    bool converged;

    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {

        // color myself with others' color
        if(!HAS_COLOR(v.get_data())) {
            converged = false;
            for(int i=0; i < v.num_edges(); i++) {
                if(HAS_COLOR(v.edge(i)->get_data())) {
                    v.set_data(v.edge(i)->get_data());
                    break;
                }
            }
        }

        if(v.get_data() == SEED) {              // If I'm a seed
            v.set_data(rand()); // seed 用自己id，随机可能重复
            converged = false;
        } else if(HAS_COLOR(v.get_data())) {    // If I've got a color
            int color = v.get_data();
            for(int i=0; i < v.num_edges(); i++) {
                if(!HAS_COLOR(v.edge(i)->get_data())) {
                    std::cout << "color uncolored edges" << std::endl;
                    v.edge(i)->set_data(color);
                    converged = false;
                }
            }
        }
    }

    void before_iteration(int iteration, graphchi_context &info) {
        converged = true;
    }
    // what if no converged?
    void after_iteration(int iteration, graphchi_context &ginfo) {
        if (converged) {
            std::cout << "Converged!" << std::endl;
            ginfo.set_last_iteration(iteration);
        }
    }

    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {}
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {}
};


int main(int argc, const char ** argv) {
    graphchi_init(argc, argv);
    metrics m("voronoi-partition");
    global_logger().set_log_level(LOG_DEBUG);
    srand(clock());

    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    double rate             = get_option_float("rate", 0.1);
    double niter            = get_option_int("niter", 1000);
    scheduler            = get_option_int("scheduler", false);

    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));

    /* Run */
    graphchi_engine<int, int> engine(filename, nshards, scheduler, m);
    engine.set_modifies_inedges(false); // Improves I/O performance.

    VoronoiSampler sampler(rate);
    engine.run(sampler, 1);

    VoronoiPartitioner partitioner;
    engine.run(partitioner, niter);

    std::cout<< "Number of colors:" << numColors << std::endl;
    metrics_report(m);
    return 0;
}


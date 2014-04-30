
/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.

 *
 * @section DESCRIPTION
 *
 * Template for GraphChi applications. To create a new application, duplicate
 * this template.
 */



#include <string>
#include <vector>
#include "graphchi_basic_includes.hpp"
#include "twodlabelanalysis.cpp"
#include <map>
#include <time.h> // for calculate time consumption
using namespace graphchi;

/**
  * Type definitions. Remember to create suitable graph shards using the
  * Sharder-program.
  */
const int MAXXSPLIT = 10000;
std::vector<int> xsplit;
std::vector< std::vector<int> > ysplits(MAXXSPLIT);

int         iterationcount = 0;// for connected component
bool        scheduler = false;// for connected component

struct point
{
	int x;
	int y;
	point(int x, int y)
	{
		this->x = x;
		this->y = y;
	}
};
struct x_less {
    bool operator()(point const& a, point const& b) const
    {
        if (a.x < b.x)
            return true;
        else
            return false;
    }
};
struct twodvertex
{
	int x;
	int y;
	vid_t twodblock;
	vid_t ccblock;
	twodvertex(){}
	twodvertex(int a,int b):x(a),y(b),twodblock(0),ccblock(0){}
	//twodvertex(int block){ this->ccblock = block;} // for twodlabelanalysis.cpp
	//void setTwodblock(int block){ this->twodblock = block;}
	//void setCCBlock(int block){ this->ccblock = block;}
	//overloaded operator to help with twodlabelanalysis.cpp
	twodvertex& operator=(const vid_t &b){
		this->ccblock = b;
		return *this;
	}
	twodvertex(vid_t b):ccblock(b){}
	friend std::ostream& operator<<(std::ostream &out, twodvertex &v)
	{
		out << v.ccblock <<" "<< v.twodblock;
		return out;
	}
};
// overloaded operator to help with twodlabelanalysis.cpp
/*
twodvertex& twodvertex::operator=(const vid_t &b)
{
	this->ccblock = b;
	return *this;
}*/
bool operator==(const twodvertex &a, const vid_t &b)
{
	return  a.ccblock == b;
}
bool operator==(const twodvertex &a, const twodvertex &b);
bool operator==(const twodvertex &a, const twodvertex &b)
{
	return a.ccblock == b.ccblock;
}
bool operator!=(const twodvertex &a, const twodvertex &b);
bool operator!=(const twodvertex &a, const twodvertex &b)
{
	return a.ccblock != b.ccblock;
}
bool operator<(const twodvertex &a, const twodvertex &b);
bool operator<(const twodvertex &a, const twodvertex &b)
{
	return a.ccblock < b.ccblock;
}
struct twodedge
{
	vid_t twodblock;
	vid_t label;
	int value;
	twodedge(){ twodblock = 0; label = 0; value = 0;}
	twodedge(int x){twodblock = 0; label = 0; value = x;}  // at least one argument to parse original edgedata
	twodedge(vid_t a, vid_t b,int c){twodblock = a; label = b; value = c;}
	//void setBlockLabel(int b,int x){ this->label = b; this->label = x;}
	//void setLabel(int l){this->label = l;}
};

FILE * inf;
char s[1024];
typedef twodvertex VertexDataType;
typedef twodedge EdgeDataType;

struct TwoDSetCoordinateProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {


    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &ginfo) {


        if (ginfo.iteration == 0 && vertex.id() != 0) {
            /* On first iteration, initialize vertex (and its edges). This is usually required, because
               on each run, GraphChi will modify the data files. To start from scratch, it is easiest
               do initialize the program in code. Alternatively, you can keep a copy of initial data files. */
            // vertex.set_data(init_value);

        	while(fgets(s,1024,inf) && s[0] != 'v');

        	//logstream(LOG_ERROR) << "Current line: \"" << s << "\"\n";
            int len = (int) strlen(s)-1;
            if(s[len] == '\n') s[len] = 0;

        	char delimits[] = " ";
        	char *t = NULL;
        	t = strtok(s,delimits);    // ignore v

            if (t == NULL) {
                logstream(LOG_ERROR) << "Input file is not in right format. "
                << "Expecting \"<from>\t<to>\". "
                << "Current line: \"" << s << "\"\n";
                assert(false);
            }

            t = strtok(NULL,delimits);

            vid_t vertexid = atoi(t);
            t = strtok(NULL,delimits);

            int xcoordinate = atoi(t);
            t = strtok(NULL,delimits);
            int ycoordinate = atoi(t);

            if( vertex.id() == vertexid )
            {
            	vertex.set_data(twodvertex(xcoordinate,ycoordinate));
            }
            else{
                logstream(LOG_ERROR) << "Vertex can not fit:  "<<"\n"
                << "vertexid in update is: "<<vertex.id()<<"\n"
                << "Current line: \"" << s << "\"\n"
                << "vertex id input: "<<vertexid<<"\n"
                << "x coordinate:"<<xcoordinate<<"\n"
                << "y coordinate:"<<ycoordinate<<"\n";
                char tmp;
                std::cin>>tmp;
            }

        } else { // do nothing
            /* Do computation */

            /* Loop over in-edges (example) */
            for(int i=0; i < vertex.num_inedges(); i++) {
                // Do something
            //    value += vertex.inedge(i).get_data();
            }

            /* Loop over out-edges (example) */
            for(int i=0; i < vertex.num_outedges(); i++) {
                // Do something
                // vertex.outedge(i).set_data(x)
            }

            /* Loop over all edges (ignore direction) */
            for(int i=0; i < vertex.num_edges(); i++) {
                // vertex.edge(i).get_data()
            }

            // v.set_data(new_value);
        }
    }

    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
    }

    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &gcontext) {
    }

    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {
    }

    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {
    }

};

/**
 *
 *
 */


struct TwoDSampleProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {


    /**
     *  Vertex update function.
     */
	int xnum;
	int ynum;
	float sampleRate;
	std::vector<point> samples;

	TwoDSampleProgram(int xnum, int ynum, float sampleRate)
	{
		srand( (unsigned) time(NULL));
		this->xnum = xnum;
		this->ynum = ynum;
		this->sampleRate = sampleRate;
		if(xnum > MAXXSPLIT)
		{
			logstream(LOG_ERROR) << "xnum is larger than MAXXSPPLIT, please reset the max x split in"
					             << " twopartition.cpp"<<"\n";
            assert(false);
		}
	}

    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &ginfo) {

        if (ginfo.iteration == 0 && vertex.id() != 0) {
            /* On first iteration, initialize vertex (and its edges). This is usually required, because
               on each run, GraphChi will modify the data files. To start from scratch, it is easiest
               do initialize the program in code. Alternatively, you can keep a copy of initial data files. */
            // vertex.set_data(init_value);
        	float samp = ((float) rand()) / RAND_MAX;
        	/*
        	std::cout<<"samp: "<<samp<<" sampleRate:"<<sampleRate<<"\n";
        	char tmp;
        	std::cin>>tmp;*/
        	if( samp <= sampleRate)
        	{
        		point p(vertex.get_data().x, vertex.get_data().y);
        		samples.push_back(p);
        	}
        	/* for debug only
			logstream(LOG_ERROR) << "vertex id:"<<vertex.id()<<" and x : "<<vertex.get_data().x
					             << " and y: "<<vertex.get_data().y<<"\n";*/

        } else {
            /* Do computation */

            /* Loop over in-edges (example) */
            for(int i=0; i < vertex.num_inedges(); i++) {
                // Do something
            //    value += vertex.inedge(i).get_data();
            }

            /* Loop over out-edges (example) */
            for(int i=0; i < vertex.num_outedges(); i++) {
                // Do something
                // vertex.outedge(i).set_data(x)
            }

            /* Loop over all edges (ignore direction) */
            for(int i=0; i < vertex.num_edges(); i++) {
                // vertex.edge(i).get_data()
            }

            // v.set_data(new_value);
        }
    }

	void getsplits()
	{
	            //sort by x-coord
	            sort(samples.begin(), samples.end(), x_less());
	            //get splits here
	            int size = samples.size();
	            int residual = size % xnum;
	            int step;
	            if (residual == 0)
	                step = size / xnum;
	            else
	                step = size / xnum + 1;
	            for (int pos = step - 1; pos < size; pos += step)
	                xsplit.push_back(samples[pos].x);
	            //------
	            std::vector< std::vector<int> > subSamps(xnum);
	            for (int i = 0; i < xnum; i++) {
	                int start = step * i;
	                int end = start + step;
	                if (end > size)
	                    end = size;
	                for (int j = start; j < end; j++) {
	                    subSamps[i].push_back(samples[j].y);
	                }
	                sort(subSamps[i].begin(), subSamps[i].end());
	                int sz = subSamps[i].size();
	                int res = sz % ynum;
	                int ystep;
	                if (res == 0)
	                    ystep = sz / ynum;
	                else
	                    ystep = sz / ynum + 1;
	                for (int pos = ystep - 1; pos < sz; pos += ystep)
	                    ysplits[i].push_back(subSamps[i][pos]);
	            }

	}

    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &ginfo) {
    }

    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &ginfo) {
    	if(ginfo.iteration == 0)
    	{
    		//std::cout<<"sample size:"<<samples.size()<<"\n";
    		//char tmp;
    		//std::cin>>tmp;
            getsplits();
    	}
    }

    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {
    }

    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {
    }

};


/**
  * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type>
  * class. The main logic is usually in the update function.
  */
struct TwoDPartitionProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {



    /**
     *  Vertex update function.
     */
	int xnum;
	int ynum;
	TwoDPartitionProgram(int xnum, int ynum)
	{
		this->xnum = xnum;
		this->ynum = ynum;
	}
    int getXid(double x)
    {
        int size = xsplit.size();
        if (size > xnum)
            size = xnum;
        for (int i = 0; i < size; i++) {
            if (x <= xsplit[i])
                return i;
        }
        return xnum - 1;
    } //may improve by binary search

    int getYid(int xid, int y)
    {
        int size = ysplits[xid].size();
        if (size > ynum)
            size = ynum;
        for (int i = 0; i < size; i++) {
            if (y <= ysplits[xid][i])
                return i;
        }
        return ynum - 1;
    } //may improve by binary search

    int getBlkID(int x, int y)
    {
        int xid = getXid(x);
        int yid = getYid(xid,y);
        int ynum = xsplit.size();
        return xid * ynum + yid;
    }
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &ginfo) {

        if (ginfo.iteration == 0 && vertex.id()!= 0) {
            /* On first iteration, initialize vertex (and its edges). This is usually required, because
               on each run, GraphChi will modify the data files. To start from scratch, it is easiest
               do initialize the program in code. Alternatively, you can keep a copy of initial data files. */
            // vertex.set_data(init_value);
        	VertexDataType v = vertex.get_data();
        	//int blocknum = getBlkID(vertex.get_data().x,vertex.get_data().y);
        	v.twodblock = getBlkID(vertex.get_data().x,vertex.get_data().y);
        	vertex.set_data( v );

        	/*

        	std::cout<<" twodblock of vertex"<<vertex.id()<<" is "<<vertex.get_data().twodblock
        			<<" and his blocknum is"<<blocknum<<"\n";

        	if(vertex.id() % 1000 == 0){ char tmp; std::cin>>tmp; std::cout<<"\n";}
        	*/
        	//store vertexid blockid into file


        } else {
            /* Do computation */

            /* Loop over in-edges (example) */
            for(int i=0; i < vertex.num_inedges(); i++) {
                // Do something
            //    value += vertex.inedge(i).get_data();
            }

            /* Loop over out-edges (example) */
            for(int i=0; i < vertex.num_outedges(); i++) {
                // Do something
                // vertex.outedge(i).set_data(x)
            }

            /* Loop over all edges (ignore direction) */
            for(int i=0; i < vertex.num_edges(); i++) {
                // vertex.edge(i).get_data()
            }

            // v.set_data(new_value);
        }
    }

    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
    }

    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &gcontext) {
    }

    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {
    }

    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {
    }

};


struct ConnectedComponentsProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {

    bool converged;

    /**
     *  Vertex update function.
     *  On first iteration ,each vertex chooses a label = the vertex id.
     *  On subsequent iterations, each vertex chooses the minimum of the neighbor's
     *  label (and itself).
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

    	if( vertex.id() ==0 ) return;

    	VertexDataType vdata;
    	EdgeDataType edata;

        if (scheduler) gcontext.scheduler->remove_tasks(vertex.id(), vertex.id());

        if (gcontext.iteration == 0) {
        	vdata = vertex.get_data();
        	vdata.ccblock = vertex.id();
            vertex.set_data(vdata);
            if (scheduler)  gcontext.scheduler->add_task(vertex.id());
        }

        /* On subsequent iterations, find the minimum label of my neighbors */
        vid_t curmin = vertex.get_data().ccblock;
        for(int i=0; i < vertex.num_edges(); i++) {
        	if( gcontext.iteration != 0 && vertex.edge(i)->get_data().twodblock != vertex.get_data().twodblock) continue;
            vid_t nblabel = vertex.edge(i)->get_data().label;
            if (gcontext.iteration == 0) nblabel = vertex.edge(i)->vertex_id();  // Note!
            curmin = std::min(nblabel, curmin);
        }

        /* Set my label */
        vdata = vertex.get_data();
        vdata.ccblock = curmin;
        vertex.set_data(vdata);

        /**
         * Broadcast new label to neighbors by writing the value
         * to the incident edges.
         * Note: on first iteration, write only to out-edges to avoid
         * overwriting data (this is kind of a subtle point)
         */
        vid_t bigblock = vertex.get_data().twodblock;
        vid_t label = vertex.get_data().ccblock;

        if (gcontext.iteration > 0) {
            for(int i=0; i < vertex.num_edges(); i++) {
            	if( vertex.edge(i)->get_data().twodblock != bigblock) continue;
                if (label < vertex.edge(i)->get_data().label) {
                	edata = vertex.edge(i)->get_data();
                	edata.label = label;
                    vertex.edge(i)->set_data(edata);
                    /* Schedule neighbor for update */
                    if (scheduler) gcontext.scheduler->add_task(vertex.edge(i)->vertex_id(), true);
                    converged = false;
                }
            }
        } else if (gcontext.iteration == 0) {
            for(int i=0; i < vertex.num_outedges(); i++) {
            	edata = vertex.edge(i)->get_data();
            	edata.twodblock = bigblock;
            	edata.label = label;
                vertex.outedge(i)->set_data(edata);
            }
        }
    }
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &info) {
        iterationcount++;
        converged = iteration > 0;
    }

    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &ginfo) {
        if (converged) {
            std::cout << "Converged!" << std::endl;
            ginfo.set_last_iteration(iteration);
        }
    }

    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {
    }

    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {
    }

};

//This is for block connected component
std::map< vid_t, vid_t> count; // label,number of vertices. will merge all the cc with big id into cc with small id
std::map< vid_t, vid_t> father; // label, father of label;
struct BlockConnectedComponent : public GraphChiProgram<VertexDataType, EdgeDataType> {

	bool converged;
	vid_t findfather(vid_t label)
	{

		if( father[label] == label){
			return label;
		}else{
			return father[label] = findfather(father[label]) ;
		}
	}

    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &ginfo) {

        if (ginfo.iteration == 0 && vertex.id() != 0) {

            // vertex.set_data(init_value);
        	vid_t mylabel = vertex.get_data().ccblock;
            for(int i=0; i < vertex.num_outedges(); i++) {
                // Do something
            	EdgeDataType e = vertex.outedge(i)->get_data();
            	e.label = mylabel;
                vertex.outedge(i)->set_data(e);
            }
            if( count.find(mylabel) == count.end() )
            {
            	count[mylabel] = 1;
            	father[mylabel] = mylabel;
            }
            else count[mylabel]++;



        } else {
            // Do computation

        	vid_t mylabel = vertex.get_data().ccblock;

            vid_t edgelabel ;
            //Loop over in-edges (example)
            for(int i=0; i < vertex.num_edges(); i++) {
               // minlabel =  std::min(minlabel,vertex.edge(i)->vertex_id());
                edgelabel = vertex.edge(i)->get_data().label;

                if( edgelabel != mylabel) // There is label smaller than me
                {
                	vid_t edgefather = findfather(edgelabel);
                	vid_t myfather = findfather(mylabel);
                	if(edgefather != myfather ) // minfather wmay be larger than my father
                	{
                		if(edgefather < myfather)
                		{
                			father[myfather] =  edgefather;
                			count[edgefather] += count[myfather];
                		}
                		else
                		{
                			father[edgefather] = myfather;
                			count[myfather] += count[edgefather];
                		}
                		converged = false;
                	}

                }
            }
            // v.set_data(new_value);


        }

    }


    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &info) {
        iterationcount++;
        converged = iteration > 0;
    }

    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &ginfo) {
        if (converged) {
            std::cout << "Block Converged!" << std::endl;
            ginfo.set_last_iteration(iteration);
        }
    }
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {
    }


    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {
    }

};
// for Block connected component;
struct BlockLabel
{
	int count;
	int label;
	BlockLabel(int a, int b):count(a),label(b){};
	BlockLabel(){}
};
bool cmp(BlockLabel a,BlockLabel b){
	return a.count > b.count;
}
void analyze_block_labels(){
	std::map<vid_t,vid_t>::iterator it = count.begin();
	std::vector<BlockLabel> result;
	for(; it != count.end() ;++it)
	{
		if(father[  it->first ] == it->first) result.push_back( BlockLabel(it->second,it->first) );
	}
	sort(result.begin(),result.end(),cmp);
	std::cout<<" Label , "<<" size"<<"\n";
	for(unsigned int i = 0; i < 20 && i < result.size(); ++i)
	{
		std::cout<<result[i].label<<" , "<<result[i].count<<"\n";
	}
}
int main(int argc, const char ** argv) {
    /* GraphChi initialization will read the command line
       arguments and the configuration file. */
    time_t begin = clock();
    graphchi_init(argc, argv);

    /* Metrics object for keeping track of performance counters
       and other information. Currently required. */
    metrics m("my-application-name");

    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename

    scheduler       = get_option_int("scheduler", 0); // Whether to use selective scheduling

    int xnum = get_option_int("xnum",20);
    int ynum = get_option_int("ynum",20);
    if(ynum <= MAXXSPLIT) ysplits.resize(ynum);
    float sampleRate = get_option_float("sampleRate",0.01);

    char tmp;
    /* Detect the number of shards or preprocess an input to create them */

    int nshards          = convert_if_notexists<int, EdgeDataType>(filename,
                        		   	   	   	   	   	   	   	get_option_string("nshards", "0"));
    time_t verybegin = clock();
    time_t end = clock();
    std::cout<<"time consumption for sharding:"<< double(end -begin)/CLOCKS_PER_SEC<<"\n";
    std::cin>>tmp;
    begin = clock();

    std::string coordinatefile = get_option_string("coordinate");
    inf = fopen(coordinatefile.c_str(), "r");
    if(inf == NULL)
    {
		logstream(LOG_ERROR) << "Can not open file"<<coordinatefile.c_str()<<"\n";
		assert(false);
    }
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m);
    graphchi_engine<VertexDataType, EdgeDataType> setcoordinateengine(filename, nshards, scheduler, m);
    /* Run for Setting Coordinate*/
    TwoDSetCoordinateProgram setcoordinate;
    setcoordinateengine.set_exec_threads(1);
    setcoordinateengine.run(setcoordinate,1);
    fclose(inf);
    /* Run for Sampling*/
    TwoDSampleProgram sample(xnum,ynum,sampleRate);
    engine.run(sample,1);
    /* Run for Partition*/
    TwoDPartitionProgram partition(xnum,ynum);
    engine.run(partition, 1);
    /* Report execution metrics */
    std::cout<<" sampleRate is :"<<sampleRate<<std::endl;
    // debug only
    std::cout << "xsplit = [";
    for (unsigned int i = 0; i < xsplit.size(); i++)
        std::cout << xsplit[i] << " ";
    std::cout << "]" << "\n";
    for (unsigned int i = 0; i < ysplits.size(); i++) {
        std::cout << "ysplit[" << i << "] = [";
        for (unsigned int j = 0; j < ysplits[i].size(); j++)
            std::cout << ysplits[i][j] << " ";
        std::cout << "]" << "\n";

    }

    end = clock();
    std::cout<<"time consumption for getting splits:"<<double(end -begin)/CLOCKS_PER_SEC<<"\n";
    std::cin>>tmp;
    begin = clock();

    int niters           = get_option_int("niters", 1000); // Number of iterations

    ConnectedComponentsProgram ccprogram;
    engine.run(ccprogram, niters);
    analyze_labels<VertexDataType>(filename);



    end = clock();
    std::cout<<"time consumption for cc in different block:"<<double(end -begin)/CLOCKS_PER_SEC<<"\n";
    std::cin>>tmp;
    begin = clock();
// you have to disable scheduler first
    scheduler = false;
    graphchi_engine<VertexDataType, EdgeDataType> enginewithoutscheduler(filename, nshards, scheduler, m);
    BlockConnectedComponent blockprogram;
    enginewithoutscheduler.run(blockprogram,10);

    end = clock();
    std::cout<<"time consumption for block connected component:"<<double(end -begin)/CLOCKS_PER_SEC<<"\n";
    std::cin>>tmp;
    begin = clock();

    analyze_block_labels();
    time_t veryend = clock();
    std::cout<<"time consumption for BCC algorithm:"<<double(veryend - verybegin)/CLOCKS_PER_SEC<<"\n";
    std::cin>>tmp;
    begin = clock();
    metrics_report(m);

    return 0;
}

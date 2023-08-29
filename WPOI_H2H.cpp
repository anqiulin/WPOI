#include<queue>
#include<unordered_set>
#include<cfloat>
#include<climits>
#include<metis.h>
#include<memory.h>
#include<unordered_map>
#include<map>
#include<set>
#include<deque>
#include<stack>
#include<algorithm>
#include<sys/time.h>
#include<string.h>
#include<iostream>
#include<cstdio>
#include<cstring>
#include<iostream>
#include<cstdlib>
#include<vector>
#include <xmmintrin.h>
#include<cmath>
#include<bitset>
#include<limits.h>
#include<fstream>
#include <chrono>
#include <time.h>
using namespace std;

const char *FILE_INDEX = "./map/data_NY/h2h/NY.index";

double GetTime(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

 //NY-Density: 2%-158 4%-317 8%-634	16%-1268 32%-2537
// NY-KIND: 3%-634 6%-1268 9%-1903 12%-2537 15%-3172
static const int KIND = 634;	//the number of keywords kind
const int infinity = 999999999;
const int SIZEOFINT = 4;

// #define FILE_NODE "./map/data_cal/sub/cal_sub_keys_5_100.node"
#define FILE_EDGE "./map/data_NY/tol/NY_deal.edge"

typedef struct{
	double x,y;
	vector<int> adjnodes;
	vector<int> adjweight;
	bitset<KIND> keyword;
}Node;

int *toRMQ, *height, **RMQIndex;
int *belong;
int root, TreeSize;
int **rootToRoot, *rootSite;
int **dis, **pos, **pos2;
int *posSize, *pos2Size;
int *chSize;
int ** ch;
int *LOG2, *LOGD; 
int rootSize;
int *DFSList, *toDFS;
int ***BS;
int *Degree;
int **Neighbor, **Weight;
inline	int LCAQuery(int _p, int _q){
		int p = toRMQ[_p], q = toRMQ[_q];
		
		if (p > q){
			int x = p;
			p = q;
			q = x;
		}
		int len = q - p + 1;
		
		int i = LOGD[len], k = LOG2[len];
		
		q = q - i + 1;
		if (height[RMQIndex[k][p]] < height[RMQIndex[k][q]])
			return RMQIndex[k][p];
		else return RMQIndex[k][q]; 
	}
	
long long queryCnt;	
//long long aCnt;
inline	int distanceQuery(int p, int q){
		if (p == q) return 0;
		int x = belong[p], y = belong[q];	
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y){
			queryCnt++;
	//		aCnt++;
			if (lca == y){
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			return dis[y][pos[x][posSize[x] - 1]];
		}
		else {
			int res = infinity;
			int *dx = dis[x], *dy = dis[y],*p2 = pos2[lca];
			_mm_prefetch(dx, _MM_HINT_T0);
			_mm_prefetch(dy, _MM_HINT_T0);
			_mm_prefetch(p2, _MM_HINT_T0);
			int ps = pos2Size[lca];
			for (int i = 0; i < ps; i++){
				queryCnt ++;
				int tmp = dx[p2[i]] + dy[p2[i]];
				if (res > tmp)
					res = tmp;
			}
			return res;
		}
		
	}
	
int H2HPath(int ID1, int ID2, vector<int>& vPath)//, vector<int>& vPathEdge)
{

	int p = ID1+1, q = ID2+1;
	int res = distanceQuery(p, q);
	vPath.push_back(p-1);
	int tn = p;
        while (p != q){ 
        	int pq = distanceQuery(p, q); 
                for (int i = 0; i < Degree[p]; i++){
  			int x = Neighbor[p][i];
                        int xq = distanceQuery(x, q); 
                        if (xq + Weight[p][i] == pq){
 	                       p = x;
                               vPath.push_back(p-1);
		       	       //vPathEdge.push_back(adjListEdge[tn-1][i].second);
			    // cout<<tn-1<<": "<<p-1<<" = "<<adjListEdge[tn-1][i].second<<endl;
			       tn = x; 
                               break;
                         }
                 }
        }
        return res;
}
	
	void readGraph(char *filename){
		FILE * file = fopen(filename, "r");
		int n, m;
		//cout << filename <<endl;
		//fstream c("testData/graph.txt",ios::out);
		//c <<"111";
		//cout << file <<endl;
		fscanf(file, "%d %d", &n, &m);
		Degree = (int*)malloc(sizeof(int) * (n + 1));
		vector< vector<pair<int, int> > > nb;
		vector<pair<int, int> > v;
		v.clear();
		for (int i = 0; i <= n; i++){
		//	Degree[i] = 0;
			nb.push_back(v);	
		}
	//	cout << n << " " << m << endl;
		for (int i = 0; i < m; i++){
			int x, y, z;
			fscanf(file, "%d %d %d", &x, &y, &z);
	//		Degree[x]++;
	//		cout << x << " " << y << " " << z << endl;
			nb[x].push_back(make_pair(y, z));
		}
		Neighbor = (int**)malloc(sizeof(int*) * (n + 1));
		Weight = (int**)malloc(sizeof(int*) * (n + 1));
		for (int i = 1; i <= n; i++){
			Degree[i] = nb[i].size();
			Neighbor[i] = (int*)malloc(sizeof(int) * nb[i].size());
			Weight[i] = (int*)malloc(sizeof(int) * nb[i].size());
			for (int j = 0; j < nb[i].size(); j++){
				Neighbor[i][j] = nb[i][j].first;
				Weight[i][j] = nb[i][j].second;
			}
		}
	}
	inline int shortestPathQuery(int p, int q){
		int res = 0;
		while (p != q){
			res++;
			int pq = distanceQuery(p, q);
			for (int i = 0; i < Degree[p]; i++){
				int x = Neighbor[p][i];
			//	int y = Weight[p][i];
				int xq = distanceQuery(x, q);
				if (xq + Weight[p][i] == pq){
					p = x;
					break;
				}
			} 
		}
		return res;
	}

	FILE *fin;
	string TT = "";
	void scanIntArray(int *a, int n){
		fread(a, SIZEOFINT, n, fin);
	}
	int* scanIntVector(int *a){
		int _n;
		fread(&_n, SIZEOFINT, 1, fin);
		a = (int*)malloc(sizeof(int) * _n);
		scanIntArray(a, _n);
		return a;
	}

	int n;
	int *EulerSeq;
	void readIndex(){
		double _time = GetTime();
		int tree_height = 0, tree_width = 0, most_sp = 0;
		fin = fopen(FILE_INDEX, "rb");
		fread(&n, SIZEOFINT, 1, fin);
		int ts;
		fread(&ts, SIZEOFINT, 1, fin);
		TreeSize = ts;
		height = (int*)malloc(sizeof(int) * (ts + 1));
		for (int i = 0; i < ts; i++){
			fread(&height[i], SIZEOFINT, 1, fin);
		}
		belong = (int*)malloc(sizeof(int) * (n + 1));
	  	fread(belong, SIZEOFINT, n + 1, fin);
		toRMQ = (int*)malloc(sizeof(int) * (n + 1));
	  	fread(toRMQ, SIZEOFINT, n + 1, fin);
		int ris;
		fread(&ris, SIZEOFINT, 1, fin);
		fread(&ts, SIZEOFINT, 1, fin);
		EulerSeq = (int*)malloc(sizeof(int) * (ts + 1));
		RMQIndex = (int**)malloc(sizeof(int*) * (ris + 1));
		for (int i = 0; i < ris; i++){
			RMQIndex[i] = scanIntVector(RMQIndex[i]);
		}
		fread(&root, SIZEOFINT, 1, fin);
		cout << "root: " << root << endl;
		
		posSize = (int*)malloc(sizeof(int) * (n + 1));
		pos2Size = (int*)malloc(sizeof(int) * (n + 1));
		pos = (int**)malloc(sizeof(int*) * (TreeSize));
		pos2 = (int**)malloc(sizeof(int*) * (TreeSize));
		dis = (int**)malloc(sizeof(int*) * (TreeSize));
		chSize = (int*)malloc(sizeof(int) * (TreeSize));
		ch = (int**)malloc(sizeof(int*) * (TreeSize));
			
		for (int i = 0; i < TreeSize; i++){
			fread(&chSize[i], SIZEOFINT, 1, fin);
			ch[i] = (int*)malloc(sizeof(int) * chSize[i]);
			for (int j = 0; j < chSize[i]; j++){
				int x;
				fread(&x, SIZEOFINT, 1, fin);
				ch[i][j] = x;
			}
		}
		for (int i = 0; i < TreeSize; i++){
			int x;
			fread(&x, SIZEOFINT, 1, fin);
			fread(&posSize[x], SIZEOFINT, 1, fin);
			pos[x] = (int*)malloc(sizeof(int) * (posSize[x] + 1));
			fread(pos[x], SIZEOFINT, posSize[x], fin);
			if (posSize[x] > tree_width)
				tree_width = posSize[x];
			int _n;
			fread(&_n, SIZEOFINT, 1, fin);
			dis[x] = (int*)malloc(sizeof(int) * _n);
			fread(dis[x], SIZEOFINT, _n, fin);
			if (_n > tree_height)
				tree_height = _n;
		}
		printf("dis read finished!\n");
		for (int i = 0; i < TreeSize; i++){
			int x;
			fread(&x, SIZEOFINT, 1, fin);
			fread(&pos2Size[x], SIZEOFINT, 1, fin);
			pos2[x] = (int*)malloc(sizeof(int) * (pos2Size[x] + 1));
			fread(pos2[x], SIZEOFINT, pos2Size[x], fin);
			if (pos2Size[x] > most_sp)
				most_sp = pos2Size[x];
		}
		
		fclose(fin);
		//
		printf("Load Index Time : %lf sec\n", (GetTime() - _time));
		printf("tree height: %d\n", tree_height);
		printf("tree width: %d\n", tree_width);
		printf("most search space: %d\n", most_sp);
		//
	}

int noe;
vector<vector<int>> RSP;
vector<Node> Nodes;

bitset<KIND> set_str_to_binary(string ss){
	bitset<KIND>binary;
		char str[200]; //according to the keywordstr.length
			strcpy(str,ss.c_str());
			const char * split = ",";
    		char * p;
    		p = strtok (str,split);
    		while(p!=NULL) {
			        if(strcmp(p,"-1") ==0)
						break;
					binary.set(atoi(p));   
        		p = strtok(NULL,split);
      		}
        return binary;
}

// set binary_keyword of node
void binaryKeyword(Node &node,string ss)
{
			char str[200];// note the size of str,maybe not save(using method strcpy()) and lead to buffer overflow
			strcpy(str,ss.c_str());
			const char * split = ",";
    		char * p;
    		p = strtok (str,split);
    		while(p!=NULL) {
			        if(strcmp(p,"-1") ==0)
						break;
					node.keyword.set(atoi(p));   
        		p = strtok(NULL,split);
        }
} 

struct PA
{
	vector<int> pathnode;
	int dis;
	int park;
	string path;
	PA(vector<int>p,int x,int y){
		pathnode = p;
		dis = x;
		park = y;
		for(int i=0;i<p.size();i++)
		{
		        path+=to_string(p[i]);
		        path+="->";
		}
	}
	bool operator<(const PA& a)const
	{
		if (park != a.park)
			return park > a.park;
		return dis > a.dis;
	}
};
priority_queue<PA>qu; 
vector<PA>R;//store result

string querystr;
bitset<KIND> querybinary;

void init(){
	FILE *fin;

	// load node
	
	// fstream fp("./map/cal/cal_keys_20_100.node");
	// fstream fp("./map/NY/NY_keys.node");
	 // fstream fp("./map/data_COL/tol/Density8/COL_keys_15_3%.node");
	 fstream fp("./map/data_NY/tol/Density8/NY_keys_15_3%.node");
	 // fstream fp("./map/data_cal/tol/Density8/cal_keys_15_3%.node");

	cout<<"LOADING NODE..."<<endl;
	int nid;
	double x,y;
	string ss;
    for(int i = 0;i < KIND; i++){
                vector< int >tmp;
                RSP.push_back(tmp);
        }
	while(fp>>nid>>x>>y>>ss){
		Node node = {x,y};
		node.keyword.reset();
		binaryKeyword(node,ss);
        for(int i = 0;i < RSP.size(); i++){
	        if(node.keyword.test(i))
        	    RSP[i].push_back(nid);
        }
		Nodes.push_back(node);
	}
	fp.close();
	printf("COMPLETE. NODE_COUNT=%d\n", (int)Nodes.size());


	// load edge

	printf("LOADING EDGE...");
	fin = fopen(FILE_EDGE, "r");
	int eid;
	int snid, enid;
	//double weight;
	int weight;
	int iweight;
	noe = 0;
	//while( fscanf(fin,"%d %d %d %lf", &eid, &snid, &enid, &weight ) == 4 ){
	while( fscanf(fin,"%d %d %d %d", &eid, &snid, &enid, &weight ) == 4 ){
		noe ++;
		// iweight = (int) (weight * WEIGHT_INFLATE_FACTOR );
		iweight = (int) (weight);
		Nodes[snid].adjnodes.push_back( enid );
		Nodes[snid].adjweight.push_back( iweight );
		Nodes[enid].adjnodes.push_back( snid );
		Nodes[enid].adjweight.push_back( iweight );
	}
	fclose(fin);
	printf("COMPLETE. EDGE_COUNT=%d\n",noe);

}
vector<int> set_str_to_vec(string ss){

	vector<int> ret;
	ret.clear();
		char str[100]; //according to the keywordstr.length
			strcpy(str,ss.c_str());
//			cout<<str<<endl;
			const char * split = ",";
    		char * p;
    		p = strtok (str,split);

    		while(p!=NULL) {
			        if(strcmp(p,"-1") ==0)
						break;
					ret.push_back(atoi(p));   
        		p = strtok(NULL,split);
        }
        return ret;
}
inline int dijkstra_p2p(int s, int t) {
    unordered_map<int, int> result;
    result[s] = 0;

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.emplace(0, s);//push 的操作可以直接用于emplace： 

    int min, minpos, adjnode, weight;
    // TIME_TICK_START

    while (!pq.empty()) {
        min = pq.top().first;
        minpos = pq.top().second;

        if (minpos == t) {
            // TIME_TICK_END
		// TIME_TICK_PRINT("P2P")
            return min;
        }

        pq.pop();

        for (int i = 0; i < Nodes[minpos].adjnodes.size(); i++) {

            adjnode = Nodes[minpos].adjnodes[i];
            weight = Nodes[minpos].adjweight[i];

            if (result.find(adjnode) == result.end() || result[adjnode] > min + weight) {
                result[adjnode] = min + weight;
                pq.emplace(min + weight, adjnode);
            }
        }
    }

    return -1;
}
vector<vector<int>> expand_wpoi(int s, int t)
{

	vector <int> qu = set_str_to_vec(querystr); 

	unordered_map<int,pair<double,double>> p;
	p.clear();
	for(auto &c : qu)
	{
        for(auto &k : RSP[c])
        {
            double pdist =  (double)(((distanceQuery(s+1, k+1) + distanceQuery(k+1, t+1)) - distanceQuery(s+1, t+1)) *1.0 / distanceQuery(s+1, t+1) *1.0);
		    bitset <KIND> tmp(querybinary & Nodes[k].keyword);
		    double pstop = (double)(((querybinary.count() - tmp.count()) *1.0) / (querybinary.count() *1.0));
		    // if(pdist == 0)
		    // {
		    // 	pdist = 0.1;
		    // }
		    p[k] = make_pair(pdist, pstop);
		    // cout << "id:" << k<<  " pdist:" << pdist << " pstop:" << pstop << endl;
        }
		
	}

	vector<double> w;
	w.clear();

	for(int i = 0; i < 100; i++)
	{
		double t = double (i*1.0 /100.0);
		w.push_back(t);
	}

	vector<vector<int>> ret;
	ret.clear();

   
	for(int wp = 0;  wp < w.size(); wp++)
	{ 
		vector <int> tmp;
		tmp.clear();
		bitset <KIND> t;
		t.reset();
		for(int i = 0; i < qu.size(); i++)
		{
			if(t.test(qu[i]))
				continue;
			double min_p= 100.0;
			int id = 0;
			for(int j = 0; j < RSP[qu[i]].size(); j++)
			{
				double pp = p[RSP[qu[i]][j]].first * w[wp] + p[RSP[qu[i]][j]].second * (1.0-w[wp]);
				if(pp < min_p)
				{
					min_p = pp;
					id = RSP[qu[i]][j];
				}
			}
			tmp.push_back(id);
			t |= Nodes[id].keyword & querybinary;
			if(t.count() == querybinary.count())
				break;
		}
		ret.push_back(tmp);
	}
	// for(auto c : ret)
	// {
	// 	for(auto k : c)
	// 	{
	// 		cout << k << " ";
	// 	}
	// 	cout << endl;
	// }
	return ret;

}
vector<int> order;
int getTSP(int s,int t,vector<int> &node){
	int dist=0;
	int index;
	order.push_back(s);
	while(node.size()> 0){
	int min=INT_MAX;
		for(int i=0;i<node.size();i++){
			if(min > distanceQuery(s+1,node[i]+1)){
				min = distanceQuery(s+1,node[i]+1);		
				index = node[i];
			}
		}
		dist+=min;
		order.push_back(index);
		for(vector<int>::iterator iter=node.begin();iter!=node.end();iter++){        //从vector中删除指定的某一个元素 
    			if(*iter==index){
        			node.erase(iter);
        			break;
    			}
		}
		s = index;	
	}
	dist+=distanceQuery(index+1,t+1);
	order.push_back(t);
	return dist;
}

int main(){
	readIndex();
	// cout << "Load Graph Finished!" << endl;
	LOG2 = (int*)malloc(sizeof(int) * (n * 2 + 10));
	LOGD = (int*)malloc(sizeof(int) * (n * 2 + 10));
	int k = 0, j = 1;
	for (int i = 0; i < n * 2 + 10; i++){
		if (i > j * 2){
			j *= 2;
			k++;
		}
		LOG2[i] = k;
		LOGD[i] = j;
	}
	// init
	init();

	// load gtree index
	// pre query init
	// pre_query();
	// cout<<"READED GTREE-INDEX.. "<<endl;	

	
	

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	std::chrono::duration<double> time_span;
	std::chrono::duration<double> IG_base_time_totol;

	
	// querystr= "28,33"; //cal
	// querystr= "10,16,23,32"; //4
	// querystr= "8,18,24,17,28,26"; //6
	// querystr= "7,10,11,24,26,31,38,45"; //8
	// querystr= "6,18,24,33,46,2,34,25,47,19"; //10

	// querystr= "156,343";//NY
	// querystr= "147,245,354,498";
	querystr= "99,154,295,334,347,456";
	// querystr= "137,247,359,468,521,570,591,621";
	// querystr= "45,89,145,169,258,290,358,402,479,587";

	// querystr= "426,635,727,730,271,556"; //6%
	// querystr= "598,676,1047,1196,1410,1411"; //9%
	// querystr= "936,1265,1417,1589,1603,1769"; //12%
	// querystr= "453,809,1505,2057,2411,2585"; //15%

	// querystr= "11,58,57,96,137,141"; //Density:2%
	// querystr= "11,95,127,216,225,272";
	// querystr= "287,369,443,576,825,975";
	// querystr= "362,428,946,1123,1411,1738";

	// querystr= "447,450"; //1000COL
	// querystr= "178,379,511,479"; //4
	// querystr= "278,384,708,930,528,875"; //6
	// querystr= "115,551,247,330,317,368,443,741"; //8
	// querystr= "978,349,407,146,421,690,643,726,829,626"; //10


	fstream infile("testData/NY_tol/query/query_1000_0-20%.txt",ios::in);
	//fstream infile("testData/res_1.txt",ios::in);
	// ofstream outfile1("testData/Cal_tol/query_res/1000_10_0-15e_IG_base.txt",ios::out);
	ofstream outfile_IG("testData/NY_tol/analysis_res/Density8/3%_6_0-20%_WPOI_H2H.txt",ios::out);
	ofstream outfile_res("testData/NY_tol/result/Density8/result.txt",ios::app);

	// fstream infile("testData/Cal_tol/query/query_1000_10-20%.txt",ios::in);

	// ofstream outfile_IG("testData/Cal_tol/analysis_res/Density8/3%_10_10-20%_IG.txt",ios::out);
	// ofstream outfile_res("testData/Cal_tol/result/Density8/result.txt",ios::app);

	// fstream infile("testData/COL_tol/query/query_1000_10-20%.txt",ios::in);

	// ofstream outfile_IG("testData/COL_tol/analysis_res/Density8/3%_6_10-20%_IG.txt",ios::out);
	// ofstream outfile_res("testData/COL_tol/result/Density8/result.txt",ios::app);

	// fstream infile("testData/NY_sub/query/query_1000_0-10%.txt",ios::in);

	// ofstream outfile_IG("testData/NY_tol/analysis_res/Density8/500_4_0-10%_IG.txt",ios::out);
	// ofstream outfile_res("testData/NY_tol/result/Density8/result.txt",ios::app);  

	string Scout = "keywords:3%  querykey:6-{99,154,295,334,347,456}  Dis:10~20%";
	// outfile_res << "********************************************************" << endl;
	
	querybinary|= set_str_to_binary(querystr);
	cout<<"querybinary:"<<querybinary.to_string()<<endl;
	// outfile1<<"querybinary:"<<querybinary.to_string()<<endl;
	int start,end;
	int leng;
	int num = 0;
	srand (time(NULL));
	int candsize = 0;
	while(infile >> start >> end >> leng)
	{	
		/*if(num < 684)
		{num++;
			continue;
		}*/
		if(num == 1000)
			break;
		num++;
		// cout << num << " strat:" << start << "   end:" << end << endl;

		t1 = std::chrono::high_resolution_clock::now();	
		vector<int>tmp;
		
		// int s = H2HPath(start,end,tmp);
		cout << num << " strat:" << start << "   end:" << end << endl;
		// outfile1 << "strat:" << start << "    end:" << end << endl;
		// outfile_IG << start << "    " << end << endl;
		// outfile_test << "strat:" << start << "    end:" << end << endl;
		//outfile_test  << start << "    " << end << endl;
		
	
		vector<vector<int>> candicate;
		vector<vector<int>>().swap(candicate); //释放内存
		

        candicate = expand_wpoi(start, end);

		
		vector<PA>RES;
		RES.clear();
		while (!qu.empty())
		{
			qu.pop();
		}

		// outfile1 << "candicate size: "<<candicate.size() << endl;
		for(int i = 0; i < candicate.size(); i++)
		{
			order.clear();
			// res_path.clear();
			//outfile1 << "111: "<<candicate.size() << endl;

			int s = getTSP(start,end,candicate[i]);
			//outfile1 << "9911: "<<candicate.size() << endl;

			// RES.push_back(PA(res_path,s,order.size()-2));
			// qu.push(PA(res_path,s,order.size()-2));
			RES.push_back(PA(order,s,order.size()-2));
			qu.push(PA(order,s,order.size()-2));

		}
		
		while (!qu.empty())
		{
			qu.pop();
		}

		for(vector<PA>::iterator it = RES.begin();it!=RES.end();it++){
		if(R.size()==0){R.push_back(*it);}
		else{
			bool flag = false;
			for(int i=0;i<R.size();i++){
				if(R[i].park<=(*it).park&&R[i].dis<=(*it).dis){
					flag = true;
					break;
				}
			}
			if(flag){
				continue;
			}else{
				vector<PA>p;
				p.clear();
				for(int i=0;i<R.size();i++){					
					if(R[i].park>=(*it).park&&R[i].dis>=(*it).dis){
						continue;
					}else{
						p.push_back(R[i]);
					}
				}
				R.clear();
				R = p;	
				R.push_back(*it);
			}
		}
		}
		for(vector<PA>::iterator it = R.begin();it!=R.end();it++){
			qu.push(*it);
		}
		R.clear();

		t2 = std::chrono::high_resolution_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
		// outfile1 << "IG_base_time: " << time_span.count() << endl;
		cout << "IG_base_time: " << time_span.count() << endl;
		IG_base_time_totol = IG_base_time_totol + time_span;
		// outfile1 << "\n" << endl;
		cout << "total_time: " << IG_base_time_totol.count() << endl;

		outfile_IG << start << " " << end << endl;
		outfile_IG << qu.size() << endl;
		while (!qu.empty())
		{
			//cout << qu.top().park << " : " << qu.top().dis<<" = "<<qu.top().path<<endl;
			// outfile1 << qu.top().park << " : " << qu.top().dis<<" = "<<qu.top().path<<endl;
			outfile_IG << qu.top().dis  << " "<<qu.top().park<<endl;
			qu.pop();
		}
		outfile_IG << endl;
	}
	// outfile1 << "IG_Maxprune : totol time = "<<IG_base_time_totol.count()<<"\t"<<"average time = "<<IG_base_time_totol.count()/1000<<endl;
	outfile_res << "IG_Maxprune : totol time = "<<IG_base_time_totol.count()<<"\t"<<"average time = "<<IG_base_time_totol.count()/num<<endl;
	 cout << "IG_Maxprune : totol time = "<<IG_base_time_totol.count()<<"\t"<<"average time = "<<IG_base_time_totol.count()/num<<endl;
	 outfile_res << "IG-Tree candidate: " << candsize / num << endl; 
	cout << Scout << endl;
	outfile_res << Scout << endl;
	outfile_res << "---------------------------------------" << endl;
	outfile_res << endl;

	return 0;
}

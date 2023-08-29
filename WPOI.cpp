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

const int KIND = 634;
#define FILE_NODE "./map/data_cal/sub/cal_sub_keys_5_100.node"
#define FILE_EDGE "./map/data_NY/tol/NY_deal.edge"

typedef struct{
	double x,y;
	vector<int> adjnodes;
	vector<int> adjweight;
	bitset<KIND> keyword;
}Node;

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
            double pdist =  (double)(((dijkstra_p2p(s, k) + dijkstra_p2p(k, t)) - dijkstra_p2p(s, t)) *1.0 / dijkstra_p2p(s, t) *1.0);
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
			if(min > dijkstra_p2p(s,node[i])){
				min = dijkstra_p2p(s,node[i]);		
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
	dist+=dijkstra_p2p(index,t);
	order.push_back(t);
	return dist;
}

int main(){
	// readIndex();
	// cout << "Load Graph Finished!" << endl;
	// LOG2 = (int*)malloc(sizeof(int) * (n * 2 + 10));
	// LOGD = (int*)malloc(sizeof(int) * (n * 2 + 10));
	// int k = 0, j = 1;
	// for (int i = 0; i < n * 2 + 10; i++){
	// 	if (i > j * 2){
	// 		j *= 2;
	// 		k++;
	// 	}
	// 	LOG2[i] = k;
	// 	LOGD[i] = j;
	// }
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

	querystr= "280,343";//NY
	// querystr= "263,285,330,489";
	// querystr= "391,501,516,222,269,393";
	// querystr= "97,234,376,326,437,385,525,603";
	// querystr= "85,182,234,142,170,189,411,432,450,525";

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


	fstream infile("testData/NY_tol/query/query_1000_10-20%.txt",ios::in);
	//fstream infile("testData/res_1.txt",ios::in);
	// ofstream outfile1("testData/Cal_tol/query_res/1000_10_0-15e_IG_base.txt",ios::out);
	ofstream outfile_IG("testData/NY_tol/analysis_res/Density8/3%_3_10-20%_IG.txt",ios::out);
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

	string Scout = "keywords:3%  querykey:6-{453,809,1505,2057,2411,2585}  Dis:10~20%";
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
		if(num == 20)
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

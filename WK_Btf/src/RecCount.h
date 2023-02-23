#ifndef RECCOUNT_H_
#define RECCOUNT_H_

#include<omp.h>
#include<string>
#include<vector>
#include<cstring>
#include<algorithm>
#include<cstdio>
#include<math.h>
#include<map>
#include<set>
#include<random>
#include<iostream>
#include<cstdlib>
#include<time.h>

#define MAXN 2000000000
#define MAXN1 1000000000

#define METHOD_NAIVE -2
#define METHOD_NAIVE_ORG -1
#define METHOD_SQR_ORG 0
#define METHOD_ARB 1
#define METHOD_SQR 2

#define METHOD_DD 0
#define METHOD_II 1
#define METHOD_DI 2

#define INC_FACTOR 1.2

using namespace std;

class IntNode {
public:
	int val, id;
	bool operator < (const IntNode& v) const;
};

class Graph {
public:
	long long m;
	long long btf_G;
	double btf_time;
	int n, n1;
	int *dat, **con, *deg;
	int *oid, *nid;
	int *cnt_dat, **cnt;
	bool is_original;

public:
	static void txt_to_bin(string path, bool gorder = 0, bool is_bipartite = false);
	static bool get_edge(char *line, int &a, int &b, int num_cnt);
	static int get_num_cnt(string path);
	static void get_order(vector<int> *con, int n, int *o);

public:
	Graph(string path, int method = METHOD_ARB, bool is_original=false);
	~Graph();
	void rec_count(int method, bool cnt_edge);
	void rec_count_arb(bool cnt_edge);
	void rec_count_sqr(bool cnt_edge);
	void rec_count_naive(bool cnt_edge);
	void init_count();
};

bool IntNode::operator < (const IntNode& v) const {
	if(val == v.val) return id < v.id; else return val < v.val;
}

int Graph::get_num_cnt(string path) {
	FILE *fin = fopen( (path + "graph.txt").c_str(), "r" );
	char line[1000], st[64];
	int cnt = 0, min_cnt = 100;

	while( fgets( line, 1000, fin ) && cnt < 10 ) {
		if( !isdigit(line[0]) ) continue;
		vector<char*> v_num;
		int len = (int) strlen(line);
		for( int i = 0; i < len; ++i )
			if( !isdigit(line[i]) && line[i] != '.' ) line[i] = '\0';
			else if(i == 0 || !line[i-1]) v_num.push_back(line+i);
		if( (int) v_num.size() < 2 ) continue;
		min_cnt = min(min_cnt, (int)v_num.size());
		++cnt;
	}
	fclose( fin );
	return min_cnt;
}

bool Graph::get_edge(char *line, int &a, int &b, int num_cnt) {

  	if( !isdigit(line[0]) ) return false; //for other datasets
	vector<char*> v_num;
	int len = (int) strlen(line);
	for( int i = 0; i < len; ++i )
		if( !isdigit(line[i]) && line[i] != '.' ) line[i] = '\0';
		else if(i == 0 || !line[i-1]) v_num.push_back(line+i);
	if( (int) v_num.size() != num_cnt ) return false;
	sscanf( v_num[0], "%d", &a );
	sscanf( v_num[1], "%d", &b );
	return true;
}

void Graph::get_order( vector<int> *con, int n, int *o) {
	IntNode *f = new IntNode[n];
	for( int i = 0; i < n; ++i )
		f[i].id = i, f[i].val = (int) con[i].size();
	sort(f, f + n);
	for(int i = 0; i < n; ++i) o[i] = f[n-1-i].id;
	delete[] f;
}

void Graph::txt_to_bin(string path, bool gorder, bool is_bipartite) {
	FILE *fin;
	if (! gorder)
		fin = fopen( (path + "graph.txt").c_str(), "r" );
	else {
		fin = fopen( (path + "graph_Gorder.txt").c_str(), "r" );
	}
	char line[1024];
	int n = 0, a, b, n1 = 0, n2 = 0, num_cnt = get_num_cnt(path);
	vector<int> *con = new vector<int>[MAXN];
	long long cnt = 0, m = 0;

	printf( "Loading text, num_cnt=%d, is_bipartite=%s...\n", num_cnt, is_bipartite?"true":"false" );
	while( fgets( line, 1024, fin ) ) {
		if( !get_edge(line, a, b, num_cnt) ) continue;
		if( a < 0 || b < 0 ) continue;
		n1 = max(n1, a+1);
		n2 = max(n2, b+1);
		if( is_bipartite ) b += MAXN1; else if( a == b ) continue;
		if( con[a].capacity() == con[a].size() ) con[a].reserve(con[a].size() * INC_FACTOR );
		if( con[b].capacity() == con[b].size() ) con[b].reserve(con[b].size() * INC_FACTOR );
		con[a].push_back(b);
		con[b].push_back(a);
		if( (++cnt) % (long long) 10000000 == 0 ) printf( "%lldM lines finished\n", cnt/1000000 );
	}
	fclose( fin );
	n = is_bipartite ? n1 + n2 : max(n1, n2);

	printf( "n1=%d,n2=%d\n", n1, n2 );

	if( is_bipartite ) {
		printf( "Re-numbering bipartite...\n" );
		for( int i = 0; i < n1; ++i )
			for( int j = 0; j < (int) con[i].size(); ++j )
				con[i][j] = con[i][j] - MAXN1 + n1;
		for( int i = 0; i < n2; ++i ) {
			con[n1+i].clear();
			con[n1+i].reserve(con[MAXN1+i].size());
			for( int j = 0; j < (int) con[MAXN1+i].size(); ++j )
				con[n1+i].push_back(con[MAXN1+i][j]);
			vector<int>().swap(con[MAXN1+i]);
		}
	}

	int maxd = 0;
	for( int i = 0; i < n; ++i )
		if( con[i].size() > 0 ){
			sort( con[i].begin(), con[i].end() );
			int p = 1;
			for( int j = 1; j < (int) con[i].size(); ++j )
				if( con[i][j-1] != con[i][j] ) con[i][p++] = con[i][j];
			con[i].resize( p ); m += p;
			maxd = max(maxd, p);
		}

	printf( "Reordering...\n" );


	int *oid = new int[n], *nid = new int[n];
	get_order(con, n, oid);
	for( int i = 0; i < n; ++i ) nid[oid[i]] = i;


	printf( "Saving binary...\n" );
	int *deg = new int[n];
	for( int i = 0; i < n; ++i) deg[i] = (int) con[oid[i]].size();

	FILE *fout = fopen( (path + "graph-sort.bin").c_str(), "wb" );
	fwrite( &n1, sizeof(int), 1, fout );
	fwrite( &n, sizeof(int), 1, fout );
	fwrite( &m, sizeof(long long), 1, fout );
	fwrite( deg, sizeof(int), n, fout );

	int *nbr = new int[maxd];
	for( int i = 0; i < n; ++i ) {
		int u = oid[i], d = (int) con[u].size();
		for( int j = 0; j < d; ++j )
			nbr[j] = nid[con[u][j]];
		sort(nbr, nbr + d);
		fwrite(nbr, sizeof(int), d, fout);
	}

	fwrite( nid, sizeof(int), n, fout );
	fwrite( oid, sizeof(int), n, fout );
	fclose( fout );
	printf( "Created binary file, n = %d, m = %lld\n", n, m );
	delete[] con; delete[] deg; delete[] oid; delete[] nid; delete[] nbr;
}

Graph::Graph(string path, int method, bool is_original) {
	{
		printf( "Loading graph...\n" );
		this->is_original = is_original;
		FILE* fin = fopen( (path+"graph-sort.bin").c_str(), "rb" );
		fread( &n1, sizeof(int), 1, fin );
		fread( &n, sizeof(int), 1, fin );
		fread( &m, sizeof(long long), 1, fin );
		deg = new int[n]; dat = new int[m]; con = new int*[n];
		nid = new int[n]; oid = new int[n];
		fread( deg, sizeof(int), n, fin );
		fread( dat, sizeof(int), m, fin );
		fread( nid, sizeof(int), n, fin );
		fread( oid, sizeof(int), n, fin );
		fclose(fin);

		cnt = NULL;
		cnt_dat = NULL;
		long long p = 0;
		for( int i = 0; i < n; ++i ) {con[i] = dat+p; p+=deg[i];}
		if( is_original ) {
			printf( "obtaining original graph...\n" );
			for( int i = 0; i < n; ++i ) {
				for( int j = 0; j < deg[i]; ++j ) con[i][j] = oid[con[i][j]];
					sort(con[i], con[i] + deg[i]);//TODO order
			}
		}
		printf( "%s: n=%d,n1=%d,m=%lld\n", path.c_str(), n, n1, m );
	}
}

Graph::~Graph() {
	delete[] deg; delete[] dat; delete[] con;
	delete[] nid; delete[] oid;
	if(cnt) delete[] cnt; if(cnt_dat) delete[] cnt_dat;
}

void Graph::init_count() {
	cnt = new int*[n]; cnt_dat = new int[m];
	memset(cnt_dat, 0, sizeof(int) * m);
	long long p = 0;
	for( int i = 0; i < n; ++i ) {cnt[i] = cnt_dat+p; p+=deg[i];}
}

void Graph::rec_count_naive(bool cnt_edge) {
	if(cnt_edge) init_count();
	double t = omp_get_wtime(), dlt = 0;
	long long cnt_rec = 0, cnt_tmp = 0, cnt_bistar = 0;
	printf( "counting rectangles, method=naive, is_original=%s, cnt_edge=%s...\n", is_original? "true":"false", cnt_edge ? "true" : "false" );

	long long sum_a = 0, sum_b = 0;
	for(int u = 0; u < n; ++u ) {
		if(oid[u] < n1) sum_a += deg[u] * (long long) deg[u];
		else {sum_b += deg[u] * (long long) deg[u];}
	}
	bool from_a = sum_a > sum_b;

	printf( "from_a=%s,sum_a=%lld,sum_b=%lld\n", from_a?"true":"false", sum_a, sum_b);
	#pragma omp parallel
	{
		int pid = omp_get_thread_num(), np = omp_get_num_threads(), u, v, w;
        	int idx = 0;
		long long local_cnt = 0, local_cnt_tmp = 0, local_cnt_bistar = 0;
		double local_dlt = 0;
		int *last_use = new int[n], *last_cnt = new int[n], pre = 0;
		memset(last_use, -1, sizeof(int) * n);

		for( int p = pid; p < n; p += np ) {
			int u = is_original ? nid[p] : p;
	 	    	for (int j = 0; j < idx; j++) {
				last_cnt[last_use[j]] = 0;
		    	}
		    	idx = 0;
			if(!is_original && ((from_a && oid[u] >= n1) || (!from_a && oid[u] < n1))) continue;
			if(is_original && ((from_a && p >= n1) || (!from_a && p < n1))) continue;
			for( int i = 0; i < deg[u]; ++i ) {
				v = is_original ? nid[con[u][i]]:con[u][i]; pre = -1;
				for( int j = 0; j < deg[v]; ++j ) {
					w = con[v][j]; ++local_cnt_tmp;  if( w >= p ) continue; 
					if(pre == -1) pre = w; else {local_dlt += abs(w-pre); pre = w;}
					//if( last_use[w] != u ) {last_use[w] = u; last_cnt[w] = 1;}
					//else {local_cnt += last_cnt[w]; ++last_cnt[w]; ++local_cnt_bistar;}
					local_cnt += last_cnt[w]; ++last_cnt[w]; ++local_cnt_bistar;
		    			if (last_cnt[w] == 1) 
						last_use[idx++] = w;
				}
			}
			if( !cnt_edge ) continue;
			for( int i = 0; i < deg[u]; ++i ) {
				v = is_original ? nid[con[u][i]]:con[u][i];
				for( int j = 0; j < deg[v]; ++j ) {
					w = con[v][j]; if( w >= p ) continue;
					int dt = last_cnt[w] - 1;
					if( dt ) {cnt[u][i] += dt; 
					#pragma omp atomic						
						cnt[v][j] += dt;}
				}
			}
		}
		if(pid==0) printf( "num_thread=%d,t_local=%0.3lfsec,", np, omp_get_wtime()-t );
		#pragma omp critical
		{
			cnt_rec += local_cnt; cnt_tmp += local_cnt_tmp; cnt_bistar += local_cnt_bistar; dlt += local_dlt;
		}
	}
	printf( "total_time=%0.3lfsec,n_rectangles=%lld,cnt_tmp=%lld,cnt_bistar=%lld,dlt=%0.0lf,ajd=%0.3lf\n", omp_get_wtime()-t, cnt_rec, cnt_tmp, cnt_bistar, dlt, dlt/cnt_tmp );
	btf_time = omp_get_wtime()-t;
	btf_G = cnt_rec;
}

void Graph::rec_count_sqr(bool cnt_edge) {
	if(cnt_edge) init_count();
	double t = omp_get_wtime(), dlt = 0;
	long long cnt_rec = 0, cnt_tmp = 0, cnt_bistar = 0;
	printf( "counting rectangles, method=SQR, is_original=%s, cnt_edge=%s...\n", is_original? "true":"false", cnt_edge ? "true" : "false" );
	
	#pragma omp parallel
	{
		int pid = omp_get_thread_num(), np = omp_get_num_threads(), u, v, w;

		long long local_cnt = 0, local_cnt_tmp = 0, local_cnt_bistar = 0;
		double local_dlt = 0;
		int *last_use = new int[n], *last_cnt = new int[n], pre = 0;
		memset(last_use, -1, sizeof(int) * n);
		int idx = 0;
		for( int p = pid; p < n; p += np ) {
			int u = is_original ? nid[p] : p;
			for (int j = 0; j < idx; j++) {
				last_cnt[last_use[j]] = 0;
			}
			idx = 0;
			//start[deg[u]]++;////count degree///////
			for( int i = 0; i < deg[u]; ++i ) {
				v = is_original ? nid[con[u][i]]:con[u][i]; pre = -1;
				//central[deg[v]]++;////count degree///////
				if( v <= u ) continue;
				for( int j = 0; j < deg[v]; ++j ) {
					w = con[v][j]; 
					if( (is_original && nid[w] <= u) || (!is_original && w <= u) ) continue;
					++local_cnt_tmp;
					//end[deg[nid[w]]]++;////count degree///////
					 if(pre == -1) pre = w; else {local_dlt += abs(w-pre); pre = w;}
					//if( last_use[w] != u ) {last_use[w] = u; last_cnt[w] = 1;}
					//else {local_cnt += last_cnt[w]; ++last_cnt[w]; ++local_cnt_bistar;}
					local_cnt += last_cnt[w]; ++last_cnt[w]; ++local_cnt_bistar;
					if (last_cnt[w] == 1) 
						last_use[idx++] = w;
				}
			}
			if( !cnt_edge ) continue;
			for( int i = 0; i < deg[u]; ++i ) {
				v = is_original ? nid[con[u][i]]:con[u][i]; pre = -1;
				if( v <= u ) continue;
				for( int j = 0; j < deg[v]; ++j ) {
					w = con[v][j];
					if( (is_original && nid[w] <= u) || (!is_original && w <= u) ) continue;
					int dt = last_cnt[w] - 1;
					if( dt ) {cnt[u][i] += dt; 
				#pragma omp atomic 
					cnt[v][j] += dt;}
				}
			}
		}
		if(pid==0) printf( "num_thread=%d,t_local=%0.3lfsec,", np, omp_get_wtime()-t );
		#pragma omp critical
		{
			cnt_rec += local_cnt; cnt_tmp += local_cnt_tmp; cnt_bistar += local_cnt_bistar; dlt += local_dlt;
		}
	}
	printf( "total_time=%0.3lfsec,n_rectangles=%lld,cnt_tmp=%lld,cnt_bistar=%lld,dlt=%0.0lf,ajd=%0.3lf\n", omp_get_wtime()-t, cnt_rec, cnt_tmp, cnt_bistar, dlt, dlt/cnt_tmp );	
	btf_time = omp_get_wtime()-t;
	btf_G = cnt_rec;
}

void Graph::rec_count_arb(bool cnt_edge) {
	if(cnt_edge) init_count();
	double t = omp_get_wtime(), dlt = 0;
	long long cnt_rec = 0, cnt_tmp = 0, cnt_bistar = 0;

	printf( "counting rectangles, method=%s, cnt_edge=%s...\n", "ARB", cnt_edge ? "true" : "false" );
	#pragma omp parallel
	{
		int pid = omp_get_thread_num(), np = omp_get_num_threads(), u, v, w;

		long long local_cnt = 0, local_cnt_tmp = 0, local_cnt_bistar = 0;
		double local_dlt = 0;
		int *last_use = new int[n], *last_cnt = new int[n], pre = 0;
		memset(last_use, -1, sizeof(int) * n);
		int idx = 0;
		for( u = pid; u < n; u += np ) {
			//start[deg[u]]++;////count degree///////
			for (int j = 0; j < idx; j++) {
				last_cnt[last_use[j]] = 0;
			}
			idx = 0;
			for( int i = 0; i < deg[u]; ++i ) {
				v = con[u][i]; pre = -1;
				//central[deg[v]]++;////count degree///////
				for( int j = 0; j < deg[v]; ++j ) {
					w = con[v][j]; ++local_cnt_tmp; if(w >= u || w >= v) break;  ++local_cnt_bistar;
					//end[deg[w]]++;////count degree///////
					 if(pre == -1) pre = w; else {local_dlt += abs(w-pre); pre = w;}
					//if( last_use[w] != u ) {last_use[w] = u; last_cnt[w] = 1;}
					//else {local_cnt += last_cnt[w]; ++last_cnt[w]; ++local_cnt_bistar;}
					local_cnt += last_cnt[w]; ++last_cnt[w]; 
		    			if (last_cnt[w] == 1) 
						last_use[idx++] = w;
				}
			}
			if( !cnt_edge ) continue;
			for( int i = 0; i < deg[u]; ++i ) {
				v = con[u][i];
				for( int j = 0; j < deg[v]; ++j ) {
					w = con[v][j]; if(w >= u || w >= v) break;
					int dlt = last_cnt[w] * (last_cnt[w] - 1) / 2;
					if( dlt ) {cnt[u][i] += dlt; 
					#pragma omp atomic
					cnt[v][j] += dlt;}
				}
			}
		}
		if(pid==0) printf( "num_thread=%d,t_local=%0.3lfsec,", np, omp_get_wtime()-t );
		#pragma omp critical
		{
			cnt_rec += local_cnt; cnt_tmp += local_cnt_tmp; cnt_bistar += local_cnt_bistar; dlt += local_dlt;
		}
	}
	printf( "total_time=%0.3lfsec,n_rectangles=%lld,cnt_tmp=%lld,cnt_bistar=%lld,dlt=%0.0lf, ajd=%0.3lf\n", omp_get_wtime()-t, cnt_rec, cnt_tmp, cnt_bistar, dlt, dlt/cnt_tmp );
	btf_time = omp_get_wtime()-t;
	btf_G = cnt_rec;
}


#endif /*RECCOUNT_H_*/
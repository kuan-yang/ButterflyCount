#ifndef BTFCOUNT_H_
#define BTFCOUNT_H_

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
#include<fstream>
#include<unordered_map>
#include<limits.h>

#define MAXN 10000000
#pragma warning(disable : 4996)


using namespace std;
typedef long long ll;
// because RANDMAX only has 15 bits
inline ll rand64() {
	return  (((ll(rand())) << 60) + ((ll(rand())) << 45) + \
		((ll(rand())) << 30) + ((ll(rand())) << 15) + ll(rand()))
		& (LLONG_MAX >> 1);
}
inline int rand32() {
	return  (((ll(rand())) << 30) + ((ll(rand())) << 15) + ll(rand()))
		& (INT_MAX >> 1);
}
inline bool smaller(int a, int a_deg, int b, int b_deg) {
	return a_deg < b_deg || (a_deg == b_deg && a < b);
}


class Graph {
public:
	ll m;		// number of edge
	ll m_hat;	// estimate of m
	ll BTF_cnt;	//  number of tirangles
	ll query_cnt;	//count for queries
	double BTF_time;	// runtime of alg
	int n;	//number of vertice
	int* dat, ** con, * deg;	//data, edge, degree
	bool* dat_mark, * deg_mark;
	//int* oid, * nid;	// old_id, new_id
	//int* cnt_dat, ** cnt;
	bool is_original;


public:
	inline int query_degree(int v);
	inline bool query_pair(int v, int u);
	inline int query_neighbor(int v, int i);

	ll estimate_with_advice(ll, double epsilon);
	void btf_count_est(double epsilon);

	void btf_count_exact();
	void btf_count_hash();
	void btf_count_ESpar(double,bool);
	void rec_count_arb();

public:
	Graph(string path);
	~Graph();
};

Graph::Graph(string path) {

	printf("Loading graph...\n");
	this->is_original = is_original;
	FILE* fin = fopen((path + "graph.bin").c_str(), "rb");
	if (!fin) exit(0);
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(ll), 1, fin);
	deg = new int[n]; dat = new int[m]; con = new int* [n];

	fread(deg, sizeof(int), this->n, fin);
	fread(dat, sizeof(int), this->m, fin);

	fclose(fin);
	m /= 2; // undirect

	ll p = 0;
	for (int i = 0; i < n; ++i) { con[i] = dat + p; p += deg[i]; }
	printf("%s: n=%d,m=%lld\n", path.c_str(), n, m);

	query_cnt = -1;
	BTF_cnt = -1;
	BTF_time = -1;
}


Graph::~Graph() {
	delete[] deg; delete[] dat; delete[] con;
	//if (cnt) delete[] cnt; if (cnt_dat) delete[] cnt_dat;
}


inline int Graph::query_degree(int v) {
#pragma omp atomic
	++query_cnt;
	return deg[v];
}

inline bool Graph::query_pair(int v, int u) {
#pragma omp atomic
	++query_cnt;
	int* tmp = lower_bound(con[v], con[v] + deg[v], u);
	return tmp != con[v] + deg[v] && *tmp == u;
}

inline int Graph::query_neighbor(int v, int i) {
#pragma omp atomic
	++query_cnt;
	//assert(i < deg[v]);
	return con[v][i];
}


ll Graph::estimate_with_advice(ll b_hat, double epsilon) {
	ll s = ll(double(log(n)) * m * n * n / b_hat  / pow(epsilon, 3)); //huge number
	s = (s>0 && s<ll(n)*n)? s: ll(n) * n;
	if (s == 0) return 0;

	ll y_sum = 0;
	for (int i = 0; i < s; ++i) {
		int u = rand32() % n, v = rand32() % n;
		if (u == v) continue;
		int u_deg = query_degree(u), v_deg = query_degree(v);
		if (u_deg == 0 || v_deg == 0) continue;
		if (!smaller(u, u_deg, v, v_deg)) {
			swap(u, v); swap(u_deg, v_deg);
		}
		int r;
		if ((ll)u_deg * u_deg <= m_hat)
			r = double(rand()) / RAND_MAX <= double(u_deg ) * u_deg / m ? 1 : 0;
		else
			r = int((ll)u_deg * u_deg / m) + 1;
		ll z_sum = 0;
		for (int j = 0; j < r; ++j) {
			int w1 = query_neighbor(u, rand32() % u_deg), w2 = query_neighbor(u, rand32() % u_deg);
			if (w1 == w2) continue;
			if (!query_pair(w1, v) || !query_pair(w2, v)) continue;
			int w1_deg = query_degree(w1), w2_deg = query_degree(w2);
			if (smaller(w1, w1_deg, u, u_deg) || smaller(w2, w2_deg, u, u_deg)) continue;
			z_sum += max(ll(u_deg) * u_deg, m_hat);
		}
		if (r)
			y_sum += z_sum / (2 * r);
	}
	ll x = ll(n) * n / (2 * s) * y_sum;
	return x;
}


void Graph::btf_count_est(double epsilon) {
	//TODO Feige Alg
	//int d_hat = 0;
	query_cnt = 0;
	m_hat = m;

	double t = omp_get_wtime();
	int itr_cnt = 0;
	ll b_max = pow(LLONG_MAX, .25) < n ? LLONG_MAX : ll(n) * n * n * n;
	ll b_hat = b_max;	//	estimate of t
	while (b_hat >= 1 ) {//&& query_cnt <= m_hat
		cout << b_hat << '\n';
		for (long long iter_b = b_max; iter_b >= b_hat; iter_b /= 2) {
			long long x = b_max;
			for (long long i = 1; i <= log(log(n)) / epsilon * 1; ++i) {  // let c = 1 here
				long long tmp = estimate_with_advice(iter_b,epsilon);
				x = min(x, tmp);
				if (x < iter_b) break; // optimize
			}
			if (x >= iter_b) {
				BTF_cnt = x;
				BTF_time = omp_get_wtime() - t;
				printf("total_time=%0.3lfsec, number of Btfs=%lld, number of queries=%lld\n"\
					, BTF_time, BTF_cnt, query_cnt);
				return;
			}
		}
		b_hat /= 2;
	}
	printf("triangles number may be less than sqrt(m)");
}


void Graph::btf_count_exact() {
	double t = omp_get_wtime();
	BTF_cnt = 0;
#pragma omp parallel
	{
		int pid = omp_get_thread_num(), np = omp_get_num_threads();
		for (int i = pid; i < n; i += np) {
			for (int j = i + 1; j < n; ++j) {
				int wedge_cnt = 0;
				for (int k = 0; k < deg[i]; ++k) {
					if (con[i][k] <= i) continue;
					int w = con[i][k];
					if (query_pair(j, w))
						++wedge_cnt;
				}
#pragma omp atomic
				BTF_cnt += wedge_cnt * (wedge_cnt - 1) / 2;
			}
			if (i % 1000 == 0)
				printf("%dk lines finished\n", i / 1000);
		}
   	}

	BTF_time = omp_get_wtime() - t;
	printf("total_time=%0.3lfsec, number of butterflies=%lld\n"\
		, BTF_time, BTF_cnt);
}


void Graph::btf_count_hash() {
	double t = omp_get_wtime();
	BTF_cnt = 0;
	unordered_map<ll, int> wedge_map;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < deg[i]; ++j) {
			for (int k = j + 1; k < deg[i]; ++k) {
				ll u = ll(query_neighbor(i, j)), v = ll(query_neighbor(i, k));
				if (u > v) swap(u, v);
				ll tmp_p = u * n + v;
				auto pnt = wedge_map.find(tmp_p);
				if (pnt != wedge_map.end())
					++pnt->second;
				else
					wedge_map[tmp_p] = 1;
				size_t a = wedge_map.size();
			}
		}
		if (i % 1000 == 0)
			printf("%dk lines finished\n", i / 1000);
		if (omp_get_wtime() - t > 300){
			cout<<"~~~~~~~~~~!! Timeout !!~~~~~~~~~~~\n";
			return;
		}
	}

	for (auto i = wedge_map.begin(); i != wedge_map.end(); ++i)
		BTF_cnt += i->second * (i->second - 1) / 2;

	BTF_cnt /= 2;

	BTF_time = omp_get_wtime() - t;
	printf("total_time=%0.3lfsec, number of butterflies=%lld\n"\
		, BTF_time, BTF_cnt);
}


void Graph::btf_count_ESpar(double p, bool fast_cnt) {
	double t = omp_get_wtime();

	m = 0;

	vector<int>* new_con = new vector<int>[MAXN];
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < deg[i]; ++j) {
			if (con[i][j] < i)
				continue;
			if (p * RAND_MAX > rand()) {
				new_con[i].emplace_back(con[i][j]);
				new_con[con[i][j]].emplace_back(i);
				m += 2;
			}
		}
	}
	delete dat;
	dat = new int[m];
	int pnt = 0;
	for (int i = 0; i < n; ++i) {
		con[i] = dat + pnt;
		deg[i] = new_con[i].size();
		for (int j = 0; j < deg[i]; ++j)
			dat[pnt++] = new_con[i][j];
	}
	m /= 2;
	cout << "Begin exact count under sparse\n";
	if (fast_cnt)
		this->btf_count_hash();
	else
		this->btf_count_exact();
	if(BTF_cnt==0) return;
	BTF_cnt = BTF_cnt / pow(p, 4);
	BTF_time = omp_get_wtime() - t;
	printf("total_time=%0.3lfsec, number of butterflies=%lld\n"\
		, BTF_time, BTF_cnt);

}

void Graph::rec_count_arb() {
	double t = omp_get_wtime(), dlt = 0;
	long long cnt_rec = 0, cnt_tmp = 0, cnt_bistar = 0;
#pragma omp parallel
	{
		int pid = omp_get_thread_num(), np = omp_get_num_threads(), u, v, w;

		long long local_cnt = 0, local_cnt_tmp = 0, local_cnt_bistar = 0;
		double local_dlt = 0;
		int* last_use = new int[n], * last_cnt = new int[n], pre = 0;
		memset(last_use, -1, sizeof(int) * n);
		int idx = 0;
		for (u = pid; u < n; u += np) {
			//start[deg[u]]++;////count degree///////
			for (int j = 0; j < idx; j++) {
				last_cnt[last_use[j]] = 0;
			}
			idx = 0;
			for (int i = 0; i < deg[u]; ++i) {
				v = con[u][i]; pre = -1;
				//central[deg[v]]++;////count degree///////
				for (int j = 0; j < deg[v]; ++j) {
					w = con[v][j]; ++local_cnt_tmp; if (w >= u || w >= v) break;  ++local_cnt_bistar;
					//end[deg[w]]++;////count degree///////
					if (pre == -1) pre = w; else { local_dlt += abs(w - pre); pre = w; }
					//if( last_use[w] != u ) {last_use[w] = u; last_cnt[w] = 1;}
					//else {local_cnt += last_cnt[w]; ++last_cnt[w]; ++local_cnt_bistar;}
					local_cnt += last_cnt[w]; ++last_cnt[w];
					if (last_cnt[w] == 1)
						last_use[idx++] = w;
				}
			}
		}
		if (pid == 0) printf("num_thread=%d,t_local=%0.3lfsec,", np, omp_get_wtime() - t);
#pragma omp critical
		{
			cnt_rec += local_cnt; cnt_tmp += local_cnt_tmp; cnt_bistar += local_cnt_bistar; dlt += local_dlt;
		}
	}
	BTF_time = omp_get_wtime() - t;
	BTF_cnt = cnt_rec;
	printf("total_time=%0.3lfsec, number of butterflies=%lld\n"\
		, BTF_time, BTF_cnt);
}

#endif /*BTFCOUNT_H_*/

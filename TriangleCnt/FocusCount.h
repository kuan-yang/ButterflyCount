#ifndef FOCUSCOUNT_H_
#define FOCUSCOUNT_H_

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

#define MAXN 10000000
#pragma warning(disable : 4996)

using namespace std;
// because RANDMAX only has 15 bits
inline long long rand64() {
	return  (((long long(rand())) << 60) + ((long long(rand())) << 45) + \
		((long long(rand())) << 30) + ((long long(rand())) << 15) + long long(rand()))
		& (LLONG_MAX >> 1);
}
inline int rand32() {
	return  (((long long(rand())) << 30) + ((long long(rand())) << 15) + long long(rand()))
		& (INT_MAX >> 1);
}
inline bool smaller(int a, int a_deg, int b, int b_deg) {
	return a_deg < b_deg || (a_deg == b_deg && a < b);
}

//get mid num
int PartSort(long long* arr, int start, int end)
{
	int left = start;
	int right = end;
	long long key = arr[end]; 
	while (left < right)
	{
		while (left < right && arr[left] <= key)  
		{
			++left;
		}
		while (left < right && arr[right] >= key)  
		{
			--right;
		}
		if (left < right)
		{
			swap(arr[left], arr[right]);  
		}
	}
	swap(arr[right], arr[end]);
	return left;
}

long long GetMidNumNoSort(long long* arr, int size)
{
	int start = 0;
	int end = size - 1;
	int mid = (size - 1) / 2;
	int div = PartSort(arr, start, end);
	while (div != mid)
	{
		if (mid < div)
			div = PartSort(arr, start, div - 1);
		else
			div = PartSort(arr, div + 1, end);
	}
	return arr[mid];
}

class Graph {
public:
	long long m;		// number of edge
	long long m_hat;	// estimate of m
	long long tri_cnt;	//  number of tirangles
	long long query_cnt;	//count for queries
	double tri_time;	// runtime of alg
	double epsilon;		// presion of alg
	double epsilon_hat;	//	epsilon` = epsilon / 3c_H
	double sqrt_m;
	int n;	//number of vertice
	int* dat, ** con, * deg;	//data, edge, degree
	bool* dat_mark, * deg_mark;
	//int* oid, * nid;	// old_id, new_id
	int* cnt_dat, ** cnt;
	bool is_original;
	

public:
	inline int query_degree(int v);
	inline bool query_pair(int v, int u);
	inline int query_neighbor(int v, int i);

	bool heavy(int v, int v_deg, long long t_hat);	// return (v_is_heavy) ? 1 : 0
	long long estimate_with_advice(long long);
	void estimate();

	void exact_cnt();

public:
	Graph(string path, int method, double ep);
	~Graph();
};

Graph::Graph(string path, int method, double ep) {

	printf("Loading graph...\n");
	this->is_original = is_original;
	FILE* fin = fopen((path + "graph.bin").c_str(), "rb");
	if (!fin) exit(0);
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(long long), 1, fin);
	deg = new int[n]; dat = new int[m]; con = new int* [n]; 
	deg_mark = new bool[n]; dat_mark = new bool[m];
	//nid = new int[n]; oid = new int[n];
	fread(deg, sizeof(int), n, fin);
	fread(dat, sizeof(int), m, fin);
	/*fread(nid, sizeof(int), n, fin);
	fread(oid, sizeof(int), n, fin);*/
	fclose(fin);

	cnt = NULL;
	cnt_dat = NULL;
	long long p = 0;
	for (int i = 0; i < n; ++i) { con[i] = dat + p; p += deg[i]; }
	printf("%s: n=%d,m=%lld\n", path.c_str(), n, m);

	memset(deg_mark, false, n*sizeof(bool));
	memset(dat_mark, false, m* sizeof(bool));
	epsilon = ep;
	const double const_H = 2000;
	epsilon_hat = epsilon / (3 * const_H);
	query_cnt = -1;
	tri_cnt = -1;
	tri_time = -1;
}


Graph::~Graph() {
	delete[] deg; delete[] dat; delete[] con;
	//delete[] nid; delete[] oid;
	if (cnt) delete[] cnt; if (cnt_dat) delete[] cnt_dat;
}


//inline int Graph::query_degree(int v) {
//	if (!deg_mark[v]) {
//#pragma omp atomic
//		++query_cnt;
//		deg_mark[v] = true;
//	}
//	return deg[v];
//}
//
//inline bool Graph::query_pair(int v, int u) {
//	int* tmp = lower_bound(con[v], con[v] + deg[v], u);
//	if (tmp != con[v] + deg[v] && *tmp == u) {
//		if (!dat_mark[tmp - dat]) {
//#pragma omp atomic
//			++query_cnt;
//			dat_mark[tmp - dat] = true;
//		}
//		return true;
//	}
//	return false;
//	//return tmp!=con[v]+deg[v] && *tmp == u;
//}
//inline int Graph::query_neighbor(int v, int i) {
//	if (!dat_mark[(con[v] - dat) + i]) {
//#pragma omp atomic
//		++query_cnt;
//		dat_mark[(con[v] - dat) + i] = true;
//	}
//	//assert(i < deg[v]);
//	return con[v][i];
//}
inline int Graph::query_degree(int v) {
#pragma omp atomic
	++query_cnt;
	return deg[v];
}

inline bool Graph::query_pair(int v, int u) {
#pragma omp atomic
			++query_cnt;
	int* pnt = lower_bound(con[v], con[v] + deg[v], u);
	return pnt != con[v] + deg[v] && *pnt == u;
}
inline int Graph::query_neighbor(int v, int i) {
#pragma omp atomic
	++query_cnt;
	//assert(i < deg[v]);
	return con[v][i];
}

bool Graph::heavy(int v, int v_deg, long long t_hat) {
	if (v_deg > 2 * m_hat / pow(epsilon_hat * t_hat, 1. / 3))
		return true;
	int len = int(10 * log(n)),
		//s = int(20 * pow(m_hat, 1.5) / (epsilon_hat * epsilon_hat * t_hat));
		s = int(20 * pow(m_hat, 1.5) / (epsilon*epsilon  * t_hat));
	s = max(s, 1);
	long long* x_arr = new long long[len];
	memset(x_arr, 0, len*sizeof(long long));
	for (int i = 0; i < len; ++i) {
#pragma omp parallel
		{
			int pid = omp_get_thread_num(), np = omp_get_num_threads();
			for (int j = pid; j < s; j += np) {
				int x = query_neighbor(v, rand32() % v_deg);
				int x_deg = query_degree(x);
				int u, u_deg;
				if (smaller(x, x_deg, v, v_deg)) {
					u = x; u_deg = x_deg;
				}
				else {
					u = v; u_deg = v_deg;
				}
				double avg_z = 0;
				int r = int(u_deg / sqrt_m) + 1;
				for (int k = 0; k < r; ++k) {
					int w = query_neighbor(u, rand32() % u_deg);
					if ((u == v && !query_pair(x, w)) ||
						(u == x && !query_pair(v, w))) 	//not form a tri
						continue;
					int w_deg = query_degree(w);
					if (smaller(w, w_deg, x, x_deg))
						continue;
					avg_z += double(u_deg) / r;
				}
#pragma omp atomic
				x_arr[i] += long long(avg_z);
			}
		}
		x_arr[i] = x_arr[i] * v_deg / s;
	}
	long long test = GetMidNumNoSort(x_arr, len);
	bool is_heavy = GetMidNumNoSort(x_arr, len) > pow(t_hat, 2. / 3) / pow(epsilon_hat, 1. / 3);
	delete[] x_arr;
	return is_heavy;
}

long long Graph::estimate_with_advice(long long t_hat) {
	const long long c_1 = 20; //10;
	int s_1 = int(c_1 * log(n / epsilon) * n / (pow(epsilon, 3) * pow(t_hat, 1. / 3)));
	//int s_1 = int(c_1 * log(n / epsilon_hat) * n / (pow(epsilon_hat, 3) * pow(t_hat, 1. / 3)));
	s_1 = min(s_1, n);
	int* sample_v = new int[s_1];
	long long* sample_deg = new long long[s_1];	//	use for sample edge by degree
	long long deg_cnt;
	for (int i = 0; i < s_1; ++i) {
		sample_v[i] = rand32() % n;
		long long v_deg = query_degree(sample_v[i]);
		sample_deg[i] = i == 0 ? v_deg - 1 : v_deg + sample_deg[i-1];
	}
	deg_cnt = sample_deg[s_1 - 1] + 1;

	const long long c_2 = 10; // 6;
	long long s_2 = long long(c_2 * pow(log(n), 2) * pow(m_hat, 1.5) / (pow(epsilon, 4) * t_hat));
	//int s_2 = int(c_2 * pow(log(n), 2) * pow(m_hat, 1.5) / (pow(epsilon_hat, 4) * t_hat));
	
	double avg_y = 0;
	for (int i = 0; i < s_2; ++i) {
		// sample v portional to degree
		int v = int(lower_bound(sample_deg, sample_deg + long long(s_1), rand64() % deg_cnt) - sample_deg);
		v = sample_v[v];
		int v_deg = deg[v];// deg of v is already required
		if (v_deg == 0) continue;
		// sample e in N(v) uniformly
		int x = query_neighbor(v, rand32() % v_deg);
		int x_deg = query_degree(x);
		int u, u_deg;
		if (smaller(x, x_deg, v, v_deg)) {
			u = x; u_deg = x_deg;
		}
		else {
			u = v; u_deg = v_deg;
		}
		int r;
		if (u_deg <= sqrt_m)
			r = double(rand()) / RAND_MAX <= double(u_deg) / sqrt_m ? 1 : 0;
		else
			r = int(u_deg / sqrt_m) + 1;
		double avg_z = 0;
		for (int j = 0; j < r; ++j) {
			int w = query_neighbor(u, rand32() % u_deg); // w may be x\v
			if ((u == v && !query_pair(x, w)) ||
				(u == x && !query_pair(v, w))) 	//not form a tri
				continue;
			int w_deg = query_degree(w);
			if (smaller(w, w_deg, x, x_deg))
				continue;
			if (heavy(v, v_deg, t_hat))
				continue;
			int l_cnt = 1 + int(!heavy(x, x_deg, t_hat)) + int(!heavy(w, w_deg, t_hat));
			avg_z += (max(double(u_deg), sqrt_m) / l_cnt) / r;
		}
		avg_y += avg_z / s_2;
	}
	delete[] sample_v;
	delete[] sample_deg;
	return n * deg_cnt * avg_y / s_1;
}

void Graph::estimate() {
	//TODO Feige Alg
	int d_hat = 0;
	query_cnt = 0;
	m_hat = m;
	sqrt_m = sqrt(m_hat);
	double t = omp_get_wtime();

	int tip = 0;
	long long t_max = pow(LLONG_MAX, 1. / 3) < n ? LLONG_MAX : long long(n) * n * n;
	long long t_hat = t_max;	//	estimate of t
	while (t_hat >= 1 && query_cnt <= m_hat) {
		cout << '\n' << ++tip << '\t';
		for (long long iter_t = t_max; iter_t >= t_hat; iter_t /= 2) {
			//cout << '#';
			long long x = t_max;
			for (long long i = 1; i <= log(log(n)) / epsilon * 1; ++i) {  // let c = 1 here
				long long tmp = estimate_with_advice(iter_t);
				x = min(x, tmp);
				if (x < iter_t) break; // optimize
			}			
			if (x >= iter_t) {
				tri_cnt = x;
				tri_time = omp_get_wtime() - t;
				printf("total_time=%0.3lfsec, number of triangles=%lld, number of queries=%lld\n"\
					, tri_time, tri_cnt, query_cnt);
				return;
			}
		}
		t_hat /= 2;
	}
	printf("triangles number may be less than sqrt(m)");
}

void Graph::exact_cnt() {
	double t = omp_get_wtime();
	tri_cnt = 0;
#pragma omp parallel
	{
		int pid = omp_get_thread_num(), np = omp_get_num_threads();
		for (int i = pid; i < n; i+=np) {
			for (int j = 0; j < deg[i]; ++j) {
				int v = con[i][j];
				if (con[i][j] <= i) continue;
				for (int k = j + 1; k < deg[i]; ++k) {
					int w = con[i][k];
					if (query_pair(v, w))
#pragma omp atomic
						++tri_cnt;
				}
			}
		}
	}
	tri_time = omp_get_wtime() - t;
	printf("total_time=%0.3lfsec, number of triangles=%lld\n"\
		, tri_time, tri_cnt);
}
#endif /*FOCUSCOUNT_H_*/
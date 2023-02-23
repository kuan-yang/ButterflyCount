#include "BtfCount.h"

bool get_edge(char* line, int& a, int& b) {

	if (!isdigit(line[0])) return false; //for other datasets
	vector<char*> v_num;
	int len = (int)strlen(line);
	for (int i = 0; i < len; ++i)
		if (!isdigit(line[i]) && line[i] != '.') line[i] = '\0';
		else if (i == 0 || !line[i - 1]) v_num.push_back(line + i);
	//if ((int)v_num.size() != 2) return false;
	sscanf(v_num[0], "%d", &a);
	sscanf(v_num[1], "%d", &b);
	return true;
}

void show_info(string path) {
	Graph* g = new Graph(path);
	printf("n=%d,m=%lld\n", g->n , g->m);
	int maxd = 0;
	double sum_sqr = 0;
	for (int i = 0; i < g->n; ++i) {
		maxd = max(maxd, g->deg[i]);
		sum_sqr += g->deg[i] * 1.0 * g->deg[i];
	}

	printf("maxd = %d, sum_sqr = %0.3lfGB\n", maxd, sum_sqr / (1000000000.0));

	delete g;
}

void txt_to_bin(string path) {
	FILE* fin;
	char line[1024];
	int n = 0, a, b;
	vector<int>* con = new vector<int>[MAXN];
	long long cnt = 0, m = 0;
	fin = fopen((path + "graph.txt").c_str(), "r");
	if (!fin) exit(0);
	printf("Loading text...\n");
	while (fgets(line, 1024, fin)) {
		if (!get_edge(line, a, b)) continue;
		if (a < 0 || b < 0) continue;
		n = max(n, a + 1);
		n = max(n, b + 1);
		if (a == b) continue;
		//if (con[a].capacity() == con[a].size()) con[a].reserve(con[a].size() * INC_FACTOR);
		//if (con[b].capacity() == con[b].size()) con[b].reserve(con[b].size() * INC_FACTOR);
		con[a].emplace_back(b);
		con[b].emplace_back(a);
		if ((++cnt) % (long long)1000000 == 0) printf("%lldM lines finished\n", cnt / 1000000);
	}
	fclose(fin);

	printf("n=%d\n", n);

	int maxd = 0;
	for (int i = 0; i < n; ++i)  //排序去重
		if (con[i].size() > 0) {
			sort(con[i].begin(), con[i].end());
			int p = 1;
			for (int j = 1; j < (int)con[i].size(); ++j)
				if (con[i][j - 1] != con[i][j]) con[i][p++] = con[i][j];
			con[i].resize(p); m += p;
			maxd = max(maxd, p);
		}

	printf("Saving binary...\n");
	int* deg = new int[n];
	for (int i = 0; i < n; ++i) deg[i] = (int)con[i].size();

	FILE* fout = fopen((path + "graph.bin").c_str(), "wb");
	if (!fout) exit(0);
	fwrite(&n, sizeof(int), 1, fout);
	fwrite(&m, sizeof(long long), 1, fout);
	fwrite(deg, sizeof(int), n, fout);

	int* nbr = new int[maxd];
	for (int i = 0; i < n; ++i) {
		int d = (int)con[i].size();
		for (int j = 0; j < d; ++j)
			nbr[j] = con[i][j];
		sort(nbr, nbr + d);
		fwrite(nbr, sizeof(int), d, fout);
	}

	fclose(fout);
	printf("Created binary file, n = %d, m = %lld\n", n, m);
	delete[] con; delete[] deg;  delete[] nbr;
}

void generate_graph(string path, int n, ll m) {
	cout << "Generating graph\n";
	ofstream fout((path + "graph.txt").c_str());
	if (!fout) { cout << "File Open Err!\n"; exit(0); }
	for (ll i = 0; i < m; ++i) {
		int u = rand32() % n, v = rand32() % n;
		if (u == v) { --i; continue; }
		fout << u << " " << v << "\n";
		if (i % (long long)1000000 == 0) printf("%lldM lines finished\n", i / 1000000);
	}
}
void generate_comgraph(string path, int n) {

	cout << "Generating complete graph\n";
	ofstream fout((path + "graph.txt").c_str());
	if (!fout) { cout << "File Open Err!\n"; exit(0); }
	for (int i = 0; i < n; ++i) {
		for (int j = i+1; j < n; ++j) {
			fout << i<< " " << j << "\n";
		}
	}	
}

void btf_count_exact(string path, int n_thread) {
	if (n_thread > 0) omp_set_num_threads(n_thread);
	Graph* g = new Graph(path);
	g->btf_count_exact();
	ofstream fout((path + "config.txt").c_str());
	fout << "n=" << g->n << ", m=" << g->m << ", Btf=" << g->BTF_cnt \
		<< ", total_time=" << g->BTF_time << "\n";
	delete g;
}

void btf_count_hash(string path, int n_thread) {
	if (n_thread > 0) omp_set_num_threads(n_thread);
	Graph* g = new Graph(path);
	g->btf_count_hash();
	ofstream fout((path + "fast_config.txt").c_str());
	fout << "n=" << g->n << ", m=" << g->m << ", Btf=" << g->BTF_cnt \
		<< ", total_time=" << g->BTF_time << "\n";
	delete g;
}

void btf_count_est(string path, double epsilon, int n_thread) {
	if (n_thread > 0) omp_set_num_threads(n_thread);
	Graph* g = new Graph(path);
	g->btf_count_est(epsilon);
	ofstream fout((path + "estimate_config.txt").c_str(),ofstream::app);
	fout << "n=" << g->n << ", m=" <<  g->m  << ", epsilon="<< epsilon << ", estimte_Btf=" << g->BTF_cnt \
		<< ", query_cnt=" << g->query_cnt << ", total_time=" << g->BTF_time << "\n";
	delete g;
}

void btf_count_Espar(string path, double p, bool fast_cnt, int n_thread) {
	if (n_thread > 0) omp_set_num_threads(n_thread);
	Graph* g = new Graph(path);
	g->btf_count_ESpar(p, fast_cnt);

	if (fast_cnt)
		path = path + "fast_Espar_config.txt";
	else
		path = path + "Espar_config.txt";
	ofstream fout(path.c_str(),ofstream::app);
	fout << "n=" << g->n << ", m=" << g->m << ", p=" << p <<", estimte_Btf=" << g->BTF_cnt \
		<< ", query_cnt=" << g->query_cnt << ", total_time=" << g->BTF_time << "\n";
	delete g;
}

void rec_count_arb(string path, int n_thread) {
	if (n_thread > 0) omp_set_num_threads(n_thread);
	Graph* g = new Graph(path);
	g->rec_count_arb();

	ofstream fout((path + "arb_config.txt").c_str(), ofstream::app);
	fout << "n=" << g->n << ", m=" << g->m  << ", estimte_Btf=" << g->BTF_cnt \
		<< ", total_time=" << g->BTF_time << "\n";
	delete g;
}

int main(int argc, char* argv[]) {
	printf("argc=%d\n", argc);
	for (int i = 0; i < argc; ++i)
		printf("argv[%d]=%s\n", i, argv[i]);

	if (setvbuf(stdout, NULL, _IONBF, 0) != 0 || setvbuf(stderr, NULL, _IONBF, 0) != 0)
		printf("failed to set up buffer for output file\n");
	else
		printf("successed to set up buffer for output file\n");

	if (argc > 1) {
		if (strcmp(argv[1], "txt-to-bin") == 0)   // 将txt文件转换为二进制
			txt_to_bin(argv[2]);
		//else if (strcmp(argv[1], "btf-count-exact") == 0)  // 精确统计图中的butterfly算法，
		//	btf_count_exact(argv[2], argc > 3 ? atoi(argv[3]) : -1);
		else if (strcmp(argv[1], "btf-count-est") == 0)  // 我们之前的一个算法
			btf_count_est(argv[2], argc > 3 ? atof(argv[3]) : 0.5, argc > 4 ? atoi(argv[4]) : -1);
		//else if (strcmp(argv[1], "btf-count-hash") == 0)  //   精确统计图中的butterfly算法，
		//	btf_count_hash(argv[2], argc > 3 ? atoi(argv[3]) : -1);
		else if (strcmp(argv[1], "btf-count-arb") == 0)  // 精确确统计图中的butterfly算法
			rec_count_arb(argv[2], argc > 3 ? atoi(argv[3]) : -1);
		else if (strcmp(argv[1], "btf-count-Espar") == 0)  // Espar 算法 -- 用作baseline
			btf_count_Espar(argv[2], argc > 3 ? atof(argv[3]) : 0.2, \
				argc > 4 ? atoi(argv[4]) : 1, argc > 5 ? atoi(argv[5]) : -1);
		else if (strcmp(argv[1], "show-ifo") == 0)
			show_info(argv[2]);
	
		return 0;
	}

	if (strcmp(argv[1], "txt-to-bin") == 0){
		txt_to_bin(argv[2]);
		return 0; 
	}
	/*
	string path = "C:/Users/admin/source/repos/BtfCount/data/wikipedia_link_bug/";
	int np = -1;
	cout << "-------------Espar------------\n";
	Espar(path, 0.2, 0, np);		
	vector<double> p_list={0.1,0.2,0.3,0.5,0.8};
	for(auto i:p_list){
		cout << "-------------fast-Espar------------\n";
		Espar(path, i, 1, np);
	}
	vector<double> acc={0.8,0.5,0.3,0.2,0.05,0.01};
	//vector<double> acc={0.05,0.01}; 
	for(auto i:acc){
		cout << "-------------estimate------------\n"; 
		estimate(path, i, np);
	}
	cout << "-------------fast-exact------------\n";
	fast_exact_count(path, np);
	cout << "-------------exact------------\n";
	exact_count(path, np);
	*/
	return 0;
}
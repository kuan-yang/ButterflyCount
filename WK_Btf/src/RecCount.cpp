#include "RecCount.h"

int get_id(map<string,int> &m_dat, int &n, string st) {
	if(m_dat.find(st) != m_dat.end()) return m_dat[st];
	m_dat[st] = n; return n++;
}

void process_amazon_rating(string path) {
	string files[] = {"ratings_Amazon_Instant_Video.csv", "ratings_Apps_for_Android.csv", "ratings_Automotive.csv", "ratings_Baby.csv",
		"ratings_Beauty.csv", "ratings_Books.csv", "ratings_CDs_and_Vinyl.csv", "ratings_Cell_Phones_and_Accessories.csv", "ratings_Clothing_Shoes_and_Jewelry.csv",
		"ratings_Digital_Music.csv", "ratings_Electronics.csv", "ratings_Grocery_and_Gourmet_Food.csv", "ratings_Health_and_Personal_Care.csv", "ratings_Home_and_Kitchen.csv",
		"ratings_Kindle_Store.csv", "ratings_Movies_and_TV.csv", "ratings_Musical_Instruments.csv", "ratings_Office_Products.csv", "ratings_Patio_Lawn_and_Garden.csv",
		"ratings_Pet_Supplies.csv", "ratings_Sports_and_Outdoors.csv", "ratings_Tools_and_Home_Improvement.csv", "ratings_Toys_and_Games.csv", "ratings_Video_Games.csv"};
	int n_file = 24;
	FILE *fout = fopen((path + string("graph.txt")).c_str(), "w");
	map<string,int> m_user, m_item;
	int n_user = 0, n_item = 0, m = 0;
	char line[1024], st_user[64], st_item[64];
	for( int i = 0; i < n_file; ++i ) {
		printf( "Processing %s\n", files[i].c_str() );
		FILE *fin = fopen( (path + files[i]).c_str(), "r" );
		bool start = true;
		while( fgets(line, 1024, fin) ) {
			if( line[0] == '\n' ) continue;
			for( int j = 0; line[j]; ++j ) if(line[j] == ',') line[j] = ' ';
			sscanf( line, "%s %s", st_user, st_item );
			int now_user = get_id(m_user, n_user, st_user);
			int now_item = get_id(m_item, n_item, st_item);
			if( start ) {printf( "%s %s %d %d\n", st_user, st_item, now_user, now_item); start = false;}
			fprintf( fout, "%d %d\n", now_user, now_item );
			++m;
		}
		fclose(fin);
		printf( "n_user = %d, n_item = %d, m = %d\n", n_user, n_item, m );
	}
	fclose(fout);
	printf( "finish processing amazon rating!\n" );
}

void gen_bipartite(string path_from, string path_to) {
	FILE *fin = fopen( (path_from + "graph.txt").c_str(), "r" );
	FILE *fout = fopen( (path_to + "graph.txt").c_str(), "w" );
	char line[1024];
	int a, b, num_cnt = Graph::get_num_cnt(path_from);
	long long cnt = 0;

	printf( "Loading text, num_cnt=%d...\n", num_cnt );
	while( fgets( line, 1024, fin ) ) {
		if( !Graph::get_edge(line, a, b, num_cnt) ) continue;
		if( a < 0 || b < 0 || a % 2 == b % 2 ) continue;
		if( a % 2 ) swap(a, b);
		a/=2; b/=2;
		fprintf( fout, "%d %d\n", a, b );
		if( (++cnt) % (long long) 10000000 == 0 ) printf( "%lldM lines finished\n", cnt/1000000 );
	}
	fclose( fin );
	fclose( fout );
}

void gen_bipartite_syn(string path_from, string path_to) {
	FILE *fin = fopen( (path_from + "graph.txt").c_str(), "r" );
	FILE *fout = fopen( (path_to + "graph.txt").c_str(), "w" );
	char line[1024];
	int a, b, num_cnt = Graph::get_num_cnt(path_from);
	long long cnt = 0;
	printf( "Loading text, num_cnt=%d...\n", num_cnt );
	while( fgets( line, 1024, fin ) ) {
		//if( !Graph::get_edge(line, a, b, num_cnt) ) continue;
		if (line[0] != 'a') continue;
		vector<char*> v_num;
		int len = (int) strlen(line);
		for( int i = 2; i < len; ++i )
			if( !isdigit(line[i]) && line[i] != '.' ) line[i] = '\0';
			else if(i == 2 || !line[i-1]) v_num.push_back(line+i);
		sscanf( v_num[0], "%d", &a );
		sscanf( v_num[1], "%d", &b );
		if( a < 0 || b < 0) continue;
		fprintf( fout, "%d %d\n", a, b );
		if( (++cnt) % (long long) 10000000 == 0 ) printf( "%lldM lines finished\n", cnt/1000000 );
	}
	fclose( fin );
	fclose( fout );
}

void txt_to_bin(string path, int gorder, int is_bipartite) {
	Graph::txt_to_bin(path, gorder != 0, is_bipartite != 0);
}

void rec_count(string path, int n_thread, int method, int cnt_edge) {
	if( n_thread > 0 ) omp_set_num_threads(n_thread);
	bool is_original = false;
	if( method == METHOD_NAIVE_ORG || method == METHOD_SQR_ORG) is_original = true;
	Graph* g = new Graph(path, method, is_original);
	if( method == METHOD_NAIVE || method == METHOD_NAIVE_ORG )
		g->rec_count_naive(cnt_edge);
	else if( method == METHOD_SQR || method == METHOD_SQR_ORG )
		g->rec_count_sqr(cnt_edge);
	else g->rec_count_arb(cnt_edge != 0);
	delete g;
}

void rec_count_app(string path, int n_thread, int method, int cnt_edge = 0) {
	if( n_thread > 0 ) omp_set_num_threads(n_thread);
	Graph* g = new Graph(path, method, false);
	g->rec_count_arb(cnt_edge != 0);
	long long exactBTF_G = g->btf_G;
	vector<double> prob_vector;
	if (path.find("wikipedia") != std::string::npos)
		prob_vector.push_back(0.002);
	if (path.find("delicious") != std::string::npos)
		prob_vector.push_back(0.003);
	if (path.find("trackers") != std::string::npos)
		prob_vector.push_back(0.023);
	if (path.find("twitter") != std::string::npos)
		prob_vector.push_back(0.008);
	double ttime = 0;
	int MAX_ITERATION = 920;
	vector < pair < double, pair <double, double> > > aux_res;
	vector<int> deg_tmp;
	vector<vector<int>> con_tmp;
	con_tmp.resize(g->n);
	for (int i = 0; i < g->n; i++) {
		deg_tmp.push_back(g->deg[i]);
	}
	for (int i = 0; i < g->n; i++) {
		for( int j = 0; j < g->deg[i]; ++j ) {
			con_tmp[i].push_back(g->con[i][j]);
		}
	}	
	for (int ii = 0; ii < prob_vector.size(); ii++) {
		double prob = prob_vector[ii];
		int iteration = 0;
		ttime = 0;
		aux_res.clear();
		for (iteration = 0; iteration < MAX_ITERATION; iteration++) {
				for( int i = 0; i < g->n; ++i ) {
					g->deg[i] = 0;
				}
				double beg_clock1 = omp_get_wtime();
				prob = prob > 1.0 ? 1.0 : prob;
				random_device rdev_edge_sprs;
				mt19937 eng_edg_sprs(rdev_edge_sprs());
				std::uniform_real_distribution<> dis(0.0, 1.0);
				//srand((unsigned)time(NULL));
				for( int i = 0; i < g->n; ++i ) {
					for( int j = 0; j < deg_tmp[i]; ++j ) {
						if (i < con_tmp[i][j]) {
							int A = i;
							int B = con_tmp[i][j];
							double coin = dis(eng_edg_sprs);
							//double coin = rand() / double(RAND_MAX);
							if (coin <= prob || abs(coin - prob) <= 1e-11) {
								g->con[A][g->deg[A]++] = B;
								g->con[B][g->deg[B]++] = A;
							}
						}
					}
				}
				double end_clock1 = omp_get_wtime();
				double time1 = end_clock1-beg_clock1;
				if( method == METHOD_NAIVE_ORG || method == METHOD_SQR_ORG) g->is_original = true;
				if( g->is_original ) {
					//printf( "obtaining original graph...\n" );
					for( int i = 0; i < g->n; ++i ) {
						for( int j = 0; j < g->deg[i]; ++j ) g->con[i][j] = g->oid[g->con[i][j]];
						sort(g->con[i], g->con[i] + g->deg[i]);//TODO order
					}
				}
				double beg_clock2 = clock();
				if( method == METHOD_NAIVE || method == METHOD_NAIVE_ORG )
					g->rec_count_naive(cnt_edge);
				else if( method == METHOD_SQR || method == METHOD_SQR_ORG )
					g->rec_count_sqr(cnt_edge);
				else g->rec_count_arb(cnt_edge != 0);			
				double end_clock2 = clock();
				double elpased_time = g->btf_time + time1;
				ttime += elpased_time;
				double res = (double)g->btf_G / (prob * prob * prob * prob);
				double error = (res - exactBTF_G) / exactBTF_G * 100.0;
				if (error < 0.0) error *= -1.0;
				aux_res.push_back(make_pair(error, make_pair(elpased_time, res)));
				sort(aux_res.begin(), aux_res.end());
				double Er = aux_res[aux_res.size() / 2].first;
				if (iteration % 10 == 0)
					cout<<"prob: "<<prob<<", method: "<<method<<", iteration: "<<iteration<<",time: "<<ttime<<", current error: "<<error<<", average error: "<<Er<<endl;
		}
		double Er = aux_res[iteration / 2].first;
		double etime = aux_res[iteration / 2].second.first;
		cout<<"prob: "<<prob<<", method: "<<method<<", iteration: "<<iteration<<",time: "<<etime<<",error: "<<Er<<",total time: "<<ttime<<endl;
	}
	delete g;
}

void tri_count(string path, int n_thread, int method, int cnt_edge) {
	if( n_thread > 0 ) omp_set_num_threads(n_thread);
	Graph* g = new Graph(path);
	//g->tri_count(method, cnt_edge != 0);
	delete g;
}

void show_info(string path) {
	Graph *g = new Graph(path);
	printf( "n=%d,n1=%d,m=%lld\n", g->n, g->n1, g->m );
	int maxd = 0;
	double sum_sqr = 0;
	for( int i = 0; i < g->n; ++i ) {
		maxd = max(maxd, g->deg[i]);
		sum_sqr += g->deg[i] * 1.0 * g->deg[i];
	}

	printf( "maxd = %d, sum_sqr = %0.3lfGB\n", maxd, sum_sqr/(1000000000.0) );

	delete g;
}

void test_rec_txt_to_bin() {
	string path = "/Users/lqin/Documents/workspace/data/bipartite/test/";
	txt_to_bin(path, 0, true);
}

void test_rec_count() {
	string path = "/Users/lqin/Documents/workspace/data/bipartite/test/";
	rec_count(path, -1, 1, 0);
}

void test_txt_to_bin() {
	string path = "/Users/lqin/Documents/workspace/data/test_triangle/";
	txt_to_bin(path, 0, false);
}

void test_tri_count() {
	string path = "/Users/lqin/Documents/workspace/data/test_triangle/";
	//tri_count(path, -1, 1, 0);
}

int main(int argc, char *argv[]) {
	printf( "argc=%d\n", argc );
	for( int i = 0; i < argc; ++i )
		printf( "argv[%d]=%s\n", i, argv[i] );

	setvbuf(stdout, NULL, _IONBF, 0);
	setvbuf(stderr, NULL, _IONBF, 0);

	if( argc > 1 ) {
		if(strcmp( argv[1], "txt-to-bin" ) == 0)
			txt_to_bin( argv[2], argc>3?atoi(argv[3]):0, argc>4?atoi(argv[4]):0 );
		else if( strcmp(argv[1], "rec-count" ) == 0 )
			rec_count( argv[2], argc>3?atoi(argv[3]):-1, argc>4?atoi(argv[4]):METHOD_ARB, argc>5?atoi(argv[5]):0 );
		else if( strcmp(argv[1], "rec-count-app" ) == 0 )
			rec_count_app( argv[2], argc>3?atoi(argv[3]):-1, argc>4?atoi(argv[4]):METHOD_ARB, argc>5?atoi(argv[5]):0 );
		else if( strcmp(argv[1], "tri-count" ) == 0 )
			tri_count( argv[2], argc>3?atoi(argv[3]):-1, argc>4?atoi(argv[4]):METHOD_DI, argc>5?atoi(argv[5]):0 );
		else if( strcmp(argv[1], "amazon") == 0 )
			process_amazon_rating(argv[2]);
		else if( strcmp(argv[1], "gen-bi") == 0 )
			gen_bipartite(argv[2], argv[3]);
		else if( strcmp(argv[1], "gen-bi-syn") == 0 )
			gen_bipartite_syn(argv[2], argv[3]);
		else if( strcmp(argv[1], "show-info" ) == 0 )
			show_info(argv[2]);
	}

	if( argc <= 1) {
		test_rec_txt_to_bin();
		test_rec_count();

		//test_txt_to_bin();
		//test_tri_count();
	}

	return 0;
}


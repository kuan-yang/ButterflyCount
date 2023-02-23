#include "FocusCount.h"

bool get_edge(char* line, int& a, int& b) {

	if (!isdigit(line[0])) return false; //for other datasets
	vector<char*> v_num;
	int len = (int)strlen(line);
	for (int i = 0; i < len; ++i)
		if (!isdigit(line[i]) && line[i] != '.') line[i] = '\0';
		else if (i == 0 || !line[i - 1]) v_num.push_back(line + i);
	if ((int)v_num.size() != 2) return false;
	sscanf(v_num[0], "%d", &a);
	sscanf(v_num[1], "%d", &b);
	return true;
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
		n = max(n, a+1);
		n = max(n, b + 1);
		if (a == b) continue;
		//if (con[a].capacity() == con[a].size()) con[a].reserve(con[a].size() * INC_FACTOR);
		//if (con[b].capacity() == con[b].size()) con[b].reserve(con[b].size() * INC_FACTOR);
		con[a].emplace_back(b);
		con[b].emplace_back(a);
		if ((++cnt) % (long long)100000 == 0) printf("%lld00k lines finished\n", cnt / 1000);
		//if ((++cnt) % (long long)1000000 == 0) printf("%lldM lines finished\n", cnt / 1000000);
	}
	fclose(fin);

	printf("n=%d\n", n);

	int maxd = 0;
	for (int i = 0; i < n; ++i)  //ÅÅÐòÈ¥ÖØ
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


int main(int argc, char* argv[]) {
	//txt_to_bin("C:/Users/admin/source/repos/FOCUS2015/data/com-lj.ungraph/");
	Graph g("C:/Users/admin/source/repos/FOCUS2015/data/com-lj.ungraph/", 0, 0.5);
	//g.exact_cnt();
	g.estimate();

	return 0;
}
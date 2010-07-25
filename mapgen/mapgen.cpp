#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <vector>
#include <algorithm>
#include <stack>
#include <queue>
#include <deque>
#include <map>
#include <set>
#include <iostream>
#include <sstream>

#define int64 long long

using namespace std;

string taskname = "map";

struct point {
	int id;
	int x, y;

	point(int tid = 0, int tx = 0, int ty = 0){
		id = tid;
		x = tx;
		y = ty;
	}

	void print() {
		printf("%d %d", x, y);
	}
};

#define MAXN 25000000
#define INF 2000000000

point a[MAXN];
int n;

inline bool inside(int x, int y, int lx, int ly, int rx, int ry) {
	return (lx <= x && x <= rx && ly <= y && y <= ry);
}

int count(int lx, int ly, int rx, int ry) {
	int ret = 0;
	for (int i = 0; i < n; i++) {
		if (inside(a[i].x, a[i].y, lx, ly, rx, ry)) {
			ret++;
		}
	}
	return ret;
}

#define MAXV 5000000 

int oldindex[MAXV];
int ind[MAXN];
int newindex[MAXV];
vector<pair<int, int> > edges[MAXV];
vector<pair<int, int> > newedges[MAXV];
int vis[MAXV];
int needv;
char * graphFileName;

void dfs(int k, int color) {
	vis[k] = color;
	for (int i = 0; i < edges[k].size(); i++) {
		int next = edges[k][i].first;
		if (!vis[next]) {
			dfs(next, color);
		}
	}
	
}

void gen(int dx, int dy) {
	int x, y, lx, ly, rx, ry;
	while (1) {
		x = rand() % dx;
		y = rand() % dy;
		if (count(x, y, dx, dy) < needv) {
			continue;
		}
		lx = x;
		ly = y;
		rx = dx;
		ry = dy;
		for (int i = 0; i < 100; i++) {
			int mx = (lx + rx) / 2;
			int my = (ly + ry) / 2;
			if (count(x, y, mx, my) < needv) {
				lx = mx;
				ly = my;
			} else {
				rx = mx;
				ry = my;
			}
		}
		break;
	}
/*	x = 0;
	y = 0;
	lx = dx;
	ly = dy;
*/
	cerr << "Done." << endl;
	cerr << "Finding maximum connected component... 	";
	int cnt = 1;
	for (int i = 0; i < n; i++) {
		if (inside(a[i].x, a[i].y, x, y, lx, ly)) {
			ind[i] = cnt;
			cnt++;
		} else {
			ind[i] = 0;
		} 
	}
	freopen(graphFileName, "r", stdin);
	char buf[1000];
	while (gets(buf)) {
		if (buf[0] == 'a') {
			int v1, v2, dist;
			sscanf(buf + 2, "%d %d %d", &v1, &v2, &dist);
			if (ind[v1] && ind[v2]) {
				edges[ind[v1]].push_back(make_pair(ind[v2], dist));	
			}
		}
	}
	cnt--;
	int color = 1;
	memset(vis, 0, sizeof(vis));
	assert(cnt <= MAXV);
	for (int i = 1; i <= cnt; i++) {
		if (!vis[i]) {
			dfs(i, color++);
		}
	}
	int * colors = new int[color + 1];
	
	for (int i = 1; i <= cnt; i++) {
		colors[vis[i]]++;
	}
	int maxcolor = 1;
	for (int i = 2; i <= color - 1; i++) {
		if (colors[maxcolor] < colors[i]) {
			maxcolor = i;
		}
	}
	int numv = 0;
	int nume = 0;
	for (int i = 1; i <= cnt; i++) {
		if (vis[i] == maxcolor) {
			numv++;
			newindex[i] = numv;
			oldindex[numv] = i;
			for (int j = 0; j < edges[i].size(); j++) {
				int next = edges[i][j].first;
				if (vis[next] == maxcolor) {
					newedges[numv].push_back(make_pair(next, edges[i][j].second));
					nume++;
				}
			}
		}
	}
	cerr << "Done." << endl;
	cerr << "Number of vertices = " << numv << endl; 
	cerr << "Number of edges = " << nume << endl;
	cerr << "Writing graph to file... 	";
	cout << "p sp " << numv << " " << nume << endl;
	for (int i = 1; i <= numv; i++) {
		for (int j = 0; j < newedges[i].size(); j++) {
			printf("a %d %d %d\n", i, newindex[newedges[i][j].first], newedges[i][j].second);
		}
	}
	cout << "p aux sp co " << numv << endl; 
	for (int i = 1; i <= numv; i++) {
		printf("v %d ", i);
		a[oldindex[i]].print();
		cout << endl;
	}
	cerr << "Done." << endl;
	
	delete[] colors;
}

int main(int argc, char * argv[]) {	
//	freopen((taskname + ".out").c_str(), "w", stdout);

	if (argc != 4) {
		cerr << "Usage:" << endl;
		cerr << "./mapgen COORDINATE_FILE GRAPH_FILE NUMBER_OF_VERTICES" << endl;
		return 0;
	}

	freopen(argv[1], "r", stdin);
	sscanf(argv[3], "%d", &needv);
	graphFileName = argv[2];

	int cnt = 0;
	char buf[1000];
	while (gets(buf)) {
		if (buf[0] == 'p') {
			sscanf(buf, "p aux sp co %d", &n);
			cerr << "Number of vertices on the map = " << n << endl;
		} else if (buf[0] == 'v') {
			int id, x, y;
			sscanf(buf + 2, "%d %d %d", &id, &x, &y);
			a[cnt] = point(id, x, y);
			cnt++;
		}
	}

	cerr << "File with coordinates read." << endl; 
	cerr << "Generating rectangle...	";

	assert(n == cnt);
	int minx = INF;
	int miny = INF;
	int maxx = -INF;
	int maxy = -INF;
	for (int i = 0; i < n; i++) {
		minx = min(minx, a[i].x);
		maxx = max(maxx, a[i].x);
		miny = min(miny, a[i].y);
		maxy = max(maxy, a[i].y);
	}
	//cerr << minx << " " << miny << " " << maxx << " " << maxy << endl;
	int dx = maxx - minx;
	int dy = maxy - miny;
	for (int i = 0; i < n; i++) {
		a[i].x -= minx;
		a[i].y -= miny;
	}
	gen(dx, dy);

	return 0;
}

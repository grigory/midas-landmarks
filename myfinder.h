// (c) Microsoft Corporation. All rights reserved.

#ifndef MYFINDER_H
#define MYFINDER_H

#include <stdlib.h>
#include <stdio.h>
#include "midas_timer.h"
#include "typmidas.h"
#include "graph.h"
#include "landfind.h"
#include "evaluator.h"
#include <vector>
#include <string>
#include <assert.h>
#include <map>
#include <set>
#include "binheap.h"
//#include <omp.h

#define int64 long long

#define GREEDYINIT true

int UPPER = 10;

int NEIGHBOURHOODRADIUS = 1000;
int MYTRY = 20000;
int MYTRY2 = 1000;
int RANDC = 10;
int POINTS = 2048; 
int FIRSTPHASE = 800; 
int LASTPHASE = 990; 

using namespace std;

class MyDijkstra {
private:
	Graph *g;
	BinaryHeap <TCost> *heap;
	TCost *dist;
	VertexId *parent;
	VertexId *scanlist;
	VertexId nscans;

public:
	MyDijkstra(Graph *_g) {
		g = _g;
		VertexId n = g->VertexCount();
		heap = new BinaryHeap<TCost>(n+1);
		parent = new VertexId[n+1];
		dist = new TCost[n+1];
		scanlist = new VertexId[n];
		for (VertexId v=1; v<=n; v++) {dist[v] = MIDAS_INFINITY;}
		nscans = 0;
	}

	~MyDijkstra() {
		delete [] scanlist;
		delete [] dist;
		delete [] parent; 
		delete heap;
	}

	//parent of vertex v in the shortest path tree from the root
	inline VertexId GetParent(VertexId v) const {return parent[v];}

	//distance from the root to vertex v
	inline TCost GetDistance(VertexId v) const {return dist[v];}

	//return the i-th vertex scanned (with i=0...n-1)
	inline VertexId GetScannedVertex(VertexId i) const {return scanlist[i];}
	
	//Run Dijkstra's algorithm from vertex r. Updates the internal state of this
	//object (with parent, distance, and scan order information). Use the functions
	//above to query the internal state.
	pair<TCost, vector<VertexId> > RunDijkstra(VertexId r, int eps = 0) {
		bool verbose = false;
		if (verbose) fprintf (stderr, "Running Dijkstra from %d... ", r);
		VertexId v, n;
		n = g->VertexCount(); //vertices have labels from 1 to n
		if (eps == 0) {
			eps = n;
		}
		if (r<1 || r>n) {
			fprintf (stderr, "ERROR: Dijkstra's starting vertex %d not in [1,%d].\n", r, n);
			exit(-1);
		}

		for (int i = 0; i < nscans; i++) {
			dist[scanlist[i]] = MIDAS_INFINITY;
		}
		
/*		for(int i = 1; i <= n; i++) {
			dist[i] = MIDAS_INFINITY;
		} 
*/		
		vector<VertexId> vertices;
		vertices.reserve(eps);
		nscans = 0;
		heap->Reset();
		
		TCost ret = 0;

		dist[r] = 0;
		parent[r] = r;
		heap->Insert(r,0);
		while (!heap->IsEmpty()) {
			TCost vdist;
			unsigned int v;
			heap->RemoveFirst(v,vdist);

			scanlist[nscans] = v;
			nscans ++;

			if (nscans > eps) {
				continue;
			}
			vertices.push_back(v);	
			ret = max(ret, dist[v]);

			Arc *a, *end;
			for (g->GetBounds(v,a,end); a<end; a++) {
				VertexId w = a->head;
				TCost wdist = vdist + a->cost;
				if (wdist < dist[w]) {
					dist[w] = wdist;
					heap->FixKey(w,wdist);
					parent[w] = v;
				}
			}
		}
		if (verbose) fprintf (stderr, "done (%d nodes visited).\n", nscans);

		return make_pair(ret, vertices);

/*		if (nscans != g->VertexCount()) {
			fprintf (stderr, "ERROR: Dijkstra visited only %d/%d vertices.\n", nscans, g->VertexCount());
			fprintf (stderr, "(Maybe graph is disconnected?)\n");
			exit(-1);
		}*/
	}
};







#define LOWERCANDIDATES 1
#define LOWERITERATIONS 1000000
#define UPPERCANDIDATES 10
#define UPPERITERATIONS 1000000 

#define LOCALOPTIMIZEUPPER false 
#define CORE false 
#define COREDEPTH 1

#define SUBMISSION true 
#define MUL ((double)1.2)

int n;
MyDijkstra * dij;
MIDASTimer * timer;
vector<VertexId> points;
Graph * g;
vector<VertexId> coreVertices;
vector<TCost> radius;
LandmarkEvaluator * le;
vector<VertexId> position;
vector<VertexId> sortedByRadius;

void myEvaluate(VertexId * list) {
	le->EvaluateLandmarks(g, 20, list, 1000000, true);
}

vector<TCost> getRadius(int eps) {
	vector<TCost> ret(n + 1);
	double dt;
	for (VertexId v = 1; v <= n; v++) {
		double t1 = timer->GetElapsedTime();
		pair<TCost, vector<VertexId> > dijres = dij->RunDijkstra(v, eps);
		ret[v] = dijres.first;
		double t2 = timer->GetElapsedTime();
		dt = t2 - t1;
	}
/*	sort(ret.begin(), ret.end());
	for (int i = 0; i < 10; i++) {
		cerr << ret[i] << endl;
	}
	for (int i = ret.size() - 10; i < ret.size(); i++) {
		cerr << ret[i] << endl;
	}*/
	cerr << "Eps Dijkstra Time = " << dt << endl;
	return ret;
}

vector<int> getDegrees(vector<int> & v) {
	vector<int> ret(n + 1, 0);
	for (int i = 1; i <= n; i++) {
		if (!v[i]) {
			continue;
		}
		Arc * start, * end;	
		g->GetBounds(i, start, end);
		for (Arc * a = start; a < end; a++) {
			VertexId w = a->head;
			if (v[w]) {
				ret[w]++;	
			}
		}
	}
	return ret;
}

void deleteDegOne(vector<int> & v) {
	vector<int> deg = getDegrees(v);
	for (int i = 1; i <= n; i++) {
		if (deg[i] == 1) {
			v[i] = 0;
		}
	}
}

vector<VertexId> getCoreVertices(int depth) {
	vector<int> v(n + 1);
	for (int i = 1; i <= n; i++) {
		v[i] = 1;
	}
	for (int i = 0; i < depth; i++) {
		deleteDegOne(v);	
	}
	vector<VertexId> ret;
	ret.reserve(n);
	for (int i = 1; i <= n; i++) {
		if (v[i]) {
			ret.push_back(i);
		}
	}
	return ret;
}

vector<VertexId> genRandomCheckpoints(int k) {
	int * perm = new int[n];
	for (int i = 0; i < n; i++) {
		perm[i] = i + 1;
	}
	random_shuffle(perm, perm + n);
	vector<VertexId> ret;
	ret.reserve(k);
	for (int i = 0; i < k; i++) {
		ret.push_back(perm[i]);
	}
	delete[] perm;
	return ret;
}

map<VertexId, vector<TCost> > dict;

void getDist(VertexId k, vector<VertexId> & points, vector<TCost> & v) {
	if (dict.count(k)) {
		v = dict[k];
		return;
	}
	dij->RunDijkstra(k);
	v.clear();
	v.reserve(points.size());
	for (int i = 0; i < (int)points.size(); i++) {
		v.push_back(dij->GetDistance(points[i]));
	}
	dict[k] = v;
}

void dfs(VertexId k, vector<VertexId> * adj, int * l, int64 * s, int * v, TCost * w) {
	v[k] = 1;
	s[k] = 0;
	for (VertexId i = 0; i < (int)adj[k].size(); i++) {
		VertexId next = adj[k][i];
		if (!v[next]) {
			dfs(next, adj, l, s, v, w);
			if (l[next]) {
				l[k] = 1;
			}
			s[k] += s[next];
		}
	}
	if (l[k]) {
		s[k] = 0;
	} else {
		s[k] += w[k];
	}
}

VertexId avoidleaf(VertexId k, int64 * s, vector<VertexId> *& adj) {
	if ((int)adj[k].size() == 0) {
		return k;
	}
	VertexId bestv = adj[k][0];
	for (int i = 1; i < (int)adj[k].size(); i++) {
		VertexId next = adj[k][i];
		if (s[bestv] < s[next]) {
			bestv = next;
		}
	}
	return avoidleaf(bestv, s, adj);
}

void generateAdjacentLists(MyDijkstra * dij, vector<VertexId> *& adj) {
	adj = new vector<VertexId>[n + 1];
	for (VertexId i = 1; i <= n; i++) {
		VertexId p = dij->GetParent(i);
		adj[p].push_back(i);
	}
}

VertexId findAvoidVertex(const vector<VertexId> & have, const vector<vector<TCost> > & distances) {
	VertexId r = MTRandom::GetInteger(1, n);
	TCost * w = new TCost[n + 1];
	for (VertexId i = 1; i <= n; i++) {
		w[i] = 0;
	}
	for (int i = 0; i < (int)have.size(); i++) {
		if (distances.empty()) {
			dij->RunDijkstra(have[i]);
		}
		for (VertexId j = 1; j <= n; j++) {
			TCost d1;
			TCost d2;
			if (distances.empty()) {
				d1 = dij->GetDistance(j);
				d2 = dij->GetDistance(r);
			} else {
				d1 = distances[i][j];
				d2 = distances[i][r];
			}
			TCost lb = (d1 > d2) ? d1 - d2 : d2 - d1;
			w[j] = max(w[j], lb);
		}
	}
	
	dij->RunDijkstra(r);
	for (VertexId i = 1; i <= n; i++) {
		w[i] = dij->GetDistance(i) - w[i]; 
	}
	vector<VertexId> * adj;
	generateAdjacentLists(dij, adj);


	int * l = new int[n + 1];
	int64 * s = new int64[n + 1];
	int * v = new int[n + 1];

	for (int i = 1; i <= n; i++) {
		l[i] = 0;
		s[i] = 0;
		v[i] = 0;
	} 

	for (int i = 0; i < (int)have.size(); i++) {
		l[have[i]] = 1;	
	}

	dfs(r, adj, l, s, v, w);

	VertexId bestv = 1;
	for (VertexId i = 2; i <= n; i++) {
		if (s[bestv] < s[i]) {
			bestv = i;
		}
	}

	VertexId result = avoidleaf(bestv, s, adj);
	
	delete[] l;
	delete[] s;
	delete[] v;
	delete[] adj;
	delete[] w;
	return result;
}

vector<TCost> getDistancesVector(MyDijkstra * dij) {
	vector<TCost> v;
	v.resize(n + 1);
	for (int j = 1; j <= n; j++) {
		v[j] = dij->GetDistance(j);
	}
	return v;
}

class SubstitutionOptimizer {
protected:
	VertexId * myList;
	int myShift, myK, myCheckpoints, myCandidates;
	string myName;
	vector<VertexId> myPoints;
	vector<TCost> * myMindist;
	vector<vector<TCost> >myDistances;
	vector <int**> myCache;
	int myLen;
	vector<bool> myOutdated;
 
public:
	SubstitutionOptimizer(VertexId * _myList, int _myShift, int _myK, string _myName, vector<VertexId> _myPoints) {
		myList = _myList;
		myShift = _myShift;
		myK = _myK;
		myCheckpoints = _myPoints.size();
		myName = _myName;
		myMindist = new vector<TCost>[myK];
		myPoints = _myPoints;
		myDistances.reserve(myK);
		for (int i = 0; i < myK; i++) {
			getDist(myList[myShift + i], myPoints, myMindist[i]);
			dij->RunDijkstra(myList[myShift + i]);
			vector<TCost> v = getDistancesVector(dij);
			myDistances.push_back(v);
		}
		myLen = myCheckpoints;
		myCache.reserve(myK);
		for (int l = 0; l < myK; l++) {
			int ** cur = new int* [myLen];
			for (int i = 0; i < myLen; i++) {
				cur[i] = new int[myLen];
			}
			myCache.push_back(cur);
		}
		myOutdated.resize(myK);
		for (int i = 0; i < myK; i++) {
			myOutdated[i] = true;
		}
	}
	
	~SubstitutionOptimizer() {
		for (int i = 0; i < myK; i++) {
			for (int j = 0; j < myLen; j++) {
				delete[] myCache[i][j];
			}
			delete[] myCache[i];
		}
		delete[] myMindist;
	}

	virtual VertexId genCandidate() = 0;
	virtual bool iteration() = 0;
	virtual int initIntValue() = 0;
	virtual int64 initLongValue() = 0;
	virtual void updateBest(int64 & best, int64 cur, int & bestpos, int c) = 0;
	inline virtual int getNewBound(int old, int d1, int d2) = 0;

	void optimize(int iterations, int candidates) {
		myCandidates = candidates;	
		int t1 = timer->GetElapsedTime();
		fprintf(stderr, ("Running " + myName + "...\n").c_str());
		for (int i = 0; i < iterations; i++) {
			fprintf(stderr, "Iteration %d...\n", i + 1);
			if (!iteration()) {
				break;
			}
		}
		fprintf(stderr, "Done.\n");
		int t2 = timer->GetElapsedTime();
		cerr << "Time elapsed: " << t2 - t1 << endl;
	}

	vector<VertexId> genCandidates() {
		vector<VertexId> ret;
		ret.reserve(myCandidates);
		for (int i = 0; i < myCandidates; i++) {
			ret.push_back(genCandidate());
		}
		return ret;
	}
	
	bool tryCandidates(int pos, vector<VertexId> & candidates) {
		double t1 = timer->GetElapsedTime();
		int len = myCheckpoints;
		int ** dist;
//		omp_set_num_threads(2);
		if (myOutdated[pos - myShift]) {
			dist = new int * [len];
			for (int i = 0; i < len; i++) {
				dist[i] = new int [len];
			}
	//		vector<vector<int> > dist(len, vector<int> (len));
			for (int i = 0; i < len; i++) {
				for (int j = 0; j < len; j++) {
					int & cur = dist[i][j];
					cur = initIntValue();
//					#pragma omp parallel for
					for (int l = 0; l < pos - myShift; l++) {
						cur = getNewBound(cur, myMindist[l][i], myMindist[l][j]);
					}
//					#pragma omp parallel for
					for (int l = pos - myShift + 1; l < myK; l++) {
						cur = getNewBound(cur, myMindist[l][i], myMindist[l][j]);
					}
				}
			}
			for (int i = 0; i < len; i++) {
				delete[] myCache[pos - myShift][i];
			}
			delete[] myCache[pos - myShift];
			myCache[pos - myShift] = dist;
			myOutdated[pos - myShift] = false;
		} else {
			dist = myCache[pos - myShift];
		}
		double t15 = timer->GetElapsedTime();
		candidates.push_back(myList[pos]);
		int64 best = initLongValue();
		int bestpos = -1;
		for (int c = 0; c < (int)candidates.size(); c++) { 
			VertexId curv = candidates[c];
			vector<TCost> curDist;
			getDist(curv, myPoints, curDist);
			int64 cur = 0;
			for (int i = 0; i < len; i++) {
//				omp_set_num_threads(2);
//				#pragma omp parallel for
				for (int j = 0; j < i; j++) {
					cur += getNewBound(dist[i][j], curDist[i], curDist[j]);
				}
			}
			updateBest(best, cur, bestpos, c);
		}
		bool ret = false;
		double t2 = timer->GetElapsedTime();
		assert(bestpos != -1);
		if (myList[pos] != candidates[bestpos]) {
			ret = true;
			cerr << "Optimized" << endl;	
			VertexId bestv = candidates[bestpos];
			vector<TCost> bestDist;
			getDist(bestv, myPoints, bestDist);
			myMindist[pos - myShift] = bestDist;
			myDistances[pos - myShift] = getDistancesVector(dij);
			myList[pos] = bestv;
//			myEvaluate(myList);
//			for (int i = 0; i < 20; i++) {
//				cerr << position[myList[i]] << " ";
//			}
//			cerr << endl;
			for (int i = 0; i < myK; i++) {
				myOutdated[i] = true;
			}
		}
		double t3 = timer->GetElapsedTime();
		cerr << "Precalc time = " << t15 - t1 << "Check time = " << t2 - t15 << " Dijkstra time = " << t3 - t2 << endl;
		return ret;
	}
};

class LowerSubstitutionOptimizer: public SubstitutionOptimizer {
private: 
	vector<VertexId> myLandmarks;
	vector<vector<TCost> > myLandmarkDistances;
public:
	LowerSubstitutionOptimizer(VertexId * _myList, int _myShift, int _myK, string _myName, vector<VertexId> _myPoints) : SubstitutionOptimizer(_myList, _myShift, _myK, _myName, _myPoints){
	}
	
	virtual VertexId genCandidate() {
		if (MTRandom::GetInteger(1, RANDC) == 1) {
			return MTRandom::GetInteger(1, n);	
		} else {
			return findAvoidVertex(myLandmarks, myLandmarkDistances);
		}
	}

	virtual bool iteration() {
		myLandmarks.clear();
		myLandmarkDistances.clear();
		myLandmarkDistances.reserve(myK);
		int pos = MTRandom::GetInteger(1, myK) - 1;
		for (int i = 0; i < myK; i++) {
			if (i != pos) {
				myLandmarks.push_back(myList[myShift + i]);
				myLandmarkDistances.push_back(myDistances[i]);
			}
		}
		double t1 = timer->GetElapsedTime();
		vector<VertexId> cand = genCandidates();
		double t2 = timer->GetElapsedTime();
		tryCandidates(myShift + pos, cand);
		double t3 = timer->GetElapsedTime();
		cerr << "Generate candidates time = " << t2 - t1 << " Try candidates time = " << t3 - t2 << endl;
		if (SUBMISSION && timer->GetElapsedTime() > LASTPHASE) {
			return false;
		}
		return true;
	}

	virtual int initIntValue() {
		return 0;
	}

	virtual int64 initLongValue() {
		return 0;
	}

	virtual int getNewBound(int old, int d1, int d2) {
		return max(old, abs(d1 - d2));
	}

	virtual void updateBest(int64 & best, int64 cur, int & bestpos, int c) {
		if (best < cur) {
			best = cur;
			bestpos = c;
		}
	}

};

class UpperSubstitutionOptimizer: public SubstitutionOptimizer {
public: 
	UpperSubstitutionOptimizer(VertexId * _myList, int _myShift, int _myK, string _myName, vector<VertexId> _myPoints) : SubstitutionOptimizer(_myList, _myShift, _myK, _myName, _myPoints){
	}
	
	virtual VertexId genCandidate() {
		int r = MTRandom::GetInteger(1, RANDC);
		if (r <= 1) {
			return MTRandom::GetInteger(1, n);	
		} else {
			return sortedByRadius[MTRandom::GetInteger(1, MYTRY)]; 
		}
	}
	
	virtual bool iteration() {
		for (int i = 0; i < myK; i++) {
			vector<VertexId> cand = genCandidates();
			tryCandidates(myShift + i, cand);
			if (SUBMISSION && timer->GetElapsedTime() > FIRSTPHASE) {
				return false;
			}
		}
		return true;
	}
	
	virtual int initIntValue() {
		return MIDAS_INFINITY;
	}
	
	virtual int64 initLongValue() {
		return (int64)MIDAS_INFINITY * (int64)MIDAS_INFINITY;
	}

	virtual int getNewBound(int old, int d1, int d2) {
		return min(old, d1 + d2);
	}

	virtual void updateBest(int64 & best, int64 cur, int & bestpos, int c) {
		if (best > cur) {
			best = cur;
			bestpos = c;
		}
	}
};

class MyLandmarkFinder : public LandmarkFinder {
	public: 
	
		
	void myGenerateAvoid(VertexId * list, VertexId pos) {
		fprintf(stderr, "Generating avoid vertex %d...\n", pos + 1);
		vector<VertexId> have;
		have.reserve(pos);
		for (VertexId i = 0; i < pos; i++) {
			have.push_back(list[i]);
		}
		list[pos] = findAvoidVertex(have, vector<vector<TCost> > ());
		fprintf(stderr, "Done.\n");
	}

	void myGenerateAvoid(VertexId * list, VertexId pos, VertexId k) {
		for (int i = 0; i < k; i++) {
			myGenerateAvoid(list, pos + i);
		}
	}


	void myGenerateRandom(VertexId * list, VertexId k) { 
	       fprintf (stderr, "Generating %d random landmarks... \n", k); 
	       for (VertexId i = 0; i < k; i++) { 
		       list[i] = MTRandom::GetInteger(1, n); 
	       } 
	       fprintf(stderr, "Done.\n"); 
	} 

	void myGenerateGreedilyWithNeighbourhoods(VertexId * list, VertexId k, int radius) {
		vector<pair<TCost, VertexId> > v;
		v.reserve(n);
		for (VertexId i = 1; i <= n; i++) {
			pair<TCost, vector<VertexId> > dijres = dij->RunDijkstra(i, radius);
			v.push_back(make_pair(dijres.first, i));	
		}
		sort(v.begin(), v.end());
		position.resize(n + 1);
		sortedByRadius.resize(n + 1);
		for (VertexId i = 1; i <= n; i++) {
			position[v[i].second] = i;
			sortedByRadius[i] = v[i].second;
		}
		set<VertexId> vis;
		int pos = 0;
		for (int i = 0; i < k; i++) {
			while (vis.count(v[pos].second)) {
				pos++;
			}
			cerr << pos << endl;
			list[i] = v[pos].second;
			pair<TCost, vector<VertexId> > dijres = dij->RunDijkstra(v[pos].second, radius);
			for (int j = 0; j < dijres.second.size(); j++) {
				vis.insert(dijres.second[j]);
			}
		}
		assert(pos < n);
	}

	void myOptimizeUpper(VertexId * list, int shift, int k, int iterations, int candidates) {
		UpperSubstitutionOptimizer uso(list, shift, k, "upper substitution optimizer", points);
		uso.optimize(iterations, candidates);
		
		if (LOCALOPTIMIZEUPPER) {
			while (true) {
				bool optimized = false;
				for (int i = 0; i < k; i++) {
					vector<VertexId> candidates;
					Arc * start, * end;
					g->GetBounds(i, start, end);
					for (Arc * a = start; a < end; a++) {
						VertexId w = a->head;
						candidates.push_back(w);
					}
					bool opt = uso.tryCandidates(shift + i, candidates);
					cerr << list[shift + i] << endl;
					optimized |= opt;
					if (opt) {
						i--;
					}
				}
				if (!optimized) {
					break;
				}
			}
		}
	}
	
	void myOptimizeLower(VertexId * list, int shift, int k, int iterations, int candidates) {
		LowerSubstitutionOptimizer lso(list, shift, k, "lower substitution optimizer", points);
		lso.optimize(iterations, candidates);
	}
	
	void GenerateMyLandmarks(Graph * _g, VertexId k, VertexId *list, MIDASTimer * _timer/*, int radius, int t1, int t2, int _points, int upper, int pool, int randc*/) {
/*		NEIGHBOURHOODRADIUS = radius;
		FIRSTPHASE = t1;
		LASTPHASE = t2;
		POINTS = _points;
		UPPER = upper;
		MYTRY = pool;
		RANDC = randc;
*/

		le = new LandmarkEvaluator();
		g = _g;
		timer = _timer;
		n = g->VertexCount();
		dij = new MyDijkstra(g);
		points = genRandomCheckpoints(POINTS);	
//		coreVertices = getCoreVertices(COREDEPTH);
//		cerr << "CORE SIZE = " << coreVertices.size() << endl;
//		radius = getRadius(100);
		
		if (!GREEDYINIT) {
			this->myGenerateRandom(list, UPPER);
		} else {
			this->myGenerateGreedilyWithNeighbourhoods(list, UPPER, NEIGHBOURHOODRADIUS);
		}
	
		this->myGenerateAvoid(list, UPPER, k - UPPER);
		
//		myEvaluate(list);
		
		this->myOptimizeUpper(list, 0, UPPER, UPPERITERATIONS, UPPERCANDIDATES);
		
//		myEvaluate(list);
		
		this->myOptimizeLower(list, UPPER, k - UPPER, LOWERITERATIONS, LOWERCANDIDATES);

//		myEvaluate(list);
	}
};

#endif

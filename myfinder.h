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
#include "hhfinder.h"
#include "vpair.h"

#include <vector>
#include <string>
#include <assert.h>
#include <map>
#include <set>

#define int64 long long

#define HH false 
#define R false
#define RAD true

#define HHRADIUS 10
#define HHCUTOFF 90000

#define CHECKPOINTEVALUATION true

#define EVALUATIONNUMBER 1000000

#define POOLSIZE 50 
#define POOLGENERATIONSIZE 10

#define GREEDYINIT true 
#define NEIGHBOURHOODRADIUS 1000

#define FACTOR 4
#define MAXCOVERITERATIONS 10

#define MYTRY 20000
#define UPPER 10 
#define POINTS 2048 
#define LOWERCANDIDATES 1
#define LOWERITERATIONS 10000000
#define UPPERCANDIDATES 100
#define UPPERITERATIONS 100000000 

#define LOCALOPTIMIZEUPPER false 
#define CORE false 
#define COREDEPTH 1

#define SUBMISSION true 
#define MUL ((double)1.2)
#define FIRSTPHASE 800 
#define LASTPHASE 1000

int n;
Dijkstra * dij;
MIDASTimer * timer;
vector<VertexId> points;
Graph * g;
vector<VertexId> coreVertices;
vector<TCost> radius;
//MyLandmarkEvaluator * le;
vector<VertexId> position;
vector<VertexId> sortedByRadius;
vector<VertexId> highwayPoints;

map<pair<VertexId, VertexId>, int > edgeNumbers;
int * edgeCover;
int edges;
map<int, vector<int> > coveredEdges;


vector<VertexId> lowerBoundsPool;

class MyLandmarkEvaluator {
	public:
		// compute a lower bound on the distance from s to t
		// using distances to a single landmark
		inline TCost LowerBound(VertexId s, VertexId t, TCost *landdist) const {
			TCost ds, dt;
			ds = landdist[s];
			dt = landdist[t];
			TCost answer = (ds > dt) ? ds - dt : dt - ds;
			return answer;
		}

		// compute the upper bound on the distance from s to t using distances
		// to a single landmark
		inline TCost UpperBound(VertexId s, VertexId t, TCost *landdist) const {
			return (landdist[s] + landdist[t]);
		}

		// compute the ratio between the lower and upper bounds
		// (with a special case when both are zero)
		inline double ComputeRatio(TCost lower, TCost upper) const {
			double ratio = 1.0; //ratio equals 1 if they are equal (including both zero)
			if (upper != lower) { 
				if (lower < 0) {cerr << "ERROR: " << lower << "is an invalid lower bound.\n"; exit(-1);}
				ratio = (double)lower / (double)upper;
			}
			return ratio;
		}

		//Compute the score of a particular query, defined as
		//score = (101^(lower/upper)) - 1
		inline double ComputeScore(TCost lower, TCost upper) const {
			double ratio = ComputeRatio(lower, upper);
			double score = exp (ratio * log (101.0)) - 1;
			return score;
		}


		//------------------------------------------------------------------------
		// Evaluate a list of k landmarks (list[0],...,list[k-1]) on npairs random
		// pairs of distinct vertices. Returns the score for these landmarks.
		// If 'verbose' is true, outputs the score and other statistics to stdout.
		// Runs in O(k*m*log n) time (m=number of edges, n=number of vertices).
		//-----------------------------------------------------------------------*/
		double EvaluateLandmarks(Graph *g, VertexId k, VertexId *list, double npairs, bool verbose) {
			VertexId i, n = g->VertexCount();

			//compute distances from each landmark to all vertices
			TCost **landdist = new TCost*[k];
			Dijkstra *dijkstra = new Dijkstra(g);
			for (i=0; i<k; i++) {
				VertexId v = list[i]; //a landmark
				dijkstra->RunDijkstra(v);
				landdist[i] = new TCost[n+1];
				for (v=1; v<=n; v++) landdist[i][v] = dijkstra->GetDistance(v);
			}
			delete dijkstra;

			double score = EvaluateLandmarks(g, k, landdist, npairs, verbose);

			for (i=0; i<k; i++) delete [] landdist[i];
			delete [] landdist;

			return score;
		}


		//-----------------------------------------------------------------------------------
		// Evaluates a set of k landmarks on a set of npairs pairs of distinct vertices.
		// Variable 'landdist' is a matrix representing the distances between all landmarks
		// and all vertices. "landist[i][v]" must contain the distance between the i-th
		// landmark (with i=0...k-1) and vertex v (with v=1...n). Runs in O(kn) time.
		// Returns the score of these landmarks. If verbose=true, outputs the score and
		// other statistics to stdout.
		//-----------------------------------------------------------------------------------
		double EvaluateLandmarks(Graph *g, VertexId k, TCost **landdist, double npairs, bool verbose) {
			VertexId n = g->VertexCount();
			double sumratio = 0; //sum of all ratios (lower/upper)
			double sumdiff = 0;  //sum of all differences (upper - lower)
			double sumscore = 0; //sum of all scores (101^(lower/upper) - 1)

			for (double p=0; p<npairs; p++) {
				VertexId s = MTRandom::GetInteger(1,n);
				VertexId t = MTRandom::GetInteger(1,n);
				if (s==t) {p--; continue;} //only allow distinct pairs

				//compute current upper and lower bounds
				TCost upper = MIDAS_INFINITY;
				TCost lower = 0;
				for (VertexId i=0; i<k; i++) {
					TCost curlower = LowerBound(s,t,landdist[i]);
					TCost curupper = UpperBound(s,t,landdist[i]);
//					if (curupper < curlower) {
//						fprintf (stderr, "%d: s:%d %d t:%d %d\n", i,s, landdist[i][s], t, landdist[i][t]);
//						exit(-1);
//					}
					lower = max(lower, curlower);
					upper = min(upper, curupper);
				}
				
				TCost diff = upper - lower;
				double ratio = ComputeRatio(lower,upper); //ratio equals 1 if they are equal (including zero)
				double score = ComputeScore(lower, upper); //exp (ratio * log (101.0)) - 1;
				
				sumscore += score;
				sumratio += ratio;
				sumdiff += (double)diff;
			}

			double score = (sumscore / npairs);
			if (verbose) {
				fprintf (stderr, "npairs %.0f\n", npairs);
				fprintf (stderr, "avgdiff %.6f\n", sumdiff /npairs);
				fprintf (stderr, "avgratio %.12f\n", sumratio / npairs);
				fprintf (stderr, "avgscore %.10f\n", score);
			}
			return score;

		}
};

vector<int> getCoveredEdges(VertexId k) {
	if (coveredEdges.count(k)) {
		return coveredEdges[k];
	}
	dij->RunDijkstra(k);
	vector<int> ret;
	ret.reserve(edges);
	for (VertexId i = 1; i <= n; i++) {
		Arc * start, * end;	
		g->GetBounds(i, start, end);
		for (Arc * a = start; a < end; a++) {
			VertexId w = a->head;
			TCost d1 = dij->GetDistance(i);
			TCost d2 = dij->GetDistance(w);
			if (d1 > d2) {
				continue;
			}
			if (d2 - d1 == a->cost) {
				ret.push_back(edgeNumbers[make_pair(i, w)]);		
			}
		}	
	}
	coveredEdges[k] = ret;
	return ret;
}

void addCoveredEdges(int pos, vector<int> addEdges) {
	for (int i = 0; i < addEdges.size(); i++) {
		edgeCover[addEdges[i]] |= (1 << pos);
	}
}

void addCoveredEdges(int pos, VertexId k) {
	vector<VertexId> v = getCoveredEdges(k);
	addCoveredEdges(pos, v);
}

void removeCoveredEdges(int pos) {
	for (int i = 1; i <= edges; i++) {
		edgeCover[i] &= ~(1 << pos);
	}
}

int getCoveredEdges() {
	int ret = 0;
	for (int i = 1; i <= edges; i++) {
		if (edgeCover[i]) {
			ret++;
		}
	}
	return ret;
}

/*void myEvaluate(VertexId * list) {
	le->EvaluateLandmarks(g, 20, list, EVALUATIONNUMBER, true);
}*/

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

void getEdgeNumbers() {
	int cnt = 1;
	for (VertexId i = 1; i <= n; i++) {
		Arc * start, * end;	
		g->GetBounds(i, start, end);
		for (Arc * a = start; a < end; a++) {
			VertexId w = a->head;
			pair<VertexId, VertexId> edge1 = make_pair(i, w);
			pair<VertexId, VertexId> edge2 = make_pair(w, i);
			if (!edgeNumbers.count(edge1)) {
				edgeNumbers[edge1] = cnt;
				edgeNumbers[edge2] = cnt;
				cnt++;
			}
		}	
	}
	cerr << "Edges = " << cnt << endl;
	edges = cnt;
	edgeCover = new int[cnt + 1];
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

void generateAdjacentLists(Dijkstra * dij, vector<VertexId> *& adj) {
	adj = new vector<VertexId>[n + 1];
	for (VertexId i = 1; i <= n; i++) {
		VertexId p = dij->GetParent(i);
		adj[p].push_back(i);
	}
}

VertexId findAvoidVertex(vector<VertexId> have, vector<vector<TCost> > distances) {
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

void generateLowerBoundsPool(int poolSize) {
	cerr << "Generating pool..." << endl;
	set<VertexId> pool;
	while (pool.size() < poolSize) {
		VertexId k = MTRandom::GetInteger(1, n);
		vector<VertexId> have;
		have.push_back(k);
		k = findAvoidVertex(have, vector<vector<TCost> > ());
		have.pop_back();
		have.push_back(k);
		while (have.size() != POOLGENERATIONSIZE) {
			have.push_back(findAvoidVertex(have, vector<vector<TCost> > ()));
		}
		for (int i = 0; i < have.size(); i++) {
			pool.insert(have[i]);
		}
		cerr << pool.size() << endl;
	}
	lowerBoundsPool = vector<VertexId>(pool.begin(), pool.end());
	cerr << "Done..." << endl;
}

vector<TCost> getDistancesVector(Dijkstra * dij) {
	vector<TCost> v;
	v.resize(n + 1);
	for (int j = 1; j <= n; j++) {
		v[j] = dij->GetDistance(j);
	}
	return v;
}

double myEvaluation(VertexId * list, int k) {
	int len = POINTS;
	vector<TCost> * dist = new vector<TCost>[k];
	assert(points.size() == len);
	cerr << "Here 1" << endl;
	for (int i = 0; i < k; i++) {
		dij->RunDijkstra(list[i]);
		dist[i].reserve(len);
		for (int j = 0; j < len; j++) {
			dist[i].push_back(dij->GetDistance(points[j]));
		}
	}
	cerr << "Here 2" << endl;
	double ret;
	int total = 0;
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < i; j++) {
			int bestlb = 0;
			int bestub = MIDAS_INFINITY;
			for (int l = 0; l < k; l++) {
				bestlb = max(bestlb, abs(dist[l][i] - dist[l][j]));
				bestub = min(bestub, dist[l][i] + dist[l][j]); 
			}
			ret += exp(log(101.0) * bestlb / bestub) - 1;
			total++;
		}
	}
	cerr << "Here 3" << endl;
	delete[] dist;
	return ret / total;
}

class SubstitutionOptimizer {
protected:
	VertexId * myList;
	int myShift, myK, myCheckpoints, myCandidates;
	string myName;
	vector<VertexId> myPoints;
	vector<TCost> * myMindist;
	vector<vector<TCost> >myDistances;
	vector <int **> myCache;
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
			int ** cur = new int * [myLen];
			for (int i = 0; i < myLen; i++) {
				cur[i] = new int [myLen];
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

/*	bool tryCandidatesWithNaiveEvaluation(int pos, vector<VertexId> & candidates) {
		int len = EVALUATIONNUMBER;
		int * dist = new int[len];
		bool ret = false;
		for (int i = 0; i < len; i++) {
			int & cur = dist[i];
			VertexId s = eval[i].first;
			VertexId t = eval[i].second;
			for (int j = 0; j < pos - myShift; j++) {
				cur = getNewBound(cur, myDijkstra[j][s], myDijkstra[j][t]);
			}
			for (int j = pos - myShift + 1; j < myK; j++) {
				cur = getNewBound(cur, myDijkstra[j][s], myDijkstra[j][t]);
			}
		}
		double best = 0;
		int bestpos = -1;

		delete[] dist;
		return ret;
	}*/

	bool tryCandidatesWithCheckpointEvaluation(int pos, vector<VertexId> & candidates) {
		double t1 = timer->GetElapsedTime();
		int len = myCheckpoints;
		int ** dist;
//		cerr << "Here 0" << endl;
		if (myOutdated[pos - myShift]) {
			dist = new int * [len];
			for (int i = 0; i < len; i++) {
				dist[i] = new int [len];
			}
	//		vector<vector<int> > dist(len, vector<int> (len));
			for (int i = 0; i < len; i++) {
				for (int j = 0; j < i; j++) {
					int & cur = dist[i][j];
					cur = initIntValue();
					for (int l = 0; l < pos - myShift; l++) {
						cur = getNewBound(cur, myMindist[l][i], myMindist[l][j]);
					}
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
//		cerr << "Here 1" << endl;
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
				for (int j = 0; j < i; j++) {
					cur += getNewBound(dist[i][j], curDist[i], curDist[j]);
				}
			}
			updateBest(best, cur, bestpos, c);
		}
//		cerr << "Here 2" << endl;
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
//			cerr << "My evaluation = " << myEvaluation(myList, 20) << endl;
//			for (int i = 0; i < 20; i++) {
//				cerr << position[myList[i]] << " ";
//			}
//			cerr << endl;
			for (int i = 0; i < myK; i++) {
				myOutdated[i] = true;
			}
		}
//		cerr << "Here 3" << endl;
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
//		return MTRandom::GetInteger(1, n);
		return findAvoidVertex(myLandmarks, myLandmarkDistances);
//		VertexId v = lowerBoundsPool[MTRandom::GetInteger(1, lowerBoundsPool.size()) - 1];
//		cerr << v << endl;
//		return v;
	}

	virtual bool iteration() {
/*		for (int i = 0; i < myK; i++) {
			vector<VertexId> cand = genCandidates();
			tryCandidatesWithCheckpointEvaluation(myShift + i, cand);
			if (SUBMISSION && timer->GetElapsedTime() > LASTPHASE) {
				return false;
			}
		}
		return true;*/
		myLandmarks.clear();
		myLandmarkDistances.clear();
		myLandmarkDistances.reserve(myK);
		for (int i = 0; i < myK; i++) {
			int rand = MTRandom::GetInteger(1, n);
			if (rand > n / 2) {
				myLandmarks.push_back(myList[myShift + i]);
				myLandmarkDistances.push_back(myDistances[i]);
			}
		}
		double t1 = timer->GetElapsedTime();
		vector<VertexId> cand = genCandidates();
		double t2 = timer->GetElapsedTime();
		for (int i = 0; i < myK; i++) {
			if (CHECKPOINTEVALUATION) {
				tryCandidatesWithCheckpointEvaluation(myShift + i, cand);
			} else {
//				tryCandidatesWithNaiveEvaluation(myShift + i, cand);
			}
			if (SUBMISSION && timer->GetElapsedTime() > LASTPHASE) {
				return false;
			}
		}
		double t3 = timer->GetElapsedTime();
		cerr << "Generate candidates time = " << t2 - t1 << " Try candidates time = " << t3 - t2 << endl;
		return true;
	}

	virtual int initIntValue() {
		return 0;
	
	}

	virtual int64 initLongValue() {
		return 0;
	}

	virtual int getNewBound(int old, int d1, int d2) {
		return max(old, (abs(d1 - d2)));
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
		if (HH) {
			return highwayPoints[MTRandom::GetInteger(1, highwayPoints.size()) - 1];
		} else if (R) {
			return MTRandom::GetInteger(1, n);
		} else if (RAD) {
			return sortedByRadius[MTRandom::GetInteger(1, MYTRY)]; 
		} else {
			assert(false);
		}
		if (!CORE) {
			return MTRandom::GetInteger(1, n);
		} else { 
			return coreVertices[MTRandom::GetInteger(0, coreVertices.size() - 1)];
		}
	}
	
	virtual bool iteration() {
		for (int i = 0; i < myK; i++) {
			vector<VertexId> cand = genCandidates();
			if (CHECKPOINTEVALUATION) {
				tryCandidatesWithCheckpointEvaluation(myShift + i, cand);
			} else {	
//				tryCandidatesWithNaiveEvaluation(myShift + i, cand);
			}
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
		
/*		if (LOCALOPTIMIZEUPPER) {
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
					bool opt = uso.tryCandidatesWithCheckpointEvaluation(shift + i, candidates);
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
		}*/
	}
	
	void myOptimizeLower(VertexId * list, int shift, int k, int iterations, int candidates) {
		LowerSubstitutionOptimizer lso(list, shift, k, "lower substitution optimizer", points);
		lso.optimize(iterations, candidates);
	}

	pair<int, int> getOptimalReplacement(vector<VertexId> v, vector<VertexId> c) {
		pair<int, int> ret = make_pair(-1, -1);
		int best = getCoveredEdges();
		for (int i = 0; i < v.size(); i++) {
			removeCoveredEdges(i);
			for (int j = 0; j < c.size(); j++) {
				addCoveredEdges(i, c[j]);
				int cur = getCoveredEdges();		
				if (best < cur) {
					best = cur;
					ret = make_pair(i, j);
				}
				removeCoveredEdges(i);
			}
			addCoveredEdges(i, v[i]);
		}
		return ret;
	}

	void myGenerateMaxCover(VertexId * list, int shift, int k, int factor) {
		set<VertexId> candidates;
		vector<VertexId> have;
		for (int i = 0; i < shift; i++) {
			have.push_back(list[i]);
		}
		for (int i = 0; i < k; i++) {
			have.push_back(findAvoidVertex(have, vector<vector<TCost> > ()));
		}
		for (int i = 0; i < have.size(); i++) {
			candidates.insert(have[i]);
		}
		vector<VertexId> newhave;
		for (int i = 0; i < k; i++) {
			newhave.push_back(have[shift + i]);
		}
		int iterations = 0;
		while (1) {
			if (candidates.size() > factor * k) {
				break;
			}
			vector<VertexId> oldhave = newhave;
			newhave.clear();
			for (int i = 0; i < k; i++) {
				if (MTRandom::GetInteger(1, n) > n / 2) {
					newhave.push_back(oldhave[i]);
				}
			}
			while (newhave.size() != k) {
				newhave.push_back(findAvoidVertex(newhave, vector<vector<TCost> > ()));
				iterations++;
			}
			for (int i = 0; i < newhave.size(); i++) {
				candidates.insert(newhave[i]);
			}
			if (iterations > factor * k) {
				break;
			}
		}
		cerr << "Candidates = " << candidates.size() << endl;
		vector<VertexId> vcandidates(candidates.begin(), candidates.end());
		getEdgeNumbers();
		vector<VertexId> best;
		int bestcover = 0;
		for (int it = 0; it < MAXCOVERITERATIONS; it++) {
			cerr << "MAXCOVER ITERATION = " << it << endl;
			for (int i = 1; i <= edges; i++) {
				edgeCover[i] = 0;
			}
			set<VertexId> init;
			while (init.size() != k) {
				init.insert(vcandidates[MTRandom::GetInteger(1, vcandidates.size()) - 1]);
			}
			vector<VertexId> cur(init.begin(), init.end());
			for (int i = 0; i < cur.size(); i++) {
				addCoveredEdges(i, cur[i]);	
			}
			cerr << getCoveredEdges() << endl;
			int mycnt = 0;
			while (1) {
				cerr << "Replacement iteration = " << mycnt << endl;
				pair<int, int> p = getOptimalReplacement(cur, vcandidates);
				if (p == make_pair(-1, -1)) {
					break;
				}
				int pos = p.first;
				removeCoveredEdges(pos);
				addCoveredEdges(pos, vcandidates[p.second]);
				cur[pos] = vcandidates[p.second];
				cerr << getCoveredEdges() << endl;
				mycnt++;
			}
			int curcover = getCoveredEdges();
			if (curcover > bestcover) {
				bestcover = curcover;
				best = cur;
			}
		}
		assert(best.size() == k);
		for (int i = 0; i < k; i++) {
			list[shift + i] = best[i];
		}
	}

	void GenerateMyLandmarks(Graph * _g, VertexId k, VertexId *list, MIDASTimer * _timer) {
		MyLandmarkEvaluator * le = new MyLandmarkEvaluator();
		g = _g;
		timer = _timer;
		n = g->VertexCount();
		dij = new Dijkstra(g);
		points = genRandomCheckpoints(POINTS);	
//		coreVertices = getCoreVertices(COREDEPTH);
//		cerr << "CORE SIZE = " << coreVertices.size() << endl;
//		radius = getRadius(100);

/*		Graph * tmp = g;
		int oldedges = -1;
		while (1) {
			if (tmp->EdgeCount() == oldedges) {
				set<VertexId> vertices;
				for (int i = 1; i <= n; i++) {
					Arc * start, * end;	
					tmp->GetBounds(i, start, end);
					for (Arc * a = start; a < end; a++) {
						VertexId w = a->head;
						vertices.insert(i);
						vertices.insert(w);
					}						
				}
				cerr << "Vertices = " << vertices.size() << endl;
				highwayPoints = vector<VertexId>(vertices.begin(), vertices.end());
				freopen("points.out", "w", stdout);
				cout << highwayPoints.size() << endl;
				for (int i = 0; i < highwayPoints.size(); i++) {
					cout << highwayPoints[i] << " ";
				}
				fclose(stdout);
				break;
			
			}
			oldedges = tmp->EdgeCount();
			HHFinder * hhf = new HHFinder(tmp, HHRADIUS);
			set<VertexPair> edges;
			for (int i = n; i >= 1; i--) {
				hhf->RunHHFinder(i, edges);	
			}
			delete hhf;
			if (tmp != g) {
				delete tmp;
			}
			VertexPair * pairs = new VertexPair[edges.size()];
			int i = 0;
			for(set<VertexPair>::iterator it = edges.begin(); it != edges.end(); it++) {
				pairs[i++] = *it;
			}
			int * vis = new int[n + 1];
			for (int i = 1; i <= n; i++) {
				vis[i] = 0;
			}

			while (1) {
				bool updated = false;
				
				int * cnt = new int[n + 1];
				for (int i = 1; i <= n; i++) {
					cnt[i] = 0;
				}

				for (int i = 0; i < edges.size(); i++) {
					int v1 = pairs[i].v;
					int v2 = pairs[i].w;
					if (!vis[v1] && !vis[v2]) {
						cnt[v1]++;
						cnt[v2]++;
					}
				}
			
				for (int i = 1; i <= n; i++) {
					if (cnt[i] < 2 && !vis[i]) {
						vis[i] = 1;
						updated = true;
					}
				}

				delete[] cnt;

				if (!updated) {
					break;
				}
			}

			VertexPair * newedges = new VertexPair[edges.size()];
			int count = 0;
			for (int i = 0; i < edges.size(); i++) {
				int v1 = pairs[i].v;
				int v2 = pairs[i].w;
				if (!vis[v1] && !vis[v2]) {
					newedges[count++] = pairs[i];
				}
			}

			tmp = new Graph(count, n);
			tmp->ProcessPairs(newedges, count);
			delete[] pairs;
			delete[] vis;
			delete[] newedges;
			cerr << "New edges = " << count << " Edges = " << edges.size() << endl;
		}


*/

		freopen("points.out", "r", stdin);
		int num;
		cin >> num;
		for (int i = 0; i < num; i++) {
			int tmp;
			cin >> tmp;
			highwayPoints.push_back(tmp);
		}

		if (!GREEDYINIT) {
			this->myGenerateRandom(list, UPPER);
		} else {
			this->myGenerateGreedilyWithNeighbourhoods(list, UPPER, NEIGHBOURHOODRADIUS);
		}
	
		this->myGenerateAvoid(list, UPPER, k - UPPER);
		
	//	this->myGenerateMaxCover(list, UPPER, k - UPPER, FACTOR);
	//	return;

//		myEvaluate(list);
//		cerr << "My evaluation = " << myEvaluation(list, k) << endl;

		this->myOptimizeUpper(list, 0, UPPER, UPPERITERATIONS, UPPERCANDIDATES);
		
//		myEvaluate(list);
//		cerr << "My evaluation = " << myEvaluation(list, k) << endl;
		
//		generateLowerBoundsPool(POOLSIZE);
		this->myOptimizeLower(list, UPPER, k - UPPER, LOWERITERATIONS, LOWERCANDIDATES);

//		myEvaluate(list);
//		cerr << "My evaluation = " << myEvaluation(list, k) << endl;
	}
};

#endif

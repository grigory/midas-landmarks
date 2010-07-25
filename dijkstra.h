// (c) Microsoft Corporation. All rights reserved.

#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include "typmidas.h"
#include "binheap.h"
#include <assert.h>
#include <vector>

using namespace std;

class Dijkstra {
private:
	Graph *g;
	BinaryHeap <TCost> *heap;
	TCost *dist;
	VertexId *parent;
	VertexId *scanlist;
	VertexId nscans;

public:
	Dijkstra(Graph *_g) {
		g = _g;
		VertexId n = g->VertexCount();
		heap = new BinaryHeap<TCost>(n+1);
		parent = new VertexId[n+1];
		dist = new TCost[n+1];
		scanlist = new VertexId[n];
		for (VertexId v=1; v<=n; v++) {dist[v] = MIDAS_INFINITY;}
		nscans = 0;
	}

	~Dijkstra() {
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


#endif

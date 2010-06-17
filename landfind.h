// (c) Microsoft Corporation. All rights reserved.

#ifndef LANDFIND_H
#define LANDFIND_H

#include "typmidas.h"
#include "dijkstra.h"
#include "graph.h"
#include "mt_random.h"
#include "midas_timer.h"
#include <stdio.h>
#include <stdlib.h>

class LandmarkFinder {
	public:
		LandmarkFinder() {}
		~LandmarkFinder() {}

		void GenerateLandmarks(Graph *g, VertexId k, VertexId *list, MIDASTimer *timer) {
			GenerateRandom(g, k, list, timer);
		}

		//-------------------------------------
		// Pick k vertices uniformly at random
		//-------------------------------------
		void GenerateRandom(Graph *g, VertexId k, VertexId *list, MIDASTimer *timer) {
			VertexId n = g->VertexCount();
			for (VertexId i=0; i<k; i++) {
				list[i] = MTRandom::GetInteger(1,n);
			}
		}

		//----------------------------------------------------------------------------
		// Generate a set of k 'farthest' landmarks:
		// - Landmark 0 is the farthest vertex from a random vertex.
		// - Landmark i (0<i<k) is the farthest vertex from the previous i landmarks.
		//----------------------------------------------------------------------------

		void GenerateFarthest(Graph *g, VertexId k, VertexId *list, MIDASTimer *timer) {
			fprintf (stderr, "Generating %d farthest landmarks... ", k);
			Dijkstra *dij = new Dijkstra(g);
			VertexId n = g->VertexCount();

			//fardist[v]: distance from v to the closest landmark
			TCost *fardist = new TCost[n+1]; 
			VertexId v, bestv;
			for (v=1; v<=n; v++) fardist[v] = MIDAS_INFINITY;

			//select the first landmark: farthest vertex from a random vertex r
			VertexId r = MTRandom::GetInteger(1,n); //pick a random reference point (not a landmark)
			dij->RunDijkstra(r); //run Dijkstra from there
			bestv = 1;
			for (v=2; v<=n; v++) {
				if (dij->GetDistance(v) > dij->GetDistance(bestv)) bestv = v;
			}
			list[0] = bestv; //farthest vertex from r is the landmark

			//remaining landmarks: pick vertex that is farthest from the existing landmarks
			for (VertexId i=1; i<k; i++) {
				//update farthest information (to include latest landmark)
				r = list[i-1]; 
				dij->RunDijkstra(r);
				for (v=1; v<=n; v++) {
					TCost vdist = dij->GetDistance(v);
					if (vdist<fardist[v]) {fardist[v] = vdist;}
				}

				//pick farthest vertex as next landmark
				bestv = 1;
				for (v=2; v<=n; v++) {if (fardist[v]>fardist[bestv]) bestv = v;}	
				list[i] = bestv;
			}

			delete [] fardist;
			delete dij;
			fprintf (stderr, "done.\n");
		}
};


#endif

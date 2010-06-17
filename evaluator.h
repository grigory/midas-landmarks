// (c) Microsoft Corporation. All rights reserved.

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <math.h>
#include "typmidas.h"
#include "dijkstra.h"

class LandmarkEvaluator {
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
					if (curupper < curlower) {
						fprintf (stderr, "%d: s:%d %d t:%d %d\n", i,s, landdist[i][s], t, landdist[i][t]);
						exit(-1);
					}
					if (curlower > lower) lower = curlower;
					if (curupper < upper) upper = curupper;
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
				fprintf (stdout, "npairs %.0f\n", npairs);
				fprintf (stdout, "avgdiff %.6f\n", sumdiff /npairs);
				fprintf (stdout, "avgratio %.12f\n", sumratio / npairs);
				fprintf (stdout, "avgscore %.10f\n", score);
			}
			return score;

		}
};

#endif

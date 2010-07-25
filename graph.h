// (c) Microsoft Corporation. All rights reserved.

#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

#include "vpair.h"

using namespace std;

struct Arc {
	int head; //the 'other' endpoint
	int cost; //cost of the arc
};

class Graph {

	private:
		int *first;
		Arc *alist;
		int nedges; //number of (undirected) edges, with ids in [0,nedges-1]
		int nvertices; //number of vertices, with ids in [1,nvertices]

		//----------------------------------------------------------
		// Create adjacency lists using the 'acount' arcs in vpair. 
		// Assumes every arc in vpair is undirected and eliminates 
		// parallel edges by keeping only the lowest-cost version
		//----------------------------------------------------------
	public:
		void ProcessPairs(VertexPair *vpair, int acount) {
			cerr << acount << endl;
//			cerr << "Here 0" << endl;
			std::sort(&vpair[0], &vpair[acount]);
			bool verbose = false;

			//eliminate parallel arcs; keep only the cheapest
			int i,v;
			int lastgood = 0;
//			cerr << "Here 1" << endl;
			for (i=1; i<acount; i++) {
				if (vpair[i]==vpair[lastgood]) { //arcs have same endpoints
					if (vpair[lastgood].cost > vpair[i].cost) {
						vpair[lastgood].cost = vpair[i].cost;
					}
				} else { //different endpoints: new arc
					lastgood ++;
					vpair[lastgood] = vpair[i];
				}
			}
//			cerr << "Here 2" << endl;
			nedges = lastgood + 1;
			if (verbose) fprintf (stderr, "There are %d edges and %d vertices.\n", nedges, nvertices);

			//compute the degrees of all vertices in the graph
			int *degree = new int [nvertices+1];
			for (v=0; v<=nvertices; v++) {degree[v] = 0;}
//			cerr << "Here 3" << endl;
			for (i=0; i<nedges; i++) {
				assert(1 <= vpair[i].v && vpair[i].v <= nvertices);
				assert(1 <= vpair[i].w && vpair[i].w <= nvertices);
				degree[vpair[i].v] ++;
				degree[vpair[i].w] ++;
			}
//			cerr << "Here 4" << endl;

			//initialize adjacency lists
			first = new int [nvertices+2]; //skip 0, need first[n+1]
			alist = new Arc [2*nedges];

			//set first values
			first[0] = 0;
			for (v=1; v<=nvertices+1; v++) {first[v] = first[v-1] + degree[v-1];}
//			cerr << "Here 5" << endl;

			//actually add arcs
			for (i=0; i<nedges; i++) {
				v = vpair[i].v;
				int w = vpair[i].w;
				alist[first[v]].cost = vpair[i].cost;
				alist[first[w]].cost = vpair[i].cost;
				alist[first[v]++].head = w;
				alist[first[w]++].head = v;
			}
//			cerr << "Here 6" << endl;

			//reset first values
			first[0] = 0;
			for (v=1; v<=nvertices+1; v++) {first[v] = first[v-1] + degree[v-1];}
//			cerr << "Here 7" << endl;

			delete [] degree;
		}


	public:
		//constructor (empty graph)
		Graph(int te = 0, int tv = 0) {
			first = NULL;
			alist = NULL;
			nedges = te;
			nvertices = tv;
		}

		//number of vertices in the graph
		inline int VertexCount() const {return nvertices;}

		//number of edges in the graph
		inline int EdgeCount() const {return nedges;}

		//degree (number of incident edges) for vertex v
		inline int GetDegree (int v) const {return (first[v+1] - first[v]);}


		//initialize pointers (start, end) to the adjacency list of v
		//- start: first arc in the list
		//- end: element *after* the last arc in the list
		inline void GetBounds (int v, Arc *&start, Arc *&end) const {
			start = &(alist[first[v]]);
			end = &(alist[first[v+1]]);
		}

		//output adjacency lists of all vertices in the graph
		void OutputGraph (FILE *output) {
			for (int v=1; v<nvertices; v++) {
				fprintf (stderr, "%d:", v);
				Arc *a, *end;
				for (GetBounds(v, a, end); a<end; a++) {
					fprintf (stderr, " %d", a->head);
				}
				fprintf (stderr, "\n");
			}
		}


		//read graph in DIMACS format
		void ReadDimacs(const char *filename) {
			bool verbose = false;
			VertexPair *vpair = NULL;

			if (verbose) fprintf (stderr, "Reading graph %s.\n", filename); 

			FILE *input = fopen (filename, "r");
			if (input == NULL) {
				fprintf (stderr, "Could not open file <%s>.\n", filename);
				exit(-1);
			}

			const int bufsize = 1024;
			char buffer[bufsize + 1];

			int narcs;
			nvertices = narcs = 0;

			int a, b, c;
			int line = 0;
			int acount = 0;
			while (fgets(buffer,bufsize,input)) {
				line ++;
				switch (buffer[0]) {
					case 'p':
						if (sscanf(buffer, "p sp %d %d", &a, &b) == 2) {
							nvertices = a;
							narcs = b;
							vpair = new VertexPair[narcs];
						}
						break;

					case 'a':
						if (sscanf(buffer, "a %d %d %d", &a, &b, &c)==3) {
							vpair[acount++].SetValues(a,b,c);
						}
						break;
				}
			}
			fclose(input);
			if (verbose) fprintf (stderr, "Read %d arcs from the input.\n", acount);

			ProcessPairs(vpair, acount);

			delete [] vpair;
		}

		~Graph() {
			if (first) delete [] first;
			if (alist) delete [] alist;
		}
};

#endif

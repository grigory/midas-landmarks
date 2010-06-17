// (c) Microsoft Corporation. All rights reserved.

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <string>
#include "graph.h"
#include "typmidas.h"
#include "mt_random.h"
#include "evaluator.h"
#include "myfinder.h"
#include "midas_timer.h"
using namespace std;


// Generate 'nl' landmarks for graph 'g' using method 'method'.
// The maximum allowed time (in seconds) for landmark generation is 'maxsec'.
void GenerateLandmarks(Graph *g, int nl, int method, double maxsec) {
	VertexId *list = new VertexId[nl];
	MyLandmarkFinder *lf = new MyLandmarkFinder();
	MIDASTimer *timer = new MIDASTimer (maxsec);

	switch (method) {
		case 0: lf->GenerateMyLandmarks(g, nl, list, timer); break; //your method
		case 1: lf->GenerateRandom(g, nl, list, timer); break;      //random
		case 2: lf->GenerateFarthest(g, nl, list, timer); break;    //farthest
		default:
			fprintf (stderr, "ERROR: %d is not a valid landmark generation method.\n", method);
			exit(-1);
	}
	double t = timer->GetElapsedTime();
	if (t < 0) t = 0;
	fprintf (stderr, "%d landmarks generated in %.4f seconds.\n", nl, t);
	

	//output landmark file
	fprintf (stdout, "p land %d\n", nl);
	fprintf (stdout, "t %.4f\n", t);
	for (int i=0; i<nl; i++) {fprintf (stdout, "l %d\n", list[i]);}

	delete timer;
	delete lf;
	delete [] list;
}

//------------------------------------------------------------
// Read landmarks from file and add them to list in positions
// 0 to k-1, where k is the number of landmarks. Returns k.
//------------------------------------------------------------
int ReadLandmarks(VertexId *&list, Graph *g, FILE *file) {
	int nl = 0;
	double gentime;

	fscanf (file, "p land %d\n", &nl);
	fscanf (file, "t %lg\n", &gentime);
	list = new VertexId [nl];
	int i = 0;
	VertexId v, n = g->VertexCount();
	while (fscanf(file, "l %u\n", &v)==1) {
		if (i>=nl) {fprintf (stderr, "WARNING: too many landmarks; keeping only the first %d.\n", nl); break;}
		if (v<1 || v>n) {
			fprintf (stderr, "ERROR: landmark %d is not in [1,%d].\n", v, n);
			exit(-1);
		}
		list[i++] = v;
	}
	if (i<nl) {
		fprintf (stderr, "ERROR: too few landmarks (%d/%d).\n", i, nl);
		exit(-1);
	}
	fprintf (stderr, "Read landmark file with %d landmarks.\n", nl);
	fprintf (stdout, "landmarks %d\n", nl);
	fprintf (stdout, "gentime %.6f\n", gentime);
	return nl;
}


void EvaluateLandmarks(Graph *g, double npairs) {
	VertexId *list = NULL;
	int nl = ReadLandmarks(list, g, stdin); //read landmarks into 'list'
	LandmarkEvaluator *le = new LandmarkEvaluator();
	le->EvaluateLandmarks(g, nl, list, npairs, true);
	delete le;
	delete [] list;
}


void ShowUsage(char *prg) {
	printf (" Generation: midas -generate <graph> [-method <nn>] [-nl <nn>] [-seed <nn>]\n");
	printf ("    Landmark file will be output to stdout.\n");
	printf ("    <graph>: input graph in DIMACS format\n");
	printf ("    -nl <nn>: number of landmarks to generate [default=20]\n");
	printf ("    -method <nn>: landmark generation method. Valid options:\n");
	printf ("      0: your method [default]\n");
	printf ("      1: random\n");
	printf ("      2: farthest\n");
	printf ("    -seed <nn>: seed for the random number generator [default=1]\n");
	printf ("    -timebound <nn>: maximum allowed time (in seconds) for landmark generation [default=300]\n\n");
	printf (" Evaluation: midas -evaluate <graph> [-np <nn>] [-seed <nn>]\n");
	printf ("    Landmark file is read from stdin.\n");
	printf ("    -np <nn>: number of random pairs to evaluate [default=1000000]\n");
	printf ("    -seed <nn>: seed for the random number generator [default=1]\n\n");
	printf (" Typical use: \"midas -generate test.gr | midas -evaluate test.gr\"\n");
	exit(-1);
}

int main(int argc, char** argv) {
	if (argc < 3) {ShowUsage(argv[0]);}
	int seed = 1; //random seed

	int landmethod = 0; 
	int nl = 20; //number of landmarks to generate
	double npairs = 1000000; //number of pairs to test (double to allow large values)
	double timebound = 1000; // time allowed for landmark generation

	char *filename = argv[2];
	Graph *g = new Graph();
	g->ReadDimacs(filename);

	if (argc % 2 == 0) {
		fprintf (stderr, "Invalid input parameters.\n");
		ShowUsage(argv[0]);
	}

	for (int i=3; i<argc; i+=2) {
		if (strcmp(argv[i], "-seed")==0) {seed = atoi(argv[i+1]); continue;}
		if (strcmp(argv[i], "-nl")==0) {nl = atoi(argv[i+1]); continue;}
		if (strcmp(argv[i], "-np")==0) {npairs = atof(argv[i+1]); continue;}
		if (strcmp(argv[i], "-method")==0) {landmethod = atoi(argv[i+1]); continue;}
		if (strcmp(argv[i], "-timebound")==0) {timebound = atof(argv[i+1]); continue;}
		fprintf (stderr, "Unrecognized input parameter (%s).\n", argv[i]);
		ShowUsage(argv[0]);
	}

	MTRandom::Randomize(seed);
	if (strcmp(argv[1], "-evaluate")==0) {
		EvaluateLandmarks(g, npairs);
	} else if (strcmp(argv[1], "-generate")==0) {
		GenerateLandmarks(g, nl, landmethod, timebound);
	} else ShowUsage(argv[0]);
	delete g;

	return 0;
}

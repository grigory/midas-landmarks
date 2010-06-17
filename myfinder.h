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

class MyLandmarkFinder : public LandmarkFinder {
	public: 
		void GenerateMyLandmarks(Graph *g, VertexId k, VertexId *list, MIDASTimer *timer) {
			fprintf (stderr, "\nMyLandmarkFinder::GenerateMyLandmark called (myfinder.h).\n");
			fprintf (stderr, "You should modify this function to generate your own landmarks.\n");
			fprintf (stderr, "Since it's not implemented, generating random landmarks for now.\n\n");
			this->GenerateRandom(g, k, list, timer);
		}
};

#endif

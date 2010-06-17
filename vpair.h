// (c) Microsoft Corporation. All rights reserved.

#ifndef VPAIR_H
#define VPAIR_H

struct VertexPair {
	int v;
	int w;
	int cost;

	void SetValues(int a, int b) {
		if (a<b) {v=a; w=b;}
		else {v=b; w=a;}
		cost = 0;
	}

	void SetValues(int a, int b, int c) {
		SetValues(a,b);
		cost = c;
	}

	const bool operator<(const VertexPair &other) const {
		if (v < other.v) return true;
		if (v > other.v) return false;
		return w<other.w;
	}

	const bool operator==(const VertexPair &other) const {
		return (v==other.v && w==other.w);
	}

};

#endif

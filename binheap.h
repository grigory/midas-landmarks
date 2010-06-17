// (c) Microsoft Corporation. All rights reserved.

#ifndef BINHEAP_H_
#define BINHEAP_H_


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

/// <summary>
/// Binary heap holding integral keys between 0 and maxid, with
/// values of type T. Lower values have higher priority. T must 
/// define comparision operators (at least "<").
/// </summary>

template <class T, bool DEBUGMODE = false> class BinaryHeap {

	private:
		/// <summary>
		/// Basic heap element structure supporting generic types
		/// </summary
		typedef struct {
			unsigned int label; //key of the element being stored
			T value; //its value
		} HeapUnit;

		// Heap elements
		unsigned int maxid;
		unsigned int lastpos;  //last position in which insertion is allowed
		HeapUnit *h;  //heap elements (elements will be stored in positions h[1] to h[lastpos])
		unsigned int nextpos;  //position of next element to be inserted
		unsigned int *heappos; //position of each element in the heap (0 is not yet inserted)

		/// <summary>
		/// Swaps elements in positions a and b.
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		inline void Swap(int a, int b) {
			HeapUnit tmp = h[a];
			h[a] = h[b];
			h[b] = tmp;
			heappos[h[a].label] = a;
			heappos[h[b].label] = b;
		}

		/// <summary>
		/// Standard 'up' operation performed on the element in the n-th position
		/// (the element is promoted until it's worse than its parent or reaches the top)
		/// </summary>
		/// <param name="n">Position.</param>
		inline void Up(int n) {
			while ((n>1) && IsBetter(GetValue(n), GetValue(n/2))) {
				Swap(n, n/2);
				n = n/2;
			}
		}

		inline T GetValue(int pos) const {return h[pos].value;}
		inline int GetLabel(int pos) const {return (h[pos].label);}

		/// <summary>
		/// Move the element currently in position n down to an appropriate position.
		/// </summary>
		/// <param name="n">Position.</param>
		void Down(int n) {
			HeapUnit tmp;
			T minchild;
			unsigned int two_n;

			if ((two_n = 2*n) < nextpos) {
				tmp = h[n]; //Save the original value
				do {
					bool first = true;
					minchild = GetValue(two_n); //guess first child is the best
					if ((two_n+1 < nextpos) && IsBetter(GetValue(two_n+1), minchild)) { //unless the second is better
						minchild = GetValue(two_n + 1);
						first = false;
					}
					if (IsBetter(minchild, tmp.value)) { //is the child better than the current value?
						if (first) { //minchild was the first
							h[n] = h[two_n];
							heappos[h[n].label] = n;
							n = two_n;
						} else { //minchild is the second
							h[n] = h[two_n + 1];
							heappos[h[n].label] = n;
							n = two_n + 1;
						}
						two_n = 2 * n;
					}
					else {
						break;
					}
				} while (two_n < nextpos);
				h[n] = tmp;
				heappos[h[n].label] = n;
			}
		}


		/// <summary>
		/// Find out if 'a' is strictly better than 'b'
		/// </summary>
		/// <param name="a">first value</param>
		/// <param name="b">second value</param>
		/// <returns>true iff a is strictly better (smaller) than b</returns>
		inline bool IsBetter(T a, T b) const {return (a<b);}

		inline void CheckRange(unsigned int v) {
			if (DEBUGMODE && v>maxid) {
				fprintf (stderr, "BinaryHeap: element %u is out of range (0,%u).\n", v, maxid);
				exit(-1);
			}
		}

	public:
		/// <summary>
		/// Creates a new heap that can hold elements from 0 to _maxid.
		/// </summary>
		/// <param name="_maxid">Maximum id (minimum is 0).</param>
		BinaryHeap(unsigned int _maxid) {
			maxid = _maxid;
			nextpos = 1; //first insertion will be in position 1
			lastpos = maxid+1; //There can be up to maxid+1 elements (0...maxid)

			h = new HeapUnit[lastpos+1]; //h[i]: i-th element in the heap (starting at h[1]).
			heappos = new unsigned int [maxid+1]; //h[i]: position of element i (between 0...maxid)
			for (unsigned int v=0; v<=maxid; v++) heappos[v] = 0;
		}

		/// <summary>
		/// Destructor
		/// </summary>
		virtual ~BinaryHeap() {
			delete [] h;
			delete [] heappos;
		};


		/// <summary>
		/// Remove the element with the highest priority (lowest value) from the heap. 
		/// Returns its label and value are returned.
		/// </summary>
		/// <param name="label"></param>
		/// <param name="value"></param>
		void RemoveFirst(unsigned &label, T &value) {
			label = h[1].label;  
			value = GetValue(1);

			h[1] = h[nextpos-1];	
			heappos[h[1].label] = 1;
			
			heappos[label] = 0;

			nextpos --;					
			Down(1);					
		}

		T MinValue() const {return GetValue(1);}

		/// <summary>
		/// Moves some elements in the heap to ensure head order
		/// </summary>
		void Heapify() {
			for (int i=(nextpos-1)/2; i>=1; i--) {Down(i);}
		}


		/// <summary>
		/// Removes the element with the specified label from the heap
		/// </summary>
		/// <param name="label"></param>
		void RemoveElement(unsigned label) {
			CheckRange(label);
			//get original position, check if element is actually in the heap
			int pos = heappos[label];
			if (pos!=0) {
				heappos[label] = 0; //mark element as no longer in the heap
				nextpos--;  // one fewer element

				// If the heap is not empty and this is not the last element, update the heap
				if ((nextpos != 1) && (nextpos != pos)) { 
					h[pos] = h[nextpos];  // Last element replaces the removed one
					heappos[h[pos].label] = pos;  // Let the element know its new position

					//unless parent is strictly better, go up
					if ((pos>1) && !IsBetter(h[pos/2].value, h[pos].value)) {Up(pos);}
					else {Down(pos);}
				}
			}
		}


		/// <summary>
		/// Insert and/or decrease the value of an element in the heap. Returns
		/// true iff successful (false if the new value is worse).
		/// </summary>
		/// <param name="label">Label of the new element.</param>
		/// <param name="value">Value of the new element.</param>
		/// <returns>True iff the element is inserted or updated.</returns>
		bool Insert(unsigned label, T value) {
			CheckRange(label);
			int prevpos = heappos[label];
			if (prevpos != 0) {  // already in the heap
				if (IsBetter(h[prevpos].value, value)) { //do not allow demotions
					return false;
				} else { //positive updates
					h[prevpos].value = value;
					Up(prevpos);
					return true;
				}
			}

			//not in the heap: insert it
			h[nextpos].label = label;
			h[nextpos].value = value;
			heappos[label] = nextpos;
			Up(nextpos);
			nextpos++;
			return true;
		}

		/// <summary>
		/// Change the key of element 'label' to value. If the element is
		/// not already in the heap, insert it; otherwise, just update it.
		/// The key will be changed regardless of whether the new key has
		/// higher or lower priority.
		/// </summary>
		/// <param name="label"></param>
		/// <param name="value"></param>
		void FixKey (unsigned label, T value) {
			CheckRange(label);
			int current_pos = heappos[label];
			if (current_pos == 0) { //element is not there: insert it
				h[nextpos].label = label;
				h[nextpos].value = value;
				heappos[label] = nextpos;
				Up(nextpos);
				nextpos++;
			} else { //element already in the heap: move it up or down as needed
				T old_value = h[current_pos].value;
				h[current_pos].value = value;

				if (IsBetter(old_value, value)) Down(current_pos); //element gets worse
				else Up(current_pos); //element gets better
			}
		}

		/// <summary>
		/// Add (label,value) as the last element of the current heap;
		/// MAY VIOLATE HEAP ORDER. The heap will be inconsistent
		/// until 'heapify' is called. This operation takes constant time.
		/// </summary>
		/// <param name="label"></param>
		/// <param name="value"></param>
		/// <returns></returns>
		bool PushBack(unsigned label, T value) {
			CheckRange(label);
			if (Contains(label)) return false;
			h[nextpos].label = label;
			h[nextpos].value = value;
			heappos[label] = nextpos;
			nextpos++;
			return true;
		}

		inline T GetElementValue(int v) { CheckRange(v); return (GetValue(heappos[v])); }
		inline bool Contains(int v) { CheckRange(v); return (heappos[v] != 0); }
		inline bool IsEmpty() const {return (nextpos == 1);}
		inline int GetSize() const {return (nextpos-1);}
		inline int GetMaxSize() const {return (lastpos);}

		/// <summary>
		/// Remove all elements from the heap.
		/// (Time proportional to the current number of elements.)
		/// </summary>
		void Reset() {
			bool full = false;
			if (full) {
				nextpos = 1;
				for (unsigned int i = 0; i <= lastpos; i++) {heappos[i] = 0;}
			} else {
				for (unsigned int p = 1; p < nextpos; p++) heappos[h[p].label] = 0;
				nextpos = 1;
			}
		}




		/// <summary>
		/// Sanity checks the heap invariants and internal mappings.
		/// DEBUG FUNCTION ONLY.
		/// </summary>
		/// <returns>true iff the data structure passes all tests</returns>
		bool Check() {
			fprintf (stderr, "Checking heap (%d elements)... ", nextpos);
			for (int i = nextpos - 1; i > 1; i--) {
				if (IsBetter(h[i].value, h[i/2].value)) {
					cerr << h[i].value << " (" << i << ") is child of " << h[i/2].value << " (" << i/2 << "). There are " << nextpos-1 << "elements.\n";
					return false;
				}
			}

			for (int i = 1; i < nextpos; i++) {
				if (heappos[h[i].label] != i) {
					cerr << "Element " << h[i].label << " things it's in position " << i << ", but it isn't.\n";
					return false;
				}
			}

			for (int i = 1; i <= lastpos; i++) {
				int pos = heappos[i];
				if ((pos > 0) && (pos <= lastpos)) {
					if (pos >= nextpos) {
						cerr << "Element " << i << "'s reported position (" << pos << ") is out of bounds (there are only " << nextpos-1 << " elements).\n";
						return false;
					}
					if (h[pos].label != i) {
						cerr << "Element " << i << " is not in the position it things it is (" << pos << "); " << h[pos].label << " is there.\n";
						return false;
					}
				}
			}

			fprintf (stderr, "passed.\n");
			return true;
		}
};

#endif

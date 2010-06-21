Compilation:
	make

Using mapgen:

1) Download file with coordinates in DIMACS format (for example, named "map.co").
2) Download file with graph description in DIMACS format (for example, named "map.gr"). 
3) Increase the stack size, running ulimit like this:
	ulimit -s 512000
4) Choose the necessary number of vertices in the resulting graph and execute mapgen:
	./mapgen FILE_WITH_COORDINATES FILE_WITH_GRAPH_DESCRIPTION NUMBER_OF_VERTICES

Note:
0) Don't forget to increase stack size! 
1) The output file with graph in DIMACS format is written to stdout.
2) The program will not generate a subgraph with exactly NUMBER_OF_VERTICES vertices, there will be slightly less vertices.

Example:
	./mapgen map.co map.gr 2000000 > map.out
	



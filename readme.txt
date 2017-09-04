
luby 1.0

This code reads in a given text file list of matrix entries (corresponding to an edge list of the underlying directed graph) and uses a version of Luby's algorithm to produce a matching for the underlying undirected network.  Edges are sorted by 1/(number of neighboring edges+1) so that leaf edges are given priority.  The resulting matching is then used to form a quotient graph with matched nodes identified.  

The program repeats this process recursively until a given condition is satisfied.  Currently it is set to loop until the returned matching has size smaller than 100.  These parameters are set in main.cc.

The end result is a partition of the original network into communities.  The program outputs a text file with lines "n c" where the n are the original nodes of the network and c is the community that n belongs to.

Other possible stopping criteria are (for positive integer k):

“match.size > k” (this says that each application of luby's algorithm joins at least k pairs of nodes.)
“proj.cols > k”  (this says that the graph has been partitioned into k communities)

 
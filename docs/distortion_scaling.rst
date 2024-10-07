Distortion Scaling
==================


When comparing the shapes of molecules, we are generally more interested in atoms and substructures which are (geometrically) nearby one another as those are more likely to interact chemically. The GW formula overweights residues which are further apart. To counteract this, given a metric (measure) space X', we can replace it with a related metric (measure) space X' with the same set of points, in a way that increases the relative distance between nearby point.
	- If we then apply GW to modified spaces X' and Y', the resulting GW distance will better capture pairs of points in X or Y which are nearby.


	- We call a function f a scaling if: 
		- f(0) = 0
		- f is strictly monotonic increasing
		- f is concave down
	- Theorem - if f is a scaling function, then given a metric space X, X' defined as the same points as X and d_X'(x1,x2) = f(d_X(x1,x2)) is also a metric space and is homeomorphic to X
	- Theorem - let GW_f(X,Y) be GW(X',Y'), then GW_f defines a metric on isomorphism classes of metric measure spaces

	- We find that this scaling generally improve the performance of GW. In our tests the square root function usually works best and has the advantage of not requiring any user-determined parameters.
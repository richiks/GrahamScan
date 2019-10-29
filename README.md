# GrahamScan

 * An algorithm for finding the convex hull of a set of points in the 2D plane
 * using the Graham Scan.  This algorithm has very good practical runtime
 * (O(n lg n)) and is fairly simple to understand.  The intuition behind the
 * algorithm is to find some point that is known to be on the convex hull of
 * the list of points, then to compute the angles between that point and every
 * other point in the set.  From there, we can sort the points in O(n lg n)
 * time.  The algorithm concludes by marching around these points in sorted
 * order, checking whether each is on the convex hull and adding it if it is.
 *
 * More specifically, the algorithm begins by picking the point with the
 * smallest Y coordinate, which must be on the convex hull.  It then computes
 * the angle made by each other point in the set and the X axis, then sorts
 * these points in ascending order.  Next, we maintain a stack of our current
 * guess of what the convex hull is, which initially is the point with the
 * lowest Y value and the point with the lowest (signed) angle with the x
 * axis.  From there, we iterate across the points in the convex hull
 * expanding our guess.  In particular, let v0 and v1 be the last two points
 * in our convex hull estimate, and let v2 be the next test point.  We then
 * compute the angle between the vector v1 - v0 and the vector v2 - v1.  If
 * this angle is in [pi, 2pi), then we are in a situation like this:
 *
 *                                             /
 *                                v1          /
 *                                *----------*
 *                               /          v0
 *                              /
 *                             /
 *                            *
 *                           v2
 *
 * This means that the point v1 is not actually on the convex hull; it's
 * contained in the hull between v0 and v2.  Of course, depending on what the
 * point in our convex hull estimate is that comes before v0, it's possible
 * that v0 itself isn't on the convex hull either.  In this case, we continue
 * removing nodes from our convex hull estimate until we find that the angle
 * between the last two nodes in the estimate and the new node is in the
 * range [0, pi).  We then update the convex hull by appending this new
 * point.
 *
 * We can now argue the correctness of the algorithm, along with the O(n lg n)
 * runtime.  We can easily see that the algorithm produces a convex polygon,
 * since if we follow the edges from the starting vertex v0 around, each edge
 * turns inward toward the center of the polygon. (A more rigorous proof of
 * this fact exists, but I think this previous one is more intuitive).  To see
 * that the algorithm produces a convex polygon containing all the points in
 * the initial input set, suppose for the sake of contradiction that this
 * isn't true; that there is some node v that isn't in the hull.  It can't be
 * lower than v0, since v0 has the lowest y coordinate, and so it must have
 * been considered by the algorithm at some point.  This means that it was
 * added to the convex hull candidate, and because it wasn't returned it must
 * have been removed at some point.  This means that there was some new point
 * in which the angle between the edge ending at v and the edge containing the
 * new point was in [0, pi).  But by removing this node from consideration,
 * the new edge added between the predecessor of v and the new point has point
 * v in its negative half-space, and so v is in the hull, a contradiction.
 *
 * To argue that the runtime is O(n lg n), we note that we can find the point
 * with the smallest y coordinate in O(n) time, compute the angle each point
 * makes with the x axis in constant time per node (taking O(n) net time), and
 * can then sort them in O(n lg n).  If we can show that the step of growing
 * the hull takes O(n lg n) time, then the overall runtime will be O(n lg n).
 * Initially, it might not seem that this step of the algorithm runs in this
 * time.  Each time we add a node we may have to backtrack all the way through
 * the potentially O(n) nodes on the convex hull, and since we're considering
 * O(n) nodes, this takes O(n^2) time.  However, this analysis is not tight.
 * Notice that once we remove a node from the stack, no future backtracking
 * can ever visit this node again.  This means that in all of the backtracking
 * operations, we can remove at most O(n) nodes from the stack.  Since each
 * backtracking operation takes time linear in the number of nodes removed,
 * this means that the total runtime for all of the backtracking operation is
 * O(n), giving us an overall runtime of O(n lg n).

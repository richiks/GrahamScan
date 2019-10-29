/**
 * @headerfile grahamscan.h
 * @author: Richik Vivek Sen (rsen9@gatech.edu)
 * @date 10/28/2019
 * @brief Header file for grahamscan.cpp
 */

#ifndef GRAHAMSCAN_H
#define GRAHAMSCAN_H

#include "matrix.h" // For Vector, DotProduct, Norm, NormSquared, CrossProduct
#include <iterator>  // For distance
#include <algorithm> // For copy, min_element, sort
#include <vector>

/**
 * Function: GrahamScan(ForwardIterator begin, ForwardIterator end,
 *                      OutputIterator out);
 * Usage: GrahamScan(pts.begin(), pts.end(), back_inserter(convexHull);
 * ---------------------------------------------------------------------------
 * Given a range of iterators [begin, end) spanning a range of Vector<2>s that
 * encode points in space, produces the convex hull of those points using the
 * Graham Scan algorithm.  The points are stored in counter-clockwise order
 * around the convex hull, so the resulting hull is the intersections of the
 * positive half-spaces of all the edges.  The result is stored in the
 * sequence beginning with out, which is assumed to have enough space to hold
 * the resulting sequence.  The return value is an iterator one past the
 * location of the last value written.
 */
template <typename ForwardIterator, typename OutputIterator>
OutputIterator GrahamScan(ForwardIterator begin, ForwardIterator end,
                          OutputIterator out);

/* * * * * Implementation Below This Point * * * * */
namespace grahamscan_detail {
  /**
   * Function: CompareYCoordinates(const Vector<2>& lhs, const Vector<2>& rhs)
   * Usage: if (CompareYCoordinates(pt1, pt2)) { ... }
   * -------------------------------------------------------------------------
   * Compares two vectors by their y coordinates, returning whether the lhs
   * has a y coordinate strictly less than the rhs's y coordinate.  In the
   * event of a tie, they are then compared by their x coordinate.
   */
  bool CompareYCoordinates(const Vector<2>& lhs, const Vector<2>& rhs);

  /**
   * Functor: CompareByAngle(const Vector<2>& lhs, const Vector<2>& rhs);
   * Usage: if (CompareByAngle(lhs, rhs)) { ... }
   * -------------------------------------------------------------------------
   * A functor class that sorts points according to the angle that they make
   * with the X axis.  All comparisons are made assuming that the points are
   * in a coordinate system that is a translated to some new origin.
   */
  class CompareByAngle {
  public:
    /**
     * Constructor: CompareByAngle(const Vector<2>& origin)
     * Usage: std::sort(begin, end, CompareByAngle(origin))
     * -----------------------------------------------------------------------
     * Constructs a new comparator with the indicated point as the origin.
     */
    explicit CompareByAngle(const Vector<2>& origin);

    /**
     * bool operator() (const Vector<2>& lhs, const Vector<2>& rhs) const
     * Usage: myComp(pt1, pt2);
     * -----------------------------------------------------------------------
     * Compares lhs and rhs, returning whether lhs makes a smaller angle with
     * the origin than rhs.  If lhs and rhs make the same angle, the distance
     * from the origin to these points is used as a tiebreaker, with the
     * closer point winning the tiebreaker.
     */
    bool operator() (const Vector<2>& lhs, const Vector<2>& rhs) const;

  private:
    const Vector<2> origin;
  };
}

/* Actual implementation of the Graham scan. */
template <typename ForwardIterator, typename OutputIterator>
OutputIterator GrahamScan(ForwardIterator begin, ForwardIterator end,
                          OutputIterator out) {
  /* Grant access to all of the utility functions and classes from above. */
  using namespace grahamscan_detail;

  /* Edge cases - if the range has fewer than three elements, the convex hull
   * is just those points.
   */
  if (size_t(std::distance(begin, end)) < 3)
    return std::copy(begin, end, out);

  /* Locate the element with the smallest y value, breaking ties by choosing
   * coordinates as far to the right (-x) as possible.
   */
  ForwardIterator minY = std::min_element(begin, end, CompareYCoordinates);

  /* Get an iterator one step past minY; it's the start of the sequence of
   * values that come after it in the input.
   */
  ForwardIterator next = minY; ++next;

  /* We now need to sort the points by their angle with the X axis.  Because
   * we aren't allowed to rearrange the input sequence, we'll make a local
   * copy of the sequence, then will sort that.  We'll leave the lowest point
   * out of the copy so that we don't end up including it in the result.
   */
  std::vector< Vector<2> > points;
  points.insert(points.end(), begin, minY); // First portion of the points.
  points.insert(points.end(), next, end);   // Remainder of the points.

  /* Sort by angle with the X axis.  To avoid issues where two adjacent points
   * in the sequence have an 180 degree angle between them, break ties by
   * choosing the point closest to the bottommost point.
   */
  std::sort(points.begin(), points.end(), CompareByAngle(*minY));

  /* For simplicity, add the minimum point onto the end of the ordering.  This
   * allows us to confirm that the last point we add in the sweep is correct
   * without having to special-case it.
   */
  points.push_back(*minY);

  /* Now, start building up the list of the points in the convex hull.
   * Initially this is the lowest point and the point with the lowest angle,
   * which happens to be the first element of the sorted sequence.
   */
  std::vector< Vector<2> > result;
  result.push_back(*minY);
  result.push_back(points[0]);

  /* Now, continuously refine the convex hull until we end up coming back
   * around to the beginning.
   */
  for (size_t i = 1; i < points.size(); ++i) {
    /* Expand the convex hull by factoring in this next point.  This may
     * entail removing some of our previous points, but it always ends by
     * adding this new point.
     */
    while (true) {
      /* Compute two vectors - one consisting of the last two points of the
       * candidate hull, and one consisting of of the last point and the next
       * point in the list.
       */
      const Vector<2> last = result[result.size() - 1] - result[result.size() - 2];
      const Vector<2> curr = points[i] - result[result.size() - 1];

      /* Check whether the angle between these vectors is in the range [0, pi)
       * or [pi, 2*pi).  If it's in the first group, we can add it.  Otherwise
       * we need to remove the last point from the hull and try again.
       *
       * Rather than directly computing the angle between the two vectors, we
       * can instead compute the sine of the angle.  If it's between [0, pi)
       * this will be nonnegative, and if it's between [pi, 2*pi) this would
       * be negative.
       *
       * We can compute the sine of the angle between the vectors by using the
       * 2D cross-product:
       *
       *             |   1   1   1 |
       *   |A x B| = | A.x A.y   0 | = A.x B.y - A.y B.x = |A| |B| sin(theta)
       *             | B.x B.y   0 |
       *
       * Since |A| |B| >= 0, this quantity is positive iff sin(theta) is
       * positive.
       */
      if (last[0] * curr[1] - last[1] * curr[0] >= 0) break;

      /* If we're here, it means that this angle was negative and so our last
       * point isn't going to work.  Undo it.
       */
      result.pop_back();
    }

    /* Finally, add the point. */
    result.push_back(points[i]);
  }

  /* At the very end, we now have our convex hull, with the lowest point added
   * twice.  We'll get rid of this point, then return the hull we found.
   */
  result.pop_back();

  /* Move the hull into the output range, then return the endpoint. */
  return std::copy(result.begin(), result.end(), out);
}

#endif

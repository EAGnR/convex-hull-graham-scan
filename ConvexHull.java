import java.util.ArrayList;
import java.util.Stack;

public class ConvexHull
{
    private ConvexHull(){}

    /**
     * Uses Graham Scan algorithm to compute the convex hull of a given set of 
     * points, and returns the vertices of the convex hull in clockwise order.
     *
     * @param points Unordered finite set of points.
     * @return The vertices of the convex hull as a {@code SimplePolygon}.
     */
    public static SimplePolygon grahamScan(Point[] points)
    {
        heapSort(points);
        // TODO: Implement rest of Graham Scan, point array should be sorted
        // after this line.

        return null;
    }

    /**
     * Heap Sort is used to order the given finite set of points, based on the 
     * size of their angle from the x-axis, the points are sorted in a
     * counter-clockwise circular manner. This is needed for the Graham Scan
     * algorithm to iterate through the points in the correct order.
     * Note that the point array passed to this method will be modified.
     *
     * @param points Unordered finite set of points, the array will be sorted.
     */
    private static void heapSort(Point[] points) 
    { 
        int n = points.length; 
        
        // Find the point with the lowest y value,to use as a reference point.
        Point lowestPoint = points[0];
        for (Point point : points)
        {
            if (point.getY() < lowestPoint.getY())
            {
                lowestPoint = point;
            }
            else if (Math.abs(point.getY() - lowestPoint.getY()) < Globals.POINT_EPSILON)
            {
                // If two points have equal y values, then get the one with
                // the smallest x value.
                if (point.getX() < lowestPoint.getX())
                {
                    lowestPoint = point;
                }
            }
        }
  
        // Build heap (rearranges array)
        // i starts at last non-leaf node, heapfiying by sift-down technique.
        for (int i = n/2 - 1; i >= 0; i--) 
            heapify(points, n, i, lowestPoint); 
  
        // One by one extract an element from heap, and sort array.
        for (int i = n-1; i >= 0; i--) 
        { 
            // Move current root to end.
            Point swap = points[0]; 
            points[0] = points[i]; 
            points[i] = swap; 
  
            // Call max heapify on the reduced heap.
            heapify(points, i, 0, lowestPoint); 
        } 
    } 
  
    /**
     * Heapifies a subtree rooted with node index i, producing a max heap
     * through Floyd's method which utilizes the sift-down technique.
     * The building of the heap is done in an optimal manner, taking O(n) time
     * overall for all points.
     * Note that the point array passed to this method will be modified.
     *
     * @param points Finite set of points, the array will be heapified.
     * @param size Given size of the heap.
     * @param i The root index of the subtree.
     * @param p The lowest point from the set of points, used as reference point.
     */
    private static void heapify(Point[] points, int size, int i, Point p) 
    { 
        int largest = i; // Initialize largest as root of subtree.
        int left = 2*i + 1; // left child = 2*i + 1 
        int right = 2*i + 2; // right child = 2*i + 2

        // Instead of computing the angle between the x-axis and a given point,
        // we calculate the Cosine of the angle as it is monotic in [0,pi],
        // this is more efficient to compute.
        double largestCos = getCos(p, points[largest]);
        double leftCos = getCos(p, points[left]);
        double rightCos = getCos(p, points[right]);
  
        // If left child is larger than parent.
        if (left < size && leftCos * -1 > largestCos * -1)
            largest = left; 
  
        // If right child is larger than parent.
        if (right < size && rightCos * -1 > largestCos * -1) 
            largest = right; 
  
        // If largest is not root.
        if (largest != i) 
        { 
            Point swap = points[i]; 
            points[i] = points[largest]; 
            points[largest] = swap; 
  
            // Recursively heapify the affected sub-tree.
            heapify(points, size, largest, p); 
        } 
    }

    /**
     * Returns the Cosine of the angle formed by a vector pq, and a unit vector
     * in the direction of the x-axis. i.e., the angle between pq and the x-axis.
     * We calculate the Cosine using the dot product of the vectors
     *
     * @param p Reference point, should be lowest point from set of points.
     * @param q Given point to calculate the angle with.
     * @return The Cosine of the angle between vector pq and x-axis.
     */
    private static double getCos(Point p, Point q)
    {
        return (q.getX() - p.getX()) / p.distance(q);
    }
}
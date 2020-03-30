import java.util.ArrayList;
import java.util.Arrays;
import java.util.Stack;

public class ConvexHull
{
    private ConvexHull(){}

    /**
     * Uses Graham Scan algorithm to compute the convex hull of a given set of 
     * points, and returns the vertices of the convex hull in clockwise order.
     *
     * @param points Unordered finite set of points in ArrayList form.
     * @return The vertices of the convex hull as a {@code SimplePolygon}.
     */
    public static SimplePolygon grahamScan(ArrayList<Point> points)
    {
        return grahamScan((Point[]) points.toArray());
    }

    /**
     * Uses Graham Scan algorithm to compute the convex hull of a given set of 
     * points, and returns the vertices of the convex hull in clockwise order.
     * Note that the point array passed to this method will be modified.
     *
     * @param points Unordered finite set of points.
     * @return The vertices of the convex hull as a {@code SimplePolygon}.
     */
    public static SimplePolygon grahamScan(Point[] points)
    {
        // Find the point with the lowest y value,to use as a reference point.
        Point lowestPoint = points[0];
        int lowestPointIndex = 0;
        for (int i = 0; i < points.length; i++)
        {
            if (points[i].getY() < lowestPoint.getY())
            {
                lowestPoint = points[i];
                lowestPointIndex = i;
            }
            else if (Math.abs(points[i].getY() - lowestPoint.getY()) < Globals.POINT_EPSILON)
            {
                // If two points have equal y values, then get the one with
                // the smallest x value.
                if (points[i].getX() < lowestPoint.getX())
                {
                    lowestPoint = points[i];
                    lowestPointIndex = i;
                }
            }
        }

        // Swap the first element of the array with the lowest point
        swap(points, 0, lowestPointIndex);
        heapSort(points);

        int size = 1;

        // If there are multiple points that give the same angle, we only keep the 
        // one farthest from the lowest point and discard the rest. Since we are
        // possibly discarding points, we will refer to the number of points by
        // the new variable size, and no longer by the length of the array.
        for(int i = 1; i < points.length; i++)
        {
            // Our heapsort implementation ensures that the farthest point will always 
            // be at the end of a group of points with the same angle.
            while(i < points.length - 1 
                && Math.abs(orientation(lowestPoint, points[i], points[i + 1])) < Globals.POINT_EPSILON)
            {
                i++;
            }

            points[size] = points[i];
            size++;
        }

        Stack<Point> stack = new Stack<Point>();
        for(int i = 0; i < size; i++)
        {
            // Remove points on stack which make a right or clockwise turn.
            // We keep points which make a left or counter-clockwise turn,
            // as well as collinear points.
            while(stack.size() > 1 
            && orientation(nextToTop(stack), stack.peek(), points[i]) < 0.0)
            {
                stack.pop();
            }


            stack.push(points[i]);
        }

        SimplePolygon polygon = new SimplePolygon();

        // Create a simple polygon with the points of the convex hull.
        // When popping the points out of the stack, we get the vertices in
        // clockwise order.
        while(!stack.empty())
        {
            polygon.addVertex(stack.pop());
        }

        // Return the convex hull as a simple polygon.
        return polygon;
    }

    /**
     * Heap Sort is used to order the given finite set of points, based on the 
     * size of their angle from the x-axis, the points are sorted in a
     * counter-clockwise circular manner. This is needed for the Graham Scan
     * algorithm to iterate through the points in the correct order. Since the
     * lowest point is set as the first point, we will only be sorting the rest
     * of the points.
     * Note that the point array passed to this method will be modified.
     * 
     * @param points Unordered finite set of points, the array will be sorted.
     */
    private static void heapSort(Point[] points) 
    { 
        Point[] subPoints = Arrays.copyOfRange(points, 1, points.length);
        int n = subPoints.length;
  
        // Build heap (rearranges array)
        // i starts at last non-leaf node, heapfiying by sift-down technique.
        for (int i = n/2 - 1; i >= 0; i--) 
            heapify(subPoints, n, i, points[0]); 
  
        // One by one extract an element from heap, and sort array.
        for (int i = n - 1; i >= 0; i--) 
        { 
            // Move current root to end.
            swap(subPoints, 0, i);
            // Call max heapify on the reduced heap.
            heapify(subPoints, i, 0, points[0]); 
        }

        for (int i = 0; i < subPoints.length; i++)
        {
            points[i + 1] = subPoints[i];
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
        // By largest we refer to the points that are considered greater by the 
        // max heap, and thus should be higher on the binary heap.
        int largest = i; // Initialize largest as root of subtree.
        int left = 2*i + 1; // left child = 2*i + 1
        int right = 2*i + 2; // right child = 2*i + 2

        // Instead of computing the angle between the x-axis and a given point,
        // we calculate the Cosine of the angle as it is monotic in [0,pi],
        // this is more efficient to compute and can be used for sorting the points
        // just the same.
        double largestCos = getCos(p, points[largest]);
        
        if(left < size)
        {
            double leftCos = getCos(p, points[left]);

            // If the two points make the same angle, the largest point in the 
            // max heap will be the one with the greatest distance from the
            // reference point.
            if (Math.abs(leftCos - largestCos) < Globals.POINT_EPSILON)
            {
                if(p.distance(points[left]) > p.distance(points[largest]))
                    largest = left;
            }
            // Otherwise the largest point is the one with the greatest Cosine.
            else if (leftCos * -1  > largestCos * -1)
            {
                largest = left; 
                largestCos = getCos(p, points[largest]);
            }
        }

        if(right < size)
        {
            double rightCos = getCos(p, points[right]);

            // If the two points make the same angle, the largest point in the 
            // max heap will be the one with the greatest distance from the
            // reference point.
            if (Math.abs(rightCos - largestCos) < Globals.POINT_EPSILON)
            {
                if(p.distance(points[right]) > p.distance(points[largest]))
                    largest = right;
            }
            // Otherwise the largest point is the one with the greatest Cosine.
            else if (rightCos * -1 > largestCos * -1)
            {
                largest = right;
                largestCos = getCos(p, points[largest]);
            }
        }
        
        // If the largest point is not the root of the subtree.
        if (largest != i) 
        { 
            // Swap parent and left or right child.
            swap(points, i, largest);
  
            // Recursively heapify the affected sub-tree.
            heapify(points, size, largest, p); 
        } 
    }

    /**
     * Returns the Cosine of the angle formed by a vector pq, and a unit vector
     * in the direction of the x-axis. i.e., the angle between pq and the x-axis.
     * We calculate the Cosine using the dot product of the vectors.
     *
     * @param p Reference point, should be lowest point from set of points.
     * @param q Given point to calculate the angle with.
     * @return The Cosine of the angle between vector pq and x-axis.
     */
    private static double getCos(Point p, Point q)
    {
        if(Math.abs(q.getX() - p.getX()) < Globals.POINT_EPSILON) return 0;

        return (q.getX() - p.getX()) / p.distance(q);
    }

    /**
     * Determines whether three points p1, p2, and p3 make a counter-clockwise 
     * turn, a clockwise turn, or if the points are collinear. To do this we 
     * calculate the z-coordinate of the cross product of the vectors p1p2 and p1p3.
     * (for counter-clockwise numbered points).
     * @param p1 The first point.
     * @param p2 The second point.
     * @param p3 The third point.
     * @return A positive number if the points make a left turn (counter-clockwise), 
     * zero if the points are collinear, a negative number otherwise (clockwise).
     */
    private static double orientation(Point p1, Point p2, Point p3)
    {

        return (p2.getX() - p1.getX()) * (p3.getY() - p1.getY()) 
             - (p2.getY() - p1.getY()) * (p3.getX() - p1.getX());
    }
    
    /**
     * Returns the element right below the top element of the stack.
     * @param stack the stack of points
     * @return the element right below the top of the stack
     */
    private static Point nextToTop(Stack<Point> stack)
    {
        Point top = stack.pop();

        // If top was the only element in the stack.
        if(stack.empty()) 
        {
            stack.push(top);
            return null;
        }
        Point nextToTop = stack.peek();
        stack.push(top);
        return nextToTop;
    }

    /**
     * Swaps two elements of an array.
     * @param points the array
     * @param x the index of the first element to be swapped.
     * @param y the index of the second element to be swapped.
     */
    private static void swap(Point[] points, int x, int y)
    {
        Point swap = points[x];
        points[x] = points[y];
        points[y] = swap;
    }
}
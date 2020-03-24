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
        System.out.println(lowestPoint);

        //swap the first element of the array with the lowest point
        swap(points, 0, lowestPointIndex);

        heapSort(points, lowestPoint);

        for(int i = 0; i < points.length; i++)
            System.out.print(points[i] + " ");
        System.out.println();

        Stack<Point> stack = new Stack<Point>();

        for(Point point : points)
        {
            while(stack.size() > 1 
            && counterClockwise(nextToTop(stack), stack.peek(), point) < 0)
            {
                stack.pop();
            }
            stack.push(point);
        }

        SimplePolygon polygon = new SimplePolygon();

        //add all vertices from the stack into the polygon
        while(!stack.empty())
        {
            polygon.addVertex(stack.pop());
        }

        return polygon;
    }

    /**
     * Returns the element right below the top element of the stack.
     * @param stack the stack of points
     * @return the element right below the top of the stack
     */
    private static Point nextToTop(Stack<Point> stack)
    {
        Point top = stack.pop();

        //if top was the only element in the stack
        if(stack.empty()) return null;

        Point nextToTop = stack.peek();
        stack.push(top);
        return nextToTop;
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
    private static void heapSort(Point[] points, Point lowestPoint) 
    { 
        int n = points.length; 
  
        // Build heap (rearranges array)
        // i starts at last non-leaf node, heapfiying by sift-down technique.
        for (int i = n/2 - 1; i >= 0; i--) 
            heapify(points, n, i, lowestPoint); 
  
        // One by one extract an element from heap, and sort array.
        for (int i = n - 1; i >= 0; i--) 
        { 
            // Move current root to end.
            swap(points, 0, i);
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
        int left = 2*i + 1; // left child = 2*i 
        int right = 2*i + 2; // right child = 2*i + 1

        // Instead of computing the angle between the x-axis and a given point,
        // we calculate the Cosine of the angle as it is monotic in [0,pi],
        // this is more efficient to compute.
        double largestCos = getCos(p, points[largest]);
        
        //if left is a valid index in the heap (as a node may not have a left
        //child)
        if(left <= size)
        {
            double leftCos = getCos(p, points[left]);

            // If left child is larger than parent.
            if (leftCos * -1 > largestCos * -1)
                largest = left; 
        }

        //if right is a valid index in the heap (as a node may not have a right
        //child)
        if(right <= size)
        {
            double rightCos = getCos(p, points[right]);

            // If right child is larger than parent.
            if (rightCos * -1 > largestCos * -1) 
                largest = right; 
        }
        
        // If largest is not root.
        if (largest != i) 
        { 
            swap(points, i, largest);
  
            // Recursively heapify the affected sub-tree.
            heapify(points, size, largest, p); 
        } 
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

    /**
     * Determines whether three points p1, p2, and p3 make a counterclockwise 
     * turn, a clockwise turn, or if the points are collinear. To do this we 
     * calculate the z-coordinate of the vectors p1p2 and p1p3.
     * @param p1 The first point.
     * @param p2 The second point.
     * @param p3 The third point.
     * @return a positive number if the points make a left turn, zero if the 
     * points are collinear, a negative number otherwise
     */
    private static double counterClockwise(Point p1, Point p2, Point p3)
    {

        return (p2.getX() - p1.getX()) * (p3.getY() - p1.getY()) 
             - (p2.getY() - p1.getY()) * (p3.getX() - p1.getX());
    }
}
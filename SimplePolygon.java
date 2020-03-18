
import java.util.ArrayList;

public class SimplePolygon
{

    ArrayList<Point> vertexList; //list of vertices in clockwise order

    /**
     * Instantiates an empty SimplePolygon object.
     */
    public SimplePolygon()
    {
        vertexList = new ArrayList<>();
    }

    /**
     * Adds a vertex as the last one in the sequence of vertices.
     *
     * @param v the new vertex
     */
    public void addVertex(Point v)
    {
        vertexList.add(v);
    }

    /**
     * Adds (i.e. inserts) a vertex at the specified location.
     *
     * @param v the new vertex
     * @param location position in the sequence of vertices to add the new
     * vertex in; locations are 0-indexed
     * @throws IndexOutOfBoundsException if <code>location</code> is out of
     * range
     */
    public void addVertex(Point v, int location) throws IndexOutOfBoundsException
    {
        if (location < 0 || location >= vertexList.size())
            throw new IndexOutOfBoundsException();
        else
            vertexList.add(location, v);
    }

    /**
     * Computes the area of this polygon.
     *
     * @return area of this polygon
     */
    public double area()
    {
        // No such value if the list is empty..
        if (vertexList.isEmpty()) return Double.NaN;

        double sum1 = 0.0, sum2 = 0.0;
        Point vertI, vertIplus1;

        for (int i = 0; i < vertexList.size(); i++) 
        {
            vertI = vertexList.get(i);
            vertIplus1 = vertexList.get((i+1) % vertexList.size());

            sum1 += vertI.getX() * vertIplus1.getY();
            sum2 += vertIplus1.getX() * vertI.getY();
        }

        // This is using the Shoelace Formula to calculate the area.
        // Only works for simple polygons.
        double area = 0.5 * Math.abs(sum1 - sum2);

        return area;
    }

    /**
     * Returns vertex at the specified location.
     *
     * @param location position of vertex in the sequence of vertices; locations
     * are 0-indexed
     * @return vertex at <code>location</code>
     * @throws IndexOutOfBoundsException if <code>location</code> is out of
     * range
     */
    public Point getVertex(int location) throws IndexOutOfBoundsException
    {
        if (location < 0 || location >= vertexList.size())
            throw new IndexOutOfBoundsException();
        else
            return vertexList.get(location);
    }

    /**
     * Determines the (minimum) bounding rectangle of this polygon.
     *
     * @return bounding rectangle
     */
    public Rectangle boundingRectangle()
    {
        Point upperLeft = new Point(smallestX(), greatestY());
        Point lowerRight = new Point(greatestX(), smallestY());
        Rectangle boundingRect = new Rectangle(upperLeft, lowerRight);

        return boundingRect;
    }

    /**
     * Determines the greatest x.
     *
     * @return the greatest x
     */
    public double greatestX()
    {
        // No such value if the list is empty.
        if (vertexList.isEmpty()) return Double.NaN;

        double maxX = vertexList.get(0).getX();
        for(int i = 1; i < vertexList.size(); i++)
        {
            double xVal = vertexList.get(i).getX();

            if(xVal > maxX) maxX = xVal;
        }

        return maxX;
    }

    /**
     * Determines the greatest y.
     *
     * @return the greatest y
     */
    public double greatestY()
    {
        // No such value if the list is empty.
        if (vertexList.isEmpty()) return Double.NaN;

        double maxY = vertexList.get(0).getY();
        for(int i = 1; i < vertexList.size(); i++)
        {
            double yVal = vertexList.get(i).getY();

            if(yVal > maxY) maxY = yVal;
        }

        return maxY;
    }

    /**
     * Computes the height of this object
     *
     * @return height of this object.
     */
    public double height()
    {
        return Math.abs(greatestY() - smallestY());
    }

    /**
     * Determines if given point is on the given line segment, this serves
     * as a helper method for isPointOnPolygon().
     *
     * @param p given point
     * @param p1 first endpoint of given line segment.
     * @param p2 second endpoint of given line segment.
     * @return true if given point is on the given line segment, false otherwise.
     */
    private static boolean isPointOnLineSegment(Point p, Point p1, Point p2)
    {
        // True if point is equal to an endpoint.
        if (p.equals(p1) || p.equals(p2)) 
        {
            return true;
        }

        // (Vertical Line Subcase)
        if (Math.abs(p1.getX() - p2.getX()) < Globals.POINT_EPSILON)
        {
            // True if point is collinear and within the vertical line segment.
            return Math.abs(p.getX() - p1.getX()) < Globals.POINT_EPSILON
                && (p.getY() > p1.getY() && p.getY() < p2.getY() 
                    || p.getY() < p1.getY() && p.getY() > p2.getY());
        }
        // (Horizontal Line Subcase)
        else if (Math.abs(p1.getY() - p2.getY()) < Globals.POINT_EPSILON)
        {
            // True if point is collinear and within the horizontal line segment.
            return Math.abs(p.getY() - p1.getY()) < Globals.POINT_EPSILON
                && (p.getX() > p1.getX() && p.getX() < p2.getX() 
                    || p.getX() < p1.getX() && p.getX() > p2.getX());
        }
        else // (Slanted Line Subcase)
        {
            // True if point is on line, and within the line segment.
            if (p.getX() > p1.getX() && p.getX() < p2.getX() 
                || p.getX() < p1.getX() && p.getX() > p2.getX()
                || Math.abs(p.getX() - p1.getX()) < Globals.POINT_EPSILON
                || Math.abs(p.getX() - p2.getX()) < Globals.POINT_EPSILON)
            {
                // Check if point is on line using point slope formula.
                double slope = (p2.getY() - p1.getY()) / (p2.getX() - p1.getX());
                double expectedY = slope * (p.getX() - p2.getX()) + p2.getY();

                return Math.abs(p.getY() - expectedY) < Globals.POINT_EPSILON;
            }

            return false;
        }
    }

    /**
     * Determines if given point is on the boundary or the interior of this
     * simple polygon.
     *
     * @param p given point
     * @return true if given point is on the boundary or the interior of this
     * polygon, false otherwise
     */
    public boolean isPointOnPolygon(Point p)
    {
        if (vertexList.size() == 1) // (Point Case)
        {
            // True if both points are equal.
            return p.equals(vertexList.get(0));
        }
        else if (vertexList.size() == 2) // (Line Case)
        {
            // True if point is within line segment.
            return isPointOnLineSegment(p, vertexList.get(0), vertexList.get(1));
        }
        else // (General Polygon Case)
        {
            // Using horizontal right-moving ray from point towards polygon to 
            // check for intersections.
            int crossings = 0;
            Point p1, p2;

            for (int i = 0; i < vertexList.size(); i++)
            {
                p1 = vertexList.get(i);
                p2 = vertexList.get((i+1) % vertexList.size());

                // True if point is equal to any other point, or is on the 
                // boundary of an edge.
                if (isPointOnLineSegment(p, p1, p2)) 
                {
                    return true;
                }

                // The Crossing Number Method Algorithm ([Franklin, 2000], [O'Rourke, 1998])
                // Rules:
                // 1) An upward edge includes its starting endpoint, and excludes 
                // its final endpoint.
                // 2) A downward edge excludes its starting endpoint, and includes 
                // its final endpoint.
                // 3) Horizontal edges are excluded.
                // 4) The edge-ray intersection point must be strictly right of 
                // the point P.

                // By using these rules we can ensure that the crossing number of 
                // the ray will always be correct even in special cases where a 
                // vertex or a horizontal edge is crossed. Relying on the condition
                // that the polygon is simple and not complex, so that using 
                // in-out parity remains reliable, by the Jordan Curve Theorem.
                // 
                // Note that applying rules 1 and 2 leads to a point on the 
                // top boundary being determined as out, and applying rule 4 leads 
                // to a point on the right boundary being determined as out. 
                // So only a point being on the left and bottom boundaries would 
                // be determined as in, which is a standard convention in computer 
                // graphics. But we correct for this by checking for boundary cases 
                // beforehand, so the overall method should return true if a point 
                // lies on any of the boundaries of the polygon.

                    // Rule 1: Upward edge crossing.
                if (( (p.getY() > p1.getY() || Math.abs(p.getY() - p1.getY()) < Globals.POINT_EPSILON) 
                    && p.getY() < p2.getY() )
                    // Rule 2: Downward edge crossing.
                    || ( p.getY() < p1.getY() 
                    && (p.getY() > p2.getY() || Math.abs(p.getY() - p2.getY()) < Globals.POINT_EPSILON) ))
                    // Rule 3: We ignore horizontal edges.
                {
                    // Calculate the x-coordinate of the edge-ray crossing 
                    // intersection. Using manipulated point-slope formula:
                    // x = ((y - y_1) / m) + x_1
                    // Since a vertical edge will lead to division by zero,
                    // seeing that m = (y_2 - y_1) / (x_2 - x_1), then we get:
                    // x = ((y - y_1) * (x_2 - x_1) / (y_2 - y_1)) + x_1
                    // Which will prevent division by zero for vertical edges,
                    // and we don't have to worry about horizontal edges as we 
                    // ignore them by Rule 3.
                    double x = ((p.getY() - p1.getY()) * (p2.getX() - p1.getX()) 
                        / (p2.getY() - p1.getY())) + p1.getX();

                    // Rule 4: The edge-ray intersection must be to the right 
                    // of point.
                    if (p.getX() < x)
                    {
                        crossings++; // Valid crossings are counted.
                    }
                }
            }

            // We return true if the number of crossings is odd,
            // otherwise (if even) we return false.
            return crossings % 2 != 0;
        }
    }

    /**
     * Computes the perimeter of this polygon.
     *
     * @return perimeter of this polygon
     */
    public double perimeter()
    {
        double perimeter = 0.0;
        
        for (int i = 0; i < vertexList.size(); i++)
        {
            perimeter += vertexList.get(i).distance(vertexList.get((i + 1) % vertexList.size()));
        }

        return perimeter;
    }

    /**
     * Replaces vertex at the specified location.
     *
     * @param v new vertex
     * @param location position in the sequence of vertices of the vertex to be
     * replaced
     * @throws IndexOutOfBoundsException if <code>location</code> is out of
     * range
     */
    public void replaceVertex(Point v, int location) throws IndexOutOfBoundsException
    {
        if (location < 0 || location >= vertexList.size())
            throw new IndexOutOfBoundsException();
        else
            vertexList.set(location, v);
    }

    /**
     * Determines the smallest x.
     *
     * @return the smallest x
     */
    public double smallestX()
    {
        // No such value if the list is empty.
        if (vertexList.isEmpty()) return Double.NaN;

        double minX = vertexList.get(0).getX();
        for(int i = 1; i < vertexList.size(); i++)
        {
            double xVal = vertexList.get(i).getX();

            if(xVal < minX) minX = xVal;
        }

        return minX;
    }

    /**
     * Determines the smallest y.
     *
     * @return the smallest y
     */
    public double smallestY()
    {
        // No such value if the list is empty.
        if (vertexList.isEmpty()) return Double.NaN;

        double minY = vertexList.get(0).getY();
        for(int i = 1; i < vertexList.size(); i++)
        {
            double yVal = vertexList.get(i).getY();

            if(yVal < minY) minY = yVal;
        }

        return minY;
    }

    /**
     * Returns an array with the vertices of the simple polygon, listed in
     * clockwise order.
     *
     * @return array containing the vertices of this simple polygon
     */
    Point[] toArray()
    {
        Point[] list = new Point[vertexList.size()];

        for (int i = 0; i < vertexList.size(); i++)
        {
            list[i] = new Point(vertexList.get(i).getX(), vertexList.get(i).getY());
        }

        return list;
    }

    /**
     * Constructs a String description of this simple polygon.
     *
     * @return String containing vertices of this polygon.
     */
    public String toString()
    {
        String str = "";

        int length = vertexList.size();

        for (int i = 0; i < length - 1; i++)
        {
            str += vertexList.get(i) + ", ";
        }

        if (length != 0)
            return "[" + str + vertexList.get(length - 1) + "]";
        else
            return "[]";
    }

    /**
     * Computes the width of this object
     *
     * @return width of this object.
     */
    public double width()
    {
        return Math.abs(greatestX() - smallestX());
    }
}


public class Point
{

    private double x;
    private double y;

    public Point()
    {
        x = y = 0;
    }

    public Point(double x1, double y1)
    {
        x = x1;
        y = y1;
    }
    /**
     * Calculates distance between this point and point <code>p</code>.
     *
     * @param p Point object
     *
     * @return distance between this point and point p
     */
    public double distance(Point p)
    {
        return Math.sqrt(Math.pow(this.x - p.x, 2) + Math.pow(this.y - p.y, 2));
    }

    public double getX()
    {
        return x;
    }

    public void setX(double d)
    {
        x = d;
    }

    public double getY()
    {
        return y;
    }

    public void setY(double d)
    {
        y = d;
    }

    public String toString()
    {
        return "(" + x + ", " + y + ")";
    }

    @Override
    public boolean equals(Object other)
    {
    	if(!other.getClass().getName().equals("Point")) return false;
    	
    	Point otherPoint = (Point) other;
        return Math.abs(this.x - otherPoint.x) < Globals.POINT_EPSILON 
            && Math.abs(this.y - otherPoint.y) < Globals.POINT_EPSILON;
    }
}

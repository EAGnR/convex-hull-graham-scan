public class Rectangle
{

    //begin and end are opposite corners; no guaranties which one
    private Point begin;
    private Point end;

    public Rectangle()
    {
        begin = new Point(0, 0);  //begin-upper corner = (0,0)
        end = new Point(1, 1);    //end-lower corner = (1,1)   
    }

    public Rectangle(Point l, Point r)
    {
        begin = l;
        end = r;
    }
    
    public Point getBegin()
    {
        return begin;
    }

    public Point getEnd()
    {
        return end;
    }

    public void setBegin(Point l)
    {
        begin = l;
    }

    public void setEnd(Point r)
    {
        end = r;
    }

    public String toString()
    {
        return "[" + begin.toString() + ", " + end.toString() + "]";
    }
}
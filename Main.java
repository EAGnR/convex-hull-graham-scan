
public class Main
{

    public static void main(String[] args)
    {
        new Main();
    }

    public Main()
    {
        
        Point[] points1 = {new Point(0,3), new Point(1,1), new Point(2,2),
                        new Point(4,4), new Point(0,0), new Point(1,2), 
                        new Point(3,1), new Point(3,3)};

        Point[] points2 = {new Point(-7,8), new Point(-4,6), new Point(2,6), 
                        new Point(6,4), new Point(8,6), new Point(7,-2), 
                        new Point(4,-6), new Point(8,-7), new Point(0,0), 
                        new Point(3,-2), new Point(6,-10), new Point(0,-6), 
                        new Point(-9,-5), new Point(-8,-2), new Point(-8,0), 
                        new Point(-10,3), new Point(-2,2), new Point(-10,4)};

        SimplePolygon polygon1 = ConvexHull.grahamScan(points1);
        System.out.println("Convex Hull 1: " + polygon1);

        SimplePolygon polygon2 = ConvexHull.grahamScan(points2);
        System.out.println("Convex Hull 2: " + polygon2);

        /*
        SimplePolygon polygon1 = new SimplePolygon();
        SimplePolygon triangle1 = new SimplePolygon();
        SimplePolygon triangle2 = new SimplePolygon();
        SimplePolygon triangle3 = new SimplePolygon();
        SimplePolygon line1 = new SimplePolygon();
        SimplePolygon line2 = new SimplePolygon();
        SimplePolygon line3 = new SimplePolygon();
        SimplePolygon hexagon = new SimplePolygon();

        polygon1.addVertex(new Point(3,4));
        polygon1.addVertex(new Point(5,11));
        polygon1.addVertex(new Point(12,8));
        polygon1.addVertex(new Point(9,5));
        polygon1.addVertex(new Point(5,6));

        triangle1.addVertex(new Point(2,4));
        triangle1.addVertex(new Point(3,-8));
        triangle1.addVertex(new Point(1,2));

        triangle2.addVertex(new Point(2,1));
        triangle2.addVertex(new Point(4,5));
        triangle2.addVertex(new Point(7,8));

        triangle3.addVertex(new Point(3,4));
        triangle3.addVertex(new Point(3,11));
        triangle3.addVertex(new Point(9,4));

        line1.addVertex(new Point(1,1));
        line1.addVertex(new Point(3,3));

        line2.addVertex(new Point(1,3));
        line2.addVertex(new Point(5,3));

        line3.addVertex(new Point(3,1));
        line3.addVertex(new Point(3,5));

        hexagon.addVertex(new Point(1,5));
        hexagon.addVertex(new Point(3,3));
        hexagon.addVertex(new Point(4,5));
        hexagon.addVertex(new Point(5,2));
        hexagon.addVertex(new Point(4,0));
        hexagon.addVertex(new Point(1,0));
        hexagon.addVertex(new Point(-5,2));

        System.out.printf("Polygon 1 greatest {x: %f, y: %f}, smallest {x: %f, y: %f}%n", 
            polygon1.greatestX(), polygon1.greatestY(), polygon1.smallestX(), polygon1.smallestY());
        System.out.printf("Polygon 1 width: %f, height %f%n", polygon1.width(), polygon1.height());
        System.out.printf("Polygon 1 perimeter: %f%n", polygon1.perimeter());
        System.out.printf("Polygon 1 bounding rectangle: %s%n", polygon1.boundingRectangle());
        System.out.printf("Polygon 1 area: %f, expected: 30.0%n", polygon1.area());
        System.out.printf("Polygon 1: %s, isPointOnPolygon(4.0,5.0): %b, expected: true%n", polygon1, polygon1.isPointOnPolygon(new Point(4.0,5.0)));
        System.out.printf("Polygon 1: %s, isPointOnPolygon(9.0,6.0): %b, expected: true%n", polygon1, polygon1.isPointOnPolygon(new Point(9.0,6.0)));
        System.out.printf("Polygon 1: %s, isPointOnPolygon(2.0,6.0): %b, expected: false%n%n", polygon1, polygon1.isPointOnPolygon(new Point(2.0,6.0)));

        System.out.printf("Triangle 1 greatest {x: %f, y: %f}, smallest {x: %f, y: %f}%n", 
            triangle1.greatestX(), triangle1.greatestY(), triangle1.smallestX(), triangle1.smallestY());
        System.out.printf("Triangle 1 width: %f, height %f%n", triangle1.width(), triangle1.height());
        System.out.printf("Triangle 1 perimeter: %f%n", triangle1.perimeter());
        System.out.printf("Triangle 1 bounding rectangle: %s%n", triangle1.boundingRectangle());
        System.out.printf("Triangle 1 area: %f, expected: 7.0%n%n", triangle1.area());

        System.out.printf("Triangle 2 greatest {x: %f, y: %f}, smallest {x: %f, y: %f}%n", 
            triangle2.greatestX(), triangle2.greatestY(), triangle2.smallestX(), triangle2.smallestY());
        System.out.printf("Triangle 2 width: %f, height %f%n", triangle2.width(), triangle2.height());
        System.out.printf("Triangle 2 perimeter: %f%n", triangle2.perimeter());
        System.out.printf("Triangle 2 bounding rectangle: %s%n", triangle2.boundingRectangle());
        System.out.printf("Triangle 2 area: %f, expected: 3.0%n%n", triangle2.area());

        System.out.printf("Triangle 3 greatest {x: %f, y: %f}, smallest {x: %f, y: %f}%n", 
            triangle3.greatestX(), triangle3.greatestY(), triangle3.smallestX(), triangle3.smallestY());
        System.out.printf("Triangle 3 width: %f, height %f%n", triangle3.width(), triangle3.height());
        System.out.printf("Triangle 3 perimeter: %f%n", triangle3.perimeter());
        System.out.printf("Triangle 3 bounding rectangle: %s%n", triangle3.boundingRectangle());
        System.out.printf("Triangle 3 area: %f, expected: 21.0%n", triangle3.area());
        System.out.printf("Triangle 3: %s, isPointOnPolygon(4.0,5.0): %b, expected: true%n", triangle3, triangle3.isPointOnPolygon(new Point(4.0,5.0)));
        System.out.printf("Triangle 3: %s, isPointOnPolygon(9.0,6.0): %b, expected: false%n", triangle3, triangle3.isPointOnPolygon(new Point(9.0,6.0)));
        System.out.printf("Triangle 3: %s, isPointOnPolygon(2.0,6.0): %b, expected: false%n%n", triangle3, triangle3.isPointOnPolygon(new Point(2.0,6.0)));

        System.out.printf("Line 1new Point(slanted): %s, isPointOnPolygon(2.0,2.0): %b, expected: true%n", line1, line1.isPointOnPolygon(new Point(2.0,2.0)));
        System.out.printf("Line 1new Point(slanted): %s, isPointOnPolygon(3.01,3.01): %b, expected: false%n", line1, line1.isPointOnPolygon(new Point(3.01,3.01)));
        System.out.printf("Line 1new Point(slanted): %s, isPointOnPolygon(1.0,2.0): %b, expected: false%n%n", line1, line1.isPointOnPolygon(new Point(1.0,2.0)));
        
        System.out.printf("Line 2new Point(horizontal): %s, isPointOnPolygon(2.5,3.0): %b, expected: true%n", line2, line2.isPointOnPolygon(new Point(2.5,3.0)));
        System.out.printf("Line 2new Point(horizontal): %s, isPointOnPolygon(6.0,3.0): %b, expected: false%n", line2, line2.isPointOnPolygon(new Point(6.0,3.0)));
        System.out.printf("Line 2new Point(horizontal): %s, isPointOnPolygon(2.5,2.99): %b, expected: false%n%n", line2, line2.isPointOnPolygon(new Point(2.5,2.99)));

        System.out.printf("Line 3new Point(vertical): %s, isPointOnPolygon(3.0,2.5): %b, expected: true%n", line3, line3.isPointOnPolygon(new Point(3.0,2.5)));
        System.out.printf("Line 3new Point(vertical): %s, isPointOnPolygon(3.0,6.0): %b, expected: false%n", line3, line3.isPointOnPolygon(new Point(3.0,6.0)));
        System.out.printf("Line 3new Point(vertical): %s, isPointOnPolygon(2.99,2.5): %b, expected: false%n%n", line3, line3.isPointOnPolygon(new Point(2.99,2.5)));
   
        System.out.printf("Hexagon: %s,isPointOnPolygon(1.0,3.0): %b, expected: true%n", hexagon, hexagon.isPointOnPolygon(new Point(1.0,3.0)));
        System.out.printf("Hexagon: %s,isPointOnPolygon(5.0,3.0): %b, expected: false%n", hexagon, hexagon.isPointOnPolygon(new Point(5.0,3.0)));
        System.out.printf("Hexagon: %s, perimeter: %.3f, expected: 26.496%n", hexagon, hexagon.perimeter());
        System.out.printf("Hexagon: %s, area: %.3f, expected: 29.500%n", hexagon, hexagon.area());
        */
    }
}

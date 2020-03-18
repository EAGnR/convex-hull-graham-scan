/**
 * This class contains declarations that are used among several other classes.
 */
public class Globals
{
    // Max precision epsilon, for floating point number comparison.
    public static final double POINT_EPSILON = 1.0 / Double.MAX_VALUE;
    public static enum GeometricObjectType
    {
        LINE_SEGMENT, NO_GEOMETRIC_OBJECT, POINT
    }
}
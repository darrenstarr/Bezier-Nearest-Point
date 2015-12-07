namespace Bezier_Nearest_Point
{
    public class Point2
    {
        public double X { get; set; } = 0.0;
        public double Y { get; set; } = 0.0;

        public double SquaredLength
        {
            get
            {
                return (X * X) + (Y * Y);
            }
        }

        public static Point2 operator +(Point2 left, Point2 right)
        {
            return new Point2 { X = left.X + right.X, Y = left.Y + right.Y };
        }

        public static Point2 operator -(Point2 left, Point2 right)
        {
            return new Point2 { X = left.X - right.X, Y = left.Y - right.Y };
        }

        public static Point2 operator *(Point2 left, Point2 right)
        {
            return new Point2 { X = left.X * right.X, Y = left.Y * right.Y };
        }

        public Point2 Scale(double scale)
        {
            return new Point2 { X = X * scale, Y = Y * scale };
        }

        public double DotProduct(Point2 against)
        {
            var product = this * against;
            return product.X + product.Y;
        }

        public override string ToString()
        {
            return "(" + X.ToString() + ", " + Y.ToString() + ")";
        }
    }
}

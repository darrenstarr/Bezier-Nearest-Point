using System;
using System.Linq;

// Solving the Nearest Point-on-Curve Problem
// and
// A Bezier Curve-Based Root-Finder
// by Philip J.Schneider
// from "Graphics Gems", Academic Press, 1990
// 
// Quick and dirty port to C# and .NET by 
//  Darren R. Starr darren at nocturnal.no 

namespace Bezier_Nearest_Point
{
    public class Bezier
    {
        // Ugly return value for Curve function
        internal class CurveResult
        {
            public Point2[] Left = null;
            public Point2[] Right = null;
            public Point2 Point = null;
        }

        // Maximum depth for recursion
        private static int MaxDepth = 64;
        // Flatness control value 
        private static double Epsilon { get; } = Math.Pow(2, -MaxDepth);
        // Cubic Bezier curve
        private static int Degree = 3;
        // Degree of eqn to find roots of
        private static int WDegree = 5;

        ///
        /// Bezier : 
        ///	Evaluate a Bezier curve at a particular parameter value
        ///     Fill in control points for resulting sub-curves if "Left" and
        ///	"Right" are non-null.
        ///
        ///     - int degree;      Degree of bezier curve	
        ///     - Point2* V;       Control pts			
        ///     - double t;        Parameter value		
        private static CurveResult Curve(Point2[] V, int degree, double t)
        {
            Point2[,] Vtemp = new Point2[WDegree + 1, WDegree + 1];

            // Copy control points	
            for (var j = 0; j <= degree; j++)
                Vtemp[0, j] = V[j];

            // Triangle computation	
            for (var i = 1; i <= degree; i++)
            {
                for (var j = 0; j <= degree - i; j++)
                {
                    Vtemp[i, j] = new Point2
                    {
                        X = (1.0 - t) * Vtemp[i - 1, j].X + t * Vtemp[i - 1, j + 1].X,
                        Y = (1.0 - t) * Vtemp[i - 1, j].Y + t * Vtemp[i - 1, j + 1].Y
                    };
                }
            }

            var result = new CurveResult
            {
                Point = Vtemp[degree, 0],
                Left = new Point2[WDegree + 1],
                Right = new Point2[WDegree + 1]
            };

            for (var i = 0; i <= degree; i++)
            {
                result.Left[i] = Vtemp[i, 0];
                result.Right[i] = Vtemp[degree - i, i];
            }

            return result;
        }

        /// ConvertToBezierForm :
        ///		Given a point and a Bezier curve, generate a 5th-degree
        ///		Bezier-format equation whose solution finds the point on the
        ///     curve nearest the user-defined point.
        ///         - Point2 P;      The point to find t for	
        ///         - Point2* V;	 The control points		
        static Point2[] ConvertToBezierForm(Point2 P, Point2[] V)
        {
            // TODO : static 
            double[,] z = new double[3, 4]
            {	
                // Precomputed "z" for cubics
	            {1.0, 0.6, 0.3, 0.1},
                {0.4, 0.6, 0.6, 0.4},
                {0.1, 0.3, 0.6, 1.0},
            };

            // Determine the c's -- these are vectors created by subtracting
            // point P from each of the control points				
            var c = new Point2[Degree + 1];      // V(i)'s - P
            for (var i = 0; i <= Degree; i++)
                c[i] = V[i] - P;

            // Determine the d's -- these are vectors created by subtracting
            // each control point from the next					
            var d = new Point2[Degree];    // V(i+1) - V(i)
            for (var i = 0; i <= Degree - 1; i++)
                d[i] = (V[i + 1] - V[i]).Scale(3.0);

            // Create the c,d table -- this is a table of dot products of the 
            // c's and d's							
            var cdTable = new double[3, 4]; // Dot product of c, d
            for (var row = 0; row <= Degree - 1; row++)
            {
                for (var column = 0; column <= Degree; column++)
                    cdTable[row, column] = d[row].DotProduct(c[column]);
            }

            // Now, apply the z's to the dot products, on the skew diagonal
            // Also, set up the x-values, making these "points"		
            var w = new Point2[WDegree + 1];
            for (var i = 0; i <= WDegree; i++)
            {
                w[i] = new Point2
                {
                    X = Convert.ToDouble(i) / Convert.ToDouble(WDegree),
                    Y = 0
                };
            }

            var n = Degree;
            var m = Degree - 1;
            for (var k = 0; k <= n + m; k++)
            {
                var lb = Math.Max(0, k - m);
                var ub = Math.Min(k, n);
                for (var i = lb; i <= ub; i++)
                {
                    var j = k - i;
                    w[i + j].Y += cdTable[j, i] * z[j, i];
                }
            }

            return (w);
        }

        /// 
        /// CrossingCount :
        /// Count the number of times a Bezier control polygon 
        ///  crosses the 0-axis. This number is >= the number of roots.
        /// 
        ///    Point2* V;            Control pts of Bezier curve	
        ///    int degree;			Degreee of Bezier curve 	
        static int CrossingCount(Point2[] V, int degree)
        {
            int result = 0;    // Number of zero-crossings	

            //  Sign of coefficients	
            var sign = Math.Sign(V[0].Y);
            var old_sign = sign;

            for (var i = 1; i <= degree; i++)
            {
                sign = Math.Sign(V[i].Y);
                if (sign != old_sign)
                    result++;

                old_sign = sign;
            }

            return result;
        }

        /// ControlPolygonFlatEnough :
        ///   Check if the control polygon of a Bezier curve is flat enough
        ///    for recursive subdivision to bottom out.
        /// 
        /// Corrections by James Walker, jw@jwwalker.com, as follows:
        /// 
        /// There seem to be errors in the ControlPolygonFlatEnough function in the
        /// Graphics Gems book and the repository (NearestPoint.c). This function
        /// is briefly described on p. 413 of the text, and appears on pages 793-794.
        /// I see two main problems with it.
        /// 
        /// The idea is to find an upper bound for the error of approximating the x
        /// intercept of the Bezier curve by the x intercept of the line through the
        /// first and last control points. It is claimed on p. 413 that this error is
        /// bounded by half of the difference between the intercepts of the bounding
        /// box. I don't see why that should be true. The line joining the first and
        /// last control points can be on one side of the bounding box, and the actual
        /// curve can be near the opposite side, so the bound should be the difference
        /// of the bounding box intercepts, not half of it.
        /// 
        /// Second, we come to the implementation. The values distance[i] computed in
        /// the first loop are not actual distances, but squares of distances. I
        /// realize that minimizing or maximizing the squares is equivalent to
        /// minimizing or maximizing the distances.  But when the code claims that
        /// one of the sides of the bounding box has equation
        /// a * x + b * y + c + max_distance_above, where max_distance_above is one of
        /// those squared distances, that makes no sense to me.
        /// 
        /// I have appended my version of the function. If you apply my code to the
        /// cubic Bezier curve used to test NearestPoint.c,
        /// 
        /// static Point2 bezCurve[4] = {    /  A cubic Bezier curve    /
        ///   { 0.0, 0.0 },
        ///   { 1.0, 2.0 },
        ///   { 3.0, 3.0 },
        ///   { 4.0, 2.0 },
        /// };
        /// 
        /// my code computes left_intercept = -3.0 and right_intercept = 0.0, which you
        /// can verify by sketching a graph. The original code computes
        /// left_intercept = 0.0 and right_intercept = 0.9
        static bool ControlPolygonFlatEnough(Point2[] V, int degree)
        {
            // Derive the implicit equation for line connecting first 
            //  and last control points 
            var a = V[0].Y - V[degree].Y;
            var b = V[degree].X - V[0].X;
            var c = V[0].X * V[degree].Y - V[degree].X * V[0].Y;

            var max_distance_above = 0.0;
            var max_distance_below = 0.0;

            for (var i = 1; i < degree; i++)
            {
                var value = a * V[i].X + b * V[i].Y + c;

                if (value > max_distance_above)
                    max_distance_above = value;
                else if (value < max_distance_below)
                    max_distance_below = value;
            }

            //  Implicit equation for zero line 
            var a1 = 0.0;
            var b1 = 1.0;
            var c1 = 0.0;

            //  Implicit equation for "above" line 
            var a2 = a;
            var b2 = b;
            var c2 = c - max_distance_above;

            var det = a1 * b2 - a2 * b1;
            var dInv = 1.0 / det;

            var intercept_1 = (b1 * c2 - b2 * c1) * dInv;

            //  Implicit equation for "below" line 
            a2 = a;
            b2 = b;
            c2 = c - max_distance_below;

            det = a1 * b2 - a2 * b1;
            dInv = 1.0 / det;

            var intercept_2 = (b1 * c2 - b2 * c1) * dInv;

            // Compute intercepts of bounding box   
            var left_intercept = Math.Min(intercept_1, intercept_2);
            var right_intercept = Math.Max(intercept_1, intercept_2);

            //Precision of root
            var error = right_intercept - left_intercept;

            return error < Epsilon;
        }

        /// ComputeXIntercept :
        /// Compute intersection of chord from first control point to last
        ///  	with 0-axis.
        /// 
        /// NOTE: "T" and "Y" do not have to be computed, and there are many useless
        ///   * operations in the following(e.g. "0.0 - 0.0").
        /// Point2* V;          - Control points	
        /// int degree; 		- Degree of curve	
        static double ComputeXIntercept(Point2[] V, int degree)
        {
            var XLK = 1.0 - 0.0;
            var YLK = 0.0 - 0.0;
            var XNM = V[degree].X - V[0].X;
            var YNM = V[degree].Y - V[0].Y;
            var XMK = V[0].X - 0.0;
            var YMK = V[0].Y - 0.0;

            var det = XNM * YLK - YNM * XLK;
            var detInv = 1.0 / det;

            var S = (XNM * YMK - YNM * XMK) * detInv;
            //  T = (XLK*YMK - YLK*XMK) * detInv; */

            var X = 0.0 + XLK * S;
            //  Y = 0.0 + YLK * S; */

            return X;
        }


        /// FindRoots :
        ///	Given a 5th-degree equation in Bernstein-Bezier form, find
        ///	all of the roots in the interval [0, 1].  Return the number
        ///	of roots found.
        ///     - Point2* w;           The control points		
        ///     - int degree;          The degree of the polynomial	
        ///     - double* t;           RETURN candidate t-values	
        ///     - int depth;		      The depth of the recursion	
        private static double[] FindRoots(Point2[] w, int degree, int depth)
        {
            int crossings = CrossingCount(w, degree);

            if (crossings == 0) // No solutions here	
                return new double[0];

            if (crossings == 1) // Unique solution
            {
                // Stop recursion when the tree is deep enough	
                // if deep enough, return 1 solution at midpoint 	
                if (depth >= MaxDepth)
                {
                    return new double[] {
                        (w[0].X + w[WDegree].X) / 2.0
                    };
                }

                if (ControlPolygonFlatEnough(w, degree))
                {
                    return new double[] {
                        ComputeXIntercept(w, degree)
                    };
                }
            }

            // Otherwise, solve recursively after	
            // subdividing control polygon		
            var curve = Curve(w, degree, 0.5);

            var leftT = FindRoots(curve.Left, degree, depth + 1);
            var rightT = FindRoots(curve.Right, degree, depth + 1);

            // Solutions from kids
            var result = new double[leftT.Count() + rightT.Count()];
            leftT.CopyTo(result, 0);
            rightT.CopyTo(result, leftT.Count());

            return result;
        }

        /*
         *  NearestPointOnCurve :
         *  	Compute the parameter value of the point on a Bezier
         *		curve segment closest to some arbtitrary, user-input point.
         *		Return the point on the curve at that parameter value.
         *
         */
        public static Point2 NearestPointOnCurve(Point2 P, Point2[] V)
        {
            //  Convert problem to 5th-degree Bezier form
            var w = ConvertToBezierForm(P, V);

            // Find all possible roots of 5th-degree equation 
            var t_candidate = FindRoots(w, WDegree, 0);

            /* Compare distances of P to all candidates, and to t=0, and t=1 */

            double new_dist;

            // Check distance to beginning of curve, where t = 0	
            var v = P - V[0];
            var dist = v.SquaredLength;

            //Parameter value of closest pt
            var t = 0.0;

            /* Find distances for candidate points	*/
            for (var i = 0; i < t_candidate.Count(); i++)
            {
                var p = Curve(V, Degree, t_candidate[i]).Point;

                //p = Curve(V, Degree, t_candidate[i], null, null);
                new_dist = (P - p).SquaredLength;
                if (new_dist < dist)
                {
                    dist = new_dist;
                    t = t_candidate[i];
                }
            }

            // Finally, look at distance to end point, where t = 1.0
            new_dist = (P - V[Degree]).SquaredLength;
            if (new_dist < dist)
            {
                dist = new_dist;
                t = 1.0;
            }

            // Return the point on the curve at parameter value t 
            //System.Diagnostics.Debug.WriteLine("t : " + t.ToString("N4.12"));
            return Curve(V, Degree, t).Point;
        }
    }
}

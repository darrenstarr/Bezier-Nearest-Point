using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using Windows.Foundation;
using Windows.Foundation.Collections;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;
using Windows.UI.Xaml.Controls.Primitives;
using Windows.UI.Xaml.Data;
using Windows.UI.Xaml.Input;
using Windows.UI.Xaml.Media;
using Windows.UI.Xaml.Navigation;

// The Blank Page item template is documented at http://go.microsoft.com/fwlink/?LinkId=402352&clcid=0x409

namespace Bezier_Nearest_Point
{
    /// <summary>
    /// An empty page that can be used on its own or navigated to within a Frame.
    /// </summary>
    public sealed partial class MainPage : Page
    {
        Point2 [] bezCurve = new Point2[] {
                new Point2 { X = 20, Y = 20 },
                new Point2 { X = 870, Y = 530 },
                new Point2 { X = 920, Y = 900 },
                new Point2 { X = 1000, Y = 160 }
            };

        public MainPage()
        {
            this.InitializeComponent();
            testCurve();
        }

        //	Given a cubic Bezier curve (i.e., its control points), and some
        //	arbitrary point in the plane, find the point on the curve
        //	closest to that arbitrary point.
        public void testCurve()
        {
            // Some arbitrary point
            var arbPoint = new Point2 { X=3.5, Y=2.0 }; 
            
            //    /*  Find the closest point */
            var pointOnCurve = Bezier.NearestPointOnCurve(arbPoint, bezCurve);
            System.Diagnostics.Debug.WriteLine("pointOnCurve : " + pointOnCurve.ToString());
        }

        private void Canvas_PointerMoved(object sender, PointerRoutedEventArgs e)
        {
            Windows.UI.Xaml.Input.Pointer ptr = e.Pointer;
            Windows.UI.Input.PointerPoint pointer = e.GetCurrentPoint(Target);

            var position = pointer.Position;

            var bPoint = new Point2 { X = position.X, Y = position.Y };
            var pointOnCurve = Bezier.NearestPointOnCurve(bPoint, bezCurve);

            box.SetValue(Canvas.LeftProperty, pointOnCurve.X - box.Width / 2);
            box.SetValue(Canvas.TopProperty, pointOnCurve.Y - box.Height / 2);
            System.Diagnostics.Debug.WriteLine("pointOnCurve : " + pointOnCurve.ToString());
        }
    }
}

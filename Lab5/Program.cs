using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Lab5
{
    class Program
    {
        static void Main(string[] args)
        {

            double a = 1, b = 10;
            double eps = 1e-2;
            int points = 100;
            LeastSquaresAppr appr = new LeastSquaresAppr(a, b, eps,points);
            appr.Out();
        }
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Lab5
{
    class Integral
    {

        static int N, N2;
        static double IN, I2N;

        public delegate double Funct(double x,int k, int t);

        static void FindN1(double eps, double a, double b)
        {
            int n = (int)(1 / Math.Sqrt(eps)) + 1;
            N = n;
        }

        static double TrapeziumCalc(double a, double b,Funct function1,int n, int k,int t)
        {
            double Summ = 0, h = (double)(b-a)/n;
            for (int i = 1; i < n; i++)
            {
                Summ += function1(a + i * h,k,t); 
            }
            return (h * (function1(a,k,t) / 2 + function1(b,k,t) / 2 + Summ));
        }

        static void CompTrapeziumRule(double eps, double a, double b, Funct function1,int k,int t)
        {
            int n = N;
            double result = TrapeziumCalc(a, b, function1,n,k,t);
            IN = result; 
        }

        static void RefinedCalc(double eps, double a, double b, Funct function1,int k,int t)
        {
            int n = N << 2;
            double In = IN;
            double I2n = TrapeziumCalc(a, b, function1,n, k,t), del = Math.Abs(In - I2n);
            while (del > 3*eps)
            {
                In = I2n;
                n = n << 2;
                I2n = TrapeziumCalc(a, b, function1,n,k,t);
                del = Math.Abs((In - I2n)/I2n);
            }
            N2 = n;
            I2N = I2n;
        }

        static public double DefineIntegral(double eps, double a, double b,Funct function1, int k,int t)
        {
            FindN1(eps, a, b);
            CompTrapeziumRule(eps, a, b, function1,k,t);
            RefinedCalc(eps, a, b, function1,k,t);
            double defintegral = I2N;
            return defintegral;
        }

    }
}

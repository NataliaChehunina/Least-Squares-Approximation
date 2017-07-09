using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading;



namespace Lab5
{
    public class NoRootsException : System.ApplicationException
    {
        public NoRootsException() { }
        public NoRootsException(string message) : base(message) { }
    }

    class LeastSquaresAppr
    {
        double interval;
        static int points;
        double x1, x2, epsilon, epsilonInt;
        double[,] M;
        double[] B;
        int bufferSize;

        double[] Ak;
        List<double> ptList = new List<double>();

        public LeastSquaresAppr(double A, double B, double eps,int P)
        {
            x1 = A;
            x2 = B;
            epsilon = eps;
            interval = Math.Abs(x1-x2);
            points = P;
            epsilonInt = eps*1e-7;
            bufferSize = 0;
        }

        double F(double x, int k, int t)//function for interpolation
        {
            return (7.5 * Math.Log10(x) * Math.Sin(x));
        }

        double FL(double x,int k, int t)
        {
            return LegandrePoly(x, k) * F(x, k, t);
        }



        double SPhi(double x,int k,int t)
        {
            return LegandrePoly(x, k) * LegandrePoly(x, t);
        }

        double LegandrePoly(double x, int n)
        {
            double P = 0, P0 = 1, P1 = x;
            if (n == 0) { return P0; }
            if (n == 1) { return P1; }
            for (int i = 2; i <= n; i++)
            {
                P = ((2 * i - 1) * x * P1 - (i - 1) * P0) / i;
                P0 = P1;
                P1 = P;
            }
            return P;
        }

        void FillM(int len, Integral.Funct f1, Integral.Funct f2)
        {
            if (len <= bufferSize)
                return;
            double[,] Mb = new double[len, len];
            double[] Bb = new double[len];

            for (int i = 0; i < bufferSize; i++)
            {
                for (int j = 0; j < bufferSize; j++)
                {
                    Mb[i, j] = M[i, j];
                }
                Bb[i] = B[i];
            }
            for (int i = bufferSize; i < len; i++)
            {
                for (int j = 0; j < i+1; j++)
                {
                    Mb[i, j] = IntegralF(i, j, f1);
                    Mb[j, i] = Mb[i, j];
                }
                Bb[i] = IntegralF(i, i, f2);
            }
            B = new double[len];
            M = new double[len, len];
            bufferSize = len;
            Bb.CopyTo(B, 0);
            Array.Copy(Mb, M, Mb.LongLength);
        }

        void GElimination(int n, double[,] _Ac, double[] _Bc)//Ak
        {
            double[,] Ac = new double[n, n];
            double[] Bc = new double[n];
            Array.Copy(_Ac, Ac, _Ac.LongLength);
            _Bc.CopyTo(Bc, 0);
            double[] Buf = new double[n + 1];
            for (int k = 0; k < n; k++)
            {
                if (Ac[k, k] == 0) { throw new NoRootsException("There is no roots in this system of equations !"); }
                else
                {
                    double elem = Ac[k, k];
                    for (int j = 0; j < n; j++)
                    {
                        Ac[k, j] /= elem;
                        Buf[j] = Ac[k, j];
                    }
                    Bc[k] /= elem;
                    Buf[n] = Bc[k];
                    for (int i = k + 1; i < n; i++)
                    {
                        double coef = Ac[i, k];
                        for (int j = 0; j < n; j++)
                        {
                            Ac[i, j] -= Buf[j] * coef;
                        }
                        Bc[i] -= Buf[n] * coef;
                    }
                }
            }
            BackSub(n, Ac, Bc);
        }

        void BackSub(int n, double[,] Ac, double[] Bc)
        {
            double summ;
            Ak = new double [Bc.Length];
            for (int k = n - 1; k >= 0; k--)
            {
                summ = 0;
                for (int j = k + 1; j < n; j++)
                {
                    summ += Ac[k, j] * Bc[j];
                }
                Bc[k] = Bc[k] - summ;
            }
            Ak = Bc;
        }

        void PrintMatrix(int n, double[,] Ac, double[] Bc)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write("|{0,5}", Ac[i, j]);
                }
                Console.Write(" |{0,5}", Bc[i]);
                Console.WriteLine();
            }
        }

        void PrintRoots(int n, double[] Bc)
        {
            for (int i = 0; i < n; i++)
            {
                Console.Write("x{0} = {1,10}", i + 1, Ak[i]);
                Console.WriteLine();
            }
        }

        double IntegralF(int k,int t, Integral.Funct f)
        {
            return Integral.DefineIntegral(epsilonInt, x1, x2, f, k, t);
        }

        double LPolynomials(int len, double x)
        {
            double p = 0;
            FillM(len, SPhi, FL);
            GElimination(len, M, B);
            for (int i = 0; i < len; i++)
            {
                p += LegandrePoly(x,i)*Ak[i];
            }
            return p;
        }

        int FindN(double a, double b,double eps)
        {
            int n = 0;
            double dx = interval/points, delta = 1,tmp = 0,x=a,funct;
            while (delta > eps) 
            {
                n++;
                tmp = 0;
                x = a;
                ptList.Clear();
                delta = 0;
                for (int i = 0; i < points; i++)
                {                 
                    funct = LPolynomials(n, x);
                    ptList.Add(x);
                    ptList.Add(funct);
                    double foo = F(x, 1, 1);
                    tmp += (F(x, 1, 1) - funct) * (F(x, 1, 1) - funct);
                    x += dx;
                }
                delta = Math.Sqrt(tmp / (points + 1));
                Console.WriteLine("N = {0}    LSDeviation = {1}", n-1, delta);

            } 
            return n-1;
        }

        public void Out()
        {
            int N = FindN(x1, x2,epsilon);
            Console.WriteLine("The power of Legendre polynomials is {0} .", N);
            double[] P = new double[ptList.Count];
            ptList.CopyTo(P, 0);
            StreamWriter sw = new StreamWriter("Points.txt", false, Encoding.Default);
            for (int i = 0; i < points<<1; i += 2)
            {
                Console.WriteLine("{0,8}   {1}", P[i], P[i + 1]);
                sw.WriteLine("{0,8}   {1}", P[i], P[i + 1]);
            }
            sw.Close();
         }

    }
}

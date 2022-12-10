using System;
using DelaunatorSharp;
using Numpy;
using DecimalMath;


double[][] test_points = new double[][]
{
new double[] { 5,4, 3 },
new double[] { 1, 1.1, 2 },
new double[] { 1, 3, 4 },
new double[] { 1, 2, 5 }
};
double[][] test_point = new double[][]
{
new double[] { 2.5,2 },
};
double[][] test_points2 = new double[][]
{
new double[] { 1, 2, 3},
new double[] { 3, 4, 2},
new double[] { 5, 6, 4 },
new double[] { -2, 4, 5 },
new double[] { -1, 2, 6 }
};
double[][] test_point2 = new double[][]
{
new double[] { 0,2.5 },
};
double[][] bad_point = new double[][]
{
new double[] { -2.5, 0, 3},
new double[] { 0, 5, 2},
new double[] { 0.5, 6, 4 },
new double[] { -5, -5, 5 },
};
void search_triangles(double[][] points, out double[][,] triangles, out double[][] triangles_norm, out double[][] centroid_coords)
{
    // points - массив точек для триангуляции вида [ [x1,y1,z1],[x2,y2,z2],... ]
    IPoint[] ipoints = new IPoint[points.Length];// массив для триангуляции делоне
    for (int i = 0; i < points.Length; i++)
    {
        ipoints[i] = new Point(points[i][0], points[i][1]);
    }
    Delaunator dn = new Delaunator(ipoints);// библиотека для триангуляции, выбрасывает exception, если треугольники невозможно построить
    triangles = new double[dn.GetTriangles().ToArray().Length][,];// массив массивов , который содержит кординаты трех точек каждого треугольника(покординатно)
    triangles_norm = new double[dn.GetTriangles().ToArray().Length][];//массив норм для каждого треуольника A,B,C, а также D (ур. плоскоти в общем виде Ax+By+Cz+D=0)
    centroid_coords = new double[dn.GetTriangles().ToArray().Length][];// массив координат центроида каждого треугольника(x+y+z/3)
    for (int i = 0; i < dn.GetTriangles().ToArray().Length; i++)
    {
        triangles[i] = new double[,] { { points[dn.Triangles[0 + 3 * i]][0], points[dn.Triangles[0 + 3 * i]][1], points[dn.Triangles[0 + 3 * i]][2] }, { points[dn.Triangles[1 + 3 * i]][0], points[dn.Triangles[1 + 3 * i]][1], points[dn.Triangles[1 + 3 * i]][2] }, { points[dn.Triangles[2 + 3 * i]][0], points[dn.Triangles[2 + 3 * i]][1], points[dn.Triangles[2 + 3 * i]][2] } };
    }
    for (int i = 0; i < triangles.Length; i++)
    {
        var res = np.cross(np.array(new double[] { triangles[i][0, 0] - triangles[i][1, 0], triangles[i][0, 1] - triangles[i][1, 1], triangles[i][0, 2] - triangles[i][1, 2] }), np.array(new double[] { triangles[i][0, 0] - triangles[i][2, 0], triangles[i][0, 1] - triangles[i][2, 1], triangles[i][0, 2] - triangles[i][2, 2] }));
        var qq = ((double)res[0]);
        var qq1 = ((double)res[1]);
        var qq2 = ((double)res[2]);
        var qq3 = np.dot(np.array(new double[] { ((double)res[0]), ((double)res[1]), ((double)res[2]) }), np.array(new double[] { triangles[i][2, 0], triangles[i][2, 1], triangles[i][2, 2] }));
        triangles_norm[i] = new double[] {qq, qq1, qq2, (double)qq3 };
    }
    for (int i = 0; i < triangles.Length; i++)
    {
        double h = 0;
        double h1 = 0;
        for (int j = 0; j < 3; j++)
        {
            h += triangles[i][j, 0];
            h1 += triangles[i][j, 1];
        }
        centroid_coords[i] = new double[] { h / 3.0, h1 / 3.0 };
    }
}

void triangulation(double[][] points, double[][] point)
{
    // point - точка, которую мы хотим интерполировать(найти ее координату z) вида:[x,y]
    double[][,] tr; // массив массивов точек каждого треугольника вида: [ [ [x1,y1,z1],[x2,y2,z2],[x3,y3,z3] ],...], получаем по ссылке из search_triangles
    double[][] tr_norm; // массив норм для каждого треуольника вида: [ [A1,B1,C1,D1],[A2,B2,C2,D2],... ], получаем по ссылке из search_triangles
    double[][] centroid_coords; // массив кординат центроида для каждого треугольника вида: [ [x1,y1,z1],[x2,y2,z2],... ], получаем по ссылке из search_triangles
    try 
    {
        search_triangles(points, out tr, out tr_norm, out centroid_coords);
    }
    catch
    {
        Console.WriteLine("Для данных точек невозможно построить треугольники делоне(точки лежат на одной прямой или количество точек <3)");
        return;
    }
    search_triangles(points, out tr, out tr_norm,out centroid_coords);
    double a1 = 0;
    double a2 = 0;
    double a3 = 0;
    for (int i = 0; i < tr.Length; i++)
    {
        a1 = (tr[i][0, 0] - point[0][0]) * (tr[i][1, 1] - tr[i][0, 1]) - (tr[i][1, 0] - tr[i][0, 0]) * (tr[i][0, 1] - point[0][1]);
        a2 = (tr[i][1, 0] - point[0][0]) * (tr[i][2, 1] - tr[i][1, 1]) - (tr[i][2, 0] - tr[i][1, 0]) * (tr[i][1, 1] - point[0][1]);// проверка принадлежности точки point, к i-му треугольнику (включая грани треугольника) (https://www.cyberforum.ru/algorithms/thread144722.html)
        a3 = (tr[i][2, 0] - point[0][0]) * (tr[i][1, 1] - tr[i][0, 1]) - (tr[i][0, 0] - tr[i][2, 0]) * (tr[i][2, 1] - point[0][1]);

        if ((a1 >= 0 && a2 >= 0 && a3 >= 0) || (a1 <= 0 && a2 <= 0 && a3 <= 0))
        {
            var zzz = (tr_norm[i][3] - tr_norm[i][0] * point[0][0] - tr_norm[i][1] * point[0][1]) / tr_norm[i][2];
            Console.WriteLine($"точка x:({point[0][0]}) y:({point[0][1]}) лежит в треугольнике(включая грани) с вершинами:\n{np.array(new double[,] { { tr[i][0, 0], tr[i][0, 1], tr[i][0, 2] }, { tr[i][1, 0], tr[i][1, 1], tr[i][1, 2] }, { tr[i][2, 0], tr[i][2, 1], tr[i][2, 2] } })}");
            Console.WriteLine($"точка x:({point[0][0]}) y:({point[0][1]})\nинтерполирована к точке треугольника с координатами:\nx:({point[0][0]}) y:({point[0][1]}) и значением функции:({zzz})");
            return;
        }
    }
    var measure = Math.Sqrt(Math.Pow(point[0][0] - centroid_coords[0][0], 2) + Math.Pow(point[0][1] - centroid_coords[0][1], 2));// если, точка не лежит внутри какого-нибудь из имеющихся треугольников, ищем ближайший треугольник исходя из минимального расстояние до центроида треугольников
    int j1 = 0;
    for (int i = 1; i < tr.Length; i++)
    {
        var step = Math.Sqrt(Math.Pow(point[0][0] - centroid_coords[i][0], 2) + Math.Pow(point[0][1] - centroid_coords[i][1], 2));
        if (step < measure)
        {
            measure = step;
            j1 = i;
        }
    }
    var z = (tr_norm[j1][3] - tr_norm[j1][0] * point[0][0] - tr_norm[j1][1] * point[0][1]) / tr_norm[j1][2];
    Console.WriteLine($"точка x:({point[0][0]}) y:({point[0][1]}) лежит вне треугольников, ближайший треугольник с вершинами:\n{np.array(new double[,] { { tr[j1][0,0], tr[j1][0, 1], tr[j1][0, 2] }, { tr[j1][1, 0], tr[j1][1, 1], tr[j1][1, 2] }, {tr[j1][2, 0], tr[j1][2, 1], tr[j1][2, 2] } })}");
    Console.WriteLine($"точка x:({point[0][0]}) y:({point[0][1]})\nинтерполирована к точке треугольника с координатами:\nx:({point[0][0]}) y:({point[0][1]}) и значением функции:({z})");
    return;
}
triangulation(test_points, test_point);
triangulation(test_points2, test_point2);
triangulation(bad_point, test_point);
//search_triangles(bad_point);

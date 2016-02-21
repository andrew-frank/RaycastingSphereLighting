using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace SphereLighting
{
    public partial class MainWindow : Window
    {

        private static int s1X = 3;
        private static int s1Y = 8;
        private static int s1Z = 10;



        public struct Vector3Df
        {

            public float X;

            public float Y;

            public float Z;

            public float _X
            {
                get { return this.X; }
                set { this.X = value; }
            }
            public float _Y
            {
                get { return this.Y; }
                set { this.Y = value; }
            }
            public float _Z
            {
                get { return this.Z; }
                set { this.Z = value; }
            }


            #region --- Methods ---
            #region Constructors

            public Vector3Df(float x, float y, float z)
            {
                this.X = x;
                this.Y = y;
                this.Z = z;
            }

  

            public Vector3Df(ref Vector3Df vector)
            {
                this.X = vector.X;
                this.Y = vector.Y;
                this.Z = vector.Z;
            }

            public Vector3Df(ref float[] Data)
            {
                this.X = Data[0];
                this.Y = Data[1];
                this.Z = Data[2];
            }
            #endregion Constructors

            #region Magnitude
            public float MagSq()
            {
                return (this.X * this.X + this.Y * this.Y + this.Z * this.Z);
            }
            public float Mag()
            {
                return (float)System.Math.Sqrt((double)this.MagSq());
            }
            #endregion

            #region Vector Operations
            public Vector3Df UnitVector()
            {
                float tmp = 1 / this.Mag();
                Vector3Df result = this * tmp;
                return result;
            }
            public static float DotProduct(ref Vector3Df Vector1, ref Vector3Df Vector2)
            {
                float result = Vector1.X * Vector2.X +
                               Vector1.Y * Vector2.Y +
                               Vector1.Z * Vector2.Z;
                return result;
            }
            public static void CrossProduct(ref Vector3Df VectorLeft, ref Vector3Df VectorRight, out Vector3Df Result)
            {
                Result = new Vector3Df((VectorLeft.Y * VectorRight.Z - VectorLeft.Z * VectorRight.Y),
                                    (VectorLeft.Z * VectorRight.X - VectorLeft.X * VectorRight.Z),
                                    VectorLeft.X * VectorRight.Y - VectorLeft.Y * VectorRight.X);
            }
            public static Vector3Df CrossProduct(ref Vector3Df VectorLeft, ref Vector3Df VectorRight)
            {
                return new Vector3Df((VectorLeft.Y * VectorRight.Z - VectorLeft.Z * VectorRight.Y),
                                    (VectorLeft.Z * VectorRight.X - VectorLeft.X * VectorRight.Z),
                                    VectorLeft.X * VectorRight.Y - VectorLeft.Y * VectorRight.X);
            }
            #endregion

            #region Misc.
            public override string ToString()
            {
                string result = "X = " + this.X + "     Y = " + this.Y + "     Z = " + this.Z;
                return result;
            }
            public float[] ToArray3f()
            {
                float[] result = new float[] { this.X, this.Y, this.Z };
                return result;
            }
            #endregion
            #endregion

            #region --- Operator Overloads ---
            #region AGE_Vector3Df operator +(AGE_Vector3Df a, AGE_Vector3Df b)
            /// <summary>
            /// Adds two vectors.  Result is (a.X + b.X, a.Y + b.Y, a.Z + b.Z).
            /// </summary>
            /// <param name="a">First vector.</param>
            /// <param name="b">Second vector.</param>
            /// <returns></returns>
            public static Vector3Df operator +(Vector3Df a, Vector3Df b)
            {
                return new Vector3Df(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
            }
            #endregion AGE_Vector3Df operator +(AGE_Vector3Df a, AGE_Vector3Df b)

            #region AGE_Vector3Df operator *(AGE_Vector3Df vector, float scalar)
            /// <summary>
            /// Multiply vector by a scalar.  Result is (vector.X * scalar, vector.Y * scalar, vector.Z * scalar).
            /// </summary>
            /// <param name="vector">The vector.</param>
            /// <param name="scalar">The scalar</param>
            /// <returns></returns>
            public static Vector3Df operator *(Vector3Df vector, float scalar)
            {
                return new Vector3Df(vector.X * scalar, vector.Y * scalar, vector.Z * scalar);
            }
            #endregion AGE_Vector3Df operator *(AGE_Vector3Df vector, float scalar)

            #region AGE_Vector3Df operator -(AGE_Vector3Df a, AGE_Vector3Df b)
            /// <summary>
            /// Subtract two vectors.  Result is (a.X - b.X, a.Y - b.Y, a.Z - b.Z).
            /// </summary>
            /// <param name="a">The first vector.</param>
            /// <param name="b">The second vector.</param>
            /// <returns></returns>
            public static Vector3Df operator -(Vector3Df a, Vector3Df b)
            {
                return new Vector3Df(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
            }
            #endregion AGE_Vector3Df operator -(AGE_Vector3Df a, AGE_Vector3Df b)
            #endregion
        }




        private void test()
        {
            Sphere sampleSphere = new Sphere();
            sampleSphere.Radius = 10;
            sampleSphere.SpehereColor = Colors.Red;
            sampleSphere.R = 255;
            sampleSphere.G = 0;
            sampleSphere.B = 0;
            sampleSphere.Center = new SinglePoint();
            sampleSphere.Center.x = 5170;
            sampleSphere.Center.y = 60;
            sampleSphere.Center.z = 8;
            SinglePoint pointOfTheViewer = new SinglePoint();
            pointOfTheViewer.x = 100;
            pointOfTheViewer.y = 80;
            pointOfTheViewer.z = 1;

            int canvasWidth = 500;
            int canvasHeight = 300;
            int[,] intersections = new int[canvasWidth, canvasHeight];
            SphereInterscetionCheck(sampleSphere, pointOfTheViewer, canvasWidth, canvasHeight);
        }


        private void DrawPixel(int x, int y, Color color)
        {
            Ellipse el = new Ellipse();
            el.Width = 1;
            el.Height = 1;
            el.Fill = new SolidColorBrush(color);
            Canvas.SetLeft(el, x);
            Canvas.SetTop(el, y);
            //couranvas.Children.Add(el);
        }



        /*   A sphere is given by its center (cx, cy, cz), its radius R, and its color (SR, SG, SB).
         *   line segment (ray) is given by its endpoints: P0 = (x0, y0, z0) and P1 = (x1, y1, z1).
         *   To find visible spheres, set P0 = viewer’s coordinates, VP = (VPx, VPy, VPz) and let P1 run through
         *   all the points (x1, y1, 0) where (x1, y1) is a pixel in the display area.
         */

        private void SphereInterscetionCheck(Sphere sphere, SinglePoint viewerPoint, int width, int height)
        {
            double dx, dy, dz, a, b, c, delta;
            SinglePoint lightPoint = new SinglePoint();
            lightPoint.x = 10;
            lightPoint.y = 15;
            lightPoint.z = 5;

            int[,] resultOfIntersection = new int[width, height];
            for (int i = 0; i < width; i++) {
                for (int j = 0; j < height; j++) {
                    dx = i - viewerPoint.x;
                    dy = j - viewerPoint.y;
                    dz = 0 - viewerPoint.z;

                    a = dx * dx + dy * dy + dz * dz;
                    b = 2 * dx * (viewerPoint.z - sphere.Center.x) + 2 * dy * (viewerPoint.y - sphere.Center.y) + 2 * dz * (viewerPoint.z - sphere.Center.z);
                    c = sphere.Center.x * sphere.Center.x + sphere.Center.y * sphere.Center.y + sphere.Center.z * sphere.Center.z + viewerPoint.x * viewerPoint.x + viewerPoint.y * viewerPoint.y + viewerPoint.z * viewerPoint.z - 2 * (sphere.Center.x * viewerPoint.x + sphere.Center.y * viewerPoint.y + sphere.Center.z * viewerPoint.z) - sphere.Radius * sphere.Radius;

                    delta = b * b - 4 * a * c;
                    if (delta < 0) //Findinghadowws
                        DrawPixel(i, j, Colors.Black);
                    else if (delta == 0)
                        DrawPixel(i, j, (DiffuseShading(sphere, viewerPoint, lightPoint, dx, dy, dz, a, b, c, width, height)));//resultOfIntersection[i, j] = 1;
                    else
                        DrawPixel(i, j, (DiffuseShading(sphere, viewerPoint, lightPoint, dx, dy, dz, a, b, c, width, height)));//resultOfIntersection[i, j] = 2;
                }
            }
        }
        private Color DiffuseShading(Sphere sphere, SinglePoint viewerPoint, SinglePoint lightPoint, double dx, double dy, double dz, double a, double b, double c, int width, int height)
        {
            double x, y, z, t;
            Vector3D normalVector = new Vector3D();
            Vector3D lightVector = new Vector3D();
            t = (-b - Math.Sqrt(b * b - 4 * a * c)) / (2 * a);

            x = viewerPoint.x + t * dx;
            y = viewerPoint.y + t * dy;
            z = viewerPoint.z + t * dz;

            //Unit normal vector
            normalVector.X = ((x - sphere.Center.x) / sphere.Radius);
            normalVector.Y = (y - sphere.Center.y) / sphere.Radius;
            normalVector.Z = (z - sphere.Center.z) / sphere.Radius;

            //Vector from normal to the Light
            lightVector.X = lightPoint.x - x;
            lightVector.Y = lightPoint.y - y;
            lightVector.Z = lightPoint.z - z;

            //constants
            double kd = 0.8;
            double ka = 0.2;

            double fctr = Vector3D.DotProduct(normalVector, lightVector);
            Color resultColor = new Color();
            resultColor.R = (byte)(ka * sphere.R + kd * fctr * sphere.R);
            resultColor.G = (byte)(ka * sphere.G + kd * fctr * sphere.G);
            resultColor.B = (byte)(ka * sphere.B + kd * fctr * sphere.B);

            return resultColor;
        }






        public struct CG_Matrix
        {
            #region Fields
            public float[][] col;
            public float Scale;
            public Vector3Df Pan;
            #endregion

            #region Methods

            public CG_Matrix(float dia)
            {
                col = new float[4][];
                col[0] = new float[4];
                col[1] = new float[4];
                col[2] = new float[4];
                col[3] = new float[4];
                for (int c = 0; c < 4; c++) {
                    for (int r = 0; r < 4; r++) {
                        col[c][r] = (c == r) ? dia : 0.0f;
                    }
                }
                this.Scale = dia;
                this.Pan = new Vector3Df();
            }



            public static CG_Matrix Identity()
            {
                return new CG_Matrix(1.0f);
            }
            public static CG_Matrix Zero()
            {
                return new CG_Matrix(0.0f);
            }
            public static CG_Matrix Intialize(float[] List)
            {
                if (List.Length != 16)
                    throw new ArgumentException("Initialization List != 16");
                CG_Matrix result = CG_Matrix.Zero();
                for (int i = 0; i < List.Length; i++) {
                    int c = i / 4;
                    int r = i - (4 * c);
                    result.col[c][r] = List[i];
                }
                return result;
            }
            public static CG_Matrix Intialize(double[] List)
            {
                if (List.Length != 16)
                    throw new ArgumentException("Initialization List != 16");
                CG_Matrix result = CG_Matrix.Zero();
                for (int i = 0; i < List.Length; i++) {
                    int c = i / 4;
                    int r = i - (4 * c);
                    result.col[c][r] = (float)List[i];
                }
                return result;
            }



            public static CG_Matrix PerspectiveProjection(float left,
                                            float right,
                                            float top,
                                            float bottom,
                                            float near,
                                            float far)
            {
                float invRDiff = 1 / (right - left + float.Epsilon);
                float invUDiff = 1 / (top - bottom + float.Epsilon);
                float invDDiff = 1 / (far - near + float.Epsilon);

                CG_Matrix result = CG_Matrix.Zero();
                result.col[0][0] = 2 * near * invRDiff;

                result.col[1][1] = 2 * near * invUDiff;

                result.col[2][0] = (right + left) * invRDiff;
                result.col[2][1] = (top + bottom) * invUDiff;
                result.col[2][2] = -1 * (far + near) * invDDiff;
                result.col[2][3] = -1;

                result.col[3][2] = -2 * (far * near) * invDDiff;

                return result;
            }

            public static CG_Matrix OrthographicProjection(float left, float right, float top, float bottom, float near, float far)
            {
                float invRDiff = 1 / (right - left);
                float invUDiff = 1 / (top - bottom);
                float invDDiff = 1 / (far - near);

                CG_Matrix result = CG_Matrix.Zero();
                result.col[0][0] = 2 * invRDiff;
                result.col[3][0] = -1 * (right + left) * invRDiff;
                result.col[1][1] = 2 * invUDiff;
                result.col[3][1] = -1 * (top + bottom) * invUDiff;
                result.col[2][2] = -2 * invDDiff;
                result.col[3][2] = -1 * (far + near) * invDDiff;
                result.col[3][3] = 1;
                return result;
            }




            public static CG_Matrix ViewSpace(ref Vector3Df Eye, ref Vector3Df Target, ref  Vector3Df UpVector)
            {
                CG_Matrix result = CG_Matrix.Zero();
                Vector3Df VDirection = (Target - Eye);
                VDirection = VDirection.UnitVector();


                Vector3Df RightDirection = Vector3Df.CrossProduct(ref VDirection, ref UpVector);
                RightDirection = RightDirection.UnitVector();
                Vector3Df UpDirection = Vector3Df.CrossProduct(ref RightDirection, ref VDirection);
                UpDirection = UpDirection.UnitVector();

                result.col[0][0] = RightDirection.X;
                result.col[1][0] = RightDirection.Y;
                result.col[2][0] = RightDirection.Z;

                result.col[0][1] = UpDirection.X;
                result.col[1][1] = UpDirection.Y;
                result.col[2][1] = UpDirection.Z;

                result.col[0][2] = -1 * VDirection.X;
                result.col[1][2] = -1 * VDirection.Y;
                result.col[2][2] = -1 * VDirection.Z;

                result.col[3][0] = -1 * Vector3Df.DotProduct(ref RightDirection, ref Eye);
                result.col[3][1] = -1 * Vector3Df.DotProduct(ref UpDirection, ref Eye);
                result.col[3][2] = Vector3Df.DotProduct(ref VDirection, ref Eye);

                result.col[3][3] = 1;

                return result;
            }



            public static CG_Matrix ScaleMatrix(ref Vector3Df Location, ref Vector3Df Scale)
            {
                CG_Matrix result = CG_Matrix.Identity();
                //Matrix_33 M = Matrix_33.Zero();
                result.col[0][0] = Scale.X;
                result.col[1][1] = Scale.Y;
                result.col[2][2] = Scale.Z;

                result.col[3][3] = 1;
                result.col[3][0] = Location.X;
                result.col[3][1] = Location.Y;
                result.col[3][2] = Location.Z;

                return result;
            }

            public static CG_Matrix ScaleMatrix(ref Vector3Df Location, float Scale)
            {
                CG_Matrix result = new CG_Matrix(Scale);
                result.col[3][3] = 1;
                result.col[3][0] = Location.X;
                result.col[3][1] = Location.Y;
                result.col[3][2] = Location.Z;
                return result;
            }


            public static CG_Matrix FromMattrix44(ref CG_Matrix Org)
            {
                CG_Matrix result = CG_Matrix.Zero();
                for (int c = 0; c < 4; c++) {
                    for (int r = 0; r < 4; r++) {
                        result.col[c][r] = Org.col[c][r];
                    }
                }
                return result;
            }

            public CG_Matrix Transpose()
            {
                float temp;
                for (int c = 0; c < 3; c++) {
                    for (int r = 0; r < 3; r++) {
                        if (r == c) break;
                        temp = col[c][r];
                        this.col[c][r] = this.col[r][c];
                        this.col[r][c] = temp;
                    }
                }
                return this;
            }

            public CG_Matrix Inverse()
            {
                CG_Matrix result = CG_Matrix.Identity();
                CG_Matrix Ref = CG_Matrix.FromMattrix44(ref this);
                for (int colm = 0; colm < Ref.col.Length; colm++) {
                    for (int row = colm; row < Ref.col.Length; row++) {
                        if (row == colm) {
                            #region Create Diagnal with a magnitiude of 1
                            if (Ref.col[colm][row] == 0) {
                                int i = row + 1;
                                bool Swaped = false;
                                while (i < 4) {
                                    if (Ref.col[colm][i] != 0) {
                                        Ref.SwapRows(row, i);
                                        result.SwapRows(row, i);
                                        Swaped = true;
                                        break;
                                    }
                                }
                                if (!Swaped)
                                    throw new Exception("Non Invertable Matrix");
                            }
                            if (Ref.col[colm][row] != 1) {
                                float nScale = -1 / Ref.col[colm][row];
                                Ref.ScaleRow(row, nScale);
                                result.ScaleRow(row, nScale);
                            }
                            #endregion
                        } else {
                            #region Deal with the lower triangle
                            if (Ref.col[colm][row] != 0) {
                                float nScale = 1 / Ref.col[colm][row];
                                Ref.AddRows(row, colm, nScale);
                                result.AddRows(row, colm, nScale);
                            }
                            #endregion
                        }
                    }
                    for (int row = 0; row < colm; row++) {
                        #region Deal with the Upper triangle
                        if (Ref.col[colm][row] != 0) {
                            float nScale = -1 / Ref.col[colm][row];
                            Ref.AddRows(row, colm, nScale);
                            result.AddRows(row, colm, nScale);
                        }
                        #endregion
                    }
                }
                return result;
            }

            public void AddRows(int DestRow, int RefRow, float nScale)
            {
                this.col[0][DestRow] += nScale * this.col[0][RefRow];
                this.col[1][DestRow] += nScale * this.col[1][RefRow];
                this.col[2][DestRow] += nScale * this.col[2][RefRow];
                this.col[3][DestRow] += nScale * this.col[3][RefRow];
            }

            public void ScaleRow(int DestRow, float nScale)
            {
                this.col[0][DestRow] *= nScale;
                this.col[1][DestRow] *= nScale;
                this.col[2][DestRow] *= nScale;
                this.col[3][DestRow] *= nScale;
            }

            public void SwapRows(int Row1, int Row2)
            {
                float[] tmp = new float[4];
                tmp[0] = this.col[0][Row1];
                tmp[1] = this.col[1][Row1];
                tmp[2] = this.col[2][Row1];
                tmp[3] = this.col[3][Row1];

                this.col[0][Row1] = this.col[0][Row2];
                this.col[1][Row1] = this.col[1][Row2];
                this.col[2][Row1] = this.col[2][Row2];
                this.col[3][Row1] = this.col[3][Row2];

                this.col[0][Row2] = this.col[0][Row1];
                this.col[1][Row2] = this.col[1][Row1];
                this.col[2][Row2] = this.col[2][Row1];
                this.col[3][Row2] = this.col[3][Row1];

            }

            public static CG_Matrix operator *(CG_Matrix lfs, CG_Matrix rhs)
            {
                CG_Matrix result = new CG_Matrix(0.0f);
                for (int r = 0; r < 4; r++) {
                    for (int c = 0; c < 4; c++) {
                        for (int k = 0; k < 4; k++) {
                            result.col[c][r] += lfs.col[k][r] * rhs.col[c][k];
                        }
                    }
                }
                return result;
            }

            public static CG_Matrix op_MultiplicationAssignment(CG_Matrix lfs, CG_Matrix rhs)
            {
                lfs = lfs * rhs;
                return lfs;
            }

            public static Vector3Df operator *(CG_Matrix lfs, Vector3Df V)
            {
                Vector3Df result = new Vector3Df();
                result.X = lfs.col[0][0] * V.X + lfs.col[1][0] * V.Y + lfs.col[2][0] * V.Z + lfs.col[3][0] * 1.0f;
                result.Y = lfs.col[0][1] * V.X + lfs.col[1][1] * V.Y + lfs.col[2][1] * V.Z + lfs.col[3][1] * 1.0f;
                result.Z = lfs.col[0][2] * V.X + lfs.col[1][2] * V.Y + lfs.col[2][2] * V.Z + lfs.col[3][2] * 1.0f;
                return result;
            }

            public float[] GetListValues()
            {
                float[] result = new float[16];
                result[0] = col[0][0];
                result[1] = col[0][1];
                result[2] = col[0][2];
                result[3] = col[0][3];
                result[4] = col[1][0];
                result[5] = col[1][1];
                result[6] = col[1][2];
                result[7] = col[1][3];
                result[8] = col[2][0];
                result[9] = col[2][1];
                result[10] = col[2][2];
                result[11] = col[2][3];
                result[12] = col[3][0];
                result[13] = col[3][1];
                result[14] = col[3][2];
                result[15] = col[3][3];
                return result;
            }

            public override string ToString()
            {
                string result = "";
                StringBuilder builder = new StringBuilder();
                for (int r = 0; r < 4; r++) {
                    builder.AppendFormat("\n|{0:E}\t\t{1:E}\t\t{2:E}\t\t{3:E}|", this.col[0][r], this.col[1][r], this.col[2][r], this.col[3][r]);
                }
                result = builder.ToString();
                return result;
            }

            public float SetScale
            {
                set { this.Scale = value; }
            }

            public Vector3Df SetTranslation
            {
                set { this.Pan = value; }
            }

            #endregion
        }


        public MainWindow()
        {
            InitializeComponent();
        }



        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            this.test();
            //this.test2();
        }


        //public Point3D getPositionCircle(double R, double psi, double fi)
        //{
        //    double sinPsi = Math.Sin(psi * Math.PI / 180);
        //    double cosPsi = Math.Cos(psi * Math.PI / 180);
        //    double sinFi = Math.Sin(fi * Math.PI / 180);
        //    double cosFi = Math.Cos(fi * Math.PI / 180);

        //    Point3D point = new Point3D();
        //    point.X = R * sinPsi * cosFi;
        //    point.Y = R * sinPsi * sinFi;
        //    point.Z = R * cosPsi;

        //    return point;
        //}


        //private void test2()
        //{
        //    System.Threading.Thread.Sleep(2000);
        //    PointLight light = new PointLight();
        //    light.Color = Brushes.White.Color;
        //    //Set the light position
        //    light.Position = new Point3D(0, 0, 5);
        //    viewport3D.Children.Add(new ModelVisual3D() { Content = light });

        //    //Create circles
        //    //drawCircle(new Point3D(3, 8, 10), 2, 100, Brushes.Green.Color);
        //    drawCircle(new Point3D(s1X, s1Y, s1Z), 2, 30, Brushes.Green.Color);
        //}

        //public void drawTriangle(Point3D p0, Point3D p1, Point3D p2, Color color)
        //{
        //    MeshGeometry3D mesh = new MeshGeometry3D();

        //    mesh.Positions.Add(p0);
        //    mesh.Positions.Add(p1);
        //    mesh.Positions.Add(p2);

        //    mesh.TriangleIndices.Add(0);
        //    mesh.TriangleIndices.Add(1);
        //    mesh.TriangleIndices.Add(2);

        //    SolidColorBrush brush = new SolidColorBrush();
        //    brush.Color = color;
        //    Material material = new DiffuseMaterial(brush);

        //    GeometryModel3D geometry = new GeometryModel3D(mesh, material);
        //    ModelUIElement3D model = new ModelUIElement3D();
        //    model.Model = geometry;

        //    viewport3D.Children.Add(model);
        //}

        //public void drawCircle(Point3D center, double R, int N, Color color)
        //{
        //    Model3DGroup circle = new Model3DGroup();
        //    Point3D[,] points = new Point3D[N, N];

        //    for (int i = 0; i < N; i++) {
        //        for (int j = 0; j < N; j++) {
        //            points[i, j] = getPositionCircle(R, i * 360 / (N - 1), j * 360 / (N - 1));
        //            points[i, j] += (Vector3D)center;
        //        }
        //    }

        //    Point3D[] p = new Point3D[4];
        //    for (int i = 0; i < N - 1; i++) {
        //        for (int j = 0; j < N - 1; j++) {
        //            p[0] = points[i, j];
        //            p[1] = points[i + 1, j];
        //            p[2] = points[i + 1, j + 1];
        //            p[3] = points[i, j + 1];
        //            drawTriangle(p[0], p[1], p[2], color);
        //            drawTriangle(p[2], p[3], p[0], color);
        //        }
        //    }
        //}
    }
}

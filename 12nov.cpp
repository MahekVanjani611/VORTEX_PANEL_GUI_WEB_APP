#include <iostream>
#include <string>
#include <cstdlib>
#include <map>
#include <Eigen/Eigen>
#include <vector>
#include <cmath>
#include <fstream>

#define pi acos(-1)
#define n 200 // Total number of nodes
#define d 0.001 // Increment along the chord for discretizing the chord into N number of points.
#define x_leading_edge 0.0
#define x_trailing_edge 1.0
#define X 0.0
#define Y 0.0


using namespace Eigen;
using namespace std;


// MatrixXd geometry(double x_nd,double ymc,double xmc,double tm,double trailing_edge_type);
// void plot_airfoil(VectorXd &camber_x, VectorXd &camber_y, VectorXd &xu, VectorXd &xl, VectorXd &yu, VectorXd &yl, double x_nd, int N,double ymc, double xmc, double tm, double trailing_edge_type);
// void nodal_coordinates(double x_nd, double delta_theta, VectorXd &x_pp, VectorXd &y_pp,double ymc,double xmc,double tm, double trailing_edge_type);
// void controlpoints(VectorXd &x_pp, VectorXd &y_pp, VectorXd &x_cp, VectorXd &y_cp);
// MatrixXd influence_matrix(double point1_x, double point1_y, double point2_x, double point2_y, double desired_point_x, double desired_point_y);
// void Amatrix(MatrixXd &A, VectorXd &x_cp, VectorXd &y_cp, VectorXd &x_pp, VectorXd &y_pp);
// void random_point_flowfield(VectorXd &x_pp, VectorXd &y_pp);
// void B_vector_particular_alfa(VectorXd &x_pp, VectorXd &y_pp, double Vinf, double alpha, VectorXd &B1);
// void Gamma_vector(VectorXd &G, VectorXd &B1, MatrixXd &A);
// //void vel_any_point_flowfield_due_to_panels(VectorXd &x_pp, VectorXd &y_pp,VectorXd &G);
// tuple<double, double, double> section_lift_and_moment_coefficient(VectorXd &G, VectorXd &x_pp, VectorXd &y_pp, double alfa, double Vinf);
// void panel(VectorXd &l_x, VectorXd &l_y, VectorXd &l, VectorXd &x_pp, VectorXd &y_pp);
// //void CL_vs_angle_of_attack(MatrixXd &A, double Vinf, VectorXd &l_x, VectorXd &l_y, VectorXd &l,VectorXd &x_pp, VectorXd &y_pp);
// VectorXd velocity(VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double alfa);
// double pressure_coefficient(VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double alfa);
// VectorXd normalize_2d(VectorXd v);
// //VectorXd func(VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double alfa);
// VectorXd cross_product_2d(VectorXd a);
// void normal_function_for_panels(MatrixXd &unit_normal, VectorXd &l_x, VectorXd &l_y);
// void plot_pressure(VectorXd &x_cp, VectorXd &y_cp, VectorXd &x_pp, VectorXd &y_pp, MatrixXd &unit_normal, VectorXd &G, double Vinf, double alfa);
// MatrixXd surface_normal(double x_nd);
// MatrixXd surface_tangent(double x_nd);
// VectorXd surface_tangential_velocity(double x_nd, VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double alfa);
// VectorXd normalize(VectorXd v1);
// //void plot_tangential_vel(VectorXd &l_x, VectorXd &l_y, VectorXd &x_pp, VectorXd &y_pp, VectorXd &x_cp, VectorXd &y_cp, VectorXd &G, double Vinf, double alfa);
// //void velocity_calc(int N, double x_nd, VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double alfa);
// //void plot_streamline(int N, VectorXd &xu, VectorXd &yu, double x_nd, VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double alfa);


MatrixXd geometry(double x_nd,double ymc,double xmc,double tm,double trailing_edge_type,double alfa)
{
   double p,q;
    MatrixXd airfoil_points(3, 2);

    double yc, der_yc, t;

    p = (ymc / 100); // p non dimensionalised ymc [p=ymc/c]
    q = (xmc * 0.1); // q non dimensionalised xmc [q=xmc/c]
    tm = (tm / 100);

    if (x_nd <= q)
    {
        yc = (p * (2.0 * (x_nd / q) - (x_nd / q) * (x_nd / q)));
        der_yc = ((2.0 * p / q) * (1.0 - (x_nd / q)));
    }
    else
    {
        yc = (p * (2.0 * (1 - x_nd) / (1 - q) - pow(((1 - x_nd) / (1 - q)), 2.0)));
        der_yc = ((2.0 * p / (1.0 - q)) * (((1.0 - x_nd) / (1.0 - q)) - 1.0));
    }
    if (trailing_edge_type == 1.0) // 1 indicates open trailing edge
    {
        t = (tm * (2.969 * sqrt(x_nd) - 1.260 * (x_nd)-3.516 * pow((x_nd), 2) + 2.843 * pow((x_nd), 3.0) - 1.015 * (pow((x_nd), 4.0))));
    }
    else
    {
        t = (tm * (2.980 * sqrt(x_nd) - 1.320 * (x_nd)-3.286 * pow((x_nd), 2) + 2.441 * pow((x_nd), 3.0) - 0.815 * (pow((x_nd), 4.0))));
    }
    if (p == 0.0) //[SYMMETRIC AIRFOIL]
    {
        airfoil_points(0, 0) = x_nd;     // x_upper
        airfoil_points(0, 1) = t / 2.0;  // y_upper
        airfoil_points(1, 0) = x_nd;     // x_lower
        airfoil_points(1, 1) = -t / 2.0; // y_lower
        airfoil_points(2, 0) = x_nd;
        airfoil_points(2, 1) = 0.0;
    }
    else
    {
        airfoil_points(0, 0) = x_nd - t * der_yc / (2 * (sqrt(1 + (der_yc * der_yc)))); // x_upper
        airfoil_points(0, 1) = yc + t / (2 * (sqrt(1 + (der_yc * der_yc))));            // y_upper
        airfoil_points(1, 0) = x_nd + t * der_yc / (2 * (sqrt(1 + (der_yc * der_yc)))); // x_lower
        airfoil_points(1, 1) = yc - t / (2 * (sqrt(1 + (der_yc * der_yc))));            // y_lower
        airfoil_points(2, 0) = x_nd;
        airfoil_points(2, 1) = yc;
    }
        MatrixXd rotation(2, 2);
    VectorXd upper(2);
    VectorXd lower(2);
    VectorXd camber(2);

    rotation(0, 0) = cos(alfa);
    rotation(0, 1) = sin(alfa);
    rotation(1, 0) = -sin(alfa);
    rotation(1, 1) = cos(alfa);

    upper(0) = airfoil_points(0, 0); // xu
    upper(1) = airfoil_points(0, 1); // yu
    lower(0) = airfoil_points(1, 0); // xl
    lower(1) = airfoil_points(1, 1); // yl
    camber(0) = airfoil_points(2, 0);
    camber(1) = airfoil_points(2, 1);

    upper = rotation * upper;
    lower = rotation * lower;
    camber = rotation * camber;

    airfoil_points(0, 0) = upper(0);
    airfoil_points(0, 1) = upper(1);
    airfoil_points(1, 0) = lower(0);
    airfoil_points(1, 1) = lower(1);
    airfoil_points(2, 0) = camber(0);
    airfoil_points(2, 1) = camber(1);

    
    return airfoil_points;
}

void plot_airfoil(VectorXd &camber_x, VectorXd &camber_y, VectorXd &xu, VectorXd &xl, VectorXd &yu, VectorXd &yl, double x_nd, int N,double ymc,double xmc,double tm,double trailing_edge_type,double alfa)
{
    ofstream myfile1;
    myfile1.open("airfoil_upper_lower_points.csv");

    MatrixXd airfoil_points(3, 2);
    for (int i = 0; i < N; i++)
    {
        x_nd = i * d;
        airfoil_points = geometry(x_nd,ymc,xmc,tm,trailing_edge_type,alfa);
        xu(i) = airfoil_points(0, 0);
        yu(i) = airfoil_points(0, 1);
        xl(i) = airfoil_points(1, 0);
        yl(i) = airfoil_points(1, 1);
        camber_x(i) = airfoil_points(2, 0);
        camber_y(i) = airfoil_points(2, 1);
        myfile1 << camber_x(i) << "\t" << camber_y(i) << endl;
    }
}
void nodal_coordinates(double x_nd, double delta_theta, VectorXd &x_pp, VectorXd &y_pp,double ymc,double xmc,double tm,double trailing_edge_type,double alfa)
{// cosine clustering
    ofstream myfile2;
    myfile2.open("panel_points.csv");
    double theta;
    MatrixXd panel_points(3, 2);

    for (int i = (n / 2); i >= 1; i--) //[storing lower coordinates first]
    {
        theta = (i - 0.5) * delta_theta; //cosine clustering
        x_nd = 0.5 * (1 - cos(theta));
        panel_points = geometry(x_nd,ymc,xmc,tm,trailing_edge_type,alfa);
        x_pp(n / 2 - i) = panel_points(1, 0); //[x_lower]
        y_pp(n / 2 - i) = panel_points(1, 1); //[y_lower]
      
    }
    for (int i = 1; i <= (n / 2); i++) //[storing upper coordinates]
    {
        theta = (i - 0.5) * delta_theta; // cosine clustering
        x_nd = 0.5 * (1 - cos(theta));
        panel_points = geometry(x_nd,ymc,xmc,tm,trailing_edge_type,alfa);
        x_pp(n / 2 + i - 1) = panel_points(0, 0); //[x_upper]
        y_pp(n / 2 + i - 1) = panel_points(0, 1); //[y_upper]
        // x_camber(i)=panel_points(2, 0); 
        // y_camber(i)=panel_points(2, 1); 
    }

    for (int i = 0; i < n; i++)
    {
        myfile2 << x_pp(i) << "\t" << y_pp(i) << "\t"<< endl;
        
    }
}
// ******************************************* storing  coordinates  of control points  *****************************************//
void controlpoints(VectorXd &x_pp, VectorXd &y_pp, VectorXd &x_cp, VectorXd &y_cp)
{
    for (int j = 0; j < (n - 1); j++)
    {
        x_cp(j) = x_pp(j) - (x_pp(j) - x_pp(j + 1)) / 2;
        y_cp(j) = y_pp(j) - (y_pp(j) - y_pp(j + 1)) / 2;
    }
}
MatrixXd influence_matrix(double point1_x, double point1_y, double point2_x, double point2_y, double desired_point_x, double desired_point_y)
{
    MatrixXd P1(2, 2);
    MatrixXd P2(2, 2);
    MatrixXd P(2, 2);
    VectorXd vec(2);
    VectorXd LPC(2);

    double dx, dy, li;
    double phi, psi;
    double eta, geta;

    dx = (point2_x - point1_x);
    dy = (point2_y - point1_y);

    li = sqrt(dx * dx + dy * dy);

    P1(0, 0) = dx;
    P1(0, 1) = dy;
    P1(1, 0) = -dy;
    P1(1, 1) = dx;

    vec(0) = desired_point_x - point1_x;
    vec(1) = desired_point_y - point1_y;

    LPC = (P1 * vec) / li;
    geta = LPC(0);
    eta = LPC(1);
    phi = atan2((eta * li), ((eta * eta) + (geta * geta) - (geta * (li))));
    psi = 0.5 * log(((geta * geta) + (eta * eta)) / (((geta - li) * (geta - li)) + (eta * eta)));

    P1(0, 0) = dx;
    P1(0, 1) = -dy;
    P1(1, 0) = dy;
    P1(1, 1) = dx;

    P2(0, 0) = (li - geta) * phi + (eta * psi);
    P2(0, 1) = (geta * phi) - (eta * psi);
    P2(1, 0) = (eta * phi - (li - geta) * psi - li);
    P2(1, 1) = ((-eta * phi) - (geta * psi) + li);

    P = (P1 * P2) / (2 * pi * li * li);
    return P;
}

void Amatrix(MatrixXd &A, VectorXd &x_cp, VectorXd &y_cp, VectorXd &x_pp, VectorXd &y_pp)
{
    ofstream myfile3;
    myfile3.open("A_matrix_file.csv");

    MatrixXd pcm(2, 2); // panel coefficient matrix
    double x1, x2, y1, y2, dx, dy, li;

    for (int j = 0; j < n; j++) // j index is row rows
    {
        for (int i = 0; i < n; i++) // i index is for columns
        {
            A(j, i) = 0.0;
        }
    }
    for (int j = 0; j < n - 1; j++)
    {

        x1 = x_pp(j);
        x2 = x_pp(j + 1);
        y1 = y_pp(j);
        y2 = y_pp(j + 1);
        dx = x2 - x1;
        dy = y2 - y1;

        li = sqrt((dx * dx) + (dy * dy));

        for (int i = 0; i < n - 1; i++)
        {
            pcm = influence_matrix(x_pp(i), y_pp(i), x_pp(i + 1), y_pp(i + 1), x_cp(j), y_cp(j));
            //cout << pcm << endl<< endl;
            A(j, i) = A(j, i) + dx / li * pcm(1, 0) - dy / li * pcm(0, 0);
            A(j, i + 1) = A(j, i + 1) + dx / li * pcm(1, 1) - dy / li * pcm(0, 1);
        }
    }
    A(n - 1, 0) = 1.0;
    A(n - 1, n - 1) = 1.0;
    myfile3 << A << endl;
}
void random_point_flowfield(VectorXd &x_pp, VectorXd &y_pp)
{
    MatrixXd panel_coeff(2, 2);
    for (int i = 0; i < n - 1; i++)
    {
        panel_coeff = influence_matrix(x_pp(i), y_pp(i), x_pp(i + 1), y_pp(i + 1), X, Y);
    }
}
void B_vector_particular_alfa(VectorXd &x_pp, VectorXd &y_pp, double Vinf, double alfa, VectorXd &B)
{
    ofstream myfile4, panel_length;
    myfile4.open("B1_vector_file0.csv");
    panel_length.open("panel_length0.csv");
    alfa=0.0;

    double x1, x2, y1, y2, dx, dy, li;

    for (int i = 0; i < n - 1; i++)
    {
        x1 = x_pp(i);
        x2 = x_pp(i + 1);
        y1 = y_pp(i);
        y2 = y_pp(i + 1);
        dx = x2 - x1;
        dy = y2 - y1;
        li = sqrt((dx * dx) + (dy * dy));
        panel_length << dx << "\t" << dy << "\t" << li << endl;
        B(i) = Vinf * ((cos(alfa) * dy - sin(alfa) * dx)) / li;
    }
    B(n - 1) = 0.0;
    myfile4 << B << endl;
}

void Gamma_vector(VectorXd &G, VectorXd &B, MatrixXd &A)
{
    ofstream myfile5;
    myfile5.open("Gamma_vector_file.csv");
    G = A.fullPivHouseholderQr().solve(B);
    myfile5 << G << endl;
}

tuple<double, double, double> section_lift_and_moment_coefficient(VectorXd &G, VectorXd &x_pp, VectorXd &y_pp, double alfa, double Vinf)
{
    double CL = 0.0;
    double Cmle = 0.0;
    double Cmc4;
    double x1, x2, y1, y2, dx, dy, li;

    for (int i = 0; i < n - 1; i++)
    {
        x1 = x_pp(i);
        x2 = x_pp(i + 1);
        y1 = y_pp(i);
        y2 = y_pp(i + 1);
        dx = x2 - x1;
        dy = y2 - y1;

        li = sqrt((dx * dx) + (dy * dy));

        CL = CL + (li * (G(i) + G(i + 1))) / (Vinf);

        Cmle = Cmle + ((li / Vinf) * ((2 * x_pp(i) * G(i) + x_pp(i) * G(i + 1) + x_pp(i + 1) * G(i) + 2 * x_pp(i + 1) * G(i + 1)) * cos(alfa) + (2 * y_pp(i) * G(i) + y_pp(i) * G(i + 1) + y_pp(i + 1) * G(i) + 2 * y_pp(i + 1) * G(i + 1)) * sin(alfa)));
    }
    Cmle = Cmle * (-1.0 / 3.0);
    Cmc4 = CL * 0.25*cos(alfa) + Cmle;

    cout << "section_lift_coefficient=" << CL << endl
         << "sectional_moment_about_leading_edge=" << Cmle << endl
         << "sectional_moment_about_quarter_chord=" << Cmc4 << endl;
    ofstream outfile("results.csv");
    if (outfile.is_open())
    {
        outfile << "section_lift_coefficient=" << CL << endl;
        outfile << "sectional_moment_about_leading_edge=" << Cmle << endl;
        outfile << "sectional_moment_about_quarter_chord=" << Cmc4 << endl;
        outfile.close();
    }
    else
    {
        cerr << "Error: Unable to open file for writing results." << endl;
    }
    return make_tuple(CL,Cmle,Cmc4);
}
void panel(VectorXd &l_x, VectorXd &l_y, VectorXd &l, VectorXd &x_pp, VectorXd &y_pp)
{
    for (int i = 0; i < n - 1; i++)
    {
        l_x(i) = x_pp(i + 1) - x_pp(i);
        l_y(i) = y_pp(i + 1) - y_pp(i);
        l(i) = sqrt(l_x(i) * l_x(i) + l_y(i) * l_y(i));
    }
}
VectorXd velocity(VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double alfa)
{
    VectorXd VEL(2);
    VEL(0) = 0.0;
    VEL(1) = 0.0;

    MatrixXd P(2, 2);

    VectorXd vortex_strength(2);
    // VectorXd freestream(2);

    double resultant_velocity, C_p;

    for (int i = 0; i < n - 1; i++)
    {
        P = influence_matrix(x_pp(i), y_pp(i), x_pp(i + 1), y_pp(i + 1), x, y);
        vortex_strength(0) = G(i);
        vortex_strength(1) = G(i + 1);
        VEL = VEL + (P * vortex_strength);
    }

    VEL(0) = VEL(0) + Vinf;
    VEL(1) = VEL(1);
    return (VEL);
}
double pressure_coefficient(VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double alfa)
{
    VectorXd VEL(2);
    VEL(0) = 0.0;
    VEL(1) = 0.0;

    MatrixXd P(2, 2);

    VectorXd vortex_strength(2);
    VectorXd freestream(2);

    double resultant_velocity, C_p;

    for (int i = 0; i < n - 1; i++)
    {
        P = influence_matrix(x_pp(i), y_pp(i), x_pp(i + 1), y_pp(i + 1), x, y);
        vortex_strength(0) = G(i);
        vortex_strength(1) = G(i + 1);
        VEL = VEL + (P * vortex_strength);
    }
    freestream(0) = Vinf * cos(alfa);
    freestream(1) = Vinf * sin(alfa);
    VEL = VEL + freestream;

    resultant_velocity = sqrt((VEL(0) * VEL(0) + VEL(1) * VEL(1))); // [VEL(0)=V_resultant_xcomp][VEL(1)=V_resultant_ycomp]
    C_p = 1 - ((resultant_velocity / Vinf) * (resultant_velocity / Vinf));
    return C_p;
}
VectorXd normalize_2d(VectorXd v)
{
    double mag;
    mag = sqrt(v(0) * v(0) + v(1) * v(1));
    v = v / mag;
    return v;
}
VectorXd func(VectorXd y0, VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double alfa)
{
    x = y0(0);
    y = y0(1);
    VectorXd V(2);
    VectorXd V_unit(2);
    V = velocity(x_pp, y_pp, x, y, G, Vinf, alfa);
    V_unit = normalize_2d(V);
    return (V_unit);
}
VectorXd cross_product_2d(VectorXd a) // axk [k is the unit vector along z dirn.]
{
    VectorXd b(2);
    b(0) = a(1);
    b(1) = -a(0);
    return b;
}

void normal_function_for_panels(MatrixXd &unit_normal, VectorXd &l_x, VectorXd &l_y) /*here we are finding the unit normal of n-1 panels*/
{
    VectorXd panel_length(2);
    VectorXd unit_normal_vector(2);
    for (int i = 0; i < n - 1; i++)
    {
        panel_length(0) = -l_y(i);
        panel_length(1) = l_x(i);
        unit_normal_vector = normalize_2d(panel_length);
        unit_normal(i, 0) = unit_normal_vector(0);
        unit_normal(i, 1) = unit_normal_vector(1);
    }
}

void plot_pressure(VectorXd &x_cp, VectorXd &y_cp, VectorXd &x_pp, VectorXd &y_pp, MatrixXd &unit_normal, VectorXd &G, double Vinf, double alfa)
{
    ofstream myfile7;
    myfile7.open("Cp.csv");

    VectorXd Cp(n - 1);
    double offset = 0.001;

    for (int j = 0; j < n - 1; j++)
    {
        Cp(j) = pressure_coefficient(x_pp, y_pp, x_cp(j) + unit_normal(j, 0) * offset, y_cp(j) + unit_normal(j, 1) * offset, G, Vinf, alfa);
        myfile7 << x_cp(j) << "\t" << Cp(j) << endl;
    }
}
MatrixXd surface_normal(double x_nd,double ymc,double xmc,double tm,double trailing_edge_type,double alfa)
{
    VectorXd nu(2);
    VectorXd nl(2);
    VectorXd vec(2);
    VectorXd u1(2);
    VectorXd l1(2);
    VectorXd u2(2);
    VectorXd l2(2);

    MatrixXd airfoil_points_matrix(3, 2);
    MatrixXd normal(2, 2);

    double dx, eps;
    dx = 1.e-12;
    eps = 1.e-10;

    if (fabs(x_nd - x_leading_edge) < eps)
    {
        airfoil_points_matrix = geometry(x_nd + dx,ymc,xmc,tm,trailing_edge_type,alfa);
        u1(0) = airfoil_points_matrix(0, 0);
        u1(1) = airfoil_points_matrix(0, 1);
        l1(0) = airfoil_points_matrix(1, 0);
        l1(1) = airfoil_points_matrix(1, 1);
        vec(0) = u1(0) - l1(0);
        vec(1) = u1(1) - l1(1);
        nu = cross_product_2d(vec);
        nu = normalize_2d(nu);
        nu = -nu;
        nl = nu;
    }
    else
    {
        if (fabs(x_nd - x_trailing_edge) < eps)
        {
            airfoil_points_matrix = geometry(x_nd - dx,ymc,xmc,tm,trailing_edge_type,alfa);
            u1(0) = airfoil_points_matrix(0, 0);
            u1(1) = airfoil_points_matrix(0, 1);
            l1(0) = airfoil_points_matrix(1, 0);
            l1(1) = airfoil_points_matrix(1, 1);
            vec(0) = u1(0) - l1(0);
            vec(1) = u1(1) - l1(1);
            nu = cross_product_2d(vec);
            nu = normalize_2d(nu);
            nl = nu;
        }
        else
        {
            airfoil_points_matrix = geometry(x_nd - dx,ymc,xmc,tm,trailing_edge_type,alfa);
            u1(0) = airfoil_points_matrix(0, 0);
            u1(1) = airfoil_points_matrix(0, 1);
            l1(0) = airfoil_points_matrix(1, 0);
            l1(1) = airfoil_points_matrix(1, 1);
            airfoil_points_matrix = geometry(x_nd + dx,ymc,xmc,tm,trailing_edge_type,alfa);
            u2(0) = airfoil_points_matrix(0, 0);
            u2(1) = airfoil_points_matrix(0, 1);
            l2(0) = airfoil_points_matrix(1, 0);
            l2(1) = airfoil_points_matrix(1, 1);
            vec(0) = u2(0) - u1(0);
            vec(1) = u2(1) - u2(1);
            nu = cross_product_2d(vec);
            nu = normalize_2d(nu);
            nu = -nu;
            vec(0) = l2(0) - l1(0);
            vec(1) = l2(1) - l1(1);
            nl = cross_product_2d(vec);
            nl = normalize_2d(nl);
        }
    }
    normal(0, 0) = nu(0);
    normal(0, 1) = nu(1);
    normal(1, 0) = nl(0);
    normal(1, 1) = nl(1);
    return (normal);
}
MatrixXd surface_tangent(double x_nd,double ymc,double xmc,double tm,double trailing_edge_type,double alfa)
{
    VectorXd tu(2);
    VectorXd tl(2);
    MatrixXd airfoil_points(3, 2);
    VectorXd u1(2);
    VectorXd l1(2);
    VectorXd u2(2);
    VectorXd l2(2);
    MatrixXd tangent(2, 2);

    double dx, eps;
    dx = 1.e-12;
    eps = 1.e-10;

    if (fabs(x_nd - x_leading_edge) < eps)
    {
        airfoil_points = geometry(x_nd + dx,ymc,xmc,tm,trailing_edge_type,alfa);
        u1(0) = airfoil_points(0, 0);
        u1(1) = airfoil_points(0, 1);
        l1(0) = airfoil_points(1, 0);
        l1(1) = airfoil_points(1, 1);
        // cout << x_nd+dx <<"\t" << u1(0) <<"\t"<< u1(1)<<"\t" << l1(0)<<"\t"<< l1(1) << endl;
        tu(0) = u1(0) - l1(0);
        tu(1) = u1(1) - l1(1);
        tu = normalize_2d(tu);
        tl = tu;
    }
    else
    {
        if (fabs(x_nd - x_trailing_edge) < eps)
        {
            airfoil_points = geometry(x_nd - dx,ymc,xmc,tm,trailing_edge_type,alfa);
            u1(0) = airfoil_points(0, 0);
            u1(1) = airfoil_points(0, 1);
            l1(0) = airfoil_points(1, 0);
            l1(1) = airfoil_points(1, 1);
            tu(0) = l1(0) - u1(0);
            tu(1) = l1(1) - u1(1);
            tu = normalize_2d(tu);
            tl = tu;
        }
        else
        {
            airfoil_points = geometry(x_nd - dx,ymc,xmc,tm,trailing_edge_type,alfa);
            u1(0) = airfoil_points(0, 0);
            u1(1) = airfoil_points(0, 1);
            l1(0) = airfoil_points(1, 0);
            l1(1) = airfoil_points(1, 1);
            airfoil_points = geometry(x_nd + dx,ymc,xmc,tm,trailing_edge_type,alfa);
            u2(0) = airfoil_points(0, 0);
            u2(1) = airfoil_points(0, 1);
            l2(0) = airfoil_points(1, 0);
            l2(1) = airfoil_points(1, 1);

            tu(0) = u2(0) - u1(0);
            tu(1) = u2(1) - u1(1);
            tu = normalize_2d(tu);

            tl(0) = l1(0) - l2(0);
            tl(1) = l1(1) - l2(1);
            tl = normalize_2d(tl);
        }
    }
    tangent(0, 0) = tu(0);
    tangent(0, 1) = tu(1);
    tangent(1, 0) = tl(0);
    tangent(1, 1) = tl(1);
    return (tangent);
}
VectorXd surface_tangential_velocity(double x_nd, VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf,double ymc,double xmc,double tm,double trailing_edge_type, double alfa)
{
    MatrixXd airfoil_points(3, 2);
    MatrixXd tangent(2, 2);

    VectorXd gu(2);
    VectorXd gl(2);
    VectorXd tu(2);
    VectorXd tl(2);

    VectorXd V(2);
    double Vu_t, Vl_t;
    VectorXd tangential_velocity(2);

    airfoil_points = geometry(x_nd,ymc,xmc,tm,trailing_edge_type,alfa);
    gu(0) = airfoil_points(0, 0); // xu
    gu(1) = airfoil_points(0, 1); // yu
    gl(0) = airfoil_points(1, 0); // xl
    gl(1) = airfoil_points(1, 1); // yl

    tangent = surface_tangent(x_nd,ymc,xmc,tm,trailing_edge_type,alfa);
    tu(0) = tangent(0, 0);
    tu(1) = tangent(0, 1);
    tl(0) = tangent(1, 0);
    tl(1) = tangent(1, 1);

    V = velocity(x_pp, y_pp, gu(0), gu(1), G, Vinf, alfa);
    Vu_t = tu(0) * V(0) + tu(1) * V(1);
    V = velocity(x_pp, y_pp, gl(0), gl(1), G, Vinf, alfa);
    Vl_t = tl(0) * V(0) + tl(1) * V(1);

    tangential_velocity(0) = Vu_t;
    tangential_velocity(1) = Vl_t;

    return (tangential_velocity);
}

VectorXd rk_method(VectorXd y0, double h, VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double alfa)
{
    VectorXd k1(2);
    VectorXd k2(2);
    VectorXd k3(2);
    VectorXd k4(2);
    VectorXd ynew(2);
    k1 = func(y0, x_pp, y_pp, x, y, G, Vinf, alfa);
    k2 = func(y0 + h * k1 / 2, x_pp, y_pp, x, y, G, Vinf, alfa);
    k3 = func(y0 + h * k2 / 2, x_pp, y_pp, x, y, G, Vinf, alfa);
    k4 = func(y0 + h * k3, x_pp, y_pp, x, y, G, Vinf, alfa);
    ynew = y0 + ((1.0 / 6.0) * h * (k1 + 2 * k2 + 2 * k3 + k4));
    return (ynew);
}
void plot_streamline(int N, VectorXd &xu, VectorXd &yu, double x_nd, VectorXd &x_pp, VectorXd &y_pp, double x, double y, VectorXd &G, double Vinf, double ymc,double xmc,double tm,double trailing_edge_type,double alfa)
{
    double x_nd_new, x_nd2, dx = 0.1;
    double search, upper = 1, lower = 2;
    VectorXd tangential_velocity(2);
    double Vu_t, Vl_t;
    double f, f2;
    double m; // slope
    MatrixXd airfoil(3, 2);
    MatrixXd normal(2, 2);
    double stag_point_x;
    double stag_point_y;
    double nu_x, nl_x;
    double nu_y, nl_y;
    VectorXd stag_point(2);
    VectorXd norm(2);
    VectorXd velocity_stag(2);

    x_nd = x_leading_edge;
    x_nd_new = 0.0;
    tangential_velocity = surface_tangential_velocity(x_nd, x_pp, y_pp, x, y, G, Vinf,ymc,xmc,tm,trailing_edge_type, alfa);
    // cout << tangential_velocity << endl;
    Vu_t = tangential_velocity(0);
    f = Vu_t;
    cout << "Tangential Velocity at leading edge =" << Vu_t << endl;

    if (fabs(f) < 1.e-12)
    {
        cout << "stagnation1 is at leading edge" << endl;
        airfoil = geometry(x_nd,ymc,xmc,tm,trailing_edge_type,alfa);
        stag_point_x = airfoil(2, 0);
        stag_point_y = airfoil(2, 1);
        stag_point(0) = stag_point_x;
        stag_point(1) = stag_point_y;
        nu_x = normal(0, 0);
        nu_y = normal(0, 1);
        norm(0) = nu_x;
        norm(1) = nu_y;
    }
    else
    {
        if (f > 0.0)
        {
            search = lower;
            cout << "searching lower surface for forward stagnation point" << endl;
        }
        if (f < 0.0)
        {
            search = upper;
            cout << "searching upper surface for forward stagnation point" << endl;
        }
        x_nd2 = x_nd + dx;
        while (fabs(f) > 1.0e-5 * Vinf)
        {
            tangential_velocity = surface_tangential_velocity(x_nd2, x_pp, y_pp, x, y, G, Vinf, ymc,xmc,tm,trailing_edge_type,alfa);
            Vu_t = tangential_velocity(0);
            Vl_t = tangential_velocity(1);
            if (search == upper)
            {
                f2 = Vu_t;
            }
            if (search == lower)
            {
                f2 = Vl_t;
            }
            m = (f2 - f) / (x_nd2 - x_nd);
            // cout << m << endl;
            x_nd_new = x_nd - (f / m);
            // cout << x_nd_new << endl;
            // cout << x_nd << "\t" << f << "\t" << x_nd2 << "\t" << f2 << "\t" << m << "\t" << x_nd_new << endl;
            if (x_nd_new < 0.0)
            {
                x_nd_new = 0.0;
            }
            x_nd = x_nd2;
            x_nd2 = x_nd_new;
            f = f2;
        }
        x_nd = x_nd_new;
        // cout << x_nd<< endl;
        if (search == upper)
        {
            airfoil = geometry(x_nd,ymc,xmc,tm,trailing_edge_type,alfa);
            stag_point_x = airfoil(0, 0);
            stag_point_y = airfoil(0, 1);
            stag_point(0) = stag_point_x;
            stag_point(1) = stag_point_y;
            normal = surface_normal(x_nd,ymc,xmc,tm,trailing_edge_type,alfa);
            nu_x = normal(0, 0);
            nu_y = normal(0, 1);
            norm(0) = nu_x;
            norm(1) = nu_y;
        }
        if (search == lower)
        {
            airfoil = geometry(x_nd,ymc,xmc,tm,trailing_edge_type,alfa);
            stag_point_x = airfoil(1, 0);
            stag_point_y = airfoil(1, 1);
            stag_point(0) = stag_point_x;
            stag_point(1) = stag_point_y;
            normal = surface_normal(x_nd,ymc,xmc,tm,trailing_edge_type,alfa);
            nl_x = normal(1, 0);
            nl_y = normal(1, 1);
            norm(0) = nl_x;
            norm(1) = nl_y;
        }
    }
    cout << "stagnation_point=" << stag_point << endl;
    velocity_stag = velocity(x_pp, y_pp, stag_point_x, stag_point_y, G, Vinf, alfa); //[stagnation point velocity]
    cout << "stagnation point velocity=" << velocity_stag << endl;
    tangential_velocity = surface_tangential_velocity(x_nd, x_pp, y_pp, x, y, G, Vinf,ymc,xmc,tm,trailing_edge_type, alfa); //[tangential velocity at stagnation point]
    cout << "tangential_velocity at stagnation point=" << tangential_velocity << endl;
    //velocity_stag = velocity(x_pp, y_pp, xu(N), yu(N), G, Vinf, alfa);
    // cout << "velocity at the trailing edge =" << velocity_stag << endl;

    int nlines = 7;

    double x_start = -1.0;
    double y_start; // [y start value will be the y coordinate of the forward stagnation line corresponding to x=-1.0]

    ofstream myfile9, myfile10, myfile11;
    myfile9.open("forward_stagline.csv");
    myfile10.open("aft_stagline.csv");

    VectorXd pos(2);
    double s = 0.0;
    double ds = 0.01;
    double x_lower_lim = -1.0;
    double x_upper_lim = 2.0;

    /************ plotting of forward stagnation streamline ***********/

    myfile9 << stag_point(0) << "\t" << stag_point(1) << endl;

    if (alfa == 0.0)
    {
        // stag_point = stag_point + 0.01 * norm;
        while (pos(0) >= x_lower_lim && pos(0) <= x_upper_lim && fabs(s) < 5.0)
        {
            pos = rk_method(stag_point, -ds, x_pp, y_pp, x, y, G, Vinf, alfa);
            stag_point = pos;
            s = s + ds;
            myfile9 << pos(0) << "\t" << pos(1) << endl;
        }
    }
    else
    {
        stag_point = stag_point + 0.01 * norm;
        while (pos(0) >= x_lower_lim && pos(0) <= x_upper_lim && fabs(s) < 5.0)
        {
            pos = rk_method(stag_point, -ds, x_pp, y_pp, x, y, G, Vinf, alfa);
            stag_point = pos;
            s = s + ds;
            myfile9 << pos(0) << "\t" << pos(1) << endl;
        }
    }
    y_start = pos(1);
    cout << "starting value of y" << y_start << endl;
    cout << y_start << endl;

    /**************** plotting of aftward stagnation streamline ***********/
    VectorXd stagpoint_te(2);
    // [since the trailing edge is closed xu=xl and yu=yl]
    cout << x_pp(n-1) << "\t" << y_pp(n-1) << endl;
    stagpoint_te(0) = x_pp(n-1);		// [rear stagnation point = trailing edge(kutta condition)]
    stagpoint_te(1) = y_pp(n-1);     // [rear stagnation point = trailing edge(kutta condition)]
    stagpoint_te = stagpoint_te + 0.0001 * norm;
    s = 0.0;
    ds = 0.01;
    pos(0) = 0.0;
    pos(1) = 0.0;
    pos = rk_method(stagpoint_te, ds, x_pp, y_pp, x, y, G, Vinf, alfa);
    cout << pos << endl;

    myfile10 << stagpoint_te(0) << "\t" << stagpoint_te(1) << endl;
    while (pos(0) >= x_lower_lim && pos(0) <= x_upper_lim && fabs(s) < 20.0)
    {
        pos = rk_method(stagpoint_te, ds, x_pp, y_pp, x, y, G, Vinf, alfa);
        stagpoint_te = pos;
        s = s + ds;
        myfile10 << pos(0) << "\t" << pos(1) << endl;
    }

    /**************  plotting of rest of the streamlines ****************/

    double r;
    double dy = 0.1;
    VectorXd start(2);

    for (int i = 0; i < nlines; i++)
    {
        pos(0) = 0.0;
        pos(1) = 0.0;
        s = 0.0;
        ds = 0.01;
        string name = "other_streamlines_";
        name += to_string(i);
        name += ".csv";
        myfile11.open(name.c_str());
        r = (i + 1.0) * dy;
        start(0) = x_start;
        start(1) = y_start + r;
        cout << i + 1 << endl
             << endl;
        cout << start[1] << endl;

        while (pos(0) >= x_lower_lim && pos(0) <= x_upper_lim && fabs(s) < 5.0)
        {
            pos = rk_method(start, ds, x_pp, y_pp, x, y, G, Vinf, alfa);
            start = pos;
            s = s + ds;
            myfile11 << pos(0) << "\t" << pos(1) << endl;
        }
        myfile11.close();
    }
    pos(0) = 0.0;
    pos(1) = 0.0;
    s = 0.0;
    ds = 0.01;
    for (int i = 0; i < 3; i++)
    {
        pos(0) = 0.0;
        pos(1) = 0.0;
        s = 0.0;
        ds = 0.01;
        string name = "other_streamlines_";
        name += to_string(i + 7);
        name += ".csv";
        myfile11.open(name.c_str());

        r = (i + 1) * dy;
        start(0) = x_start;
        start(1) = y_start - r;
        cout << i + 1 << endl
             << endl;
        cout << start[1] << endl;
        while (pos(0) >= x_lower_lim && pos(0) <= x_upper_lim && fabs(s) < 5.0)
        {
            pos = rk_method(start, ds, x_pp, y_pp, x, y, G, Vinf, alfa);
            start = pos;
            s = s + ds;
            myfile11 << pos(0) << "\t" << pos(1) << endl;
        }
        myfile11.close();
    }
}

int main(int argc, char* argv[]) {
    // Check if enough command-line arguments are passed
    if (argc != 6) { // Change to != 5 to ensure exactly 4 inputs
        cerr << "Usage: .\\a.exe <ymc> <xmc> <tm> <trailing_edge_type>" << endl;
        return 1; // Exit with error code if the number of arguments is incorrect
    }

    // Parse command-line arguments
    double ymc = atof(argv[1]);
    double xmc = atof(argv[2]);
    double tm = atof(argv[3]);
    double trailing_edge_type = atof(argv[4]);
    double alfa= atof(argv[5]);


    cout << "Parsed values:" << endl;
    cout << "ymc: " << ymc << endl;
    cout << "xmc: " << xmc << endl;
    cout << "tm: " << tm << endl;
    cout << "trailing_edge_type: " << trailing_edge_type << endl;
    cout << "alfa: " << alfa << endl;

   
     int N;
    double x_nd;
    double delta_theta;
    double  Vinf;

    alfa = alfa* pi / 180.0;
    Vinf=5;

    N = 1001;
    delta_theta = (2 * pi) / (n - 1);

    VectorXd xu(N + 1);
    VectorXd yu(N + 1);
    VectorXd xl(N + 1);
    VectorXd yl(N + 1);
    VectorXd camber_x(N + 1);
    VectorXd camber_y(N + 1);

    VectorXd x_pp(n);
    VectorXd y_pp(n);
    VectorXd x_camber(n/2);
    VectorXd y_camber(n/2);
    VectorXd x_cp(n - 1);
    VectorXd y_cp(n - 1);

    VectorXd l(n - 1);   // l is the panel length which will differ from panel to panel.// l=(l_x)^2+(l_y)^2
    VectorXd l_x(n - 1); // (delta_x=(x_i+1)-(x_i)
    VectorXd l_y(n - 1); // delta_y =(y_i+1)-(y_i)


    double x, y;

    MatrixXd A(n, n);
    VectorXd B(n);
    VectorXd G(n);

    VectorXd a(2);
    VectorXd v(2);
    VectorXd v1(3);

    VectorXd y0(2);

    MatrixXd unit_normal(n - 1, 2);
    double section_lift_coefficient, sectional_moment_le, sectional_moment_qc;

    plot_airfoil(camber_x, camber_y, xu, xl, yu, yl, x_nd, N,ymc,xmc,tm,trailing_edge_type,alfa);
    nodal_coordinates(x_nd, delta_theta, x_pp, y_pp,ymc,xmc,tm,trailing_edge_type,alfa);
    controlpoints(x_pp, y_pp, x_cp, y_cp);
    Amatrix(A, x_cp, y_cp, x_pp, y_pp);
    random_point_flowfield(x_pp, y_pp);
    B_vector_particular_alfa(x_pp, y_pp, Vinf, alfa, B);
    Gamma_vector(G, B, A);
   // vel_any_point_flowfield_due_to_panels(x_pp,y_pp,G);
    tie(section_lift_coefficient, sectional_moment_le, sectional_moment_qc)= section_lift_and_moment_coefficient(G, x_pp, y_pp, alfa, Vinf);
    panel(l_x, l_y, l, x_pp, y_pp);
   // CL_vs_angle_of_attack(A, Vinf, l_x, l_y, l,x_pp,y_pp);
    normal_function_for_panels(unit_normal, l_x, l_y);
    plot_pressure(x_cp, y_cp, x_pp, y_pp, unit_normal, G, Vinf, alfa);
    plot_streamline(N, xu, yu, x_nd, x_pp, y_pp, x, y, G, Vinf,ymc,xmc,tm,trailing_edge_type,alfa);

    cout << "section_lift_coefficient=" << section_lift_coefficient << endl;
    cout << "sectional_moment_about_leading_edge=" << sectional_moment_le << endl;
    cout << "sectional_moment_about_quarter_chord=" << sectional_moment_qc << endl;
    

    return 0;
}

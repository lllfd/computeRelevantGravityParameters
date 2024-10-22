/// @Author       : linfd 3039562364@qq.com
/// @Date         : 2024-10-21 08:01:58
/// @LastEditTime : 2024-10-22 11:13:01
/// @FilePath     : \computeRelevantGravityParameters\src\computeRelevantGravityParameters.cpp
/// @Description  :

#include "computeRelevantGravityParameters.h"

int main()
{
    ReferenceEllipsoid WGS84("WGS84", 6378137.0, 1.0 / 298.257223563, 3.986004418e14, 6371000);
    WGS84.outputReferenceEllipsoidParameters();
    if (WGS84.readCmn_AndSmn_("../data/EGM96.cyh"))
    {
        ofstream ofsN("../out/N.txt");
        ofstream ofsT("../out/T.txt");
        ofstream ofsDeltaG("../out/DdeltaG.txt");
        ofstream ofsdeltaG("../out/ddeltaG.txt");

        double r = 6378137.0;
        for (double B = 21; B <= 26; B += 5 / 60.0)
        {
            vector<vector<double>> Pnm_;
            double theta = 90 - B;
            WGS84.computePmn_(theta, Pnm_);
            for (double L = 119; L <= 124; L += 5 / 60.0)
            {
                auto [X, Y, Z] = BLH2XYZ(DEG2RAD(B), DEG2RAD(L), 0, WGS84.getA(), WGS84.getF());
                auto [r, theta, lambda] = XYZ2R_THETA_LAMBDA(X, Y, Z);
                double phi = 0.5 * PI - theta;
                auto [T, delta_g, Delta_g, N] = WGS84.computeGeodeticGravity(theta, lambda, r, Pnm_);
                ofsT << format("{:15.8f} {:15.8f} {:20.8E}\n", B, L, T);
                ofsDeltaG << format("{:15.8f} {:15.8f} {:20.8E}\n", B, L, Delta_g);
                ofsdeltaG << format("{:15.8f} {:15.8f} {:20.8E}\n", B, L, delta_g);
                ofsN << format("{:15.8f} {:15.8f} {:20.8E}\n", B, L, N);
            }
        }
        clog << "Compute relevant gravity parameters successfully!" << endl;
        ofsT.close();
        ofsDeltaG.close();
        ofsdeltaG.close();
        ofsN.close();
    }
    return 0;
}

tuple<double, double, double> BLH2XYZ(double B, double L, double H, double a, double f)
{
    double ee = 2.0 * f - f * f;
    double W = sqrt(1.0 - ee * sin(B) * sin(B));
    double N = a / W;
    double X = (N + H) * cos(B) * cos(L);
    double Y = (N + H) * cos(B) * sin(L);
    double Z = (N * (1.0 - ee) + H) * sin(B);
    return make_tuple(X, Y, Z);
}

tuple<double, double, double> XYZ2R_THETA_LAMBDA(double X, double Y, double Z)
{
    double r = sqrt(X * X + Y * Y + Z * Z);
    double theta = atan2(sqrt(X * X + Y * Y), Z);
    double lambda = atan2(Y, X);
    return make_tuple(r, theta, lambda);
}

ReferenceEllipsoid::ReferenceEllipsoid(const string name, double a, double f, double GM, double omega) : m_name(name), m_a(a), m_f(f), m_GM(GM), m_omega(omega)
{
    m_b = m_a * (1.0 - m_f);
    m_c = m_a * m_a / m_b;
    m_e2 = 1.0 - pow(m_b, 2) / pow(m_a, 2);
    m_ee2 = pow(m_a, 2) / pow(m_b, 2) - 1.0;
    m_E = sqrt(pow(m_a, 2) - pow(m_b, 2));
    m_m = m_omega * m_omega * pow(m_a, 2) * m_b / m_GM;
    m_U0 = m_GM / m_E * atan(m_E / m_b) + 1.0 / 3.0 * pow(m_omega, 2) * pow(m_a, 2);
    double q0 = 1.0 / 2.0 * ((1.0 + 3.0 * pow(m_b, 2) / pow(m_E, 2)) * atan2(m_E, m_b) - 3.0 * m_b / m_E);
    double qq0 = 3.0 * (1.0 + pow(m_b, 2) / pow(m_E, 2)) * (1.0 - m_b / m_E * atan2(m_E, m_b)) - 1.0;
    m_J2_ = (pow(m_E, 2) / 3.0 / pow(m_a, 2) * (1.0 - 2.0 / 15.0 * m_m * sqrt(m_ee2) / q0)) / sqrt(5);
    m_garmaA = GM / (a * m_b) * (1.0 - m_m - m_m / 6.0 * sqrt(m_ee2) * qq0 / q0);
    m_garmaB = GM / pow(m_a, 2) * (1.0 + m_m / 3.0 * sqrt(m_ee2) * qq0 / q0);
}

double ReferenceEllipsoid::computeNormalGravityOnEllipsoidalSurface(double phi)
{
    return (m_a * m_garmaA * pow(cos(phi), 2) + m_b * m_garmaB * pow(sin(phi), 2)) / (sqrt((pow(m_a, 2) * pow(cos(phi), 2) + pow(m_b, 2) * pow(sin(phi), 2))));
}

double ReferenceEllipsoid::computeNormalGravityUpEllipsoidalSurface(double phi, double h)
{
    double garma = computeNormalGravityOnEllipsoidalSurface(phi);
    return garma * (1.0 - 2.0 / m_a * (1.0 + m_f + m_m - 2.0 * m_f * pow(sin(phi), 2)) * h + 3.0 / pow(m_a, 2) * pow(h, 2));
}

void ReferenceEllipsoid::outputReferenceEllipsoidParameters()
{
    ofstream ofsRefEllPara("../out/ReferenceEllipsoidParameters.txt");
    ofsRefEllPara << fixed << setprecision(4);
    ofsRefEllPara << setw(10) << "Name: " << setw(30) << m_name << endl;
    ofsRefEllPara << setw(10) << "a: " << setw(30) << m_a << "  m" << endl;
    ofsRefEllPara << setw(10) << "GM: " << setw(30) << m_GM << endl;
    ofsRefEllPara << setw(10) << "1/f: " << setw(30) << 1 / m_f << endl;
    ofsRefEllPara << setw(10) << "omega: " << setw(30) << m_omega << "  rad/s" << endl;
    ofsRefEllPara << setw(10) << "b: " << setw(30) << m_b << "  m" << endl;
    ofsRefEllPara << setw(10) << "E: " << setw(30) << m_E << "  m^2" << endl;
    ofsRefEllPara << setw(10) << "c: " << setw(30) << m_c << "  m^2" << endl;
    ofsRefEllPara << fixed << setprecision(14);
    ofsRefEllPara << setw(10) << "e^2: " << setw(30) << m_e2 << endl;
    ofsRefEllPara << setw(10) << "e'^2: " << setw(30) << m_ee2 << endl;
    ofsRefEllPara << fixed << setprecision(4);
    ofsRefEllPara << setw(10) << "U0: " << setw(30) << m_U0 << "  m/s^2" << endl;
    ofsRefEllPara << fixed << setprecision(10);
    ofsRefEllPara << setw(10) << "garmaA: " << setw(30) << m_garmaA << endl;
    ofsRefEllPara << setw(10) << "garmaB: " << setw(30) << m_garmaB << endl;
    ofsRefEllPara << setw(10) << "m: " << setw(30) << m_m << endl;
    ofsRefEllPara << scientific << setprecision(14);
    ofsRefEllPara << setw(10) << "J2_: " << setw(30) << m_J2_ << endl;
    ofsRefEllPara.close();
}

bool ReferenceEllipsoid::readCmn_AndSmn_(const string &fileName)
{
    m_Cnm_.resize(m_order + 1, std::vector<double>(m_order + 1, 0));
    m_Snm_.resize(m_order + 1, std::vector<double>(m_order + 1, 0));
    FILE *fp = fopen(fileName.c_str(), "r");
    if (fp == NULL)
    {
        cerr << "Open file " << fileName << " failed!" << endl;
        return false;
    }
    while (!feof(fp))
    {
        char line[1000];
        fgets(line, sizeof(line), fp);
        if (((string)line).find("END_OF_HEAD") != string::npos)
        {
            int m, n;
            double c, s, dc, ds;
            while (fscanf(fp, "%d %d %lf %lf %lf %lf", &m, &n, &c, &s, &dc, &ds) == 6)
            {
                m_Cnm_[m][n] = c;
                m_Snm_[m][n] = s;
            }
            clog << "Read normal gravity parameters successfully!" << endl;
            fclose(fp);
            return true;
        }
    }
    cerr << "Not Find END_OF_HEAD" << endl;
    return false;
}

bool ReferenceEllipsoid::computePmn_(double theta, vector<vector<double>> &Pnm_)
{
    auto getDiagonalLines = [](int l, double theta, double p) -> double
    { return sqrt((2 * l * 1.0 + 1) / (2 * l * 1.0)) * sin(theta) * p; };
    auto getVerticalValues = [](int l, int m, double theta, double p1, double p2) -> double
    { return sqrt((2 * l + 1) * 1.0 / (l * l - m * m)) * (sqrt(2 * l - 1) * cos(theta) * p1 - sqrt((l - m - 1) * 1.0 * (l + m - 1) / (2 * l - 3)) * p2); };

    Pnm_.resize(m_order + 1, std::vector<double>(m_order + 1, 0));
    Pnm_[0][0] = 1;
    Pnm_[1][1] = sin(theta) * sqrt(3);
    Pnm_[1][0] = cos(theta) * sqrt(3);
    for (int i = 2; i <= m_order; i++)
        Pnm_[i][i] = getDiagonalLines(i, theta, Pnm_[i - 1][i - 1]);
    for (int m = 0; m <= m_order - 1; m++)
        for (int n = m + 1 + (m == 0); n <= m_order; n++)
            Pnm_[n][m] = getVerticalValues(n, m, theta, Pnm_[n - 1][m], Pnm_[n - 2][m]);
    return true;
}

double ReferenceEllipsoid::computeGravityDisturbance(double theta, double lambda, double m_r, vector<vector<double>> Pnm_)
{
    double delta_g = 0;
    double sum2 = 0;
    for (int n = 0; n <= m_order; n++)
    {
        double sum1 = 0;
        for (int m = 0; m <= n; m++)
            sum1 += (m_Cnm_[n][m] * cos(m * lambda) + m_Snm_[n][m] * sin(m * lambda)) * Pnm_[n][m];
        sum2 += sum1 * (n + 1) * pow(m_a / m_r, n);
    }
    delta_g = m_GM * sum2 / pow(m_r, 2);
    return delta_g;
}

double ReferenceEllipsoid::computeGravityAnomaly(double theta, double lambda, double m_r, vector<vector<double>> Pnm_)
{
    double Delta_g = 0;
    double sum2 = 0;
    for (int n = 0; n <= m_order; n++)
    {
        double sum1 = 0;
        for (int m = 0; m <= n; m++)
            sum1 += (m_Cnm_[n][m] * cos(m * lambda) + m_Snm_[n][m] * sin(m * lambda)) * Pnm_[n][m];
        sum2 += sum1 * (n - 1) * pow(m_a / m_r, n);
    }
    Delta_g = m_GM * sum2 / pow(m_r, 2);
    return Delta_g;
}

double ReferenceEllipsoid::computePerturbationPosition(double theta, double lambda, double m_r, vector<vector<double>> Pnm_)
{
    double T = 0;
    double sum2 = 0;
    for (int n = 0; n <= m_order; n++)
    {
        double sum1 = 0;
        for (int m = 0; m <= n; m++)
            sum1 += (m_Cnm_[n][m] * cos(m * lambda) + m_Snm_[n][m] * sin(m * lambda)) * Pnm_[n][m];
        sum2 += sum1 * pow(m_a / m_r, n);
    }
    T = m_GM * sum2 / m_r;
    return T;
}

double ReferenceEllipsoid::computeGeoidDifference(double theta, double lambda, double m_r, vector<vector<double>> Pnm_)
{
    double N = 0;
    double garma = computeNormalGravityOnEllipsoidalSurface(0.5 / PI - theta);
    double sum2 = 0;
    for (int n = 0; n <= m_order; n++)
    {
        double sum1 = 0;
        for (int m = 0; m <= n; m++)
            sum1 += (m_Cnm_[n][m] * cos(m * lambda) + m_Snm_[n][m] * sin(m * lambda)) * Pnm_[n][m];
        sum2 += sum1 * pow(m_a / m_r, n);
    }
    N = m_GM * sum2 / m_r / garma;
    return N;
}

tuple<double, double, double, double> ReferenceEllipsoid::computeGeodeticGravity(double theta, double lambda, double m_r, vector<vector<double>> Pnm_)
{
    double garma = computeNormalGravityOnEllipsoidalSurface(0.5 * PI - theta);
    double T = 0, delta_g = 0, Delta_g = 0, N = 0;
    double sum2 = 0, sum3 = 0, sum4 = 0;
    for (int n = 0; n <= m_order; n++)
    {
        double sum1 = 0;
        for (int m = 0; m <= n; m++)
            sum1 += (m_Cnm_[n][m] * cos(m * lambda) + m_Snm_[n][m] * sin(m * lambda)) * Pnm_[n][m];
        sum2 += sum1 * pow(m_a / m_r, n);
        sum3 += sum1 * (n + 1) * pow(m_a / m_r, n);
        sum4 += sum1 * (n - 1) * pow(m_a / m_r, n);
    }
    T = m_GM * sum2 / m_r;
    delta_g = m_GM * sum3 / pow(m_r, 2);
    Delta_g = m_GM * sum4 / pow(m_r, 2);
    N = m_GM * sum2 / m_r / garma;
    return make_tuple(T, delta_g, Delta_g, N);
}
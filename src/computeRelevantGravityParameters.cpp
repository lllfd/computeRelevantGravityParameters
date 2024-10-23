/// @Author       : linfd 3039562364@qq.com
/// @Date         : 2024-10-21 08:01:58
/// @LastEditTime : 2024-10-23 09:17:29
/// @FilePath     : /computeRelevantGravityParameters/src/computeRelevantGravityParameters.cpp
/// @Description  : 计算重力场模型的相关参数

#include "computeRelevantGravityParameters.h"

int main()
{
    clog << "\nProgram started..." << endl;
    ReferenceEllipsoid WGS84("WGS84", 6378137.0, 1.0 / 298.257223563, 3.986004418e14, 7.29211500000E-5);
    WGS84.outputReferenceEllipsoidParameters();
    if (!WGS84.readCmn_AndSmn_("../data/EGM96.cyh"))
        return -1;

    clog << "\nComputing relevant gravity parameters..." << endl;
    ofstream ofsN("../out/N.txt");
    ofstream ofsT("../out/T.txt");
    ofstream ofsDeltaG("../out/DDeltaG.txt");
    ofstream ofsdeltaG("../out/deltaG.txt");

    for (double B = 21; B <= 26; B += 5 / 60.0)
    {
        auto [r, theta] = BH2RTheta(DEG2RAD(B), 0, WGS84.getA(), WGS84.getF());
        vector<vector<double>> Pnm_;
        WGS84.computePmn_(theta, Pnm_);
        for (double L = 119; L <= 124; L += 5 / 60.0)
        {
            double lambda = DEG2RAD(L);
            double phi = 0.5 * PI - theta;
            auto [T, delta_g, Delta_g, N] = WGS84.computeGeodeticGravity(theta, lambda, r, Pnm_);
            ofsT << format("{:15.8f} {:15.8f} {:20.8E}\n", B, L, T);
            ofsDeltaG << format("{:15.8f} {:15.8f} {:20.8E}\n", B, L, Delta_g * 1000);
            ofsdeltaG << format("{:15.8f} {:15.8f} {:20.8E}\n", B, L, delta_g * 1000);
            ofsN << format("{:15.8f} {:15.8f} {:20.8E}\n", B, L, N);
        }
    }
    clog << "\nCompute relevant gravity parameters successfully!" << endl;
    ofsT.close();
    ofsDeltaG.close();
    ofsdeltaG.close();
    ofsN.close();

    clog << "\nData visualizing..." << endl;
    _putenv("PYTHONHOME=");
    Py_Initialize();
    PyObject *pName = PyUnicode_DecodeFSDefault(/* pythonFileName */ "draw");
    PyObject *pModule = PyImport_Import(pName);
    PyObject *pFunc = PyObject_GetAttrString(pModule, /* funcName */ "drawInCpp");
    PyObject_CallObject(pFunc, NULL);
    Py_DECREF(pName);
    Py_DECREF(pModule);
    Py_DECREF(pFunc);
    Py_Finalize();
    clog << "\nData visualize successfully!" << endl;
    clog << "\nProgram finished!" << endl;

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

tuple<double, double, double> XYZ2RThetaLambda(double X, double Y, double Z)
{
    double r = sqrt(X * X + Y * Y + Z * Z);
    double theta = atan2(sqrt(X * X + Y * Y), Z);
    double lambda = atan2(Y, X);
    return make_tuple(r, theta, lambda);
}

tuple<double, double> BH2RTheta(double B, double H, double a, double f)
{
    double ee = 2.0 * f - f * f;
    double W = sqrt(1.0 - ee * sin(B) * sin(B));
    double N = a / W;
    double sum_xx_yy = (pow((N + H) * cos(B), 2));
    double Z = (N * (1.0 - ee) + H) * sin(B);
    double r = sqrt(sum_xx_yy + Z * Z);
    double theta = atan2(sqrt(sum_xx_yy), Z);
    return make_tuple(r, theta);
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
    ofsRefEllPara << format("{:<10} {:<30}\n", "Name:", m_name)
                  << format("{:<10} {:<30} m\n", "a:", m_a)
                  << format("{:<10} {:<30}\n", "GM:", m_GM)
                  << format("{:<10} {:<30}\n", "1/f:", 1 / m_f)
                  << format("{:<10} {:<30} rad/s\n", "omega:", m_omega)
                  << format("{:<10} {:<30} m\n", "b:", m_b)
                  << format("{:<10} {:<30} m^2\n", "E:", m_E)
                  << format("{:<10} {:<30} m^2\n", "c:", m_c)
                  << format("{:<10} {:<30}\n", "e^2:", m_e2)
                  << format("{:<10} {:<30}\n", "e'^2:", m_ee2)
                  << format("{:<10} {:<30} m/s^2\n", "U0:", m_U0)
                  << format("{:<10} {:<30}\n", "garmaA:", m_garmaA)
                  << format("{:<10} {:<30}\n", "garmaB:", m_garmaB)
                  << format("{:<10} {:<30}\n", "m:", m_m)
                  << format("{:<10} {:<30}\n", "J2_:", m_J2_);
    ofsRefEllPara.close();
}

bool ReferenceEllipsoid::readCmn_AndSmn_(const string &fileName)
{
    m_Cnm_.resize(m_order + 1, std::vector<double>(m_order + 1, 0));
    m_Snm_.resize(m_order + 1, std::vector<double>(m_order + 1, 0));
    FILE *fp = fopen(fileName.c_str(), "r");
    if (fp == NULL)
    {
        cerr << "\nOpen file " << fileName << " failed!" << endl;
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
            clog << "\nRead normal gravity parameters successfully!" << endl;
            fclose(fp);
            return true;
        }
    }
    cerr << "\nNot Find END_OF_HEAD" << endl;
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

double ReferenceEllipsoid::computeGravityDisturbance(double theta, double lambda, double r, vector<vector<double>> Pnm_)
{
    vector<vector<double>> Cnm_T = m_Cnm_;
    vector<vector<double>> Snm_T = m_Snm_;
    for (int n = 0; n <= Cnm_T.size(); n += 2)
        Cnm_T[n][0] = 0;
    double delta_g = 0;
    double sum2 = 0;
    for (int n = 0; n <= m_order; n++)
    {
        double sum1 = 0;
        for (int m = 0; m <= n; m++)
            sum1 += (Cnm_T[n][m] * cos(m * lambda) + Snm_T[n][m] * sin(m * lambda)) * Pnm_[n][m];
        sum2 += sum1 * (n + 1) * pow(m_a / r, n);
    }
    delta_g = m_GM * sum2 / pow(r, 2);
    return delta_g;
}

double ReferenceEllipsoid::computeGravityAnomaly(double theta, double lambda, double r, vector<vector<double>> Pnm_)
{
    vector<vector<double>> Cnm_T = m_Cnm_;
    vector<vector<double>> Snm_T = m_Snm_;
    for (int n = 0; n <= Cnm_T.size(); n += 2)
        Cnm_T[n][0] = 0;
    double Delta_g = 0;
    double sum2 = 0;
    for (int n = 0; n <= m_order; n++)
    {
        double sum1 = 0;
        for (int m = 0; m <= n; m++)
            sum1 += (Cnm_T[n][m] * cos(m * lambda) + Snm_T[n][m] * sin(m * lambda)) * Pnm_[n][m];
        sum2 += sum1 * (n - 1) * pow(m_a / r, n);
    }
    Delta_g = m_GM * sum2 / pow(r, 2);
    return Delta_g;
}

double ReferenceEllipsoid::computePerturbationPosition(double theta, double lambda, double r, vector<vector<double>> Pnm_)
{
    vector<vector<double>> Cnm_T = m_Cnm_;
    vector<vector<double>> Snm_T = m_Snm_;
    for (int n = 0; n <= Cnm_T.size(); n += 2)
        Cnm_T[n][0] = 0;
    double T = 0;
    double sum2 = 0;
    for (int n = 0; n <= m_order; n++)
    {
        double sum1 = 0;
        for (int m = 0; m <= n; m++)
            sum1 += (Cnm_T[n][m] * cos(m * lambda) + Snm_T[n][m] * sin(m * lambda)) * Pnm_[n][m];
        sum2 += sum1 * pow(m_a / r, n);
    }
    T = m_GM * sum2 / r;
    return T;
}

double ReferenceEllipsoid::computeGeoidDifference(double theta, double lambda, double r, vector<vector<double>> Pnm_)
{
    vector<vector<double>> Cnm_T = m_Cnm_;
    vector<vector<double>> Snm_T = m_Snm_;
    for (int n = 0; n <= Cnm_T.size(); n += 2)
        Cnm_T[n][0] = 0;
    double N = 0;
    double garma = computeNormalGravityOnEllipsoidalSurface(0.5 / PI - theta);
    double sum2 = 0;
    for (int n = 0; n <= m_order; n++)
    {
        double sum1 = 0;
        for (int m = 0; m <= n; m++)
            sum1 += (Cnm_T[n][m] * cos(m * lambda) + Snm_T[n][m] * sin(m * lambda)) * Pnm_[n][m];
        sum2 += sum1 * pow(m_a / r, n);
    }
    N = m_GM * sum2 / r / garma;
    return N;
}

tuple<double, double, double, double> ReferenceEllipsoid::computeGeodeticGravity(double theta, double lambda, double r, vector<vector<double>> Pnm_)
{
    vector<vector<double>> Cnm_T = m_Cnm_;
    vector<vector<double>> Snm_T = m_Snm_;
    for (int n = 0; n <= Cnm_T.size(); n += 2)
        Cnm_T[n][0] = 0;
    double garma = computeNormalGravityOnEllipsoidalSurface(0.5 * PI - theta);
    double T = 0, delta_g = 0, Delta_g = 0, N = 0;
    double sum2 = 0, sum3 = 0, sum4 = 0;
    for (int n = 0; n <= m_order; n++)
    {
        double sum1 = 0;
        for (int m = 0; m <= n; m++)
            sum1 += (Cnm_T[n][m] * cos(m * lambda) + Snm_T[n][m] * sin(m * lambda)) * Pnm_[n][m];
        sum2 += sum1 * pow(m_a / r, n);
        sum3 += sum1 * (n + 1) * pow(m_a / r, n);
        sum4 += sum1 * (n - 1) * pow(m_a / r, n);
    }
    T = m_GM * sum2 / r;
    delta_g = m_GM * sum3 / pow(r, 2);
    Delta_g = m_GM * sum4 / pow(r, 2);
    N = m_GM * sum2 / r / garma;
    return make_tuple(T, delta_g, Delta_g, N);
}
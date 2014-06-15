#include <HIntLib/integrand.h>
#include <HIntLib/integrator.h>
#include <HIntLib/mcintegrator.h>
#include <HIntLib/mersennetwister.h>
#include <HIntLib/qmcintegrator.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/esterr.h>
#include <HIntLib/digitalnet2.h>
#include <HIntLib/sobolmatrix.h>
#include <HIntLib/niederreitermatrix.h>
#include <HIntLib/defaults.h>
#include <HIntLib/mcpointset.h>
#include <interceptableintegrand.h>
#include <rqmcintegrator.h>
#include <qintintegrator.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

using namespace HIntLib;
using namespace std;

class TestFunction : public Integrand
{
public:
    TestFunction(int s) : Integrand(s) {}
    virtual ~TestFunction() {}
    virtual real operator() (const real *x)
    {
      vector<real> v(x, x + this->getDimension());
      real sum = 0;
      for (auto n : v) sum += n;
      return sum;
    }
};

class TestFunctionInterceptable : public InterceptableIntegrand
{
public:
    TestFunctionInterceptable(int s) : InterceptableIntegrand(s) {}
    virtual ~TestFunctionInterceptable() {}
    virtual InterceptedSet intercept(const real *x)
    {
        vector<real> v(x, x + this->getDimension());
        real sum = 0;
        for (auto n : v) sum += n;
        InterceptedSet res {sum, v};
        return res;
    }
};

class SequenceInterceptor : public Integrand
{
private:
    ofstream f;
public:
    SequenceInterceptor(int s, const string filename = string("seq.txt")) :
        Integrand(s)
    {
        struct passwd *pw = getpwuid(getuid());
        const char *homedir = pw->pw_dir;
        string address = string(homedir) + "/" + filename;
        f.open(address);
    }
    virtual ~SequenceInterceptor()
    {
        f.close();
    }
    virtual real operator() (const real *x)
    {
        vector<real> v(x, x + this->getDimension());
        string s;
        for (auto n : v)
        {
            s += to_string(n);
            s += ", ";
        }
        s.pop_back();
        s.pop_back();
        f << s << "\n";
        return 0;
    }
};

void calculateIntegral(Integrand& f, Integrator& integrator)
{
    Hypercube h(f.getDimension());
    EstErr ee;
    auto result = integrator.integrate(f, h, 10, 0, 0, ee);
    std::cout << "Estimate: " << ee.getEstimate() << "(+/-)" << ee.getError() << " "
              << "Result: " << (result == Integrator::Status::MAX_EVAL_REACHED ? "OK" : "Not OK") << "\n";
}

int main()
{
    int s=10;
    int rc=2;
    int seed=42;
    TestFunction f(s);
    SequenceInterceptor si(s);
    TestFunctionInterceptable fi(s);

    MonteCarloPointSet<MersenneTwister> ps_mc;
    MCIntegrator integrator_mc(&ps_mc);
    integrator_mc.randomize(seed);
    integrator_mc.setMinNumEval(300);

    SobolMatrix matrix_sobol;
    DigitalSeq2PointSet<real> ps_sobol(matrix_sobol, true);
    RQMCIntegrator integrator_sobol(&ps_sobol, rc, seed);
    QINTIntegrator integrator_sobol_qint(&ps_sobol, rc, seed);

    //calculateIntegral(si, integrator_sobol);
    //calculateIntegral(f, integrator_mc);
    calculateIntegral(f, integrator_sobol);
    //calculateIntegral(f, integrator_nied);
    calculateIntegral(fi, integrator_sobol_qint);

}

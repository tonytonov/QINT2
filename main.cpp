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
#include <utils.h>
#include <testcollection.h>

using namespace HIntLib;
using namespace std;

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

void runMethodComparison(InterceptableIntegrand& f, std::vector<Integrator*>& integratorList, int maxEval)
{
    Hypercube h(f.getDimension());
    EstErr ee;
    std::cout << "Integral: " << f.getDimension() << " dimensions, " << "exact value " << f.getExactValue() << "\n";
    for (auto integrator : integratorList)
    {
        auto status = integrator->integrate(f, h, maxEval, 0, 0, ee);
        std::cout << "Estimate: " << ee.getEstimate() << "(+/-)" << ee.getError() << "; "
                  << "Result: " << (status == Integrator::Status::MAX_EVAL_REACHED ? "OK" : "Not OK") << "\n";
    }
}

int main()
{
    int s=5;
    int rc=20;
    int sparam=8;
    int seed=42;
    int maxEval=10000;

    //FI00_simpleSum f(s);
    FI01_fMorCaf f(s);

    MonteCarloPointSet<MersenneTwister> ps_mc;
    MCIntegrator integrator_mc(&ps_mc);
    integrator_mc.randomize(seed);

    SobolMatrix matrix_sobol;
    DigitalSeq2PointSet<real> ps_sobol(matrix_sobol, true);
    RQMCIntegrator integrator_sobol(&ps_sobol, rc, seed);
    QintIntegrator integrator_sobol_qint(&ps_sobol, rc, sparam, seed);

    std::vector<Integrator*> integratorList
    {
        &integrator_mc, &integrator_sobol, &integrator_sobol_qint
    };
    runMethodComparison(f, integratorList, maxEval);
}

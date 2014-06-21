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
#include <sqlite3.h>
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
    for (const auto integrator : integratorList)
    {
        auto status = integrator->integrate(f, h, maxEval, 0, 0, ee);
        std::cout << "Estimate: " << ee.getEstimate() << "(+/-)" << ee.getError() << "; "
                  << "Result: " << (status == Integrator::Status::MAX_EVAL_REACHED ? "OK" : "Not OK") << "\n";
    }
}

void writeMethodComparison(InterceptableIntegrand& f, std::vector<Integrator*>& integratorList,
                           int maxEval, int randCount, int sParamQint, int seed)
{
    Hypercube h(f.getDimension());
    EstErr ee;
    sqlite3 *db;
    //TODO : do smth with file locations and rc
    int dbStatus = sqlite3_open("../qint.sqlite", &db);
    std::string createDBQuery =
            R"sql(CREATE TABLE IF NOT EXISTS results(
            function TEXT NOT NULL,
            dim INT NOT NULL,
            exactval REAL,
            method TEXT NOT NULL,
            maxeval INT NOT NULL,
            status BOOLEAN NOT NULL,
            estimate REAL NOT NULL,
            stddev REAL NOT NULL,
            randcount INT,
            sparam INT,
            seed INT NOT NULL)
            )sql";
    for (const auto integrator : integratorList)
    {
        auto status = integrator->integrate(f, h, maxEval, 0, 0, ee);
        std::string methodName = "Undefined method";
        int rc = -1;
        int sparam = -1;
        if (dynamic_cast<MCIntegrator*>(integrator))
        {
            methodName = "MC";
        }
        if (dynamic_cast<RQMCIntegrator*>(integrator))
        {
            methodName = "RQMC";
            rc = randCount;
        }
        if (dynamic_cast<QintIntegrator*>(integrator))
        {
            methodName = "Qint";
            rc = randCount;
            sparam = sParamQint;
        }
        std::string addResultsQuery =
                "INSERT INTO results("
                "function, dim, exactval, method, "
                "maxeval, status, estimate, stddev, "
                "randcount, sparam, seed) VALUES"
                "("
                "\'" + f.name + "\'," +
                S(f.getDimension()) + "," +
                S(f.getExactValue()) + "," +
                "\'" + methodName + "\'," +
                S(maxEval) + "," +
                S(status == Integrator::Status::MAX_EVAL_REACHED) + "," +
                S(ee.getEstimate()) + "," +
                S(ee.getError()) + "," +
                S(rc) + "," +
                S(sparam) + "," +
                S(seed) + ")";
        sqlite3_exec(db, addResultsQuery.c_str(), 0, 0, 0);
    }

    sqlite3_exec(db, createDBQuery.c_str(), 0, 0, 0);
    if (db)
    {
      sqlite3_close(db);
    }
}

int main()
{
    int s=5;
    int rc=16;
    int sparam=0;
    int seed=42;
    int maxEval=std::pow(2, 13);

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
    //runMethodComparison(f, integratorList, maxEval);
    writeMethodComparison(f, integratorList, maxEval, rc, sparam, seed);
}

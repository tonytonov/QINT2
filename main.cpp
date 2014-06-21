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

void writeMethodComparison(InterceptableIntegrand& f, std::vector<Integrator*>& integratorList, int maxEval)
{
    Hypercube h(f.getDimension());
    EstErr ee;
    sqlite3 *db;
    //TODO : do smth with file locations and rc
    int rc = sqlite3_open("../qint.sqlite", &db);
    char* createDBQuery = R"sql(
            CREATE TABLE IF NOT EXISTS results(
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
        char* addResultsQuery = R"sql(
                INSERT INTO results(
                function, dim, exactval, method,
                maxeval, status, estimate, stddev,
                randcount, sparam, seed) VALUES
                (
                    'test', 99, 0.1, 'somemethod',
                    999, 1, 42, 24,
                    333, 2, 911
                )
                )sql";
        sqlite3_exec(db, addResultsQuery, 0, 0, 0);
    }

    sqlite3_exec(db, createDBQuery, 0, 0, 0);
    if (db)
    {
      sqlite3_close(db);
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
    //runMethodComparison(f, integratorList, maxEval);
    writeMethodComparison(f, integratorList, maxEval);
}

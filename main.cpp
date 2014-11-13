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

class SequenceInterceptor : public Integrand
{
private:
    std::ofstream f;
public:
    SequenceInterceptor(int s, const std::string filename = std::string("seq.txt")) :
        Integrand(s)
    {
        struct passwd *pw = getpwuid(getuid());
        const char *homedir = pw->pw_dir;
        std::string address = std::string(homedir) + "/" + filename;
        f.open(address);
    }
    virtual ~SequenceInterceptor()
    {
        f.close();
    }
    virtual real operator() (const real *x)
    {
        std::vector<real> v(x, x + this->getDimension());
        std::string s;
        for (auto n : v)
        {
            s += std::to_string(n);
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
    int dbStatus = sqlite3_open("../qint_reborn.sqlite", &db);
    std::string createDBQuery =
            R"sql(CREATE TABLE IF NOT EXISTS results(
            hash INT PRIMARY KEY,
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
            sqlite3_exec(db, createDBQuery.c_str(), 0, 0, 0);
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
        if (methodName != "Qint" && sParamQint > 0) continue;
        if (std::isnan(ee.getEstimate()) || std::isnan(ee.getError())) continue;
        std::string exportRecord =
                "\'" + f.name + "\'," +
                S(f.getDimension()) + "," +
                S(f.getExactValue()) + "," +
                "\'" + methodName + "\'," +
                S(maxEval) + "," +
                S(std::log2(maxEval) - std::log2(randCount) == sParamQint
                  ? 2 : status != Integrator::Status::MAX_EVAL_REACHED) + "," +
                S(ee.getEstimate()) + "," +
                S(ee.getError()) + "," +
                S(rc) + "," +
                S(sparam) + "," +
                S(seed);
        auto hashRecord = std::hash<std::string>()(exportRecord);
        std::string addResultsQuery =
                "INSERT INTO results("
                "hash, function, dim, exactval, "
                "method, maxeval, status, estimate, "
                "stddev, randcount, sparam, seed) VALUES"
                "(" + S(hashRecord) + "," + exportRecord + ")";
        sqlite3_exec(db, addResultsQuery.c_str(), 0, 0, 0);
    }

    if (db)
    {
        sqlite3_close(db);
    }
}

int main()
{
    std::vector<int> s {1, 2, 3, 4, 5, 10, 15, 20, 30};
    std::vector<int> sparam(8); // sparam(k) is 0:k
    std::iota(std::begin(sparam), std::end(sparam), 0);
    int seed = 42;
    int rc = std::pow(2, 4);
    int maxEval = std::pow(2, 1 + *std::max_element(sparam.begin(), sparam.end())) * rc;

    //    SequenceInterceptor x(45);
    //    SobolMatrix matrix_sobol;
    //    DigitalSeq2PointSet<real> ps_sobol(matrix_sobol, true);
    //    RQMCIntegrator integrator_sobol(&ps_sobol, rc, seed);
    //    Hypercube h(x.getDimension());
    //    EstErr ee;
    //    integrator_sobol.integrate(x, h, maxEval, 0, 0, ee);

    for (auto i_s : s)
    {
        for (auto i_sparam : sparam)
        {
            FI05_CubicPolynomial f(i_s);

            MonteCarloPointSet<MersenneTwister> ps_mc;
            MCIntegrator integrator_mc(&ps_mc);
            integrator_mc.randomize(seed);

            SobolMatrix matrix_sobol;
            DigitalSeq2PointSet<real> ps_sobol(matrix_sobol, true);
            RQMCIntegrator integrator_sobol(&ps_sobol, rc, seed);
            QintIntegrator integrator_sobol_qint(&ps_sobol, rc, i_sparam, seed, 1);
            QintIntegrator integrator_sobol_qint_mc(&ps_sobol, rc, i_sparam, seed, 2);
            QintIntegrator integrator_sobol_qint_rqmc(&ps_sobol, rc, i_sparam, seed, 3);

            std::vector<Integrator*> integratorList
            {
                &integrator_mc,
                &integrator_sobol,
                &integrator_sobol_qint,
//                &integrator_sobol_qint_mc,
//                &integrator_sobol_qint_rqmc
            };
            std::cout << std::to_string(i_s) << " " << std::flush;
            //runMethodComparison(f, integratorList, maxEval);
            writeMethodComparison(f, integratorList, maxEval, rc, i_sparam, seed);
        }
    }
}

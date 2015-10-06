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
#include <stdexcept>
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
    int dbStatus = sqlite3_open_v2("../qint.sqlite", &db, SQLITE_OPEN_CREATE | SQLITE_OPEN_READWRITE, NULL);
    if (dbStatus) throw std::runtime_error("Database cannot be opened!");
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
        auto status = integrator->integrate(f, h, maxEval, 0, 0, ee);
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
    std::vector<int> dim {1, 2, 3, 4, 5, 10, 15, 20, 30};
    std::vector<int> limit {8, 9, 10, 11, 12, 13, 14, 15, 16};
    std::vector<int> randCount {1, 16};
    std::vector<int> seed(15);
    std::iota(std::begin(seed), std::end(seed), 0);

    //    SequenceInterceptor x(45);
    //    SobolMatrix matrix_sobol;
    //    DigitalSeq2PointSet<real> ps_sobol(matrix_sobol, true);
    //    RQMCIntegrator integrator_sobol(&ps_sobol, rc, seed);
    //    Hypercube h(x.getDimension());
    //    EstErr ee;
    //    integrator_sobol.integrate(x, h, maxEval, 0, 0, ee);

    int progress = 0;
    int total = seed.size() * dim.size() * limit.size() * randCount.size();
    int barWidth = 70;

    for (auto i_rc : randCount)
    {
        for (auto i_seed : seed)
        {
            for (auto i_dim : dim)
            {
                for (auto i_lim : limit)
                {
                    std::vector<int> sparam(i_lim + 1); // sparam(k) is 0:k
                    std::iota(std::begin(sparam), std::end(sparam), 0);
                    int maxEval = std::pow(2, i_lim);
                    for (auto i_sparam : sparam)
                    {
                        FI07_Singular f(i_dim);

                        MonteCarloPointSet<MersenneTwister> ps_mc;
                        MCIntegrator integrator_mc(&ps_mc);
                        integrator_mc.randomize(i_seed);

                        SobolMatrix matrix_sobol;
                        DigitalSeq2PointSet<real> ps_sobol(matrix_sobol, true);
                        RQMCIntegrator integrator_sobol(&ps_sobol, i_rc, i_seed);
                        QintIntegrator integrator_sobol_qint(&ps_sobol, i_rc, i_sparam, i_seed, 1);
                        //QintIntegrator integrator_sobol_qint_mc(&ps_sobol, rc, i_sparam, i_seed, 2);
                        //QintIntegrator integrator_sobol_qint_rqmc(&ps_sobol, rc, i_sparam, i_seed, 3);

                        std::vector<Integrator*> integratorList
                        {
                            &integrator_mc,
                                    &integrator_sobol,
                                    &integrator_sobol_qint,
                                    //&integrator_sobol_qint_mc,
                                    //&integrator_sobol_qint_rqmc
                        };
                        //runMethodComparison(f, integratorList, maxEval);
                        writeMethodComparison(f, integratorList, maxEval, i_rc, i_sparam, i_seed);
                    }
                    //update progress bar
                    ++progress;
                    std::cout << "[";
                    int pos = barWidth * progress / total;
                    for (int i = 0; i < barWidth; ++i) {
                        if (i < pos) std::cout << "=";
                        else if (i == pos) std::cout << ">";
                        else std::cout << " ";
                    }
                    std::cout << "] " << int(1.0 * progress / total * 100) << " %\r";
                    std::cout.flush();
                }
            }
        }
    }
}

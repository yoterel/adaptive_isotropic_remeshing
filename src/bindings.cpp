#include <npe.h>
#include "adaptive_remesh_botsch.h"

npe_function(adaptive_remesh_botsch)
npe_arg(v, dense_double)
npe_arg(f, dense_int)
npe_arg(e, double)
npe_arg(i, int)
npe_arg(p, bool)
npe_arg(a, bool)

npe_doc("Performs remeshing of a mesh v,f using (e)psilon for i (i)terations, possibly (p)rojecting output to input mesh, and doing so (a)daptively.")

npe_begin_code()
    Eigen::MatrixXd V(v);
    Eigen::MatrixXi F(f);
    double epsilon(e);
    int ierations(i);
    bool project(p);
    bool adaptive(a);
    adaptive_remesh_botsch(V, F, epsilon, ierations, project, adaptive);
    return std::make_tuple(npe::move(V), npe::move(F));
npe_end_code()
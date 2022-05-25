#ifndef REMESH_BOTSCH
#define REMESH_BOTSCH



#include <Eigen/Core>

void adaptive_remesh_botsch(Eigen::MatrixXd & V, Eigen::MatrixXi & F, double epsilon_in, int iters, bool project, bool adaptive);

#endif

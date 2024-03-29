#ifndef TANGENTIAL_RELAXATION
#define TANGENTIAL_RELAXATION



#include <Eigen/Core>

void tangential_smoothing(Eigen::MatrixXd & V,Eigen::MatrixXi & F, Eigen::VectorXi & feature,
Eigen::MatrixXd & V0 ,Eigen::MatrixXi & F0, Eigen::VectorXd & lambda, Eigen::VectorXd & sizingField);


#endif

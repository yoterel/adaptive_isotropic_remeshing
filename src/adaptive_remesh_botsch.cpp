#include "equalize_valences.h"
#include "collapse_edges.h"
#include "tangential_smoothing.h"
#include <igl/is_edge_manifold.h>
#include <igl/writeOBJ.h>
#include "split_edges_until_bound.h"
#include <igl/unique_edge_map.h>
#include <igl/edge_flaps.h>
#include <igl/circulation.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/cotmatrix.h>
#include <igl/avg_edge_length.h>
#include <igl/point_mesh_squared_distance.h>
#include <iostream>

void calc_sizing_field(Eigen::MatrixXd &V, Eigen::MatrixXi &F, double epsilon, Eigen::VectorXd &sizingField, bool adaptive)
{
	if (adaptive)
	{
		Eigen::VectorXd epsilon_vec = Eigen::VectorXd::Constant(V.rows(), epsilon);
		Eigen::MatrixXd K_tot; // maximum absolute curvature
		Eigen::VectorXd K;	   // gaussian curvature
		// Compute integral of Gaussian curvature
		igl::gaussian_curvature(V, F, K);
		// Compute mass matrix
		Eigen::SparseMatrix<double> L, M, Minv;
		igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
		igl::invert_diag(M, Minv);
		// Divide by area to get integral average
		K = (Minv * K).eval();
		// mean curvature
		Eigen::MatrixXd HN, H;
		igl::cotmatrix(V, F, L);
		HN = -Minv * (L * V);
		H = HN.rowwise().norm().array() / 4; // up to sign.
		Eigen::VectorXd delta = H.array().square() - K.array();
		for (int i = 0; i < delta.size(); i++)  // clip delta to 0 to avoid nans
		{
			if (delta[i] < 0)
			{
				delta(i) = 0;
			}
		}
		K_tot = H.array() + delta.array().sqrt();
		K_tot = K_tot.matrix();
		sizingField = (6 * epsilon_vec.array() * K_tot.array().inverse() - 3 * epsilon_vec.array().square()).sqrt();  // this could potentially become nan. Authors did not discuss this.
		// std::cout << "max: " << sizingField.array().isNaN().select(0,sizingField).maxCoeff() << "\n";
		// std::cout << "min: " << sizingField.array().isNaN().select(0,sizingField).minCoeff() << "\n";
		double my_max = 4;
		double my_min = 0.01;
		for (int i = 0; i < sizingField.size(); i++)  // clip sizingfield if it is below or above values
		{
			if ((sizingField(i) > my_max)) // || std::isnan(sizingField(i))
			{
				sizingField(i) = my_max;
			}
			if ((sizingField(i) < my_min))
			{
				sizingField(i) = my_min;
			}
			// std::cout << "found nan in sizing field. i: " << i << "\n";
			// std::cout << "ktot: " << K_tot.array()(i) << "\n";
			// std::cout << "epsilon: " << epsilon.array()(i) << "\n";
			// std::cout << "H: " << H.array()(i) << "\n";
			// std::cout << "K: " << K.array()(i) << "\n";
		}
		sizingField = sizingField.matrix();
	}
	else
	{
		sizingField = Eigen::VectorXd::Constant(V.rows(), epsilon);
	}

}

void adaptive_remesh_botsch(Eigen::MatrixXd &V, Eigen::MatrixXi &F, double epsilon, int iters, bool project, bool adaptive)
{
	Eigen::MatrixXd V0;
	Eigen::MatrixXi F0;

	Eigen::VectorXd high, low, lambda, sizingField;
	// high = 1.4*target;
	// low = 0.7*target;
	F0 = F;
	V0 = V;
	Eigen::VectorXi feature;
	feature.resize(0);
	// Iterate the four steps
	for (int i = 0; i < iters; i++)
	{
		std::cout << "iter: " << i << "/" << iters << "\n";
		calc_sizing_field(V, F, epsilon, sizingField, adaptive);  // set target length per edge using sizingfield
		high = 1.4 * sizingField;
		low = 0.7 * sizingField;
		std::cout << "splitting...\n";
		split_edges_until_bound(V, F, feature, high, low); // Split
		calc_sizing_field(V, F, epsilon, sizingField, adaptive);
		high = 1.4 * sizingField;
		low = 0.7 * sizingField;
		std::cout << "collapsing...\n";
		collapse_edges(V, F, feature, high, low); // Collapse
		std::cout << "flipping...\n";
		equalize_valences(V, F, feature); // Flip
		int n = V.rows();
		lambda = Eigen::VectorXd::Constant(n, 1.0);
		std::cout << "smoothing...\n";
		tangential_smoothing(V, F, feature, V0, F0, lambda, sizingField); // Smooth
		if (project)
		{
 			std::cout << "projecting..." << "\n";
			Eigen::VectorXi sqrI;
    		Eigen::VectorXd sqrD;
    		Eigen::MatrixXd V_projected;
    		igl::point_mesh_squared_distance(V, V0, F0, sqrD, sqrI, V_projected);  // Project
    		V = V_projected;
		}
		std::cout << "finished iter, faces amount: " << F.rows() << "\n";
	}
}

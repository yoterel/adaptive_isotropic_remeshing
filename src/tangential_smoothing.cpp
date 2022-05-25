#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/avg_edge_length.h>
#include <igl/massmatrix.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/barycenter.h>
#include <igl/pinv.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>
#include <Eigen/SparseCore>
#include <igl/adjacency_list.h>
#include <igl/adjacency_matrix.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/avg_edge_length.h>
#include <igl/edge_flaps.h>
#include <igl/unique_edge_map.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/principal_curvature.h>
#include <igl/collapse_edge.h>
#include <igl/C_STR.h>
#include <igl/flip_edge.h>
#include <igl/remove_duplicate_vertices.h>


void tangential_smoothing(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &feature,
                           Eigen::MatrixXd &V0, Eigen::MatrixXi &F0, Eigen::VectorXd &lambda, Eigen::VectorXd &sizingField)
{
    using namespace Eigen;
    MatrixXd Q, P, N, V_fixed;
    VectorXd dblA;
    std::vector<std::vector<int>> A;
    Matrix3d I, NN;
    I.setIdentity();
    Eigen::MatrixXd SV;
    Eigen::MatrixXi SVI, SVJ;

    V_fixed = V;

    int n = V.rows();
    int m = F.rows();

    // igl::doublearea(V,F,dblA);

    // std::vector<double> vertex_areas;
    // vertex_areas.setZero(m);

    // for (int j = 0; j < m; j++) {
    //     vertex_areas[F(j,0)] = vertex_areas[F(j,0)] + (abs(dblA(j))/6);
    //     vertex_areas[F(j,1)] = vertex_areas[F(j,1)] + (abs(dblA(j))/6);
    //     vertex_areas[F(j,2)] = vertex_areas[F(j,2)] + (abs(dblA(j))/6);
    // }

    Eigen::MatrixXd N_before, N_after;
    igl::adjacency_list(F, A);

    int num_feat = feature.size();
    std::vector<bool> is_feature_vertex;
    is_feature_vertex.resize(n);

    for (int s = 0; s < num_feat; s++)
    {
        is_feature_vertex[feature(s)] = true;
    }

    Q.resize(n, 3);
    P.resize(n, 3);
    //           Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);
    std::vector<std::vector<int>> V2F, _1;
    igl::vertex_triangle_adjacency(n, F, V2F, _1);
    VectorXd area;
    MatrixXd BC;
    igl::doublearea(V, F, area);
    area = area.array() / 2;
    igl::barycenter(V, F, BC);
    for (int i = 0; i < n; i++)
    {
        // std::cout << "i: " << i << "/" << n << "\n";
        bool is_feature = is_feature_vertex[i];
        if (!is_feature)
        {
            Eigen::RowVector3d nom, q, p;
            double denom = 0;
            q.setZero();
            nom.setZero();
            p.setZero();
            double denominator = 0.0;
            /* mine */
            for (int j = 0; j < V2F[i].size(); j++)
            {
                double cur_denom = (area(V2F[i][j]) * (sizingField(F(V2F[i][j], 0)) + sizingField(F(V2F[i][j], 1)) + sizingField(F(V2F[i][j], 2))) / 3);
                denom = denom + cur_denom;
                nom = nom + (cur_denom * BC.row(V2F[i][j]).array()).matrix();
            }
            // std::cout << "V[0]: " << V.row(0) << "\n";
            // std::cout << "nom: " << nom << "\n";
            // std::cout << "denom: " << denom << "\n";
            q = nom.array() / denom;
            q = q.matrix();
            // std::cout << "q[0]: " << q << "\n";
            /* mine */
            /* theirs */
            // for (int j = 0; j < A[i].size(); j++)
            // {
            //     q = q + (V.row(A[i][j]) / A[i].size());
            //     // q = q + (V.row(A[i][j])*vertex_areas[A[i][j]]);
            //     // std::cout << q << std::endl;
            //     // denominator = denominator + vertex_areas[A[i][j]];
            // } // q is )( barycenter?
            /* theirs */
            // q = q/denominator;
            // N.row(i) = N.row(i)/N.row(i).norm();
            NN = lambda(i) * (Eigen::MatrixXd::Identity(3, 3) - N.row(i).transpose() * (N.row(i)));
            p = (V.row(i).transpose() - (NN * (V.row(i).transpose() - q.transpose()))).transpose();
            // p = q;
            // std::cout << N.row(i) << std::endl;

            V.row(i) = p;

            // igl::per_face_normals(V_projected,F,Eigen::Vector3d(0,0,0),N_after);
            //            for (int j = 0; j < m ; j++) {
            //                if (N_before.row(j).dot(N_after.row(j)) < 0) {
            //                    // std::cout << "Avoided face flipping, I think." << std::endl;
            //                    V.row(i) = V_fixed.row(i);
            //                }
            //            }
        }
    }
    //        igl::remove_duplicate_vertices(V,0,SV,SVI,SVJ);
    //        std::cout << V.rows()-SV.rows() << std::endl;
    //	igl::writeOBJ("pre-project.obj",V,F);
    //
    
    //	igl::writeOBJ("post-project.obj",V,F);
    //    igl::remove_duplicate_vertices(V,0,SV,SVI,SVJ);
    //    std::cout << V.rows()-SV.rows() << std::endl;
    //   std::cout << "not projecting!" << std::endl;
}
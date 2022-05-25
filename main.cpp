#include "version.h"
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include "src/adaptive_remesh_botsch.h"

const char *get_exec_name(char *argv0);
void print_usage(const char *progName);
void parse_args(int argindex, int argc, char *argv[], std::string &out, int &iterations, double &e, bool &project, bool &adaptive);

int main(int argc, char *argv[])
{
	const char *exec_name = get_exec_name(argv[0]);
	std::cout << exec_name << " version " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << std::endl;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	std::string in;
	std::string out = "output.obj";
	int iterations = 3;
	double e = 0.01;
	bool project = false;
	bool adaptive = true;
	int argindex = 0;
	if (argc > 1)
	{
		in = argv[1];
		argindex = 1;
	}
	else
	{
		print_usage(exec_name);
		return argindex;
	}
	parse_args(argindex, argc, argv, out, iterations, e, project, adaptive);
	igl::read_triangle_mesh(in, V, F);
	std::cout << "epsilon: " << e << "\n";
	std::cout << "iters: " << iterations << "\n";
	std::cout << std::boolalpha;
	std::cout << "project: " << project << "\n";
	std::cout << "input faces amount: " << F.rows() << "\n";
	adaptive_remesh_botsch(V, F, e, iterations, project, adaptive);
	std::cout << "output faces amount: " << F.rows() << "\n";
	igl::write_triangle_mesh(out, V, F);
}

const char *get_exec_name(char *argv0)
{
	const char *exename = strrchr(argv0, '/');
	if (exename)
		// skip past the last /
		++exename;
	else
		exename = argv0;
	return exename;
}

void print_usage(const char *progName)
{
	std::cout << progName << " [input_mesh.ext] [output_mesh.ext] [options] \n\n"
			  << "performs adaptive remeshing of a mesh, which must be closed and manifold, and ext can be either of obj, ply, off, stl or mesh.\n\n"
			  << "options:\n"
			  << "-i  | Number of iterations to run remeshing\n"
			  << "-e  | Epsilon (for adaptive) / Target Edge Length (not adaptive), controls the res of output\n"
			  << "-p  | Project result vertices of each iteration on original mesh\n"
			  << "-na | Not Adaptive, will use original botsch & kobbelt remeshing algorithm\n";
}

void parse_args(int argindex, int argc, char *argv[], std::string &out, int &iterations, double &e, bool &project, bool &adaptive)
{
	while (argindex + 1 < argc)
	{
		if (strncmp(argv[argindex + 1], "-i", 2) == 0)
		{
			iterations = atoi(argv[argindex + 2]);
			argindex = argindex + 2;
		}
		else
		{
			if (strncmp(argv[argindex + 1], "-e", 2) == 0)
			{
				e = atof(argv[argindex + 2]);
				argindex = argindex + 2;
			}
			else
			{
				if (strncmp(argv[argindex + 1], "-p", 2) == 0)
				{
					project = true;
					argindex = argindex + 1;
				}
				else
				{
					if (strncmp(argv[argindex + 1], "-na", 2) == 0)
					{
					adaptive = false;
					argindex = argindex + 1;
					}
					else
					{
						out = argv[argindex + 1];
						argindex = argindex + 1;
					}
				}
			}
		}
	}
}
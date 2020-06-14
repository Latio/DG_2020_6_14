#include"NdgPhysMat.h"

MeshUnion mesh;
const MeshUnion *meshunion = &mesh;
int main()
{

	NdgPhysMat Solver;
	Solver.matSolver();

	system("pause");
	return 0;
}

#include "Solver.h"
#include "Material.h"


#include <iostream>

int main() {
	Material* material = new Material();

	Quadrature* quad = new Quadrature(4);
	//Geometry* geo = new Geometry(16,35);
	std::cout<<".o.o."<<endl;
//	quad->initializeArray();
	Solver* solver = new Solver(84,143,4);
	std::cout<<".o.o.111111"<<endl;
	solver->setQuadrature(quad);
	solver->setMaterial(material);
	std::cout<<".o.o.111111"<<endl;
	solver->startIteration();

	std::cout << "Residual1 = " << solver->getResidual1() << std::endl;
	std::cout << "Residual2 = " << solver->getResidual2() << std::endl;
	std::cout << "Residual3 = " << solver->getResidual3() << std::endl;
	std::cout << "Residual4 = " << solver->getResidual4() << std::endl;
	//std::cout << "R_left_45 " << geo->get_R_left(45)<< std::endl;
	std::cout << "R_0000000 "<< std::endl;


	return 0;

}

#ifndef SOLVER_H_
#define SOLVER_H_


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <memory.h>


#include "Quadrature.h"
#include "Material.h"
using namespace std;


class Solver {

private:
    /* Geometry */
	int _Nx;
	int _Ny;
	int _Lx;
	int _Ly;
	int _deltax;
	int _deltay;
	int _idx1;
	int _idx2;
	int _idy1;
	int _idy2;
	double _S0;
	double pi;

	double** a;
	double** b;

	/* Material */	
	double _sigma_tt[4];
	double _sigma_t;
	double _sigma_s0;
	double _sigma_s1;
	double _sigma_s2;
	double _sigma_s[4][4];
	double _externalsource;
	double _nu_sigma_f[4];
	double _distribution[4];
	double _A_left[143][84];
	double _A_right[143][84];
	double _B[84];
	double _V[143][84];

	
	/* Quadrature */
	int _n; /* Number of directions */
	int _M4;
	int _M1;
	int _M2;
	int _M3;
	Quadrature* _quadrature;
	Material* _material;

    int e;//Material ID!!!!!
    int o;//nengqunshu_intal
	/* Scalar flux */
	double*** _flux_x;
	double*** _flux_y;
	double*** _flux_c;//three dimension arry
	double*** _flux_m;
	double** _flux;
	double** _flux_pre;

	double** _flux_old_1;
	double** _flux_old_2;
	double** _flux_old_3;
	double** _flux_old_4;

	double** _flux_1;
	double** _flux_2;
	double** _flux_3;
	double** _flux_4;

	
	/* Angular flux */
	//double* _angular_flux;
	
	/* Boundary flux */
	//double* _boundary_flux;

	/* Source */
	double*** _scatter_source_1;
	double*** _scatter_source_2;
	double*** _scatter_source_3;
	double*** _scatter_source_4;
	double*** _fiss_source_1;
	double*** _fiss_source_2;
	double*** _fiss_source_3;
	double*** _fiss_source_4;
	double*** _tot_source;
	//double*** _ex_source;
	//double* _fiss_source;

	double _keff;
	
	/* Iteration */
	int _max_iter;
	double _converge_thresh;
	//     double _res_keff;
	double _res_flux_1;
	double _res_flux_2;
	double _res_flux_3;
	double _res_flux_4;
	double _res_keff;

public:

	Solver(int Nx=84,int Ny=143, int n=4);
	virtual ~Solver();

	void setQuadrature(Quadrature* quad);
	void setMaterial(Material* material);

	void startIteration();//key project

	void InitializeFlux();
	void InitializeSource();
	void InitialParameter();
	void InputAnVoidFlux();
    

    void storeOldflux();
	void Internal(double*** _scatter_source,double*** _fiss_source ,int Number);
	void scattersourceIteration1();
	void scattersourceIteration2();
	void scattersourceIteration3();
	void scattersourceIteration4();

	//void computeTotalSource();
	//void computeFissionSource();


	void Quadrat_1(int o);
	void Quadrat_2(int o);
	void Quadrat_3(int o);
	void Quadrat_4(int o);
	void initializeAllangleFlux();
	void computeFlux(int Number);
	void resetAngleFlux();

	void storeFlux();
	//void zeroFlux();

	void sourceIteration();
	void externalsource();//pull-in a fixed source

	void computeResidual();
	void computeTotalSource();
	void computeKeff();
	void computeFissionSource();

	double getResidual1() { return _res_flux_1; };
	double getResidual2() { return _res_flux_2; };
	double getResidual3() { return _res_flux_3; };
	double getResidual4() { return _res_flux_4; };

};

Solver::Solver(int Nx,int Ny, int n) {

	_Nx = Nx;
	_Ny = Ny;
	e=0;
	o=0;
	_n = n;
	_M4 = n*(n+2)/2;
	_M1 = n*(n+2)/8;
	_M2 = n*(n+2)/4;
	_M3 = 3*n*(n+2)/8;
	pi = 3.14159;

	//_Lx = 100;
	//_Ly = 100;
	//_deltax = _Lx/_Nx;
	//_deltay = _Ly/_Ny;
	/*_sigma_tt[0] = 0.12935;
	_sigma_tt[1] = 0.20274;
	_sigma_tt[2] = 0.20930;
	_sigma_tt[3] = 0.21275;
	_sigma_s[0][0]=0.12083;
	_sigma_s[0][1]=0.0;
	_sigma_s[0][2]=0.0;
	_sigma_s[0][3]=0.0;
	_sigma_s[1][0]=0.19763;
	_sigma_s[1][1]=0.84760e-2;
	_sigma_s[1][2]=0.0;
	_sigma_s[1][3]=0.0;
	_sigma_s[2][0]=0.19937;
	_sigma_s[2][1]=0.49168e-2;
	_sigma_s[2][2]=0.33889e-9;
	_sigma_s[2][3]=0.0;
	_sigma_s[3][0]=0.21045;
	_sigma_s[3][1]=0.75837e-2;
	_sigma_s[3][2]=0.91769e-9;
	_sigma_s[3][3]=0.1151e-12;
	_sigma_s1 = 0.01;
	_sigma_s2 = 0.0025;
	_nu_sigma_f[0]=0.50404e-4;
	_nu_sigma_f[1]=0.66532e-4;
	_nu_sigma_f[2]=0.42091e-3;
	_nu_sigma_f[3]=0.34020e-2;*/
	_distribution[0]=0.98439;
	_distribution[1]=0.156e-1;
	_distribution[2]=0.67690e-7;
	_distribution[3]=0.0;
	_keff =1.0;
	_max_iter = 1000;
	_S0=0.0;//no ex source

	_converge_thresh = 0.00001;
	_res_keff = 1.0;
	//_res_flux = 1.0;
	_quadrature = NULL;
	//_idx1=1;
	//_idx2=25/_deltax;
	//_idy1=25/_deltay+1;
	//_idy2 = _idy1 + 25 / _deltay - 1;
	
}

Solver::~Solver() {
	/*if (_quadrature != NULL)
		delete _quadrature;

	for (int i=0;i<_Ny+1;i++){
		for (int j=0;j<_Nx;j++){
			if (_flux_x[i][j] != NULL)
			{
				delete [] _flux_x[i][j];
			}
		}
	}

	for (int i=0;i<_Ny+1;i++){
		if (_flux_x[i] != NULL)
			{
				delete [] _flux_x[i];
			}
	}

	if (_flux_x != NULL)
			{
				delete [] _flux_x;
			}


	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx+1;j++){
			if (_flux_y[i][j] != NULL)
			{
				delete [] _flux_y[i][j];
			}
		}
	}

	for (int i=0;i<_Ny;i++){
		if (_flux_y[i] != NULL)
			{
				delete [] _flux_y[i];
			}
	}

	if (_flux_y != NULL)
			{
				delete [] _flux_y;
			}

    for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			if (_flux_c[i][j] != NULL)
			{
				delete [] _flux_c[i][j];
			}
		}
	}

	for (int i=0;i<_Ny;i++){
		if (_flux_c[i] != NULL)
			{
				delete [] _flux_c[i];
			}
	}

	if (_flux_c != NULL)
			{
				delete [] _flux_c;
			}

	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			if (_tot_source[i][j] != NULL)
			{
				delete [] _tot_source[i][j];
			}
		}
	}

	for (int i=0;i<_Ny;i++){
		if (_tot_source != NULL)
			{
				delete [] _tot_source;
			}
	}

	if (_tot_source != NULL)
			{
				delete [] _tot_source;
			}
    

	for (int i=0;i<_Ny;i++){
		if (_flux[i] != NULL)
			{
				delete [] _flux[i];
			}
	}

	if (_flux != NULL)
			{
				delete [] _flux;
			}

    for (int i=0;i<_Ny;i++){
		if (_flux_pre[i] != NULL)
			{
				delete [] _flux_pre[i];
			}
	}

	if (_flux_pre != NULL)
			{
				delete [] _flux_pre;
			}*/
}//should change delete in three dimension 

void Solver::setQuadrature(Quadrature* quad) {
	_quadrature = quad;
}

void Solver::setMaterial(Material* material){
	_material = material;
}

void Solver::startIteration() {
	int num = 0;
	int i = 0;
	_keff = 1.0;
	InitializeFlux();
    std::cout << "Flux array created!!!!!--!..." << std::endl;
    InputAnVoidFlux();//Need to give a suitable Flux    
	InitializeSource();
	InitialParameter();
    std::cout << "Flux array created!!!..." << std::endl;
	//externalsource();
    //std::cout << "Flux array created!!!!!!..." << std::endl;    
	//storeFlux();//if there should be a storeFlux?????
	//std::cout << "Flux array created!!!!!22222!..." << std::endl; 
	//0();
	//computeTotalSource();
	computeFissionSource();




	for (int i=0; i<500/*_max_iter*/; i++) {

	    //computeFissionSource();
        storeOldflux();

		scattersourceIteration1();
		//actually this function only add fission source

        Internal(_scatter_source_1, _fiss_source_1 ,1);
        scattersourceIteration2();
        Internal(_scatter_source_2, _fiss_source_2 ,2);
        scattersourceIteration3();
        Internal(_scatter_source_3, _fiss_source_3 ,3);
        scattersourceIteration4();
        Internal(_scatter_source_4, _fiss_source_4 ,4);

		//sourceIteration();
		
		//externalsource();
		//storeFlux();
		
		computeKeff();//There have a computeFissionSource!!!!
		computeResidual();
		/*for (int i=0; i<_Nx; i++) {
		   for (int j=0; j<_Ny; j++){
			   std::cout << "Mesh #" <<i<<j<< " = " << _flux[i][j] << '\t';
			}
		}*/
		




		//normalize();
		//computeTotalSource();
		//for (int j=0; j<_n/2; j++)
			//sweepMuMinus(j);
		//for (int j=_n/2-1; j>=0; j--)
			//sweepMuPlus(j);	
		//computeKeff();
		//computeResidual();
		//storeFlux();
		num++;
		if (i>1 && _res_flux_1 < _converge_thresh && _res_flux_2 < _converge_thresh && _res_flux_3 < _converge_thresh && _res_flux_4 < _converge_thresh)
			break;
		std::cout << "# " << num << '\t' << "resflux_1 = " << _res_flux_1 << std:: endl;
		std::cout << "# " << num << '\t' << "resflux_2 = " << _res_flux_2 << std:: endl;
		std::cout << "# " << num << '\t' << "resflux_3 = " << _res_flux_3 << std:: endl;
		std::cout << "# " << num << '\t' << "resflux_4 = " << _res_flux_4 << std:: endl;
		std::cout << "# " << num << '\t' << "Keff = " << _keff/2 << std:: endl;
		//std::cout << "0.0.0.0...45r zuo"<< '\t' <<_material->get_R_left(83)<<std::endl;
		//std::cout << "0.0.0.0...45r you"<< '\t' <<_material->get_R_right(83)<<std::endl;
		//std::cout << "0.0.0.0...45z shang"<< '\t' <<_material->get_Z_up(142)<<std::endl;
		//std::cout << "0.0.0.0...45z xia"<< '\t' <<_material->get_Z_down(142)<<std::endl;
		//std::cout << "0.0.0.0...45jiemian"<< '\t' <<_material->get_property(0,1,2)<<std::endl;
		//std::cout << "0.0.0.0...45xiwang"<< '\t' <<_material->get_xi(142,83)<<std::endl;
		/*for(int i=0;i<84;i++){
			std::cout << "0.0.0.0...45r zuo"<<i<< '\t' <<_material->get_R_left(i)<<std::endl;
			std::cout << "0.0.0.0...45r you"<<i<< '\t' <<_material->get_R_right(i)<<std::endl;
		}*/


	}
	/*for (int i=0; i<_Nx; i++) {
		   for (int j=0; j<_Ny; j++){
			   std::cout << "Mesh #" <<i<<j<< " = " << 0.0<< '\t';
			}
		}*/


	/*std::cout << "Scalar Flux: " << std::endl;
	for (int i=0; i<_Nx; i++) {
		for (int j=0; j<_Ny; j++){
			std::cout << "Mesh #" <<i<<j<< " = " << _flux[i][j] << '\t';
		if ((i+1)%4 ==0)
			std::cout << '\n';
		}
	}*/
}

void Solver::InitialParameter(){
	for(int i=0;i<_Ny;i++){
		for(int j=0;j<_Nx;j++){
           _A_left[i][j]=2*pi*_material->get_R_left(j)*((_material->get_Z_up(i))-(_material->get_Z_down(i)));
           _A_right[i][j]=2*pi*_material->get_R_right(j)*((_material->get_Z_up(i))-(_material->get_Z_down(i)));
           _V[i][j]=pi*(((_material->get_R_right(j))*(_material->get_R_right(j)))-(_material->get_R_left(j))*(_material->get_R_left(j)))*((_material->get_Z_up(i))-(_material->get_Z_down(i)));
		}
	}
	for(int j=0;j<_Nx;j++){
		_B[j]=pi*(((_material->get_R_right(j))*(_material->get_R_right(j)))-(_material->get_R_left(j))*(_material->get_R_left(j)));
	}
}

void Solver::storeOldflux(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			
		    _flux_old_1[i][j]=_flux_1[i][j];
			_flux_old_2[i][j]=_flux_2[i][j];
			_flux_old_3[i][j]=_flux_3[i][j];
			_flux_old_4[i][j]=_flux_4[i][j];			
		}
	}
}



void Solver::Internal(double*** _scatter_source,double*** _fiss_source,int Number){
	if (Number==1)
	{
		for (int i=0;i<_Ny;i++){
			for (int j=0;j<_Nx;j++){
				for (int m=0;m<_M4;m++){
					_tot_source[i][j][m]=_scatter_source_1[i][j][m]+_fiss_source_1[i][j][m]/_keff;
					//std::cout << "Mesh #" <<i<<j<< " = " << _tot_source[i][j][m] << '\t';
				}
			}
		}
		//_sigma_t=_sigma_tt[0];
		o=0;
	}
	else if (Number==2)
	{
		for (int i=0;i<_Ny;i++){
			for (int j=0;j<_Nx;j++){
				for (int m=0;m<_M4;m++){
					_tot_source[i][j][m]=_scatter_source_2[i][j][m]+_fiss_source_2[i][j][m]/_keff;
				}
			}
		}
		//_sigma_t=_sigma_tt[1];
		o=1;
	}
	else if (Number==3)
	{
		for (int i=0;i<_Ny;i++){
			for (int j=0;j<_Nx;j++){
				for (int m=0;m<_M4;m++){
					_tot_source[i][j][m]=_scatter_source_3[i][j][m]+_fiss_source_3[i][j][m]/_keff;
				}
			}
		}
		//_sigma_t=_sigma_tt[2];
		o=2;
	}
	else
	{
		for (int i=0;i<_Ny;i++){
			for (int j=0;j<_Nx;j++){
				for (int m=0;m<_M4;m++){
					_tot_source[i][j][m]=_scatter_source_4[i][j][m]+_fiss_source_4[i][j][m]/_keff;
				}
			}
		}
		//_sigma_t=_sigma_tt[3];
		o=3;
	}
	resetAngleFlux();
	Quadrat_3(o);

	Quadrat_4(o);

    Quadrat_2(o);
		std::cout << "Flux array created!!!!!!11111..." << std::endl;
	Quadrat_1(o);
		std::cout << "Flux array created!!!!!!11111..." << std::endl;
	computeFlux(Number);
	resetAngleFlux();
}



void Solver::resetAngleFlux(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_flux_x[i][j][m]=0.0;
				_flux_y[i][j][m]=0.0;
				_flux_m[i][j][m]=0.0;
				_flux_c[i][j][m]=0.0;
			}
		}
	}
}


void Solver::computeFlux(int Number){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			if (Number==1)
			{
				_flux_1[i][j]=0.0;
			}
			if (Number==2)
			{
				_flux_2[i][j]=0.0;
			}
			if (Number==3)
			{
				_flux_3[i][j]=0.0;
			}
			if (Number==4)
			{
				_flux_4[i][j]=0.0;
			}
			for (int m=0;m<_M4;m++){
				if (Number==1)
				{
					_flux_1[i][j]+=2*pi*_flux_c[i][j][m]/_M4;
				}
				else if (Number==2)
				{
					_flux_2[i][j]+=2*pi*_flux_c[i][j][m]/_M4;
				}
				else if (Number==3)
				{
					_flux_3[i][j]+=2*pi*_flux_c[i][j][m]/_M4;
				}
				else if (Number==4)
				{
					_flux_4[i][j]+=2*pi*_flux_c[i][j][m]/_M4;
				}
				else{
					std::cout <<"computeFlux is error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
				}
			}
		}
	}

}


void Solver::scattersourceIteration1(){
	double C00=0.0;
	//int e;
    for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			e=_material->get_xi(i,j);
			for (int m=0;m<_M4;m++){
				_scatter_source_1[i][j][m]=_material->get_property(0,e,5)*_flux_1[i][j]/(2*pi);
			}
		}
	}
	e=0;
}


void Solver::scattersourceIteration2(){
	double C00=0.0;
	double C01=0.0;
	//int e;
    for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			e=_material->get_xi(i,j);
			for (int m=0;m<_M4;m++){
				_scatter_source_2[i][j][m]=_material->get_property(1,e,6)*_flux_2[i][j]/(2*pi)+ _material->get_property(1,e,5)*_flux_1[i][j]/(2*pi);
			}
		}
	}
}

void Solver::scattersourceIteration3(){
	double C00=0.0;
	double C01=0.0;
	double C02=0.0;
	//int e;
    for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			e=_material->get_xi(i,j);
			for (int m=0;m<_M4;m++){
				_scatter_source_3[i][j][m]=_material->get_property(2,e,7)*_flux_3[i][j]/(2*pi)+ _material->get_property(2,e,6)*_flux_2[i][j]/(2*pi)+_material->get_property(2,e,5)*_flux_1[i][j]/(2*pi);
			}
		}
	}
}

void Solver::scattersourceIteration4(){
	double C00=0.0;
	double C01=0.0;
	double C02=0.0;
	double C03=0.0;
	//int e;
    for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			e=_material->get_xi(i,j);
			for (int m=0;m<_M4;m++){
				_scatter_source_4[i][j][m]=_material->get_property(3,e,8)*_flux_4[i][j]/(2*pi)+ _material->get_property(3,e,7)*_flux_3[i][j]/(2*pi)+_material->get_property(3,e,6)*_flux_2[i][j]/(2*pi)+_material->get_property(3,e,5)*_flux_1[i][j]/(2*pi);
			}
		}
	}
}


void Solver::InputAnVoidFlux(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){	
		  _flux_1[i][j]=100.0;
		  _flux_2[i][j]=100.0;
		  _flux_3[i][j]=100.0;
		  _flux_4[i][j]=100.0;	
		}
	}
}

void Solver::computeFissionSource(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			e=_material->get_xi(i,j);
			for (int m=0;m<_M4;m++){
				_fiss_source_1[i][j][m]=(_distribution[0]/*/_keff*/)*(_flux_1[i][j]*_material->get_property(0,e,2)+_flux_2[i][j]*_material->get_property(1,e,2)+_flux_3[i][j]*_material->get_property(2,e,2)+_flux_4[i][j]*_material->get_property(3,e,2));
				_fiss_source_2[i][j][m]=(_distribution[1]/*/_keff*/)*(_flux_1[i][j]*_material->get_property(0,e,2)+_flux_2[i][j]*_material->get_property(1,e,2)+_flux_3[i][j]*_material->get_property(2,e,2)+_flux_4[i][j]*_material->get_property(3,e,2));
				_fiss_source_3[i][j][m]=(_distribution[2]/*/_keff*/)*(_flux_1[i][j]*_material->get_property(0,e,2)+_flux_2[i][j]*_material->get_property(1,e,2)+_flux_3[i][j]*_material->get_property(2,e,2)+_flux_4[i][j]*_material->get_property(3,e,2));
				_fiss_source_4[i][j][m]=(_distribution[3]/*/_keff*/)*(_flux_1[i][j]*_material->get_property(0,e,2)+_flux_2[i][j]*_material->get_property(1,e,2)+_flux_3[i][j]*_material->get_property(2,e,2)+_flux_4[i][j]*_material->get_property(3,e,2));
			}
			//e=0;
		}
	}
}

void Solver::computeKeff(){
	double old_keff = _keff;
	double old_fiss_source = 0.0;
	double new_fiss_source = 0.0;
	for (int i=0; i<_Ny; i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				old_fiss_source =old_fiss_source + _V[i][j]*(_fiss_source_1[i][j][m]+_fiss_source_2[i][j][m]+_fiss_source_3[i][j][m]+_fiss_source_4[i][j][m]);
			}
		}
	}
	computeFissionSource();
	for (int i=0; i<_Ny; i++){
		for (int j=0;j<_Nx; j++){
			for (int m=0;m<_M4;m++){
				new_fiss_source =new_fiss_source + _V[i][j]*(_fiss_source_1[i][j][m]+_fiss_source_2[i][j][m]+_fiss_source_3[i][j][m]+_fiss_source_4[i][j][m]);
			}
		}
	}
	_keff = _keff*new_fiss_source/old_fiss_source;
	_res_keff = fabs((_keff-old_keff)/_keff);
	std::cout << "Show old fission source" <<'\t'<<old_fiss_source<< std::endl;
	std::cout << "Show new fission source" <<'\t'<<new_fiss_source<< std::endl;
}


/*void Solver::computeTotalSource(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_tot_source[i][j][m]=_tot_source[i][j][m]+_fiss_source[i][j][m];
			}
		}
	}
}*/

void Solver::InitializeFlux(){
	std::cout << "Flux array created..." << std::endl;
	_flux_x = (double***) new double**[_Ny + 1];
	for (int i=0;i<_Ny+1;i++){
		_flux_x[i] = (double**) new double* [_Nx];
	}
	for (int i=0;i<_Ny+1;i++){
		for (int j=0;j<_Nx;j++){
			_flux_x[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny+1;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_flux_x[i][j][m]=0.0;
			}
		}
	}

	_flux_y = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_y[i] = (double**) new double* [_Nx+1];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx+1;j++){
			_flux_y[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx+1;j++){
			for (int m=0;m<_M4;m++){
				_flux_y[i][j][m]=0.0;
			}
		}
	}


	_flux_c = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_c[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_flux_c[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_flux_c[i][j][m]=0.0;
			}
		}
	}

	_flux_m = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_m[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_flux_m[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_flux_m[i][j][m]=0.0;
			}
		}
	}


	_flux_1 = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_1[i] = new double [_Nx];
	};
	_flux_2 = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_2[i] = new double [_Nx];
	};
	_flux_3 = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_3[i] = new double [_Nx];
	};
	_flux_4 = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_4[i] = new double [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
		
			_flux_1[i][j]=0.0;
			_flux_2[i][j]=0.0;
			_flux_3[i][j]=0.0;
			_flux_4[i][j]=0.0;
		
		}
	}


	_flux_old_1 = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_old_1[i] = new double [_Nx];
	};
	_flux_old_2 = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_old_2[i] = new double [_Nx];
	};
	_flux_old_3 = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_old_3[i] = new double [_Nx];
	};
	_flux_old_4 = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_old_4[i] = new double [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
		
			_flux_old_1[i][j]=0.0;
			_flux_old_2[i][j]=0.0;
			_flux_old_3[i][j]=0.0;
			_flux_old_4[i][j]=0.0;
		}
	}

	a = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		a[i] = new double [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
		
			a[i][j]=0.0;
		
		}
	}

	b = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		b[i] = new double [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
		
			b[i][j]=0.0;
		
		}
	}
}

void Solver::InitializeSource(){
	_scatter_source_1 = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_scatter_source_1[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_scatter_source_1[i][j] = new double[_M4];
		}
	}
	_scatter_source_2 = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_scatter_source_2[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_scatter_source_2[i][j] = new double[_M4];
		}
	}
	_scatter_source_3 = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_scatter_source_3[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_scatter_source_3[i][j] = new double[_M4];
		}
	}
	_scatter_source_4 = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_scatter_source_4[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_scatter_source_4[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_scatter_source_1[i][j][m]=0.0;
				_scatter_source_2[i][j][m]=0.0;
				_scatter_source_3[i][j][m]=0.0;
				_scatter_source_4[i][j][m]=0.0;
			}
		}
	}
	
	_fiss_source_1 = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_fiss_source_1[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_fiss_source_1[i][j] = new double[_M4];
		}
	}
	_fiss_source_2 = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_fiss_source_2[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_fiss_source_2[i][j] = new double[_M4];
		}
	}
	_fiss_source_3 = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_fiss_source_3[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_fiss_source_3[i][j] = new double[_M4];
		}
	}
	_fiss_source_4 = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_fiss_source_4[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_fiss_source_4[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_fiss_source_1[i][j][m]=0.0;
				_fiss_source_2[i][j][m]=0.0;
				_fiss_source_3[i][j][m]=0.0;
				_fiss_source_4[i][j][m]=0.0;
			}
		}
	}
	_tot_source = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_tot_source[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_tot_source[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_tot_source[i][j][m]=0.0;
			}
		}
	}

}

/*void Solver::externalsource(){
	for (int i=_idy1-1;i<_idy2;i++){
		for (int j=_idx1-1;j<_idx2;j++){
			for (int m=0;m<_M4;m++){
				_tot_source[i][j][m]=_tot_source[i][j][m]+_S0/(2*pi);
			}

		}
	}
}*/

void Solver::Quadrat_1(int o){
	//int e;
	for (int j=0;j<_Nx;j++){
		for (int m=_M3;m<_M4;m++){
			_flux_x[0][j][m]=0.0;//there is 0
		}
	}
	for (int i=0;i<_Ny;i++){
		_flux_y[i][0][_M3]=_flux_y[i][0][6];
		_flux_y[i][0][_M3+1]=_flux_y[i][0][8];
		_flux_y[i][0][_M3+2]=_flux_y[i][0][7];
	}
	for (int m=_M3;m<_M4;m++){
		if (m-_M3==0)
		{
			//double G = (2*pi*_deltax*_deltay)*( _quadrature->getAlpha(m) + _quadrature->getAlpha(m+1) ) / _quadrature->getWeight();
			for (int i=0;i<_Ny;i++){
			   for (int j=0;j<_Nx;j++){
			   	   e=_material->get_xi(i,j);
			   	   double G = (_A_right[i][j]-_A_left[i][j])*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
				   double E=fabs(_quadrature->getMux(m))*(_A_left[i][j]+_A_right[i][j]);
		           double F=2*fabs(_quadrature->getMuy(m))*_B[j];
		           double C=E+F+G+_material->get_property(o,e,0)*_V[i][j];
		          _flux_c[i][j][m]=(E*_flux_y[i][j][m]+F*_flux_x[i][j][m]+G*a[i][j]+_tot_source[i][j][m]*_V[i][j])/C;
				  _flux_x[i+1][j][m]=2*_flux_c[i][j][m]-_flux_x[i][j][m];
				  _flux_y[i][j+1][m]=2*_flux_c[i][j][m]-_flux_y[i][j][m];
				  _flux_m[i][j][m]=2*_flux_c[i][j][m]-a[i][j];
				   if (_flux_x[i+1][j][m]<0.0)
				    {
					  _flux_x[i+1][j][m]=0.0;
				    }
				   if (_flux_y[i][j+1][m]<0.0)
				   {
					   _flux_y[i][j+1][m] = 0.0;
				   }
				   if (_flux_m[i][j][m]<0.0)
				   {
					  _flux_m[i][j][m] = 0.0;
				    }
				   if (_flux_c[i][j][m]<0.0)
				   {
					  _flux_c[i][j][m] = 0.0;
				    }
			    }
		    }
	    }
		else if (m-_M3==1)
		{
			//double G = (2*pi*_deltax*_deltay)*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
			for (int i=0;i<_Ny;i++){
			   for (int j=0;j<_Nx;j++){
			   	   e=_material->get_xi(i,j);
			   	   double G = (_A_right[i][j]-_A_left[i][j])*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
				   double E=fabs(_quadrature->getMux(m))*(_A_left[i][j]+_A_right[i][j]);
		           double F=2*fabs(_quadrature->getMuy(m))*_B[j];
		           double C=E+F+G+_material->get_property(o,e,0)*_V[i][j];
		          _flux_c[i][j][m]=(E*_flux_y[i][j][m]+F*_flux_x[i][j][m]+G*b[i][j]+_tot_source[i][j][m]*_V[i][j])/C;
				  _flux_x[i+1][j][m]=2*_flux_c[i][j][m]-_flux_x[i][j][m];
				  _flux_y[i][j+1][m]=2*_flux_c[i][j][m]-_flux_y[i][j][m];
				  _flux_m[i][j][m]=2*_flux_c[i][j][m]-b[i][j];
				  if (_flux_x[i+1][j][m]<0.0)
				   {
					  _flux_x[i+1][j][m]=0.0;
				   }
				  if (_flux_y[i][j+1][m]<0.0)
				   {
				      _flux_y[i][j+1][m] = 0.0;
				    }
				  if (_flux_m[i][j][m]<0.0)
				   {
					  _flux_m[i][j][m] = 0.0;
				    }
				   if (_flux_c[i][j][m]<0.0)
				   {
					  _flux_c[i][j][m] = 0.0;
				    }
			    }
		  
		     }
	    }
		else{
			//double G = (2*pi*_deltax*_deltay)*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
		   for (int i=0;i<_Ny;i++){
			   for (int j=0;j<_Nx;j++){
				   e=_material->get_xi(i,j);
			   	   double G = (_A_right[i][j]-_A_left[i][j])*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
				   double E=fabs(_quadrature->getMux(m))*(_A_left[i][j]+_A_right[i][j]);
		           double F=2*fabs(_quadrature->getMuy(m))*_B[j];
		           double C=E+F+G+_material->get_property(o,e,0)*_V[i][j];
				   _flux_c[i][j][m]=(E*_flux_y[i][j][m]+F*_flux_x[i][j][m]+G*_flux_m[i][j][m-1]+_tot_source[i][j][m]*_V[i][j])/C;
				   _flux_x[i+1][j][m]=2*_flux_c[i][j][m]-_flux_x[i][j][m];
				   _flux_y[i][j+1][m]=2*_flux_c[i][j][m]-_flux_y[i][j][m];
				   _flux_m[i][j][m]=2*_flux_c[i][j][m]-_flux_m[i][j][m-1];
				   if (_flux_x[i+1][j][m]<0.0)
				    {
					   _flux_x[i+1][j][m]=0.0;
				    }    
				    if (_flux_y[i][j+1][m]<0.0)
				    {
					   _flux_y[i][j+1][m] = 0.0;
				    }
				    if (_flux_m[i][j][m]<0.0)
				    {
					   _flux_m[i][j][m] = 0.0;
				    }
				    if (_flux_c[i][j][m]<0.0)
				   {
					  _flux_c[i][j][m] = 0.0;
				    } 
			    }
		    }
	    }
    }
   for (int i=0;i<_Ny;i++){
	    for (int j=0;j<_Nx;j++){
			a[i][j]=0.0;
			b[i][j]=0.0;
		}
	}
}

void Solver::Quadrat_2(int o){
	//int e;
	for (int j=0;j<_Nx;j++){
		for (int m=_M1;m<_M2;m++){
			_flux_x[0][j][m]=0.0;
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int m=_M1;m<_M2;m++){
			_flux_y[i][_Nx][m]=0.0;
		}
	}
	for (int m=_M2;m<_M3;m++){
		if ((m-_M2==0)||(m-_M2==1))
		{
			double G = 0.0;
			for (int i=0;i<_Ny;i++){
			   for (int j=_Nx-1;j>-1;j--){
			   	  e=_material->get_xi(i,j);
				   double E=fabs(_quadrature->getMux(m))*(_A_left[i][j]+_A_right[i][j]);
		           double F=2*fabs(_quadrature->getMuy(m))*_B[j];
		           double C=E+F+G+_material->get_property(o,e,0)*_V[i][j];
		          _flux_c[i][j][m]=(F*_flux_x[i][j][m]+E*_flux_y[i][j+1][m]+_tot_source[i][j][m]*_V[i][j])/C;
				  _flux_x[i+1][j][m]=2*_flux_c[i][j][m]-_flux_x[i][j][m];
				  _flux_y[i][j][m]=2*_flux_c[i][j][m]-_flux_y[i][j+1][m];
				  _flux_m[i][j][m]=_flux_c[i][j][m];
				   if (_flux_x[i+1][j][m]<0.0)
				   {
					   _flux_x[i+1][j][m]=0.0;
				    }
				   if (_flux_y[i][j][m]<0.0)
				   {
					   _flux_y[i][j][m] = 0.0;
				    }
				   if (_flux_m[i][j][m]<0.0)
				   {
					   _flux_m[i][j][m] = 0.0;
				    }
				    if (_flux_c[i][j][m]<0.0)
				   {
					  _flux_c[i][j][m] = 0.0;
				    }
			    }
		    }
	    }
		else{
			//double G = (2*pi*_deltax*_deltay)*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
		   for (int i=0;i<_Ny;i++){
			   for (int j=_Nx-1;j>-1;j--){
				   e=_material->get_xi(i,j);
			   	   double G = (_A_right[i][j]-_A_left[i][j])*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
				   double E=fabs(_quadrature->getMux(m))*(_A_left[i][j]+_A_right[i][j]);
		           double F=2*fabs(_quadrature->getMuy(m))*_B[j];
		           double C=E+F+G+_material->get_property(o,e,0)*_V[i][j];
				   _flux_c[i][j][m]=(F*_flux_x[i][j][m]+E*_flux_y[i][j+1][m]+G*_flux_m[i][j][m]+_tot_source[i][j][m]*_V[i][j])/C;
				   _flux_x[i+1][j][m]=2*_flux_c[i][j][m]-_flux_x[i][j][m];
				   _flux_y[i][j][m]=2*_flux_c[i][j][m]-_flux_y[i][j+1][m];
				   _flux_m[i][j][m]=2*_flux_c[i][j][m]-_flux_m[i][j][m-1];
				    if (_flux_x[i+1][j][m]<0.0)
				   {
					   _flux_x[i+1][j][m]=0.0;
				    }
				    if (_flux_y[i][j][m]<0.0)
				    {
					   _flux_y[i][j][m] = 0.0;
				    }
				    if (_flux_m[i][j][m]<0.0)
				    {
					   _flux_m[i][j][m] = 0.0;
				    }
				    if (_flux_c[i][j][m]<0.0)
				   {
					  _flux_c[i][j][m] = 0.0;
				    }
			    }
		    }
	    }
    }
    for (int i=0;i<_Ny;i++){
		    for (int j=0;j<_Nx;j++){

			    a[i][j]=_flux_m[i][j][_M2];
			    b[i][j]=_flux_m[i][j][_M2+2];
		    }
	}
}

void Solver::Quadrat_3(int o){
	//int e;
	for (int j=0;j<_Nx;j++){
		for (int m=0;m<_M1;m++){
			_flux_x[_Ny][j][m]=0.0;
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int m=0;m<_M1;m++){
			_flux_y[i][_Nx][m]=0.0;
		}
	}

	for (int m=0;m<_M1;m++){
		if ((m==0)||(m==1))//there need to think!!!
		{
			double G = 0.0;
			for (int i=_Ny-1;i>-1;i--){
			   for (int j=_Nx-1;j>-1;j--){
			   	  e=_material->get_xi(i,j);
			   	  double E=fabs(_quadrature->getMux(m))*(_A_left[i][j]+_A_right[i][j]);
		          double F=2*fabs(_quadrature->getMuy(m))*_B[j];
		          double C=E+F+G+_material->get_property(o,e,0)*_V[i][j];
		          _flux_c[i][j][m]=(F*_flux_x[i+1][j][m]+E*_flux_y[i][j+1][m]+_tot_source[i][j][m]*_V[i][j])/C;
				  _flux_x[i][j][m]=2*_flux_c[i][j][m]-_flux_x[i+1][j][m];
				  _flux_y[i][j][m]=2*_flux_c[i][j][m]-_flux_y[i][j+1][m];
				  _flux_m[i][j][m]=_flux_c[i][j][m];
				  //std::cout << "Mesh #" <<i<<j<< " = " << E<< '\t'<<F<< '\t'<<C << '\t';
				if (_flux_x[i][j][m]<0.0)
				{
					_flux_x[i][j][m]=0.0;
				}
				if (_flux_y[i][j][m]<0.0)
				{
					_flux_y[i][j][m] = 0.0;
				}
				if (_flux_m[i][j][m]<0.0)
				{
					_flux_m[i][j][m] = 0.0;
				}
				if (_flux_c[i][j][m]<0.0)
				   {
					  _flux_c[i][j][m] = 0.0;
				    }
			    }
		    }
	    }
    
      else{
			//double G = (2*pi*_deltax*_deltay)*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
		   for (int i=_Ny-1;i>-1;i--){
			   for (int j=_Nx-1;j>-1;j--){
			   	   e=_material->get_xi(i,j);
			   	   double G = (_A_right[i][j]-_A_left[i][j])*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
				   double E=fabs(_quadrature->getMux(m))*(_A_left[i][j]+_A_right[i][j]);
		           double F=2*fabs(_quadrature->getMuy(m))*_B[j];
		           double C=E+F+G+_material->get_property(o,e,0)*_V[i][j];
				   _flux_c[i][j][m]=(F*_flux_x[i+1][j][m]+E*_flux_y[i][j+1][m]+G*_flux_m[i][j][m-1]+_tot_source[i][j][m]*_V[i][j])/C;
				   _flux_x[i][j][m]=2*_flux_c[i][j][m]-_flux_x[i+1][j][m];
				   _flux_y[i][j][m]=2*_flux_c[i][j][m]-_flux_y[i][j+1][m];
				   _flux_m[i][j][m]=2*_flux_c[i][j][m]-_flux_m[i][j][m-1];
				    if (_flux_x[i][j][m]<0.0)
				    {  
					   _flux_x[i][j][m]=0.0;
				    }
				    if (_flux_y[i][j][m]<0.0)
				    {
					   _flux_y[i][j][m] = 0.0;
				    }
				    if (_flux_m[i][j][m]<0.0)
				    {
					  _flux_m[i][j][m] = 0.0;
				    }
				    if (_flux_c[i][j][m]<0.0)
				   {
					  _flux_c[i][j][m] = 0.0;
				    }
			    }
			}
	    }
    }
    for (int i=0;i<_Ny;i++){
		   for (int j=0;j<_Nx;j++){
			  a[i][j]=_flux_m[i][j][0];
			  b[i][j]=_flux_m[i][j][2];
		    }
	}
}

void Solver::Quadrat_4(int o){
	//int e;
	for (int j=0;j<_Nx;j++){
		for (int m=_M1;m<_M2;m++){
			_flux_x[_Ny][j][m]=0.0;
		}
	}
	for (int i=0;i<_Ny;i++){
	    _flux_y[i][0][_M1]=_flux_y[i][0][0];
	    _flux_y[i][0][_M1+1]=_flux_y[i][0][2];
	    _flux_y[i][0][_M1+2]=_flux_y[i][0][1];
	}
	for (int m=_M1;m<_M2;m++){
		if (m-_M1==0)
		{
			//double G = (2*pi*_deltax*_deltay)*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
			for (int i=_Ny-1;i>-1;i--){
			   for (int j=0;j<_Nx;j++){
			   	   e=_material->get_xi(i,j);
			   	   double G = (_A_right[i][j]-_A_left[i][j])*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
				   double E=fabs(_quadrature->getMux(m))*(_A_left[i][j]+_A_right[i][j]);
		           double F=2*fabs(_quadrature->getMuy(m))*_B[j];
		           double C=E+F+G+_material->get_property(o,e,0)*_V[i][j];
		          _flux_c[i][j][m]=(F*_flux_x[i+1][j][m]+E*_flux_y[i][j][m]+G*a[i][j]+_tot_source[i][j][m]*_V[i][j])/C;
				  _flux_x[i][j][m]=2*_flux_c[i][j][m]-_flux_x[i+1][j][m];
				  _flux_y[i][j+1][m]=2*_flux_c[i][j][m]-_flux_y[i][j][m];
				  _flux_m[i][j][m]=2*_flux_c[i][j][m]-a[i][j];
				if (_flux_x[i][j][m]<0.0)
				{
					_flux_x[i][j][m]=0.0;
				}
				if (_flux_y[i][j+1][m]<0.0)
				{
					_flux_y[i][j+1][m] = 0.0;
				}
				if (_flux_m[i][j][m]<0.0)
				{
					_flux_m[i][j][m] = 0.0;
				}
				if (_flux_c[i][j][m]<0.0)
				   {
					  _flux_c[i][j][m] = 0.0;
				    }
			    }
		    }
	    }
		else if (m-_M1==1)
		{
			//double G = (2*pi*_deltax*_deltay)*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
			for (int i=_Ny-1;i>-1;i--){
			   for (int j=0;j<_Nx;j++){
			   	   e=_material->get_xi(i,j);
			   	   double G = (_A_right[i][j]-_A_left[i][j])*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
				   double E=fabs(_quadrature->getMux(m))*(_A_left[i][j]+_A_right[i][j]);
		           double F=2*fabs(_quadrature->getMuy(m))*_B[j];
		           double C=E+F+G+_material->get_property(o,e,0)*_V[i][j];
		          _flux_c[i][j][m]=(F*_flux_x[i+1][j][m]+E*_flux_y[i][j][m]+G*b[i][j]+_tot_source[i][j][m]*_V[i][j])/C;
				  _flux_x[i][j][m]=2*_flux_c[i][j][m]-_flux_x[i+1][j][m];
				  _flux_y[i][j+1][m]=2*_flux_c[i][j][m]-_flux_y[i][j][m];
				  _flux_m[i][j][m]=2*_flux_c[i][j][m]-b[i][j];
				if (_flux_x[i][j][m]<0.0)
				{
					_flux_x[i][j][m]=0.0;
				}
				if (_flux_y[i][j+1][m]<0.0)
				{
					_flux_y[i][j+1][m] = 0.0;
				}
				if (_flux_m[i][j][m]<0.0)
				{
					_flux_m[i][j][m] = 0.0;
				}
				if (_flux_c[i][j][m]<0.0)
				   {
					  _flux_c[i][j][m] = 0.0;
				    }
			    }
		  
		    }
	    }
		else{
			//double G = (2*pi*_deltax*_deltay)*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
		   for (int i=_Ny-1;i>-1;i--){
			   for (int j=0;j<_Nx;j++){
			   	   e=_material->get_xi(i,j);
			   	   double G = (_A_right[i][j]-_A_left[i][j])*(_quadrature->getAlpha(m)+_quadrature->getAlpha(m+1))/_quadrature->getWeight();
				   double E=fabs(_quadrature->getMux(m))*(_A_left[i][j]+_A_right[i][j]);
		           double F=2*fabs(_quadrature->getMuy(m))*_B[j];
		           double C=E+F+G+_material->get_property(o,e,0)*_V[i][j];
				   _flux_c[i][j][m]=(F*_flux_x[i+1][j][m]+E*_flux_y[i][j][m]+G*_flux_m[i][j][m-1]+_tot_source[i][j][m]*_V[i][j])/C;
				   _flux_x[i][j][m]=2*_flux_c[i][j][m]-_flux_x[i+1][j][m];
				   _flux_y[i][j+1][m]=2*_flux_c[i][j][m]-_flux_y[i][j][m];
				   _flux_m[i][j][m]=2*_flux_c[i][j][m]-_flux_m[i][j][m-1];
				if (_flux_x[i][j][m]<0.0)
				{
					_flux_x[i][j][m]=0.0;
				}
				if (_flux_y[i][j+1][m]<0.0)
				{
					_flux_y[i][j+1][m] = 0.0;
				}
				if (_flux_m[i][j][m]<0.0)
				{
					_flux_m[i][j][m] = 0.0;
				}
				if (_flux_c[i][j][m]<0.0)
				   {
					  _flux_c[i][j][m] = 0.0;
				    }
			    }
		    }
	    }

    }
    for (int i=0;i<_Ny;i++){
	    for (int j=0;j<_Nx;j++){
			  a[i][j]=0.0;
			  b[i][j]=0.0;
		}
	}
}

/*void Solver::sourceIteration(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_tot_source[i][j][m] = 0.0;
			}
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			double C00 = 0.0;
			double C1_1 = 0.0;
			double C11 = 0.0;
			double C2_2 = 0.0;
			double C20 = 0.0;
			double C22 = 0.0;
			for (int m=0;m<_M4;m++){
				C00=C00+_quadrature->getY00(m)*_flux_c[i][j][m]/_M4;
				C1_1=C1_1+_quadrature->getY1_1(m)*_flux_c[i][j][m]/_M4;
				C11=C11+_quadrature->getY11(m)*_flux_c[i][j][m]/_M4;
				C2_2=C2_2+_quadrature->getY2_2(m)*_flux_c[i][j][m]/_M4;
				C20=C20+_quadrature->getY20(m)*_flux_c[i][j][m]/_M4;
				C22=C22+_quadrature->getY22(m)*_flux_c[i][j][m]/_M4;
			}
			for (int m=0;m<_M4;m++){
				_tot_source[i][j][m]=_sigma_s0*_quadrature->getY00(m)*C00+3*_sigma_s1*(_quadrature->getY11(m)*C11+_quadrature->getY1_1(m)*C1_1)+5*_sigma_s2*(_quadrature->getY2_2(m)*C2_2/12+_quadrature->getY20(m)*C20/12+_quadrature->getY22(m)*C22/12);
				_tot_source[i][j][m]=_tot_source[i][j][m]/(2*pi);
			}
		}
	}
}*/

/*void Solver::storeFlux(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_flux_pre[i][j]=_flux[i][j];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_flux[i][j]=0.0;
		}
	}
	double XX=0.0;
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				XX=XX+4*pi*_flux_c[i][j][m]/_M4;
			}
			_flux[i][j]=XX;
			XX=0.0;
		}
	}
}*/

void Solver::computeResidual(){
	double max_scalarflux_1 = 0.0;
	double max_scalarflux_2 = 0.0;
	double max_scalarflux_3 = 0.0;
	double max_scalarflux_4 = 0.0;
	double mm1;
	double mm2;
	double mm3;
	double mm4;
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			mm1 = fabs(_flux_1[i][j]-_flux_old_1[i][j]);
			if (mm1 > max_scalarflux_1)
			{
				max_scalarflux_1 = mm1;
			}
			mm2 = fabs(_flux_2[i][j]-_flux_old_2[i][j]);
			if (mm2 > max_scalarflux_2)
			{
				max_scalarflux_2 = mm2;
			}
			mm3 = fabs(_flux_3[i][j]-_flux_old_3[i][j]);
			if (mm3 > max_scalarflux_3)
			{
				max_scalarflux_3 = mm3;
			}
			mm4 = fabs(_flux_4[i][j]-_flux_old_4[i][j]);
			if (mm4 > max_scalarflux_4)
			{
				max_scalarflux_4 = mm4;
			}
		 }
	}
	_res_flux_1 = max_scalarflux_1;
	_res_flux_2 = max_scalarflux_2;
	_res_flux_3 = max_scalarflux_3;
	_res_flux_4 = max_scalarflux_4;
}


#endif

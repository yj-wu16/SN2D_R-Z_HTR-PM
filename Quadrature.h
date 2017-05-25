#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

class Quadrature{

private:
 int _n;
 int M1;
 int M4;
 int M5;
 double* _mux;
 double* _muy;
 double* _muz;
 double _w;
 double* _Y00;
 double* _Y1_1;
 double* _Y11;
 double* _Y2_2;
 double* _Y20;
 double* _Y22;
 double* _alpha;

 

public:
 Quadrature(int n);
 virtual ~Quadrature();

 void initializeArray();//need to do more work
 void initializeSphericalHarmonics();//new pull-in
 void initialize_Mu(int M);
 void initialize_w(int M);
 void initialize_alpha();


 void setNumDirection(int num);
 int getNumDirection() { return _n; }
 double getMux(int d);
 double getMuy(int d);
 double getMuz(int d);
 double getY00(int d);
 double getY1_1(int d);
 double getY11(int d);
 double getY2_2(int d);
 double getY20(int d);
 double getY22(int d);
 double getWeight();
 double getAlpha(int d);
};

Quadrature::Quadrature(int num) {

	setNumDirection(num);
	_mux = NULL;
	_muy = NULL;
	_muz = NULL;
	_w = 1.0/3.0;
	initializeArray();
	initializeSphericalHarmonics();
    initialize_alpha();
}

Quadrature::~Quadrature() {
	if (_mux != NULL)
		delete [] _mux;
	if (_muy != NULL)
		delete [] _muy;
	if (_muz != NULL)
		delete [] _muz;
}//need to do

void Quadrature::setNumDirection(int num) {
	if (num % 2 != 0)
		std::cout << "Can not set num of Direction to " << num << " , cause it is not even." << std::endl;
	_n = num;//here i choose to substitude _n with N
}

double Quadrature::getMux(int d) {
	return _mux[d];
}

double Quadrature::getMuy(int d) {
	return _muy[d];
}

double Quadrature::getMuz(int d) {
	return _muz[d];
}

double Quadrature::getWeight() {
	return _w;
}

void Quadrature::initializeArray(){
	_mux = new double [_n*(_n+2)/2];
	_muy = new double [_n*(_n+2)/2];
	_muz = new double [_n*(_n+2)/2];
	//_w = new double [_n];
	if (_n == 0){
		std::cout << "ERROR, Quadrature can NOT creat array..." << std::endl;
	}
	if (_n == 4){
		initialize_Mu(4);
		initialize_w(4);
	}
	else if (_n == 6){
		initialize_Mu(6);
		initialize_w(6);
	}

}
void Quadrature::initialize_Mu(int M){
	switch (M){
	case 4:
		M1=M*(M+2)/8;
		double cos[2]={0.3500212,0.8688903};
		double ux1[12]={-cos[1],-cos[2],-cos[1], cos[1], cos[1], cos[2],-cos[1],-cos[2],-cos[1], cos[1], cos[1], cos[2]};
        double uy1[12]={-cos[2],-cos[1],-cos[1],-cos[2],-cos[1],-cos[1], cos[2], cos[1], cos[1], cos[2], cos[1], cos[1]};
        double uz1[3]={cos[2],cos[1],cos[1]};
        for (int i=0; i<(4*M1); i++){
        	_mux[i]=ux1[i];
        } 	
        /*for (int i=M1; i<(2*M1); i++){
        	_mux[i]=-ux1[i-M1];
        } 
        for (int i=2*M1; i<(3*M1); i++){
        	_mux[i]=-ux1[i-2*M1];
        } 
        for (int i=3*M1; i<(4*M1); i++){
        	_mux[i]=ux1[i-3*M1];
        } //x direction quadrature*/
        for (int i=0; i<(4*M1); i++){
        	_muy[i]=uy1[i];
        } 	
        /*for (int i=M1; i<(2*M1); i++){
        	_muy[i]=uy1[i-M1];
        } 
        for (int i=2*M1; i<(3*M1); i++){
        	_muy[i]=-uy1[i-2*M1];
        } 
        for (int i=3*M1; i<(4*M1); i++){
        	_muy[i]=-uy1[i-3*M1];
        } //y direction quadrature*/
        for (int i=0; i<(M1); i++){
        	_muz[i]=uz1[i];
        } 	
        for (int i=M1; i<(2*M1); i++){
        	_muz[i]=uz1[i-M1];
        } 
        for (int i=2*M1; i<(3*M1); i++){
        	_muz[i]=uz1[i-2*M1];
        } 
        for (int i=3*M1; i<(4*M1); i++){
        	_muz[i]=uz1[i-3*M1];
        } //z direction quadrature
		break;
	/*case 6:
	
		M1=M*(M+2)/8;
		double cos[4]={0.2561429,0.2663443,0.6815646,0.9320846};
		double ux1[6]={cos[1],cos[2],cos[1],cos[3],cos[3],cos[4]};
        double uy1[6]={cos[1],cos[3],cos[4],cos[2],cos[3],cos[1]};
        double uz1[6]={cos[4],cos[3],cos[1],cos[3],cos[2],cos[1]};
        for (int i=0; i<(M1); i++){
        	_mux[i]=ux1[i];
        } 	
        for (int i=M1; i<(2*M1); i++){
        	_mux[i]=-ux1[i];
        } 
        for (int i=2*M1; i<3*M1; i++){
        	_mux[i]=-ux1[i];
        } 
        for (int i=3*M1; i<(4*M1); i++){
        	_mux[i]=ux1[i];
        } //x direction quadrature
        for (int i=1; i<(M1); i++){
        	_muy[i]=uy1[i];
        } 	
        for (int i=M1; i<(2*M1); i++){
        	_muy[i]=uy1[i];
        } 
        for (int i=2*M1; i<(3*M1); i++){
        	_muy[i]=-uy1[i];
        } 
        for (int i=3*M1; i<(4*M1); i++){
        	_muy[i]=-uy1[i];
        } //y direction quadrature
        for (int i=1; i<(M1); i++){
        	_muz[i]=uz1[i];
        } 	
        for (int i=M1; i<(2*M1); i++){
        	_muz[i]=uz1[i];
        } 
        for (int i=2*M1; i<(3*M1); i++){
        	_muz[i]=uz1[i];
        } 
        for (int i=3*M1; i<(4*M1); i++){
        	_muz[i]=uz1[i];
        } //z direction quadrature
       break;*/
}
}

void Quadrature::initializeSphericalHarmonics(){
	M4 = _n*(_n+2)/2;
	_Y00 = new double [M4];
	_Y11 = new double [M4];
	_Y1_1 = new double [M4];
	_Y20 = new double [M4];
	_Y22 = new double [M4];
	_Y2_2 = new double [M4];
	for (int i=0;i<M4;i++){
		_Y00[i] = 0;
		_Y11[i] = -_mux[i];
		_Y1_1[i] = -_muy[i];
		_Y2_2[i] = 6*(_mux[i]*_muy[i]);
		_Y20[i] = (3*_muz[i]*_muz[i]-1)/2;
		_Y22[i] = 3*(_mux[i]*_mux[i]-_muy[i]*_muy[i]);
	}
}

double Quadrature::getY00(int d){
        return _Y00[d];
}

double Quadrature::getY11(int d){
        return _Y11[d];
}

double Quadrature::getY1_1(int d){
        return _Y1_1[d];
}

double Quadrature::getY22(int d){
        return _Y22[d];
}

double Quadrature::getY20(int d){
        return _Y20[d];
}

double Quadrature::getY2_2(int d){
        return _Y2_2[d];
}

void Quadrature::initialize_w(int M){
	_w=3.14159/24.0;
    std::cout << "whats wrong ????" << std::endl;
    std::cout << "w="<<_w << std::endl;
	/*case 6 : { _w = 1/6;}*/
}

void Quadrature::initialize_alpha(){
    M5 = _n*(_n+2)/2+1;
    _alpha = new double [M5];
    _alpha[0] = 0.0;
    for (int n=1;n<M5;n++){
        _alpha[n]=-_mux[n-1]*3.14159/24+_alpha[n-1];

    }
}

double Quadrature::getAlpha(int d){
    return _alpha[d];
}


#endif

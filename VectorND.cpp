#include <iostream>
#include <cstdlib>
#include "VectorND.hpp"
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>


/*Este metodo pide al sistema el tamaño y tipo de memoria*/
/*que necesitamos para un vector*/
/*Se usa como   	VectorND vecGo;      */
/*					vecGo.initVectorND(5, 1, 1, 1, 1, 1)*/
void VectorND::initVectorND(int n, ...){

	va_list list;
	this->n = n;
	this->A = (double *) malloc (n * sizeof(double));
	int i;

	va_start(list, n);
	/**/
	/**/	for (int i = 0; i < n; i++)
	/**/		this->A[i] = va_arg(list, double);
	/**/
	va_end(list);
	
	double norma = 0;
	for (i = 0; i < this->n; i++)
		norma += this->A[i] * this->A[i];

	this->norma = sqrt(norma);
}

/*O simplemente crea un espacio de memoria vacío pero ya listo para almacenar objetos*/
void VectorND::initVectorNDType2(int n){

	this->n = n;
	this->A = (double *) malloc (n * sizeof(double));
}

/*Con este metodo puedes acceder al número de filas.*/
/*int a = vecGo.dim();*/
int VectorND::dim() {
	return this->n;
}

double VectorND::norm() {
	return this->norma;
}

void VectorND::zeroVectorND(int n) {

	this->initVectorNDType2(n);
	for (int i = 0; i < n; i++)
			this->A[i] = 0;

	this->norma = 0.0;
}

/*Cuenta cuantos digitos tiene un double*/
/*int a = vecGo.digitCounter(0.3434)*/
int VectorND::digitCounter(double a) {
	
	string hello = to_string(a);
	int ret = hello.length();
	return ret;
}

/*Este metodo escribe en pantalla de manera ordenada*/
/*matGo.escVectorND();*/
void VectorND::escVectorND() {

	int k = 20;
	printf("(");
	for (int i = 0; i < this->n; i++) {
			cout << to_string(this->A[i]);
			int l = (this)->digitCounter(this->A[i]);
			for (int r = 0; r < k - l; r++)
				cout << " ";
	}
	printf(")");

}

/*La mayoría de las veces trabajaremos con Vectores en 3D entonces haremos un metodo*/
/*Que trabaja con vectores de dimensión 3*/
void VectorND::updateVector3D(double x, double y, double z) {

	if (this->n != 3) printf("\n\nVECTORND CANDIDATE NOT AVAILABLE FOR UPDATE\n\nCHECK DIMENSION\n\n");
	else {

		/*Actualizamos cada una de las nuevas entradas*/
		this->A[0] = x;
		this->A[1] = y;
		this->A[2] = z;

		/*La norma no se actualiza sola también es necesario recalcularla*/
		double norma = 0;
		for (int i = 0; i < this->n; i++)
			norma += this->A[i] * this->A[i];

		this->norma = sqrt(norma);		
	}
}

void VectorND::updateVector3DP(VectorND p) {

	this->updateVector3D(p.A[0], p.A[1], p.A[2]);
}

double VectorND::access(int i)
{	
	if(i >= this->n) 
	{	
		return printf("The dimension is wrong");
	}
	else
	{
		return this->A[i];
	}
}

void VectorND::AddVector3D(VectorND B)
{
	if(B.dim() != this->n){
		printf("\n\n THIS OPERATION CANNOT BE EXECUTED \n\n CHECK DIMENSION \n\n");
	}

	else{
	for (int i = 0; i < this->n; i++)
		this->A[i] += B.access(i);

	double norma = 0;
	for (int i = 0; i < this->n; i++)
		norma += this->A[i] * this->A[i];

	this->norma = sqrt(norma);	
	}
}

VectorND VectorND::SumVector3D(VectorND B, VectorND C)
{
	VectorND ret;
	ret.initVectorND(3, 0.0, 0.0, 0.0);

	this-> n = 3; 

	if(C.dim() != B.dim()) {
		printf("\n\n THIS OPERATION CANNOT BE EXECUTED \n\n CHECK DIMENSION \n\n");
		return ret;
	}
	
	else{
		
		for(int i = 0; i < this->n; i++)
			ret.A[i] = B.access(i) + C.access(i);

		double norma = 0;
		for (int i = 0; i < this->n; i++)
			norma += ret.A[i] * ret.A[i];

		this->norma = sqrt(norma);	

		return ret;
	}
}

void VectorND::SumVectorND(VectorND B, VectorND C)
{
	if(C.dim() != B.dim())
		printf("\n\n SumVectorND :: THIS OPERATION CANNOT BE EXECUTED \n\n CHECK DIMENSION \n\n");
	
	else{
		
		for(int i = 0; i < this->n; i++)
			this->A[i] = B.access(i) + C.access(i);

		double norma = 0;
		for (int i = 0; i < this->n; i++)
			norma += this->A[i] * this->A[i];

		this->norma = sqrt(norma);	

	}
}

void VectorND::AddVectorND(VectorND B)
{
	if(B.dim() != this->n){
		printf("\n\n AddVectorND :: THIS OPERATION CANNOT BE EXECUTED \n\n CHECK DIMENSION \n\n");
	}

	else{
	for (int i = 0; i < this->n; i++)
		this->A[i] += B.access(i);

	double norma = 0;
	for (int i = 0; i < this->n; i++)
		norma += this->A[i] * this->A[i];

	this->norma = sqrt(norma);	
	}
}

void VectorND::MinusVector3D(VectorND B){

if(B.dim() != this->n){
		printf("\n\n THIS OPERATION CANNOT BE EXECUTED \n\n CHECK DIMENSION \n\n");
	}

	else{
	for (int i = 0; i < this->n; i++)
		this->A[i] -= B.access(i);

	double norma = 0;
	for (int i = 0; i < this->n; i++)
		norma += this->A[i] * this->A[i];

	this->norma = sqrt(norma);	
	}
}	

void VectorND::SubVector3D(VectorND B, VectorND C)
{
	this-> n = 3; 

	if(C.dim() != B.dim() != this->n)
		printf("\n\n THIS OPERATION CANNOT BE EXECUTED \n\n CHECK DIMENSION \n\n");
	
	else{
		
		for(int i = 0; i < this->n; i++)
			this->A[i] = B.access(i) - C.access(i);

		double norma = 0;
		for (int i = 0; i < this->n; i++)
			norma += this->A[i] * this->A[i];

		this->norma = sqrt(norma);	

	}
}

void VectorND::subVectorND(VectorND a, VectorND b) {

	if (a.dim() != b.dim())
		printf("\n\n THIS OPERATION CANNOT BE EXECUTED \n\n CHECK DIMENSION \n\n");
	else {

		for (int i = 0; i < this->n; i++)
			this->A[i] = a.A[i] - b.A[i];

		double norma = 0;
		for (int i = 0; i < this->n; i++)
			norma += this->A[i] * this->A[i];

		this->norma = sqrt(norma);
	}
}

void VectorND::scaleVectorND(double s, VectorND a) {

	if (a.dim() != this->dim())
		printf("\n\n scaleVectorND :: THIS OPERATION CANNOT BE EXECUTED \n\n CHECK DIMENSION \n\n");
	else {

		for (int i = 0; i < this->n; i++)
			this->A[i] = s * a.A[i];

		double norma = 0;
		for (int i = 0; i < this->n; i++)
			norma += this->A[i] * this->A[i];

		this->norma = sqrt(norma);
	}
}

void VectorND::scaleND(double s) {
	
		for (int i = 0; i < this->n; i++)
			this->A[i] *= s;

		double norma = 0;
		for (int i = 0; i < this->n; i++)
			norma += this->A[i] * this->A[i];

		this->norma = sqrt(norma);
}

double VectorND::dotProd(VectorND b) {

	double norma = 0;
	for (int i = 0; i < this->n; i++)
		norma += this->A[i] * b.A[i];
	return norma;
}

void VectorND::Cross3D(VectorND B, VectorND C){
	
	this->n = 3;

	if(C.dim() != B.dim())
		printf("\n\n THIS OPERATION CANNOT BE EXECUTED \n\n CHECK DIMENSION \n\n");	
	else{
	this->A[0] = ( B.access(1) * C.access(2)) - (B.access(2) * C.access(1));
	this->A[1] = (-B.access(0) * C.access(2)) + (B.access(2) * C.access(0));
	this->A[2] = ( B.access(0) * C.access(1)) - (B.access(1) * C.access(0));	
	
	double norma = 0;	
	for (int i = 0; i < this->n; i++)
		norma += this->A[i] * this->A[i];

	this->norma = sqrt(norma);
	}
				
}

void VectorND::unit() {
	
	for (int i = 0; i < this->n; i++) {
		this->A[i] *= 1/(this->norm()); 
	}

	this->norma = 1.0;
}

double VectorND::distance(VectorND b) {

	if (this->n != b.n) {
		printf("\n\ndistance :: ERROR DIMENSION");
		return -11;
	}
	else {
		double ret = 0;
		for (int i = 0; i < this->n; i++)
			ret += (this->access(i) - b.access(i)) * (this->access(i) - b.access(i));
		ret = sqrt(ret);
		return ret;
	}
}

int VectorND::equalVector(VectorND b) {

	double epsilon = 1e-20;
	double ret = this->distance(b);
	if (ret < epsilon) return 1;
	else return 0;
}

void VectorND::dilate(VectorND head, VectorND source, double lambda) {

	this->subVectorND(head, source);
	this->scaleND(lambda);
	this->AddVectorND(source);
}

VectorND dilate3D(VectorND head, VectorND source, double lambda) {

	VectorND ret;
	ret.initVectorND(3, 0.0, 0.0, 0.0);
	ret.subVectorND(head, source);
	ret.scaleND(lambda);
	ret.AddVectorND(source);

	return ret;
}

VectorND midPoint(VectorND head, VectorND source) {

	VectorND ret, up;
	ret.initVectorND(3, 0.0, 0.0, 0.0);
	up.initVectorND(3, 0.0, 0.25, 0.0);
	ret.subVectorND(head, source);
	ret.scaleND(0.35);
	ret.AddVectorND(source);
	ret.AddVectorND(up);
	return ret;
}

double& VectorND::operator [] (int index)
{   return this->A[index%(this->n)];    }

double  VectorND::operator [] (int index) const
{   return this->A[index%(this->n)];    }













//================CLASS VEC_ARRAY========================


VecArray::VecArray(initializer_list<VectorND> List){
    this->n = List.size();
    this->A = (VectorND*) malloc (List.size() * sizeof(VectorND));
    for(int index=0; index<List.size(); ++ index)
        this->A[index] = List.begin()[index];
}

VecArray::VecArray(int dim){
    if(dim>=0){
        this->n = int(dim);
        this->A = (VectorND*) malloc (dim*sizeof(VectorND));
    }
    else
        cout << "La dimension debe ser un entero positivo.\n";
}

VecArray::~VecArray(void)
{   this->clear();  }


int  VecArray::dim(void)
{   return this->n;     }

void VecArray::clear(void){
    this->n = 0;
    free(this->A);
    this->A = NULL;
}

void VecArray::resize(int dim){
    this->clear();
    this->n = dim;
    this->A = (VectorND*) malloc (n*sizeof(VectorND));
}

bool VecArray::checkDim(void){
    bool Same_Dimension = true;
    for(int index=1; index<this->n; ++index){
        if( this->A[0].dim()!=this->A[index].dim() ){
            Same_Dimension = false;
            break;
        }
    }
    return Same_Dimension;
}

void VecArray::append(VectorND& new_vector){
    if( this->dim()==0 ){
        this->n = 1;
        this->A = (VectorND*) malloc (sizeof(VectorND));
        this->A[0] = new_vector;
    }
    else if( this->dim()>0 ){
        if( this->checkDim() ){
            this->n += 1;
            this->A = (VectorND*) realloc (this->A,this->n*sizeof(VectorND));
            this->A[this->n - 1] = new_vector;
        }
    }
    else
        cout << "Error: " << endl << "Ingresas un vector de dimension distinta" << endl;
}

void VecArray::append(initializer_list<double> List){
    VectorND new_vector (List); 
    if( this->dim()==0 ){
        this->n = 1;
        this->A = (VectorND*) malloc (sizeof(VectorND));
        this->A[0] = new_vector;
    }
    else if( this->dim()>0 ){
        if( this->checkDim() ){
            this->n += 1;
            this->A = (VectorND*) realloc (this->A,this->n*sizeof(VectorND));
            this->A[this->n - 1] = new_vector;
        }
    }
    else
        cout << "Error: " << endl << "Ingresas un vector de dimension distinta" << endl;
}

VectorND& VecArray::operator [] (int index)
{   return this->A[index%(this->n)];    }

VectorND  VecArray::operator [] (int index) const
{   return this->A[index%(this->n)];    }

VecArray& VecArray::operator =  (VecArray& copy_array){
    this->clear();
    this->n = copy_array.dim();
    this->A = (VectorND*) malloc (this->n * sizeof(VectorND));
    for(int index=0; index<copy_array.dim(); ++index)
        this->A[index] = copy_array[index];
    return *this;
}

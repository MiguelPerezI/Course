#include <stdarg.h>

#ifndef VECTORND_H
#define VECTORND_H

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

class VectorND {
    private:
        
        int n;
        double * A;
        double norma;

    public:
        
    /*El metodo genera un espacio de memoria para almacenar las 
    * entradas del vector, dado un arreglo local en el constructor.
    * Se utiliza de la siguiente forma:
    *          VectorND nombre_vector ( {entrada_1, entrada_2, ... } );*/
        VectorND(initializer_list<double> List){
            double suma = 0.0;
            this->n = List.size();
            this->A = (double*) malloc (List.size()*sizeof(double));
            for(int index=0; index<List.size(); ++index){
                this->A[index] = List.begin()[index];
                suma += List.begin()[index];
            }
            this->norma = sqrt(suma);
        };

        /*El metodo pide la dimension del arreglo y genera espacio en 
        * la memoria tal que pueda utilizarse. Este metodo ingresa ceros
        * en las entradas si la dimension es mayor a cero*/
        VectorND(int dim=0){
            this->n = dim;
            this->A = (double*) malloc (dim*sizeof(double));
            for(int index=0; index<dim; ++index)
                this->A[index] = 0.0;
            this->norma = 0.0;
        };

        void initVectorND(int, ...);
        void initVectorNDType2(int);
        int dim();
        double norm();
        void zeroVectorND(int);
        void escVectorND();
        int digitCounter(double);
        void updateVector3D(double, double, double);
        void updateVector3DP(VectorND);



        
        double access(int i);
        void AddVector3D(VectorND);
        VectorND SumVector3D(VectorND, VectorND);
        void SumVectorND(VectorND, VectorND);
        void MinusVector3D(VectorND);
        void SubVector3D(VectorND, VectorND);
        double dotProd(VectorND);

        void subVectorND(VectorND, VectorND);
        void scaleVectorND(double, VectorND);
        void Cross3D(VectorND, VectorND);
        void unit();
        double distance(VectorND);
        int equalVector(VectorND);
        void scaleND(double);
        void AddVectorND(VectorND);
        void dilate(VectorND, VectorND, double);


        double& operator [] (int);
        double  operator [] (int) const;
};

VectorND dilate3D(VectorND head, VectorND source, double lambda);
VectorND midPoint(VectorND, VectorND);

class Vector3D: public VectorND {

    int n;
    double * A;
    double norma;
}; 

class VecArray: public VectorND{
    private:
        int n;
        VectorND* A;

    public:
        VecArray(initializer_list<VectorND>);
        VecArray(int=0);
        ~VecArray(void);

        int  dim   (void);
        void clear (void);
        void resize(int=0);
        bool checkDim(void);
        void append(VectorND&);
        void append(initializer_list<double>);
        
        VectorND& operator [] (int);
        VectorND  operator [] (int) const;
        VecArray& operator =  (VecArray&);
};

#endif

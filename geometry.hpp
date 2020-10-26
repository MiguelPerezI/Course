#include <stdarg.h>
#include "VectorND.hpp"
#include "matrix.hpp"

#ifndef GEOMETRY
#define GEOMETRY

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

class Facet {
    private:
        VectorND a;
        VectorND b;
        VectorND c;
        VectorND aux1;
        VectorND aux2;
        VectorND barycenter;
        VectorND aPrime;
        VectorND bPrime;
        VectorND cPrime;
        VectorND normal;
        
    public:
        void updateNormal();
        void initFacet(VectorND a, VectorND b, VectorND c);
        void renderFacetOpenGL(int mod, RotationMats U);
        void updatePoints(VectorND a, VectorND b, VectorND c);
        void updateFacet(VectorND a, VectorND b, VectorND c);
        VectorND pointA();
        VectorND pointB();
        VectorND pointC();
        VectorND normalK();
        void updateBarycenter();
        void updatePrime();
        VectorND getBarycenter();
        VectorND getPrimeA();
        VectorND getPrimeB();
        VectorND getPrimeC();
};

class Inversion {
    private:
        double radius;
        VectorND gamma;
        double d, check, dd;
        VectorND rest, rest0, rest1;
        VectorND aux, a, b, c;
        Facet facet;
    public:
        void initInversion(double, VectorND);
        void updateInversion(double, VectorND);
        void applyInversion(VectorND);
        VectorND returnInversion(VectorND);
        void facetInversion(Facet);
};

class MatrixFacet {
    private:
        int m, n;
        Facet * * A;
    public:
        void initMatrixFacet(int, int);
        void inputFacet(int, int, VectorND, VectorND, VectorND);
        void updateFacet(int, int, VectorND, VectorND, VectorND);
        Facet getFacet(int, int);
        void renderMatrixFacet(int, RotationMats);
        int getM();
        int getN();
};

class MatrixFacetList {
    private:
        int n;
        MatrixFacet * A;
    public:
        Inversion inversion;
        void initMatrixFacetList(int);
        void inputMatrixFacet(int, MatrixFacet);
        void updateMatrixFacet(int, MatrixFacet);
        void inputMatrixFacet2(int, int, int);
        MatrixFacet getMatrixFacet(int);
        void updateMatrixFacet(int, int, int, VectorND, VectorND, VectorND);
        void updateMatrixFacetInversion(int, int, int, VectorND, VectorND, VectorND);
        int getLength();
        void renderMatrixFacetList(int, int, RotationMats);
};

class CrackFacet {
    private:
        Facet mainFacet;
        MatrixFacet mesh;
    public:
        void initCrackFacet(Facet);
        void updateCrackFacet(Facet);
        Facet getFacet(int);
        void updateSubFacet(int, VectorND, VectorND, VectorND);
};

class Square {
    public:
        VectorND a, b, c, d;
        MatrixFacet mesh;
        VectorND centerM;

        void initSquare(VectorND, VectorND, VectorND, VectorND);
        void renderSquare(int, RotationMats);
};

class MatrixSquare {
    public:
        int m, n;
        Square * * A;
        void initMatrixSquare(int, int);
};

class Sphere {
    public:
        int n;
        double R;
        VectorND center;
        VectorND aux[4];
        MatrixSquare mesh;

        void initSphere(int, double, VectorND);
        void renderSphere(int, RotationMats);
};

class MatrixSphere {
    public:
        int m, n;
        Sphere * * A;

        void initMatrixSphere(int, int);
        void renderMatrixSphere(int, RotationMats);
};

class MatrixSphereList {
    public:
        int m;
        MatrixSphere * matrixSphere;
        void initMatrixSphereList(int);
};

class MatR3 {
    public:
        int m, n;
        VectorND * * A;

        void initMatR3Space(int m, int n);
        void updateA(int i, int j,  double x, double y, double z);
        void escMatR3Space();
        VectorND getA(int, int);
        void push(VectorND);
};

class MatR3List {
    public:
        int m;
        MatR3 * matR3;

        void initMatR3List(int m);
};

class ElectricField {
    public:
        double twoPI = 2*3.141592653589793233846126;
        double ke = 40.0;
        int n;
        MatR3 points;
        MatrixSquare mesh;
        MatR3 r;
        Matrix dQ;
        Matrix dS;
        double dsigma;
        MatrixSphere normals;
        MatR3 electrons;
        MatrixSphere electronsForm;
        Matrix electronCharge;
        double q0;
        MatR3 dirR;
        int numElectrons;
        MatR3List transformAux;
        MatrixSphereList transform;

        void initElectricFieldDistribution(int n, double size);
        void renderElectricField(int mod, RotationMats U);
        void initElectrons(VectorND);
        void induceElecticField(VectorND);
        void superpositionElecticField(double charge, VectorND orig);
};

class Torus {
    public:
        int n;
        double R;
        VectorND center;
        VectorND aux[4];
        MatrixSquare mesh;

        void initTorus(int, double, double, VectorND);
        void renderTorus(int, RotationMats);
};

class Plane {
    public:
        int n;
        double R;
        VectorND center;
        VectorND aux[4];
        MatrixSquare mesh;

        void initPlane(int, double, double, VectorND);
        void renderPlane(int, RotationMats);
};

class MatrixTorus {
    public:
        int m, n;
        Torus * * A;

        void initMatrixTorus(int, int);
        void renderMatrixTorus(int, RotationMats);
};

class CubeNeighborhood {
    private:
        VectorND center;
        
        double radius;
        MatrixFacet mesh;
        int color;
        Inversion systemInversion;
    public:
        VectorND vertex[20];
        void initCubeNeighborhood(VectorND center, double, int);
        CubeNeighborhood copyCube1();
        void copyCube(CubeNeighborhood copyFrom);
        void renderCubeNeighborhood(RotationMats U);
        VectorND getVertex(int);
        VectorND getCenter();
        double getRadius();
        void updateRadius(double);
        void updateCenter(VectorND);
        void applyInversionF(double, VectorND);
        Facet getFacet(int i);
        void applyDilation(VectorND head, double lambda);
};

class MatrixCubeNeighborhood {
    public:
        int m, n;
        CubeNeighborhood * * A;

        void initMatrixCubeNeighborhood(int, int);
        void inputCubeNeighborhood(int, int, VectorND, double, int);
        void updateCubeNeighborhood(int, int, VectorND, VectorND, VectorND);
        CubeNeighborhood getCube(int, int);
        void updateRadiusMatrix(int, int, double);
        void updateInversionMatrix(int, int, double, VectorND);
        void updateInversion(double, VectorND);
        void updateDilation(int i, int j, double newR, VectorND gamma);
        void dilateMatrixCube(double newR, VectorND gamma);
        void dilateMatrixCubeInterval(int b, int f, double newR, VectorND gamma);

};

class MatrixCubeList {
    public:
        int m;
        MatrixCubeNeighborhood * B;

        void initMatrixCubeList(int);
        MatrixCubeNeighborhood getA(int);
        void initA(int i, int m, int n);
        void initInitA(int k, int i, int j, VectorND center, double raius, int mod);
        void initInitCopy(int k, int i, int j, CubeNeighborhood cube);
        void renderMatrixCubeList(RotationMats U);
};

class MengerSponge {
    public:
        int n;
        MatrixCubeList order;
        MatR3 skeleton;
        MatrixSphere v;

        void initMengerSponge(int n, double, VectorND);
        void renderMengerSponge(RotationMats U, int);
};

class MatrixMengerList {
    public:
        int m;
        MengerSponge * B;

        void initMatrixMengerList(int);
        MengerSponge getA(int);
        void initA(int i, int n, double radius, VectorND center);
        void initInitA(int k, int i, int j, VectorND center, double raius, int mod);
        void renderMatrixMengerList(RotationMats U, int);
        void applyInversion(double r, VectorND gamma);
};

class Lattice2D { 
    private:
        int m, n;
        VectorND center;
        double radius;
        VectorND pivote;
        MatrixCubeNeighborhood A;
        MatrixCubeNeighborhood A0;
        MatrixCubeNeighborhood A1;
        MatrixCubeNeighborhood A2;
    public:
        void initLattice2D(int, int, VectorND, double);
        CubeNeighborhood getCell(int, int);
        CubeNeighborhood getCell0(int, int);
        CubeNeighborhood getCell1(int, int);
        CubeNeighborhood getCell2(int, int);
        int getM();
        int getN();
};

class Arista {
    public:
        VectorND head;
        VectorND tail;
        VectorND diff;
        double length;

        MatrixFacet mesh;
        VectorND points[11];
        Facet orientationPlane;
        RotationMats U;
    
        void initArista(VectorND, VectorND);
        void initMesh(double);
        void renderMesh(int, RotationMats);
        Facet getOrientationPlane();
        VectorND getHead();
        VectorND getTail();
        double getLength();
        void escArista();
        void updateArista(VectorND, VectorND);
};

class MatrixArista {
    public: 
        int m, n;
        Arista * * A;
        void initMatrixArista(int, int);
        void inputArista(int, int, VectorND, VectorND);
        void updateArista(int, int, VectorND, VectorND);
        Arista getArista(int i, int j);
};

class MatrixAristaList {
        
    public:
        int n;
        MatrixArista * A;
        void initMatrixAristaList(int);
        void inputMatrixArista(int, MatrixArista);
};

#endif

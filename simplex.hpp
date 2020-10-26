#include <stdarg.h>
#include "VectorND.hpp"
#include "geometry.hpp"
#include "matrix.hpp"

#ifndef SIMPLEX
#define SIMPLEX

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

//////////////////////////////////////
//									//
//	 ESTRUCTURA DE ALMACENAMIENTO	//
//									//
//////////////////////////////////////

class MatrixInt {
    private:
    	int m, n;
        int * * A;
        
    public:
        void initMatrixInt(int, int);
        void zeroMatrixInt(int, int);
        int getM();
        int getN();
        int getA(int, int);
        void updateA(int, int, int);
        void makeEqualMatrixInt(MatrixInt);
        void escMatrixInt();
};

//////////////////////////////////////
//									//
//	 ESTRUCTURA DE SIMPLEJOS		//
//									//
//////////////////////////////////////

class Simplex {
	private:
		int dim;
		MatrixInt A;
	public:
		void initSimplex(int);
		int getDim();
		void escSimplex();
		void makeEqualSimplex(Simplex);
		void zeroSimplex();
		void simplexProd(Simplex, Simplex);
};

class MatrixSimplex {
	private:
		int m, n;
		Simplex * * A;
	public:
		void initMatrixSimplex(int, int);
		void zeroMatrixSimplex(int, int);
		int getM();
		int getN();
		Simplex getA(int, int);
		void updateA(int, int, Simplex);
		void makeEqualMatrixSimplex(MatrixSimplex);
		void escMatrixSimplex();
};

class VectorInt {
	private:
		int x, y;
	public:
		void initVectorInt(int, int);
		int getX();
		int getY();
		void updateVectorInt(int, int);
		void updateX(int);
		void updateY(int);
		void escVectorInt();
		void makeEqual(VectorInt);
		void zeroVectorInt();
		void upOneX();
		void upOneY();
		void upOneXY();
		void downOneX();
		void downOneY();
		void downOneXY();
		void upCero();
};

class MatrixVectorInt {
	private:
		int m, n;
	public:
		VectorInt * * A;
		void initMatrixVectorInt(int, int);
        void zeroMatrixVectorInt(int, int);
        int getM();
        int getN();
        VectorInt getA(int, int);
        void updateA(int, int, VectorInt);
        void makeEqualMatrixVectorInt(MatrixVectorInt);
        void escMatrixVectorInt();
        void updateAX(int, int, int);
        void updateAY(int, int, int);
};

class SimplexND {
	private:
		int dim, numVertex, id;
	public:
		MatrixVectorInt A;
		void initSimplexND(int, ...);
		void initSimplexNDM(MatrixVectorInt);
		void initSimplexNDP(int n);
		int getId();
		void updateId(int);
		int getDim();
		int getNumVertex();
		void escSimplexND();
		void makeEqual(SimplexND);
		void zeroSimplexND();
		void simplexNDProd(SimplexND, SimplexND);
		void updateNest(int, VectorInt);
		VectorInt getA(int);
};

class MatrixVectorND {
	
	public:
		int m, n;
		VectorND * * A;
		void initMatrixVectorND(int, int);
		void updateA(int, int, VectorND);
		void zeroMatrixVector3D();
		void vectorArray3D(int);
		void escVectorArray3D();
		VectorND getA(int, int);
		void updateVectorArray3D(int, VectorND);
};

class MatrixVectorNDList {
        
    public:
        int n;
        MatrixVectorND * A;
        void initMatrixVectorNDList(int);
        void inputMatrixVectorND(int, MatrixArista);
};

class SimplexRND {
	private:
		int dim, numVertex, id;
	public:
		MatrixVectorND A, B;
		MatrixFacet mesh; 
		void initSimplex3D(int, VectorND, VectorND, VectorND, VectorND);
		void updateSimplex3D(VectorND, VectorND, VectorND, VectorND);
		void renderSimplex3D(int, RotationMats);
		void dilateSimplexRND(VectorND, double);
		SimplexRND dilateSimplexRNDType2(VectorND, double);
};

class MatrixSimplexRND {
	public:
		int m, n;
		SimplexRND * * A;
		void initMatrixSimplexRND(int, int);
		void initA(int, int, SimplexRND);
		void updateA(int, int, SimplexRND);
		void renderMatrixSimplexRND(int, RotationMats);
};

class ListMatrixSimplexRND {
	public:
		int m;
		MatrixSimplexRND * List;
		void initListMatrixSimplexRND(int);
		void initList(int, MatrixSimplexRND);
		void initListA(int, int, int, SimplexRND);
		void renderListMatrixSimplexRND(int, RotationMats);
};

class MatrixSimplexND {
	private:
		int m, n;
	public:
		SimplexND * * A;
		void initMatrixSimplexND(int, int);
        void zeroMatrixSimplexND(int, int);
        int getM();
        int getN();
        SimplexND getA(int, int);
        void updateA(int, int, SimplexND);
        void makeEqualMatrixSimplexND(MatrixSimplexND);
        void escMatrixSimplexND();
        void updateId(int, int, int);
        void initA(int, int, SimplexND);
};

class TetahedronTree {
	public:
		int n;
		VectorND a0, a1, a2, a3;
		SimplexRND simplex;
		ListMatrixSimplexRND storage;

		void initTetahedronTree(int, SimplexRND);
		void renderTetahedronTree(int, RotationMats);
};

class TreeNode {		
	public:
		VectorInt data;
		SimplexRND simplex;
		VectorND dataPrime;
		TreeNode * branch_1;
		TreeNode * branch_2;
		TreeNode * branch_3;
		int n;
		int boundX;
		int boundY;
		MatrixVectorInt * image;
		TetahedronTree tetaTree;

		void initTreeNode(VectorInt);
		void copyTree(TreeNode*);
		void push(VectorInt);
		void pop();
		void escTreeNode();
		void initMemorySpaceTreeNode();
		void equalBranches(TreeNode*);
		void equalData(TreeNode * head);
		void simplicialGene();
		void plantSeed(VectorInt, SimplexRND);
		void growX();
		void growY();
		void grow(int, TreeNode*);
		MatrixVectorInt initPath(MatrixInt);
		TreeNode * localMapChoice(int);
};

TreeNode * auxFunctions(TreeNode*);

class Complex3x3 {
    public:
    	TreeNode * arbol;
    	MatrixVectorND keyGen;
    	MatrixVectorND skeleton0;
    	MatrixSphere skeleton0M;
    	SimplexRND simplex0;
    	SimplexRND simplex1;
    	SimplexRND simplex2;
    	SimplexRND simplex3;

    	void initComplex3x3(VectorInt);
    	void renderVerticies(int mod, RotationMats U);
};

class Point {
    public:
        void renderPoint(VectorND, RotationMats);
};

#endif
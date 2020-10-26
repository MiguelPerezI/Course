#include <iostream>
#include <cstdlib>
#include "VectorND.hpp"
#include <math.h>
#include <string>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "simplex.hpp"
#include <GL/glut.h>
#include "geometry.hpp"
#include "matrix.hpp"

//////////////////////////////////////
//									//
//	 ESTRUCTURA DE ALMACENAMIENTO	//
//									//
//////////////////////////////////////

void MatrixInt::initMatrixInt(int m, int n) {

	this->m = m;
	this->n = n;
	this->A = (int **) malloc (m * sizeof(int*));
	for (int i = 0; i < m; i++)
		this->A[i] = (int *) malloc (n * sizeof(int));

}

int MatrixInt::getM() {
	return this->m;
}

int MatrixInt::getN() {
	return this->n;
}

int MatrixInt::getA(int i, int j) {
	return this->A[i][j];
}

void MatrixInt::updateA(int i, int j, int x) {
	this->A[i][j] = x;
}

void MatrixInt::makeEqualMatrixInt(MatrixInt mat) {

	if (this->m == mat.getM() && this->n == mat.getN()) {

		for (int i = 0; i < this->m; i++)
			for (int j = 0; j < this->n; j++)
				this->updateA(i, j, mat.getA(i, j));
	} else {
		printf("\n\n--> In makeEqualMatrixInt :: unable to make equal\n--> Please check dimensions\n\n");
	}
}

void MatrixInt::zeroMatrixInt(int m, int n) {

	this->initMatrixInt(m, n);
	for (int i = 0; i < this->m; ++i)
		for (int j = 0; j < this->n; ++j)
			this->A[i][j] = 0;
}

void MatrixInt::escMatrixInt() {

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++)
			printf("   %d", this->getA(i, j));
		printf("\n\n");
	}
}

//////////////////////////////////////
//									//
//	 ESTRUCTURA DE SIMPLEJOS		//
//									//
//////////////////////////////////////

void Simplex::initSimplex(int dim) {

	this->dim = dim;
	this->A.initMatrixInt(1, dim);

	for (int i = 0; i < dim; i++)
		this->A.updateA(0, i, i);
}

int Simplex::getDim() {
	return this->dim;
}

void Simplex::escSimplex() {
	
	printf("[");
	for (int i = 0; i < this->dim; i++){
		printf("  %d", this->A.getA(0, i));
	}
	printf("  ]");
}

void Simplex::makeEqualSimplex(Simplex sim) {

	this->dim = sim.getDim();
	this->A.makeEqualMatrixInt(sim.A);
}

void Simplex::zeroSimplex() {
	this->dim = 0;
	this->A.initMatrixInt(1, 1);
	this->A.updateA(0, 0, 0);
}

void Simplex::simplexProd(Simplex sig0, Simplex sig1) {

	/*SKETCH*/

	MatrixInt prodData;
	prodData.initMatrixInt(sig0.A.getN() * sig1.A.getN(), 2);
	int count = 0;
	for (int i = 0; i < sig0.A.getN(); i++)
		for (int j = 0; j < sig1.A.getN(); j++) {
			prodData.updateA(count, 0, sig0.A.getA(0, i));
			prodData.updateA(count, 1, sig1.A.getA(0, j));
			count += 1;
		}
		printf("\n");
	prodData.escMatrixInt();

	/*SKETCH OF SIMPLEX GENERATOR*/

	//Starting point
	int x = 0;
	int y = 0;

	//Ending Point
	int xf = sig0.A.getN() - 1;
	int yf = sig1.A.getN() - 1;

}

//////////////////////////////////////
//									//
//	 ESTRUCTURA DE ALMACENAMIENTO	//
//									//
//////////////////////////////////////

void MatrixSimplex::initMatrixSimplex(int m, int n) {

	this->m = m;
	this->n = n;
	this->A = (Simplex **) malloc (m * sizeof(Simplex*));
	for (int i = 0; i < m; i++)
		this->A[i] = (Simplex *) malloc (n * sizeof(Simplex));

}

int MatrixSimplex::getM() {
	return this->m;
}

int MatrixSimplex::getN() {
	return this->n;
}

Simplex MatrixSimplex::getA(int i, int j) {
	return this->A[i][j];
}

void MatrixSimplex::updateA(int i, int j, Simplex x) {
	this->A[i][j].makeEqualSimplex(x);
}

void MatrixSimplex::makeEqualMatrixSimplex(MatrixSimplex mat) {

	if (this->m == mat.getM() && this->n == mat.getN()) {

		for (int i = 0; i < this->m; i++)
			for (int j = 0; j < this->n; j++)
				this->updateA(i, j, mat.getA(i, j));
	} else {
		printf("\n\n--> In makeEqualMatrixSimplex :: unable to make equal\n--> Please check dimensions\n\n");
	}
}

void MatrixSimplex::zeroMatrixSimplex(int m, int n) {

	this->initMatrixSimplex(m, n);
	for (int i = 0; i < this->m; ++i)
		for (int j = 0; j < this->n; ++j)
			this->A[i][j].zeroSimplex();
}

void MatrixSimplex::escMatrixSimplex() {

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++){
			printf("\n");
			this->A[i][j].escSimplex();
		}
		printf("\n\n");
	}
}

//////////////////////////////////////
//									//
//	 		VECTOR INT				//
//									//
//////////////////////////////////////

void VectorInt::initVectorInt(int x, int y) {
	this->x = x;
	this->y = y;
}

int VectorInt::getX() {
	return this->x;
}

int VectorInt::getY() {
	return this->y;
}

void VectorInt::updateVectorInt(int newX, int newY) {

	this->x = newX;
	this->y = newY;
}

void VectorInt::updateX(int newX) {
	this->x = newX;
}

void VectorInt::updateY(int newY) {
	this->y = newY;
}

void VectorInt::escVectorInt() {

	printf("(%d, %d)", this->x, this->y);
}

void VectorInt::makeEqual(VectorInt v) {
	
	this->x = v.getX();
	this->y = v.getY();
}

void VectorInt::zeroVectorInt() {
	this->initVectorInt(0, 0);
}

void VectorInt::upOneX() {
	this->updateX(this->getX() + 1);
}

void VectorInt::upOneY() {
	this->updateY(this->getY() + 1);
}

void VectorInt::upOneXY() {
	this->updateX(this->getX() + 1);
	this->updateY(this->getY() + 1);
}

void VectorInt::downOneX() {
	this->updateX(this->getX() - 1);
}

void VectorInt::downOneY() {
	this->updateY(this->getY() - 1);
}

void VectorInt::downOneXY() {
	this->updateX(this->getX() - 1);
	this->updateY(this->getY() - 1);
}

void VectorInt::upCero() {}

//////////////////////////////////////
//									//
//	 ESTRUCTURA DE ALMACENAMIENTO	//
//									//
//////////////////////////////////////

void MatrixVectorInt::initMatrixVectorInt(int m, int n) {

	this->m = m;
	this->n = n;
	this->A = (VectorInt **) malloc (m * sizeof(VectorInt*));
	for (int i = 0; i < m; i++)
		this->A[i] = (VectorInt *) malloc (n * sizeof(VectorInt));

}

int MatrixVectorInt::getM() {
	return this->m;
}

int MatrixVectorInt::getN() {
	return this->n;
}

VectorInt MatrixVectorInt::getA(int i, int j) {
	return this->A[i][j];
}

void MatrixVectorInt::updateA(int i, int j, VectorInt x) {
	this->A[i][j].makeEqual(x);
}

void MatrixVectorInt::makeEqualMatrixVectorInt(MatrixVectorInt mat) {

	if (this->m == mat.getM() && this->n == mat.getN()) {

		for (int i = 0; i < this->m; i++)
			for (int j = 0; j < this->n; j++)
				this->updateA(i, j, mat.getA(i, j));
	} else {
		printf("\n\n--> In makeEqualMatrixVectorInt :: unable to make equal\n--> Please check dimensions\n\n");
	}
}

void MatrixVectorInt::zeroMatrixVectorInt(int m, int n) {

	this->initMatrixVectorInt(m, n);
	for (int i = 0; i < this->m; ++i)
		for (int j = 0; j < this->n; ++j)
			this->A[i][j].zeroVectorInt();
}

void MatrixVectorInt::escMatrixVectorInt() {

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++) {
			printf("\n");
			this->A[i][j].escVectorInt();
		}
	}
}

void MatrixVectorInt::updateAX(int i, int j, int X) {
	this->A[i][j].updateX(X);
}

void MatrixVectorInt::updateAY(int i, int j, int X) {
	this->A[i][j].updateY(X);
}

//////////////////////////////////////
//									//
//	 ESTRUCTURA DE SIMPLEJOS_ND		//
//									//
//////////////////////////////////////

void SimplexND::initSimplexND(int n, ...) {

	va_list list;
	this->id = 0;
	this->dim = n - 1;
	this->numVertex = n;
	this->A.initMatrixVectorInt(1, n);

	va_start(list, n);

		for (int i = 0; i < n; i++)
			this->A.updateA(0, i, va_arg(list, VectorInt));
	va_end(list);

	for (int i = 0; i < this->A.getN(); i++)
		this->A.updateAX(0, i, this->id);
}

void MatrixVectorND::initMatrixVectorND(int m, int n) {

	this->m = m;
	this->n = n;

	this->A = (VectorND * *) malloc (m * sizeof(VectorND*));
	for (int i = 0; i < m; i++)
		this->A[i] = (VectorND *) malloc (n * sizeof(VectorND));
}

void MatrixVectorNDList::initMatrixVectorNDList(int n) {

	this->n = n;
	this->A = (MatrixVectorND *) malloc (n * sizeof(MatrixVectorND));	
}

void MatrixVectorND::zeroMatrixVector3D() {
	
	this->initMatrixVectorND(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			this->A[i][j].initVectorND(3, 0.0, 0.0, 0.0);
}

void MatrixVectorND::vectorArray3D(int size) {

	this->initMatrixVectorND(1, size);
	for (int i = 0; i < size; i++)
		this->A[0][i].initVectorND(3, 0.0, 0.0, 0.0);
}

void MatrixVectorND::updateVectorArray3D(int i, VectorND v) {

	this->A[0][i].updateVector3DP(v);
}

void MatrixVectorND::escVectorArray3D() {
	
	for (int i = 0; i < this->n; i++) {
		printf("\n");
		this->A[0][i].escVectorND();
	}
}

VectorND MatrixVectorND::getA(int i, int j) {
	return this->A[i][j];
}

void SimplexRND::initSimplex3D(int id, VectorND a0, VectorND a1, VectorND a2, VectorND a3) {

	va_list list;
	this->id = 0;
	this->dim = 3;
	this->numVertex = 4;
	this->A.vectorArray3D(4);
	this->B.vectorArray3D(4);
	this->A.updateVectorArray3D(0, a0);
	this->A.updateVectorArray3D(1, a1);
	this->A.updateVectorArray3D(2, a2);
	this->A.updateVectorArray3D(3, a3);

	this->B.updateVectorArray3D(0, a0);
	this->B.updateVectorArray3D(1, a1);
	this->B.updateVectorArray3D(2, a2);
	this->B.updateVectorArray3D(3, a3);

	this->mesh.initMatrixFacet(1, 4);
	this->mesh.inputFacet(0, 0, a0, a1, a2);
	this->mesh.inputFacet(0, 1, a0, a2, a3);
	this->mesh.inputFacet(0, 2, a0, a3, a1);
	this->mesh.inputFacet(0, 3, a1, a3, a2);
}

void SimplexRND::updateSimplex3D(VectorND a0, VectorND a1, VectorND a2, VectorND a3) {

	this->A.updateVectorArray3D(0, a0);
	this->A.updateVectorArray3D(1, a1);
	this->A.updateVectorArray3D(2, a2);
	this->A.updateVectorArray3D(3, a3);

	this->B.updateVectorArray3D(0, a0);
	this->B.updateVectorArray3D(1, a1);
	this->B.updateVectorArray3D(2, a2);
	this->B.updateVectorArray3D(3, a3);

	this->mesh.inputFacet(0, 0, a0, a1, a2);
	this->mesh.inputFacet(0, 1, a0, a2, a3);
	this->mesh.inputFacet(0, 2, a0, a3, a1);
	this->mesh.inputFacet(0, 3, a1, a3, a2);
}

void SimplexRND::renderSimplex3D(int mod, RotationMats U) {

	for (int i = 0; i < 4; i++)
		this->mesh.getFacet(0, i).renderFacetOpenGL(mod, U);
}

void SimplexRND::dilateSimplexRND(VectorND resp, double lambda) {

	for (int i = 0; i < 4; i++)
		this->A.A[0][i].dilate(this->A.A[0][i], resp, lambda);

	this->mesh.updateFacet(0, 0, this->A.A[0][0], this->A.A[0][1], this->A.A[0][2]);
	this->mesh.updateFacet(0, 1, this->A.A[0][0], this->A.A[0][2], this->A.A[0][3]);
	this->mesh.updateFacet(0, 2, this->A.A[0][0], this->A.A[0][3], this->A.A[0][1]);
	this->mesh.updateFacet(0, 3, this->A.A[0][1], this->A.A[0][3], this->A.A[0][2]);
}

SimplexRND SimplexRND::dilateSimplexRNDType2(VectorND resp, double lambda) {

	SimplexRND ret;
	ret.initSimplex3D(0, 
						this->A.A[0][0],
						this->A.A[0][1],
						this->A.A[0][2],
						this->A.A[0][3]);

	for (int i = 0; i < 4; i++)
		ret.A.A[0][i].dilate(ret.A.A[0][i], resp, lambda);

	ret.mesh.updateFacet(0, 0, ret.A.A[0][0], ret.A.A[0][1], ret.A.A[0][2]);
	ret.mesh.updateFacet(0, 1, ret.A.A[0][0], ret.A.A[0][2], ret.A.A[0][3]);
	ret.mesh.updateFacet(0, 2, ret.A.A[0][0], ret.A.A[0][3], ret.A.A[0][1]);
	ret.mesh.updateFacet(0, 3, ret.A.A[0][1], ret.A.A[0][3], ret.A.A[0][2]);

	return ret;
}

void SimplexND::initSimplexNDM(MatrixVectorInt mat) {

	this->dim = mat.getN() - 1;
	this->numVertex = mat.getN();
	this->id =-1;
	this->A.initMatrixVectorInt(mat.getM(), mat.getN());
	this->A.makeEqualMatrixVectorInt(mat);
}

void SimplexND::initSimplexNDP(int n) {

	this->id = 0;
	this->dim = n;
	this->numVertex = n + 1;
	this->A.initMatrixVectorInt(1, n + 1);
}

int SimplexND::getId() {
	return this->id;
}

void SimplexND::updateId(int newId) {
	this->id = newId;
	for (int i = 0; i < this->A.getN(); i++)
		this->A.updateAX(0, i, this->id);
}

int SimplexND::getDim() {
	return this->dim;
}

int SimplexND::getNumVertex() {
	return this->numVertex;
}

void SimplexND::escSimplexND() {

	printf("\nSimplexND Info:\n");
	printf("id -- > %d\n", this->id);
	printf("Dimension --> %d\n", this->dim);
	printf("Number of Vertices --> %d\n", this->numVertex);
	this->A.escMatrixVectorInt();
}

void SimplexND::makeEqual(SimplexND sigma) {

	this->dim = sigma.getDim();
	this->id = sigma.getId();
}

void SimplexND::zeroSimplexND() {
	this->dim = 0;
	this->numVertex = 1;
	this->A.initMatrixVectorInt(1, 1);
}

void SimplexND::updateNest(int i, VectorInt v) {

	this->A.updateA(0, i, v);
}

VectorInt SimplexND::getA(int i) {
	return this->A.getA(0, i);
}

//////////////////////////////////////
//									//
//	 ESTRUCTURA DE ALMACENAMIENTO	//
//									//
//////////////////////////////////////

void MatrixSimplexND::initMatrixSimplexND(int m, int n) {

	this->m = m;
	this->n = n;
	this->A = (SimplexND **) malloc (m * sizeof(SimplexND*));
	for (int i = 0; i < m; i++)
		this->A[i] = (SimplexND *) malloc (n * sizeof(SimplexND));
}

int MatrixSimplexND::getM() {
	return this->m;
}

int MatrixSimplexND::getN() {
	return this->n;
}

SimplexND MatrixSimplexND::getA(int i, int j) {
	return this->A[i][j];
}

void MatrixSimplexND::initA(int i, int j, SimplexND x) {
	
	this->A[i][j].initSimplexNDP(x.getDim());
	for (int k = 0; k < x.getNumVertex(); k++)
		this->A[i][j].updateNest(k, x.getA(k));
}

void MatrixSimplexND::updateA(int i, int j, SimplexND x) {
	this->A[i][j].makeEqual(x);
}

void MatrixSimplexND::makeEqualMatrixSimplexND(MatrixSimplexND mat) {

	if (this->m == mat.getM() && this->n == mat.getN()) {

		for (int i = 0; i < this->m; i++)
			for (int j = 0; j < this->n; j++)
				this->updateA(i, j, mat.getA(i, j));
	} else {
		printf("\n\n--> In makeEqualMatrixSimplexND :: unable to make equal\n--> Please check dimensions\n\n");
	}
}

void MatrixSimplexND::zeroMatrixSimplexND(int m, int n) {

	this->initMatrixSimplexND(m, n);
	for (int i = 0; i < this->m; ++i)
		for (int j = 0; j < this->n; ++j)
			this->A[i][j].zeroSimplexND();
}

void MatrixSimplexND::escMatrixSimplexND() {

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++) {
			printf("\n");
			this->A[i][j].escSimplexND();
		}
	}
}

void MatrixSimplexND::updateId(int i, int j, int X) {
	this->A[i][j].updateId(X);
}

void MatrixSimplexRND::initMatrixSimplexRND(int m, int n) {

	this->m = m;
	this->n = n;
	this->A = (SimplexRND **) malloc (m * sizeof(SimplexRND*));
	for (int i = 0; i < m; i++)
		this->A[i] = (SimplexRND *) malloc (n * sizeof(SimplexRND));
}

void MatrixSimplexRND::initA(int i, int j, SimplexRND simplex) {

	this->A[i][j].initSimplex3D(0, 
									simplex.A.A[0][0], 
									simplex.A.A[0][1], 
									simplex.A.A[0][2], 
									simplex.A.A[0][3]);
}

void MatrixSimplexRND::updateA(int i, int j, SimplexRND simplex) {

	this->A[i][j].updateSimplex3D(
									simplex.A.A[0][0], 
									simplex.A.A[0][1], 
									simplex.A.A[0][2], 
									simplex.A.A[0][3]);
}

void MatrixSimplexRND::renderMatrixSimplexRND(int mod, RotationMats U) {

	for (int i = 0; i < this->m; i++)
		for (int j = 0; j < this->n; j++)
			this->A[i][j].renderSimplex3D(mod, U);
}

void ListMatrixSimplexRND::initListMatrixSimplexRND(int m) {

	this->m = m;
	this->List = (MatrixSimplexRND *) malloc (m * sizeof(MatrixSimplexRND));
}

void ListMatrixSimplexRND::initList(int k, MatrixSimplexRND mat) {

	this->List[k].initMatrixSimplexRND(mat.m, mat.n);
	
	for (int i = 0; i < mat.m; i++)
		for (int j = 0; j < mat.n; j++)
			this->List[k].initA(i, j, mat.A[i][j]);
}

void ListMatrixSimplexRND::renderListMatrixSimplexRND(int mod, RotationMats U) {

	for (int i = 0; i < this->m; i++)
		this->List[i].renderMatrixSimplexRND(mod, U);	
}

void TetahedronTree::initTetahedronTree(int n, SimplexRND simplex) {

	this->simplex.initSimplex3D(0, simplex.A.A[0][0], simplex.A.A[0][1], simplex.A.A[0][2], simplex.A.A[0][3]);
	this->storage.initListMatrixSimplexRND(n + 1);
	this->a0.initVectorND(3, simplex.A.A[0][0].access(0), simplex.A.A[0][0].access(1), simplex.A.A[0][0].access(2));
	this->a1.initVectorND(3, simplex.A.A[0][1].access(0), simplex.A.A[0][1].access(1), simplex.A.A[0][1].access(2));
	this->a2.initVectorND(3, simplex.A.A[0][2].access(0), simplex.A.A[0][2].access(1), simplex.A.A[0][2].access(2));
	this->a3.initVectorND(3, simplex.A.A[0][3].access(0), simplex.A.A[0][3].access(1), simplex.A.A[0][3].access(2));

	this->simplex.dilateSimplexRND(this->a0, 0.5);
	this->storage.List[0].initMatrixSimplexRND(1, 1);
	//printf("\n\n------------------------>>>>>\nM0 contains %d simplexRND\n", 1);
	this->storage.List[0].A[0][0].initSimplex3D(0, this->simplex.A.A[0][0], this->simplex.A.A[0][1], this->simplex.A.A[0][2], this->simplex.A.A[0][3]);
	//printf("Saved simplexRND in M0\n");
	int cycle = 1;
	int pow = 1;
	if (n > 0) {

		while (cycle < n + 1) {
			pow *= 3;
			this->storage.List[cycle].initMatrixSimplexRND(1, pow);
			//printf("\nM%d contains %d simplexRND\n\n", cycle, pow);

			//printf("-->%d elements in M%d\n", this->storage.List[cycle - 1].n, cycle-1);
			int mod = 0;
			for (int i = 0; i < pow/3; i ++) {
				//printf("-------->count %d\n", i);
				this->storage.List[cycle].A[0][mod+0] = this->storage.List[cycle-1].A[0][i].dilateSimplexRNDType2(this->a1, 0.5);
				//printf("Saved simplexRND in M%d\n", cycle);
				this->storage.List[cycle].A[0][mod+1] = this->storage.List[cycle-1].A[0][i].dilateSimplexRNDType2(this->a2, 0.5);
				//printf("Saved simplexRND in M%d\n", cycle);
				this->storage.List[cycle].A[0][mod+2] = this->storage.List[cycle-1].A[0][i].dilateSimplexRNDType2(this->a3, 0.5);	
				//printf("Saved simplexRND in M%d\n\n", cycle);

				mod += 3;	
			}

			//pow *= 3;
			cycle += 1;
		}

	}
}

void TetahedronTree::renderTetahedronTree(int mod, RotationMats U) {

	for (int i = 0; i < this->storage.m; i++)
		this->storage.List[i].renderMatrixSimplexRND(mod, U);
}

//////////////////////////////////////
//									//
//	            ARBOL				//
//									//
//////////////////////////////////////

void TreeNode::equalBranches(TreeNode * head) {
	
	head->branch_1 = this->branch_1;
	head->branch_2 = this->branch_2;
}

void TreeNode::equalData(TreeNode * head) {
	head->data.initVectorInt(this->data.getX(), this->data.getY());
}

void TreeNode::simplicialGene() {
	
	this->branch_1->data.upOneX();
	this->branch_2->data.upOneY();
}

void TreeNode::growX() {
	this->data.upOneX();
	this->branch_1->data.upOneX();
	this->branch_2->data.upOneX();
}

void TreeNode::growY() {
	this->data.upOneY();
	this->branch_1->data.upOneY();
	this->branch_2->data.upOneY();
}

void TreeNode::initTreeNode(VectorInt data) {

	this->data.initVectorInt(data.getX(), data.getY());
	this->initMemorySpaceTreeNode();
	this->n = 1;
}

void TreeNode::plantSeed(VectorInt seed, SimplexRND simplex) {

	this->data.initVectorInt(seed.getX(), seed.getY());
	this->boundX = 2;
	this->boundY = 2;
	this->n = 0;
	this->initMemorySpaceTreeNode();
	this->tetaTree.initTetahedronTree(5, simplex);
	this->dataPrime.initVectorND(3, 
				this->tetaTree.storage.List[0].A[0][0].A.A[0][0].access(0),
				this->tetaTree.storage.List[0].A[0][0].A.A[0][0].access(1),
				this->tetaTree.storage.List[0].A[0][0].A.A[0][0].access(2));
	this->n += 4;
}

void TreeNode::initMemorySpaceTreeNode() {
	
	this->branch_1 = NULL;
	this->branch_1 = new TreeNode();
	if (this->data.getX() < 3)
		this->branch_1->data.initVectorInt(this->data.getX() + 1, this->data.getY());
	else
		this->branch_1->data.initVectorInt(this->data.getX() + 0, this->data.getY());


	
	this->branch_2 = NULL;
	this->branch_2 = new TreeNode();
	if (this->data.getY() < 3)
		this->branch_2->data.initVectorInt(this->data.getX(), this->data.getY() + 1);
	else
		this->branch_2->data.initVectorInt(this->data.getX(), this->data.getY() + 0);



	this->branch_3 = NULL;
	this->branch_3 = new TreeNode();
	if (this->data.getX() < 3 && this->data.getY() < 3)
		this->branch_3->data.initVectorInt(this->data.getX() + 1, this->data.getY() + 1);
	else {
		if (this->data.getX() < 3)
			this->branch_3->data.initVectorInt(this->data.getX() + 1, this->data.getY());
		else
			this->branch_3->data.initVectorInt(this->data.getX() + 0, this->data.getY());

		if (this->data.getY() < 3)
			this->branch_3->data.initVectorInt(this->data.getX() + 0, this->data.getY() + 1);
		else
			this->branch_3->data.initVectorInt(this->data.getX() + 0, this->data.getY());
	}
	this->n = this->n + 3;
}

void TreeNode::copyTree(TreeNode * head) {
	
	this->equalData(head);
	this->equalBranches(head);
	head->n = this->n;
}

void TreeNode::grow(int i, TreeNode * head) {


	TreeNode * tail1 = head->branch_1;
	TreeNode * tail2 = head->branch_2;
	TreeNode * tail3 = head->branch_3;

	if (i != 0) {
		tail1->initMemorySpaceTreeNode();
		tail2->initMemorySpaceTreeNode();
		tail3->initMemorySpaceTreeNode();

		tail1->grow(i - 1, tail1);
		tail2->grow(i - 1, tail2);
		tail3->grow(i - 1, tail3);
	}
}

TreeNode * TreeNode::localMapChoice(int i) {

	TreeNode * ret;

	if (i == 0) ret = this->branch_1;
	if (i == 1) ret = this->branch_2;
	if (i == 2) ret = this->branch_3;

	return ret;
}

TreeNode * auxFunctions(TreeNode * tree) {

	TreeNode copyTree;
	copyTree = {
		data : tree->data,
		branch_1 : tree->branch_1,
		branch_2 : tree->branch_2,
		branch_3 : tree->branch_3
	};

	TreeNode * newTree;
	newTree = &copyTree;

	return newTree;
};

MatrixVectorInt TreeNode::initPath(MatrixInt map) {

	MatrixVectorInt image;
	image.initMatrixVectorInt(1, 4);
	image.updateA(0, 0, this->data);

	TreeNode * auxTree = this;

	int counted = 0;
	while (counted < map.getN()) {

		auxTree = auxTree->localMapChoice(map.getA(0, counted));
		image.updateA(0, counted + 1, auxTree->data);
		counted += 1;
	}
	return image;
}

void TreeNode::pop() {

	if (this->n > 1) {
		this->data.makeEqual(this->branch_1->data);
		this->branch_1 = this->branch_1->branch_1;
		this->n -= 1;
	}
}

void TreeNode::escTreeNode() {
	printf("\n\nSeed.-"); this->data.escVectorInt();
		
	printf("\n\nbranch1.-"); this->branch_1->data.escVectorInt();

		printf("\n\n 	     "); this->branch_1->branch_1->data.escVectorInt();
			printf("\n 	      	     "); this->branch_1->branch_1->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_1->branch_1->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_1->branch_1->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_1->branch_1->branch_3->data.escVectorInt();
			printf("\n 	      	     "); this->branch_1->branch_1->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_1->branch_2->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_1->branch_2->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_1->branch_2->branch_3->data.escVectorInt();
			printf("\n 	      	     "); this->branch_1->branch_1->branch_3->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_1->branch_3->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_1->branch_3->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_1->branch_3->branch_3->data.escVectorInt();

		printf("\n\n 	     "); this->branch_1->branch_2->data.escVectorInt();
			printf("\n 	      	     "); this->branch_1->branch_2->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_2->branch_1->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_2->branch_1->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_2->branch_1->branch_3->data.escVectorInt();

			printf("\n 	      	     "); this->branch_1->branch_2->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_2->branch_2->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_2->branch_2->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_2->branch_2->branch_3->data.escVectorInt();

			printf("\n 	      	     "); this->branch_1->branch_2->branch_3->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_2->branch_3->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_2->branch_3->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_2->branch_3->branch_3->data.escVectorInt();

		printf("\n\n 	     "); this->branch_1->branch_3->data.escVectorInt();
			printf("\n 	      	     "); this->branch_1->branch_3->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_3->branch_1->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_3->branch_1->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_3->branch_1->branch_3->data.escVectorInt();

			printf("\n 	      	     "); this->branch_1->branch_3->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_3->branch_2->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_3->branch_2->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_3->branch_2->branch_3->data.escVectorInt();

			printf("\n 	      	     "); this->branch_1->branch_3->branch_3->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_3->branch_3->branch_1->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_3->branch_3->branch_2->data.escVectorInt();
				printf("\n 	      	     	      	     "); this->branch_1->branch_3->branch_3->branch_3->data.escVectorInt();
	


	printf("\n\nbranch2.-"); this->branch_2->data.escVectorInt();

		printf("\n\n	     "); this->branch_2->branch_1->data.escVectorInt();
			printf("\n 	      	     "); this->branch_2->branch_1->branch_1->data.escVectorInt();
			printf("\n 	      	     "); this->branch_2->branch_1->branch_2->data.escVectorInt();
			printf("\n 	      	     "); this->branch_2->branch_1->branch_3->data.escVectorInt();

		printf("\n\n 	     "); this->branch_2->branch_2->data.escVectorInt();
			printf("\n 	      	     "); this->branch_2->branch_2->branch_1->data.escVectorInt();
			printf("\n 	      	     "); this->branch_2->branch_2->branch_2->data.escVectorInt();
			printf("\n 	      	     "); this->branch_2->branch_2->branch_3->data.escVectorInt();

		printf("\n\n 	     "); this->branch_2->branch_3->data.escVectorInt();
			printf("\n 	      	     "); this->branch_2->branch_3->branch_1->data.escVectorInt();
			printf("\n 	      	     "); this->branch_2->branch_3->branch_2->data.escVectorInt();
			printf("\n 	      	     "); this->branch_2->branch_3->branch_3->data.escVectorInt();



	printf("\n\nbranch3.-"); this->branch_3->data.escVectorInt();
		printf("\n\n 	     "); this->branch_3->branch_1->data.escVectorInt();
			printf("\n 	      	     "); this->branch_3->branch_1->branch_1->data.escVectorInt();
			printf("\n 	      	     "); this->branch_3->branch_1->branch_2->data.escVectorInt();
			printf("\n 	      	     "); this->branch_3->branch_1->branch_3->data.escVectorInt();

		printf("\n\n 	     "); this->branch_3->branch_2->data.escVectorInt();
			printf("\n 	      	     "); this->branch_3->branch_2->branch_1->data.escVectorInt();
			printf("\n 	      	     "); this->branch_3->branch_2->branch_2->data.escVectorInt();
			printf("\n 	      	     "); this->branch_3->branch_2->branch_3->data.escVectorInt();

		printf("\n\n 	     "); this->branch_3->branch_3->data.escVectorInt();
			printf("\n 	      	     "); this->branch_3->branch_3->branch_1->data.escVectorInt();
			printf("\n 	      	     "); this->branch_3->branch_3->branch_2->data.escVectorInt();
			printf("\n 	      	     "); this->branch_3->branch_3->branch_3->data.escVectorInt();
	
}


void Complex3x3::initComplex3x3(VectorInt seed) {

	double ss = 1.0/sqrt(2);
	this->keyGen.vectorArray3D(4);
	this->keyGen.A[0][0].updateVector3D( 2.0, 0.0, -2.0*ss);
	this->keyGen.A[0][1].updateVector3D(-2.0, 0.0, -2.0*ss);
	this->keyGen.A[0][2].updateVector3D( 0.0, 2.0,  2.0*ss);
	this->keyGen.A[0][3].updateVector3D( 0.0,-2.0,  2.0*ss);

	SimplexRND simplex;
	simplex.initSimplex3D(0, this->keyGen.A[0][0], this->keyGen.A[0][1], this->keyGen.A[0][2], this->keyGen.A[0][3]);

	this->arbol = new TreeNode();
	this->arbol->plantSeed(seed, simplex);
	this->arbol->grow(6, this->arbol);
	printf("\nTotal Number of Elements --> %d\n", this->arbol->n);
	this->skeleton0.initMatrixVectorND(4, 4);

	simplex0.initSimplex3D(0, this->keyGen.A[0][0], this->keyGen.A[0][1], this->keyGen.A[0][2], this->keyGen.A[0][3]);
    simplex1.initSimplex3D(1, this->keyGen.A[0][0], this->keyGen.A[0][1], this->keyGen.A[0][2], this->keyGen.A[0][3]);
    simplex2.initSimplex3D(2, this->keyGen.A[0][0], this->keyGen.A[0][1], this->keyGen.A[0][2], this->keyGen.A[0][3]);
    simplex3.initSimplex3D(3, this->keyGen.A[0][0], this->keyGen.A[0][1], this->keyGen.A[0][2], this->keyGen.A[0][3]);
    simplex0.dilateSimplexRND(this->keyGen.A[0][0], 0.33);
    simplex1.dilateSimplexRND(this->keyGen.A[0][1], 0.33);
    simplex2.dilateSimplexRND(this->keyGen.A[0][2], 0.33);
    simplex3.dilateSimplexRND(this->keyGen.A[0][3], 0.33);

    for (int i = 0; i < this->skeleton0.m; i++)
    	for (int j = 0; j < this->skeleton0.n; j++)
    		this->skeleton0.A[i][j].initVectorND(3, 0.0, 0.0, 0.0);

   	this->skeleton0.A[0][0].updateVector3DP(simplex0.A.getA(0, 0));
    this->skeleton0.A[0][1].updateVector3DP(simplex0.A.getA(0, 1));
    this->skeleton0.A[0][2].updateVector3DP(simplex0.A.getA(0, 2));
    this->skeleton0.A[0][3].updateVector3DP(simplex0.A.getA(0, 3));

    this->skeleton0.A[1][0].updateVector3DP(simplex1.A.getA(0, 0));
    this->skeleton0.A[1][1].updateVector3DP(simplex1.A.getA(0, 1));
    this->skeleton0.A[1][2].updateVector3DP(simplex1.A.getA(0, 2));
    this->skeleton0.A[1][3].updateVector3DP(simplex1.A.getA(0, 3));

    this->skeleton0.A[2][0].updateVector3DP(simplex2.A.getA(0, 0));
    this->skeleton0.A[2][1].updateVector3DP(simplex2.A.getA(0, 1));
    this->skeleton0.A[2][2].updateVector3DP(simplex2.A.getA(0, 2));
    this->skeleton0.A[2][3].updateVector3DP(simplex2.A.getA(0, 3));

    this->skeleton0.A[3][0].updateVector3DP(simplex3.A.getA(0, 0));
    this->skeleton0.A[3][1].updateVector3DP(simplex3.A.getA(0, 1));
    this->skeleton0.A[3][2].updateVector3DP(simplex3.A.getA(0, 2));
    this->skeleton0.A[3][3].updateVector3DP(simplex3.A.getA(0, 3));

    this->skeleton0M.initMatrixSphere(4, 4);

    for (int i = 0; i < this->skeleton0M.m; i++)
    	for (int j = 0; j < this->skeleton0M.n; j++)
    		this->skeleton0M.A[i][j].initSphere(20, 0.1, this->skeleton0.A[i][j]);
}

void Complex3x3::renderVerticies(int mod, RotationMats U) {

	for (int i = 0; i < this->skeleton0M.m; i++)
    	for (int j = 0; j < this->skeleton0M.n; j++)
    		this->skeleton0M.A[i][j].renderSphere(mod, U);
}
/*
void TreeNode::equalBranches(TreeNode * head) {
	head->branch_1 = this->branch_1;
	head->branch_2 = this->branch_2;
}

void TreeNode::equalData(TreeNode * head) {
	head->data.initVectorInt(this->data.getX(), this->data.getY());
}

void TreeNode::initMemorySpaceTreeNode() {
	this->branch_1 = NULL;
	this->branch_1 = new TreeNode();
	this->branch_2 = NULL;
	this->branch_2 = new TreeNode();
}

void TreeNode::initTreeNode(VectorInt data) {

	this->data.initVectorInt(data.getX(), data.getY());
	this->initMemorySpaceTreeNode();
	this->n = 1;
}

void TreeNode::copyTree(TreeNode * head) {
	
	this->equalData(head);
	this->equalBranches(head);
	head->n = this->n;
}

void TreeNode::push(VectorInt data) {

	TreeNode * head = NULL;
	head = new TreeNode();
	this->copyTree(head);
	this->data.makeEqual(data);
	this->branch_1 = head;
	this->branch_2 = head;
	head = NULL;
	this->n += 1;
}

void TreeNode::pop() {

	if (this->n > 1) {
		this->data.makeEqual(this->branch_1->data);
		this->branch_1 = this->branch_1->branch_1;
		this->n -= 1;
	}
}

void TreeNode::escTreeNode() {
	printf("\n\n");
	int counted = 1;
	TreeNode * head = NULL;
	head = new TreeNode();
	this->copyTree(head);
	while (counted < this->n + 1) {
		printf("%d.- Data", counted); 
		head->data.escVectorInt();
		//head->branch_1->data.escVectorInt();
		//head->branch_2->data.escVectorInt();
		head->pop();
		printf("\n");
		counted += 1;
	}
}
*/

void Point::renderPoint(VectorND v, RotationMats U) {

	glPointSize(22.7);
   	glBegin(GL_POINTS);
   	glColor3f(0, 0, 0);
   	U.rot3D(v);
   	glVertex3f(U.getAux().access(0), U.getAux().access(1), U.getAux().access(2));
   	glEnd();
}
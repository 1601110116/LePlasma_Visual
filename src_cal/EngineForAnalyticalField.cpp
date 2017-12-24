//
// Created by ylang on 17-12-22.
//

//　This Engine requires OMP and MPI parallels to be turned off

#include "EngineForAnalyticalField.h"
#include "Grid.h"
#include "CaseOfDipole.h"
#include "LePlasma.h"

Tensor3D EngineForAnalyticalField::gradAdt = Tensor3D();
EngineForAnalyticalField::EngineForAnalyticalField(Grid *grid, double deltaT, map<string, double> &units):\
 Engine(grid, deltaT, units){

}

EngineForAnalyticalField::~EngineForAnalyticalField() {

}

void EngineForAnalyticalField::update(const Range &range){

    updateP(range);
    updateX(range);
}

void EngineForAnalyticalField::calcA(Particle *p) {
    p->A.x = -units["B0r03/c"]*p->X.y*tmp3;
    p->A.y = units["B0r03/c"]*p->X.x*tmp3;
    p->A.z = 0;

    // UniformB
//    p->A.x = -1.0e-2*units["gauss"]*y/2/units["c"];
//    p->A.y = 1.0e-2*units["gauss"]*x/2/units["c"];
}

void EngineForAnalyticalField::calcGradAdt(Particle *p) {
    gradAdt.x.x = 3*units["B0r03/c"]*p->X.x*p->X.y*tmp2*deltaT;
    gradAdt.x.y = units["B0r03/c"]*(tmp3 - 3*Square(p->X.x)*tmp2)*deltaT;
    gradAdt.x.z = 0;
    gradAdt.y.x = -units["B0r03/c"]*(tmp3 - 3*Square(p->X.y)*tmp2)*deltaT;
    gradAdt.y.y = -3*units["B0r03/c"]*p->X.x*p->X.y*tmp2*deltaT;
    gradAdt.y.z = 0;
    gradAdt.z.x = 3*units["B0r03/c"]*p->X.y*p->X.z*tmp2*deltaT;
    gradAdt.z.y = -3*units["B0r03/c"]*p->X.z*p->X.x*tmp2*deltaT;
    gradAdt.z.z = 0;


    //UniformB
//    gradAdt.x.x = 0;
//    gradAdt.x.y = 1.0e-2*units["gauss"]/2/units["c"]*deltaT;
//    gradAdt.x.z = 0;
//    gradAdt.y.x = -1.0e-2*units["gauss"]/2/units["c"]*deltaT;
//    gradAdt.y.y = 0;
//    gradAdt.y.z = 0;
//    gradAdt.z = Vector3D();
}

void EngineForAnalyticalField::updateP(const Range &range) {
    static Tensor3D coefficient;
    static Vector3D RHS;
    Particle *p;
    for_each_Particle_within(grid, p, range){

                        p->X.x = p->Position.x - grid->gridX()/2.0;
                        p->X.y = p->Position.y - grid->gridY()/2.0;
                        p->X.z = p->Position.z - grid->gridZ()/2.0;
                        tmp1 = Square(p->X.x)+Square(p->X.y)+Square(p->X.z);
                        tmp3 = 1/Cube(sqrt(tmp1));
                        tmp2 = tmp3/tmp1;
//                        tmp2 = pow(tmp1, -2.5);
//                        tmp3 = pow(tmp1, -1.5);
                        calcA(p);
                        calcGradAdt(p);
                        buildCoefficientTensor(coefficient);
                        RHS.restore();
                        RHS += p->Momentum;
                        RHS -= (gradAdt*(p->A));
//                        RHS = p->Momentum - (gradAdt*(p->A));
                        rootOfLinearEquationSet(p->Momentum, coefficient, RHS);

                    }end_for_each_Particle(p)
}

void EngineForAnalyticalField::updateX(const Range &range) {
    Particle *p;
    for_each_Particle_within(grid, p, range){
                        p->Position += (p->Momentum - p->A)*deltaT;
                    }end_for_each_Particle(p)
}

void EngineForAnalyticalField::buildCoefficientTensor(Tensor3D &coeT) {
    coeT.x.x = 1 - gradAdt.x.x;
    coeT.x.y = -gradAdt.x.y;
    coeT.x.z = -gradAdt.x.z;
    coeT.y.x = -gradAdt.y.x;
    coeT.y.y = 1 - gradAdt.y.y;
    coeT.y.z = -gradAdt.y.z;
    coeT.z.x = -gradAdt.z.x;
    coeT.z.y = -gradAdt.z.y;
    coeT.z.z = 1 - gradAdt.z.z;
}


void EngineForAnalyticalField::rootOfLinearEquationSet(Vector3D &root, const Tensor3D &coefficient,
                                                       const Vector3D &RHS) {
    static Tensor3D adjacentMatrix;
    double determinant = 0;

    adjacentMatrix.x.x = coefficient.y.y * coefficient.z.z - coefficient.y.z * coefficient.z.y;
    adjacentMatrix.x.y = -coefficient.x.y * coefficient.z.z + coefficient.x.z * coefficient.z.y;
    adjacentMatrix.x.z = coefficient.x.y * coefficient.y.z - coefficient.x.z * coefficient.y.y;
    adjacentMatrix.y.x = -coefficient.y.x * coefficient.z.z + coefficient.y.z * coefficient.z.x;
    adjacentMatrix.y.y = coefficient.x.x * coefficient.z.z - coefficient.x.z * coefficient.z.x;
    adjacentMatrix.y.z = -coefficient.x.x * coefficient.y.z + coefficient.x.z * coefficient.y.x;
    adjacentMatrix.z.x = coefficient.y.x * coefficient.z.y - coefficient.y.y * coefficient.z.x;
    adjacentMatrix.z.y = -coefficient.x.x * coefficient.z.y + coefficient.x.y * coefficient.z.x;
    adjacentMatrix.z.z = coefficient.x.x * coefficient.y.y - coefficient.x.y * coefficient.y.x;

    determinant += coefficient.x.x * adjacentMatrix.x.x;
    determinant += coefficient.y.x * adjacentMatrix.x.y;
    determinant += coefficient.z.x * adjacentMatrix.x.z;

    if (fabs(determinant) < 1e-16){
        std::cerr << "Coefficient Matrix is sigular, unable to find solution." << std::endl;
        std::cerr<<"Determinant: "<<determinant<<std::endl;
        //exit(0);
    }

#if OPTIMIZE
    adjacentMatrix/=determinant;

	root.x = adjacentMatrix.x * RHS;
	root.y = adjacentMatrix.y * RHS;
	root.z = adjacentMatrix.z * RHS;
#else
    Tensor3D inverseMatrix;

    inverseMatrix = adjacentMatrix/determinant;

    root.x = inverseMatrix.x * RHS;
    root.y = inverseMatrix.y * RHS;
    root.z = inverseMatrix.z * RHS;
#endif
}
//
// Created by ylang on 17-12-28.
//



//　This Engine requires OMP and MPI parallels to be turned off

#include "EngineForModifiedField.h"
#include "Grid.h"
#include "CaseOfDipole.h"
#include "LePlasma.h"

Tensor3D EngineForModifiedField::gradAdt = Tensor3D();
EngineForModifiedField::EngineForModifiedField(Grid *grid, double deltaT, map<string, double> &units):\
 Engine(grid, deltaT, units){

}

EngineForModifiedField::~EngineForModifiedField() {

}

void EngineForModifiedField::update(const Range &range){

    updateP(range);
    updateX(range);
}

void EngineForModifiedField::calcA(Particle *p) {
    r0_1 = 1/sqrt(Square(p->X.x)+Square(p->X.y));
    ax2 = -p->X.y-0.5*Square(p->X.y)*r0_1+p->X.x*Square(p->X.z)*Square(r0_1);
    ay2 = p->X.x+0.5*p->X.x*p->X.y*r0_1+p->X.y*Square(p->X.z)*Square(r0_1);
    p->A.x = units["B0r03/c"]*r_3*ax2;
    p->A.y = units["B0r03/c"]*r_3*ay2;
    p->A.z = -1*units["B0r03/c"]*p->X.z*r_3;

}

void EngineForModifiedField::calcGradAdt(Particle *p) {
    gradAdt.x.x = deltaT*units["B0r03/c"]*(-3*p->X.x*r_5*ax2+r_3*(0.5*p->X.x*Square(p->X.y)*Cube(r0_1) \
 +Square(p->X.z)*(Square(p->X.y)-Square(p->X.x))*Square(r0_1)*Square(r0_1)));
    gradAdt.x.y = deltaT*units["B0r03/c"]*(-3*p->X.x*r_5*ay2+r_3*(1+0.5*Cube(p->X.y)*Cube(r0_1) \
 -2*p->X.x*p->X.y*Square(p->X.z)*Square(r0_1)*Square(r0_1)));
    gradAdt.x.z = deltaT*units["B0r03/c"]*3*p->X.z*p->X.x*r_5;
    gradAdt.y.x = deltaT*units["B0r03/c"]*(-3*p->X.y*r_5*ax2+r_3*(-1+0.5*(-2*Square(p->X.x)*p->X.y-Cube(p->X.y))*Cube(r0_1) \
-2*p->X.x*p->X.y*Square(p->X.z)*Square(r0_1)*Square(r0_1)));
    gradAdt.y.y = deltaT*units["B0r03/c"]*(-3*p->X.y*r_5*ay2+r_3*(0.5*Cube(p->X.x)*Cube(r0_1) \
 +Square(p->X.z)*(Square(p->X.x)-Square(p->X.y))*Square(r0_1)*Square(r0_1)));
    gradAdt.y.z = deltaT*units["B0r03/c"]*3*p->X.y*p->X.z*r_5;
    gradAdt.z.x = deltaT*units["B0r03/c"]*(-3*p->X.z*r_5*ax2+r_3*2*p->X.z*p->X.x*Square(r0_1));
    gradAdt.z.y = deltaT*units["B0r03/c"]*(-3*p->X.z*r_5*ay2+r_3*2*p->X.y*p->X.z*Square(r0_1));
    gradAdt.z.z = deltaT*units["B0r03/c"]*(-1*r_3+3*Square(p->X.z)*r_5);

}

void EngineForModifiedField::updateP(const Range &range) {
    static Tensor3D coefficient;
    static Vector3D RHS;
    Particle *p;
    for_each_Particle_within(grid, p, range){

                        calcGradAdt(p);
                        buildCoefficientTensor(coefficient);
                        RHS.restore();
                        RHS += p->Momentum;
                        RHS -= (gradAdt*(p->A));
//                        RHS = p->Momentum - (gradAdt*(p->A));
                        rootOfLinearEquationSet(p->Momentum, coefficient, RHS);

                    }end_for_each_Particle(p)
}

void EngineForModifiedField::updateX(const Range &range) {
    Particle *p;
    for_each_Particle_within(grid, p, range){
                        p->Position += (p->Momentum - p->A)*deltaT;
                        p->X.x = p->Position.x - grid->gridX()/2.0;
                        p->X.y = p->Position.y - grid->gridY()/2.0;
                        p->X.z = p->Position.z - grid->gridZ()/2.0;
                        r2 = Square(p->X.x)+Square(p->X.y)+Square(p->X.z);
                        r_3 = 1/Cube(sqrt(r2));
                        r_5 = r_3/r2;
                        calcA(p);

                    }end_for_each_Particle(p)
}

void EngineForModifiedField::buildCoefficientTensor(Tensor3D &coeT) {
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


void EngineForModifiedField::rootOfLinearEquationSet(Vector3D &root, const Tensor3D &coefficient,
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
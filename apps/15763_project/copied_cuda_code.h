#include <uipc/uipc.h>

Matrix3x3 Ds(const Vector3& x0, const Vector3& x1, const Vector3& x2, const Vector3& x3)
{
    Matrix3x3 Ds;
    Ds.col(0) = x1 - x0;
    Ds.col(1) = x2 - x0;
    Ds.col(2) = x3 - x0;
    return Ds;
}

Matrix3x3 F(const Vector3&   x0,
                         const Vector3&   x1,
                         const Vector3&   x2,
                         const Vector3&   x3,
                         const Matrix3x3& DmInv)
{
    auto ds = Ds(x0, x1, x2, x3);
    return ds * DmInv;
}

template <typename T>
void dEdVecF(Eigen::Matrix<T, 3, 3>& PEPF, const T& mu, const T& lambda, const Eigen::Matrix<T, 3, 3>& F)
{
    auto J  = F.determinant();
    auto Ic = F.squaredNorm();
    Eigen::Matrix<T, 3, 3> pJpF;

    pJpF(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
    pJpF(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
    pJpF(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);

    pJpF(1, 0) = F(2, 1) * F(0, 2) - F(2, 2) * F(0, 1);
    pJpF(1, 1) = F(2, 2) * F(0, 0) - F(2, 0) * F(0, 2);
    pJpF(1, 2) = F(2, 0) * F(0, 1) - F(2, 1) * F(0, 0);

    pJpF(2, 0) = F(0, 1) * F(1, 2) - F(1, 1) * F(0, 2);
    pJpF(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
    pJpF(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);

    PEPF = mu * (1 - 1 / (Ic + 1)) * F + (lambda * (J - 1 - 0.75 * mu / lambda)) * pJpF;
}

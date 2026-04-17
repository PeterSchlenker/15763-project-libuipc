#include <app/asset_dir.h>
#include <uipc/constitution/elastic_moduli.h>
#include <uipc/uipc.h>
#include <15763_project/copied_cuda_code.h>

#include <iostream>

using namespace uipc::constitution;

vector<Matrix3x3> calculate_cauchy_stress(
    const vector<Vector3> &Vs,
    const vector<Vector4i> &Ts,
    const ElasticModuli& moduli,
    const vector<Matrix3x3> &Dm_invs,
    const span<const Vector3> &currentVs
) {
    // code based on stable_neo_hookean_3d.cu
    vector<Matrix3x3> sigmas(Ts.size());

    // TODO: convert to a cuda parallel for
    for (int I = 0; I < Ts.size(); I++) {
        const Vector4i &tet = Ts[I];
        const Matrix3x3 &Dm_inv = Dm_invs[I];
        
        const Vector3 &x0 = currentVs[tet(0)];
        const Vector3 &x1 = currentVs[tet(1)];
        const Vector3 &x2 = currentVs[tet(2)];
        const Vector3 &x3 = currentVs[tet(3)];

        auto Fmat = F(x0, x1, x2, x3, Dm_inv);

        auto J = Fmat.determinant();

        Matrix3x3 dEdF;
        dEdVecF(dEdF, moduli.mu(), moduli.lambda(), Fmat);

        sigmas[I] = 1 / J * dEdF * Fmat.transpose();
    }

    return sigmas;
}

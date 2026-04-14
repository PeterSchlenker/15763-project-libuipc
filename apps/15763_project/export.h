#include <app/asset_dir.h>
#include <uipc/uipc.h>

using namespace uipc;
using namespace uipc::core;
using namespace uipc::geometry;

bool neighbor(const Vector3i &face, const Vector4i &tet) {
    return (face.x() == tet.x() || face.x() == tet.y() || face.x() == tet.z() || face.x() == tet.w())
        && (face.y() == tet.x() || face.y() == tet.y() || face.y() == tet.z() || face.y() == tet.w())
        && (face.z() == tet.x() || face.z() == tet.y() || face.z() == tet.z() || face.z() == tet.w());
}

SizeT faces(const vector<Vector3> &Vs, const vector<Vector4i> &Ts, vector<Vector3i> &Fs, vector<Vector2> &uvs) {
    // step 1: map vertices to all neighboring tetrahedra
    vector<vector<SizeT>> vertex_tet_map(Vs.size());

    for (SizeT i = 0; i < Ts.size(); i++) {
        vertex_tet_map[Ts[i].x()].push_back(i);
        vertex_tet_map[Ts[i].y()].push_back(i);
        vertex_tet_map[Ts[i].z()].push_back(i);
        vertex_tet_map[Ts[i].w()].push_back(i);
    }

    // step 2: for each face of each tetrahedra, find its neighbors and insert a face if needed
    Fs = vector<Vector3i>();
    uvs = vector<Vector2>();

    for (SizeT i = 0; i < Ts.size(); i++) {
        // orient faces so normals point out
        // assume tets have positive volume, i.e. (xy cross xz) dot xw > 0
        bool external_face;

        // xzy
        external_face = true;
        for (auto potential_neighbor : vertex_tet_map[Ts[i].x()]) {
            // only insert faces for neighbors we haven't got to yet
            if (potential_neighbor > i && neighbor({Ts[i].x(), Ts[i].y(), Ts[i].z()}, Ts[potential_neighbor])) {


                external_face = false;
            }
        }

        // xyw

        // xwz

        // yzw

    }
}

void to_ply(const vector<Vector3> &Vs, const vector<Vector4i> &Ts, string filename) {

}

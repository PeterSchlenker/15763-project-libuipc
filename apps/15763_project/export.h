#include <app/asset_dir.h>
#include <uipc/uipc.h>
#include <iostream>
#include <fstream>

using namespace uipc;
using namespace uipc::core;
using namespace uipc::geometry;

bool neighbor(const Vector3i &face, const Vector4i &tet) {
    return (face.x() == tet.x() || face.x() == tet.y() || face.x() == tet.z() || face.x() == tet.w())
        && (face.y() == tet.x() || face.y() == tet.y() || face.y() == tet.z() || face.y() == tet.w())
        && (face.z() == tet.x() || face.z() == tet.y() || face.z() == tet.z() || face.z() == tet.w());
}

void faces(const span<const Vector3> &Vs, const span<const Vector4i> &Ts, vector<Vector3i> &Fs, vector<Vector2i> &sides) {
    // step 1: map vertices to all neighboring tetrahedra
    vector<vector<SizeT>> vertex_tet_map(Vs.size());

    for (SizeT i = 0; i < Ts.size(); i++) {
        vertex_tet_map[Ts[i].x()].push_back(i);
        vertex_tet_map[Ts[i].y()].push_back(i);
        vertex_tet_map[Ts[i].z()].push_back(i);
        vertex_tet_map[Ts[i].w()].push_back(i);
    }

    // step 2: for each face of each tetrahedra, find its neighbors and insert a face if needed
    SizeT num_tets = Ts.size();
    Fs = vector<Vector3i>();
    sides = vector<Vector2i>();

    for (SizeT i = 0; i < Ts.size(); i++) {
        // orient faces so normals point out
        // assume tets have positive volume, i.e. (xy cross xz) dot xw > 0
        vector<Vector3i> tet_faces = {
            {Ts[i].x(), Ts[i].z(), Ts[i].y()},
            {Ts[i].x(), Ts[i].y(), Ts[i].w()},
            {Ts[i].x(), Ts[i].w(), Ts[i].z()},
            {Ts[i].y(), Ts[i].z(), Ts[i].w()}
        };

        for (auto face : tet_faces) {
            bool external_face = true;
            for (auto potential_neighbor : vertex_tet_map[face.x()]) {
                if (potential_neighbor != i && neighbor(face, Ts[potential_neighbor])) {
                    external_face = false;

                    // only insert faces for neighbors we haven't got to yet
                    if (potential_neighbor > i) {
                        Fs.push_back(face);
                        sides.push_back({potential_neighbor, i});
                    }

                    // there can only be one neighbor, no point in looking for more
                    break;
                }
            }
            if (external_face) {
                Fs.push_back(face);
                sides.push_back({num_tets, i});
            }
        }
    }
}

void write_ply(const span<const Vector3> &Vs, const span<const Vector4i> &Ts, string filename) {
    std::ofstream file(filename);

    vector<Vector3i> Fs;
    vector<Vector2i> sides;
    faces(Vs, Ts, Fs, sides);

    // header
    file
    << "ply\nformat ascii 1.0\n"
    << "element vertex " << Fs.size() * 3 << "\n"
    << "property float x\n"
    << "property float y\n"
    << "property float z\n"
    << "property float outer_tet\n"
    << "property float inner_tet\n"
    << "element face " << Fs.size() <<"\n"
    << "property list uchar uint vertex_indices\n"
    << "comment tet " << Ts.size() << "\n"
    << "end_header" << std::endl;

    // data
    // vertices
    // because each face needs to have different uv coordinates, they can't share vertices
    for (SizeT i = 0; i < Fs.size(); i++) {
        for (auto vi : {Fs[i].x(), Fs[i].y(), Fs[i].z()}) {
            Vector3 v = Vs[vi];

            file << v.x() << " " << v.y() << " " << v.z() << " " << sides[i].x() << " " << sides[i].y() << std::endl;
        }
    }

    // faces
    for (SizeT i = 0; i < Fs.size(); i++) {
        // use the index in the file, not the index in Vs
        file << "3 " << 3 * i << " " << 3 * i + 1 << " " << 3 * i + 2 << std::endl;
    }

    file.close();
}

void write_cauchy_stress_csv(const span<const Matrix3x3> &sigmas, string filename) {
    std::ofstream file(filename);

    for (auto sigma : sigmas) {
        // row-major ordering
        file << sigma(0, 0) << "," << sigma(0, 1) << "," << sigma(0, 2)
            << "," << sigma(1, 0) << "," << sigma(1, 1) << "," << sigma(1, 2)
            << "," << sigma(2, 0) << "," << sigma(2, 1) << "," << sigma(2, 2) << std::endl;
    }
}

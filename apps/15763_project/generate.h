#include <app/asset_dir.h>
#include <uipc/uipc.h>
#include <cmath>
#include <numbers>
#include <iostream>

using namespace uipc;
using namespace uipc::core;
using namespace uipc::geometry;

Float pi = std::numbers::pi_v<Float>;

// row 0 is special because it is the center, shared across all sectors
int cylinder_vertex_index(SizeT layer, SizeT row, SizeT sector, SizeT i, SizeT num_disk_vertices) {
    SizeT layer_offset = layer * num_disk_vertices;
    SizeT row_offset = row == 0 ? 0 : 1 + 6 * row * (row - 1) / 2;
    SizeT sector_offset = sector * row;

    return layer_offset + row_offset + sector_offset + i;
}

SimplicialComplex cylinder(Float radius, Float height, SizeT radius_subdivisions, SizeT height_subdivisions, vector<Vector3> &Vs, vector<Vector4i> &Ts) {
    SizeT num_disk_vertices = 1 + 6 * radius_subdivisions * (radius_subdivisions + 1) / 2;
    SizeT num_disk_faces = 6 * radius_subdivisions * radius_subdivisions;

    vector<Vector2> disk_vertices(num_disk_vertices);

    // put the middle vertex at the start of the array
    disk_vertices[0] = Vector2{0, 0};

    // insert all vertices except the middle one
    for (SizeT row = 1; row <= radius_subdivisions; row++) {
        Float r = (Float)row / radius_subdivisions * radius;
        
        for (SizeT sector = 0; sector < 6; sector++) {
            for (SizeT i = 0; i < row; i++) {
                Float theta = (sector + (Float)i / row) * pi / 3;
                Float x = r * std::cos(theta);
                Float y = r * std::sin(theta);

                disk_vertices[cylinder_vertex_index(0, row, sector, i, num_disk_vertices)] = Vector2{x, y};
            }
        }
    }

    Vs = vector<Vector3>(num_disk_vertices * (height_subdivisions + 1));
    Ts = vector<Vector4i>(3 * num_disk_faces * height_subdivisions);

    for (SizeT layer = 0; layer <= height_subdivisions; layer++) {
        Float z = -height / 2 + layer * height / height_subdivisions;

        // populate the vertices
        for (SizeT k = 0; k < num_disk_vertices; k++) {
            Vector2 xy = disk_vertices[k];
            Vs[layer * num_disk_vertices + k] = Vector3{xy.x(), xy.y(), z};
        }

        // only populate the tets after the first layer
        if (layer > 0) {
            for (SizeT row = 0; row <= radius_subdivisions - 1; row++) {
                for (SizeT sector = 0; sector < 6; sector++) {

                    // populate the pairs of triangular prisms between each row
                    for (SizeT i = 0; i < row; i++) {
                        // controls the alternating edges in two dimensions
                        // each prism has 4 configurations, determined by these variables
                        bool lean_forward = (row % 2 == 0) ^ (layer % 2 == 0);
                        bool lean_left_top = (i % 2 == 0) ^ (row % 2 == 0 && sector % 2 != 0) ^ (layer % 2 == 0);
                        bool lean_left_bottom = (i % 2 == 0) ^ (row % 2 != 0 && sector % 2 != 0) ^ (layer % 2 == 0);

                        SizeT offset = 3 * (num_disk_faces * (layer - 1) + 6 * row * row + sector * (2 * row + 1) + 2 * i);

                        // first prism
                        // swap the order of vertices if needed to keep the volume positive
                        int bottom_index_1 = cylinder_vertex_index(layer - 1, row + !lean_forward, sector, i, num_disk_vertices);
                        int bottom_index_2 = cylinder_vertex_index(layer - 1, row + 1, sector, i + (!lean_forward || !lean_left_top), num_disk_vertices);

                        Ts[offset] = Vector4i{
                            cylinder_vertex_index(layer - 1, row, sector, i, num_disk_vertices),
                            cylinder_vertex_index(layer - 1, row + 1, sector, i, num_disk_vertices),
                            cylinder_vertex_index(layer - 1, row + 1, sector, i + 1, num_disk_vertices),
                            cylinder_vertex_index(layer, row + lean_forward, sector, i + (lean_forward && lean_left_top), num_disk_vertices)
                        };
                        Ts[offset + 1] = Vector4i{
                            !lean_forward ? bottom_index_1 : bottom_index_2,
                            !lean_forward ? bottom_index_2 : bottom_index_1,
                            cylinder_vertex_index(layer, row + lean_forward, sector, i, num_disk_vertices),
                            cylinder_vertex_index(layer, row + 1, sector, i + (lean_forward || lean_left_top), num_disk_vertices)
                        };
                        Ts[offset + 2] = Vector4i{
                            cylinder_vertex_index(layer - 1, row + !lean_forward, sector, i + (!lean_forward && !lean_left_top), num_disk_vertices),
                            cylinder_vertex_index(layer, row, sector, i, num_disk_vertices),
                            cylinder_vertex_index(layer, row + 1, sector, i, num_disk_vertices),
                            cylinder_vertex_index(layer, row + 1, sector, i + 1, num_disk_vertices)
                        };

                        // second prism
                        // if we are at the end, the corner vertex needs to wrap around
                        // get the next sector and reset i instead of incrementing i
                        // i.e. in that case we can't use (_, row, sector, i + 1), we have to wrap
                        bool corner_vertex = i == row - 1;
                        SizeT left_sector = corner_vertex ? (sector + 1) % 6 : sector;
                        SizeT left_i_bottom = corner_vertex ? 0 : i + 1;
                        int bottom_index_3 = cylinder_vertex_index(layer - 1, row, lean_forward || !lean_left_bottom ? left_sector : sector, lean_forward || !lean_left_bottom ? left_i_bottom : i, num_disk_vertices);
                        int bottom_index_4 = cylinder_vertex_index(layer - 1, row + !lean_forward, sector, i + !lean_forward, num_disk_vertices);

                        Ts[offset + 3] = Vector4i{
                            cylinder_vertex_index(layer - 1, row, sector, i, num_disk_vertices),
                            cylinder_vertex_index(layer - 1, row + 1, sector, i + 1, num_disk_vertices),
                            cylinder_vertex_index(layer - 1, row, left_sector, left_i_bottom, num_disk_vertices),
                            cylinder_vertex_index(layer, row + lean_forward, !lean_forward && lean_left_bottom ? left_sector : sector, lean_forward ? i + 1 : (lean_left_bottom ? left_i_bottom : i), num_disk_vertices)
                        };
                        Ts[offset + 4] = Vector4i{
                            !lean_forward ? bottom_index_3 : bottom_index_4,
                            !lean_forward ? bottom_index_4 : bottom_index_3,
                            cylinder_vertex_index(layer, row, !lean_forward || lean_left_bottom ? left_sector : sector, !lean_forward || lean_left_bottom ? left_i_bottom : i, num_disk_vertices),
                            cylinder_vertex_index(layer, row + lean_forward, sector, i + lean_forward, num_disk_vertices)
                        };
                        Ts[offset + 5] = Vector4i{
                            cylinder_vertex_index(layer - 1, row + !lean_forward, lean_forward && !lean_left_bottom ? left_sector : sector, !lean_forward ? i + 1 : (!lean_left_bottom ? left_i_bottom : i), num_disk_vertices),
                            cylinder_vertex_index(layer, row, left_sector, left_i_bottom, num_disk_vertices),
                            cylinder_vertex_index(layer, row, sector, i, num_disk_vertices),
                            cylinder_vertex_index(layer, row + 1, sector, i + 1, num_disk_vertices)
                        };
                    }

                    // populate the last prism in the row of the sector, using wrapping
                    bool lean_forward = (row % 2 == 0) ^ (layer % 2 == 0);
                    bool lean_left_top = (row % 2 == 0) ^ (row % 2 == 0 && sector % 2 != 0) ^ (layer % 2 == 0);
                    SizeT next_sector = (sector + 1) % 6;

                    SizeT offset = 3 * (num_disk_faces * (layer - 1) + 6 * row * row + sector * (2 * row + 1) + 2 * row);

                    // swap the order of vertices if needed to keep the volume positive
                    int bottom_index_1 = cylinder_vertex_index(layer - 1, row + !lean_forward, lean_forward ? next_sector : sector, lean_forward ? 0 : row, num_disk_vertices);
                    int bottom_index_2 = cylinder_vertex_index(layer - 1, row + 1, !lean_forward || !lean_left_top ? next_sector : sector, !lean_forward || !lean_left_top ? 0 : row, num_disk_vertices);

                    Ts[offset] = Vector4i{
                        cylinder_vertex_index(layer - 1, row, next_sector, 0, num_disk_vertices),
                        cylinder_vertex_index(layer - 1, row + 1, sector, row, num_disk_vertices),
                        cylinder_vertex_index(layer - 1, row + 1, next_sector, 0, num_disk_vertices),
                        cylinder_vertex_index(layer, row + lean_forward, !lean_forward || lean_left_top ? next_sector : sector, !lean_forward || lean_left_top ? 0 : row, num_disk_vertices)
                    };
                    Ts[offset + 1] = Vector4i{
                        !lean_forward ? bottom_index_1 : bottom_index_2,
                        !lean_forward ? bottom_index_2 : bottom_index_1,
                        cylinder_vertex_index(layer, row + lean_forward, !lean_forward ? next_sector : sector, !lean_forward ? 0 : row, num_disk_vertices),
                        cylinder_vertex_index(layer, row + 1, lean_forward || lean_left_top ? next_sector : sector, lean_forward || lean_left_top ? 0 : row, num_disk_vertices)
                    };
                    Ts[offset + 2] = Vector4i{
                        cylinder_vertex_index(layer - 1, row + !lean_forward, lean_forward || !lean_left_top ? next_sector : sector, lean_forward || !lean_left_top ? 0 : row, num_disk_vertices),
                        cylinder_vertex_index(layer, row, next_sector, 0, num_disk_vertices),
                        cylinder_vertex_index(layer, row + 1, sector, row, num_disk_vertices),
                        cylinder_vertex_index(layer, row + 1, next_sector, 0, num_disk_vertices)
                    };
                }
            }
        }
    }

    return tetmesh(Vs, Ts);
}

int box_vertex_index(SizeT i, SizeT j, SizeT k, SizeT num_x_vertices, SizeT num_y_vertices) {
    return i + j * num_x_vertices + k * num_x_vertices * num_y_vertices;
}

SimplicialComplex box(Vector3 size, Vector3i subdivisions, vector<Vector3> &Vs, vector<Vector4i> &Ts) {
    SizeT num_x_vertices = subdivisions.x() + 1;
    SizeT num_y_vertices = subdivisions.y() + 1;
    SizeT num_z_vertices = subdivisions.z() + 1;
    SizeT num_vertices = num_x_vertices * num_y_vertices * num_z_vertices;
    Vs = vector<Vector3>(num_vertices);
    
    for (SizeT k = 0; k < num_z_vertices; k++) {
        for (SizeT j = 0; j < num_y_vertices; j++) {
            for (SizeT i = 0; i < num_x_vertices; i++) {
                Vs[box_vertex_index(i, j, k, num_x_vertices, num_y_vertices)] = Vector3{
                    ((Float)i / subdivisions.x() - 0.5) * size.x(),
                    ((Float)j / subdivisions.y() - 0.5) * size.y(),
                    ((Float)k / subdivisions.z() - 0.5) * size.z()
                };
            }
        }
    }

    SizeT num_voxels = subdivisions.x() * subdivisions.y() * subdivisions.z();
    Ts = vector<Vector4i>(num_voxels * 5);

    SizeT tet_index = 0;
    for (SizeT k = 0; k < subdivisions.z(); k++) {
        for (SizeT j = 0; j < subdivisions.y(); j++) {
            for (SizeT i = 0; i < subdivisions.x(); i++) {
                if ((i + j + k) % 2 == 0) {
                    // i, j, k is a corner
                    Ts[tet_index++] = Vector4i{
                        box_vertex_index(i, j, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j + 1, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j, k + 1, num_x_vertices, num_y_vertices)
                    };
                    Ts[tet_index++] = Vector4i{
                        box_vertex_index(i, j, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j + 1, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j + 1, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j + 1, k, num_x_vertices, num_y_vertices)
                    };
                    Ts[tet_index++] = Vector4i{
                        box_vertex_index(i, j, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j + 1, k + 1, num_x_vertices, num_y_vertices)
                    };
                    Ts[tet_index++] = Vector4i{
                        box_vertex_index(i, j, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j + 1, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j + 1, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j, k + 1, num_x_vertices, num_y_vertices)
                    };
                    Ts[tet_index++] = Vector4i{
                        box_vertex_index(i + 1, j + 1, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j + 1, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j + 1, k + 1, num_x_vertices, num_y_vertices)
                    };
                }
                else {
                    // i, j, k is not a corner
                    Ts[tet_index++] = Vector4i{
                        box_vertex_index(i, j, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j + 1, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j, k + 1, num_x_vertices, num_y_vertices)
                    };
                    Ts[tet_index++] = Vector4i{
                        box_vertex_index(i + 1, j, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j + 1, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j + 1, k + 1, num_x_vertices, num_y_vertices)
                    };
                    Ts[tet_index++] = Vector4i{
                        box_vertex_index(i + 1, j, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j + 1, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j + 1, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j + 1, k + 1, num_x_vertices, num_y_vertices)
                    };
                    Ts[tet_index++] = Vector4i{
                        box_vertex_index(i, j + 1, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j + 1, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i, j, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j + 1, k + 1, num_x_vertices, num_y_vertices)
                    };
                    Ts[tet_index++] = Vector4i{
                        box_vertex_index(i, j, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j, k + 1, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j, k, num_x_vertices, num_y_vertices),
                        box_vertex_index(i + 1, j + 1, k + 1, num_x_vertices, num_y_vertices)
                    };
                }
            }
        }
    }

    std::cout << "x" << std::endl;

    return tetmesh(Vs, Ts);
}

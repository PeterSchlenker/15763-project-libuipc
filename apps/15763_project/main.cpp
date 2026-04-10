#include <app/asset_dir.h>
#include <uipc/uipc.h>
#include <uipc/constitution/affine_body_constitution.h>
#include <cmath>
#include <numbers>
#include <iostream>

using namespace uipc;
using namespace uipc::core;
using namespace uipc::geometry;
using namespace uipc::constitution;
namespace fs = std::filesystem;

Float pi = std::numbers::pi_v<Float>;

// row 0 is special because it is the center, shared across all sectors
int cylinder_vertex_index(SizeT layer, SizeT row, SizeT sector, SizeT i, SizeT num_disk_vertices) {
    SizeT layer_offset = layer * num_disk_vertices;
    SizeT row_offset = row == 0 ? 0 : 1 + 6 * row * (row - 1) / 2;
    SizeT sector_offset = sector * row;

    return layer_offset + row_offset + sector_offset + i;
}

SimplicialComplex cylinder(Float radius, Float height, SizeT radius_subdivisions, SizeT height_subdivisions) {
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

    vector<Vector3> Vs(num_disk_vertices * (height_subdivisions + 1));
    vector<Vector4i> Ts(3 * num_disk_faces * height_subdivisions);

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

                        Ts[offset] = Vector4i{0, 0, 0, 0};
                        Ts[offset + 1] = Vector4i{0, 0, 0, 0};
                        Ts[offset + 2] = Vector4i{0, 0, 0, 0};

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

                        Ts[offset + 3] = Vector4i{0, 0, 0, 0};
                        Ts[offset + 4] = Vector4i{0, 0, 0, 0};
                        Ts[offset + 5] = Vector4i{0, 0, 0, 0};
                    }

                    // populate the last prism in the row of the sector, using wrapping
                    bool lean_forward = (row % 2 == 0) ^ (layer % 2 == 0);
                    bool lean_left_top = (row % 2 == 0) ^ (row % 2 == 0 && sector % 2 != 0) ^ (layer % 2 == 0);
                    SizeT next_sector = (sector + 1) % 6;

                    SizeT offset = 3 * (num_disk_faces * (layer - 1) + 6 * row * row + sector * (2 * row + 1) + 2 * row);

                    // swap the order of vertices if needed to keep the volume positive
                    int bottom_index_1 = cylinder_vertex_index(layer - 1, row + !lean_forward, !lean_forward ? next_sector : sector, !lean_forward ? 0 : row, num_disk_vertices);
                    int bottom_index_2 = cylinder_vertex_index(layer - 1, row + 1, lean_forward && !lean_left_top ? next_sector : sector, lean_forward && !lean_left_top ? 0 : row, num_disk_vertices);

                    Ts[offset] = Vector4i{
                        cylinder_vertex_index(layer - 1, row, next_sector, 0, num_disk_vertices),
                        cylinder_vertex_index(layer - 1, row + 1, sector, row, num_disk_vertices),
                        cylinder_vertex_index(layer - 1, row + 1, next_sector, 0, num_disk_vertices),
                        cylinder_vertex_index(layer, row + lean_forward, !lean_forward || lean_left_top ? next_sector : sector, !lean_forward || lean_left_top ? 0 : row, num_disk_vertices)
                    };
                    Ts[offset + 1] = Vector4i{
                        lean_forward ? bottom_index_1 : bottom_index_2,
                        lean_forward ? bottom_index_2 : bottom_index_1,
                        cylinder_vertex_index(layer, row + lean_forward, lean_forward ? next_sector : sector, lean_forward ? 0 : row, num_disk_vertices),
                        cylinder_vertex_index(layer, row + 1, !lean_forward && lean_left_top ? next_sector : sector, !lean_forward && lean_left_top ? 0 : row, num_disk_vertices)
                    };
                    Ts[offset + 2] = Vector4i{
                        cylinder_vertex_index(layer - 1, row + !lean_forward, lean_forward || lean_left_top ? next_sector : sector, lean_forward || lean_left_top ? 0 : row, num_disk_vertices),
                        cylinder_vertex_index(layer, row, next_sector, 0, num_disk_vertices),
                        cylinder_vertex_index(layer, row + 1, sector, row, num_disk_vertices),
                        cylinder_vertex_index(layer, row + 1, next_sector, 0, num_disk_vertices)
                    };

                    Ts[offset] = Vector4i{0, 0, 0, 0};
                    //Ts[offset + 1] = Vector4i{0, 0, 0, 0};
                    Ts[offset + 2] = Vector4i{0, 0, 0, 0};
                }
            }
        }
    }

    return tetmesh(Vs, Ts);
}

int main() {
    Engine engine{"cuda"};

    World world{engine};
    auto  config = Scene::default_config();
    config["gravity"] = Vector3{0, -9.8, 0};
    config["dt"] = 0.01_s;

    Scene scene{config};
    {
        // create constitution and contact model
        AffineBodyConstitution abd;
        scene.constitution_tabular().insert(abd);

        // friction ratio and contact resistance
        scene.contact_tabular().default_model(0.5, 1.0_GPa);
        auto default_element = scene.contact_tabular().default_element();

        // create a regular tetrahedron
        vector<Vector3>  Vs = {Vector3{0, 1, 0},
                               Vector3{0, 0, 1},
                               Vector3{-std::sqrt(3) / 2, 0, -0.5},
                               Vector3{std::sqrt(3) / 2, 0, -0.5}};
        vector<Vector4i> Ts = {Vector4i{0, 1, 2, 3}};

        // setup a base mesh to reduce the later work
        SimplicialComplex base_mesh = cylinder(0.5, 0.2, 5, 2);
        // apply the constitution and contact model to the base mesh
        abd.apply_to(base_mesh, 100.0_MPa);
        // apply the default contact model to the base mesh
        default_element.apply_to(base_mesh);

        // label the surface, enable the contact
        label_surface(base_mesh);
        // label the triangle orientation to export the correct surface mesh
        label_triangle_orient(base_mesh);

        SimplicialComplex mesh1 = base_mesh;
        {
            // move the mesh1 up for 1 unit
            auto pos_view = view(mesh1.positions());
            std::ranges::transform(pos_view,
                                   pos_view.begin(),
                                   [](const Vector3& v) -> Vector3
                                   { return v + Vector3::UnitY() * 1.3; });
        }

        SimplicialComplex mesh2 = base_mesh;
        {
            // find the is_fixed attribute
            auto is_fixed = mesh2.instances().find<IndexT>(builtin::is_fixed);
            // set the first instance to be fixed
            auto is_fixed_view = view(*is_fixed);
            is_fixed_view[0] = 1;
        }


        // create object with two geometries
        auto object = scene.objects().create("tets");
        {
            object->geometries().create(mesh1);
            object->geometries().create(mesh2);
        }
    }

    world.init(scene);

    SceneIO sio{scene};

    auto this_output_path = AssetDir::output_path(UIPC_RELATIVE_SOURCE_FILE);

    sio.write_surface(fmt::format("{}scene_surface{}.obj", this_output_path, 0));

    /*
    for(int i = 1; i < 50; i++)
    {
        world.advance();
        world.sync();
        world.retrieve();
        sio.write_surface(fmt::format("{}scene_surface{}.obj", this_output_path, i));
    }*/
}

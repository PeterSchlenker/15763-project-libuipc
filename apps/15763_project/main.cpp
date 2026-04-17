#include <app/asset_dir.h>
#include <uipc/uipc.h>
#include <uipc/constitution/affine_body_constitution.h>
#include <uipc/constitution/stable_neo_hookean.h>
#include <15763_project/generate.h>
#include <15763_project/export.h>
#include <15763_project/stress.h>
#include <cmath>
#include <iostream>

using namespace uipc;
using namespace uipc::core;
using namespace uipc::geometry;
using namespace uipc::constitution;
namespace fs = std::filesystem;

int main() {
    Engine engine{"cuda"};

    World world{engine};
    auto  config = Scene::default_config();
    config["gravity"] = Vector3{0, -9.8, 0};
    config["dt"] = 0.01_s;

    vector<Vector3> Vs;
    vector<Vector4i> Ts;
    ElasticModuli moduli = ElasticModuli::youngs_poisson(120.0_kPa, 0.49);
    vector<Matrix3x3> Dm_invs;

    Scene scene{config};

    // create constitution and contact model
    AffineBodyConstitution abd;
    scene.constitution_tabular().insert(abd);
    StableNeoHookean snh;
    scene.constitution_tabular().insert(snh);

    // friction ratio and contact resistance
    scene.contact_tabular().default_model(0.5, 1.0_GPa);
    auto default_element = scene.contact_tabular().default_element();

    // setup a base mesh to reduce the later work
    SimplicialComplex base_mesh = cylinder(0.5, 1, 10, 20, Vs, Ts);
    // apply the default contact model to the base mesh
    default_element.apply_to(base_mesh);

    // precompute info for stress calculation
    Dm_invs.reserve(Ts.size());
    for (auto tet : Ts) {
        Dm_invs.push_back(Dm_inv(Vs[tet(0)], Vs[tet(1)], Vs[tet(2)], Vs[tet(3)]));
    }

    // label the surface, enable the contact
    label_surface(base_mesh);
    // label the triangle orientation to export the correct surface mesh
    label_triangle_orient(base_mesh);

    SimplicialComplex rendered_mesh = base_mesh;
    {
        // apply the constitution and contact model to the base mesh
        snh.apply_to(rendered_mesh, moduli);
        // move the mesh1 up for 1 unit
        auto pos_view = view(rendered_mesh.positions());
        std::ranges::transform(pos_view,
                                pos_view.begin(),
                                [](const Vector3& v) -> Vector3
                                { return v + Vector3::UnitY() * 1.3; });
    }

    SimplicialComplex mesh2 = base_mesh;
    {
        // apply the constitution and contact model to the base mesh
        abd.apply_to(mesh2, 100.0_MPa);
        // find the is_fixed attribute
        auto is_fixed = mesh2.instances().find<IndexT>(builtin::is_fixed);
        // set the first instance to be fixed
        auto is_fixed_view = view(*is_fixed);
        is_fixed_view[0] = 1;
    }

    // create object with two geometries
    auto object = scene.objects().create("tets");
    {
        object->geometries().create(rendered_mesh);
        object->geometries().create(mesh2);
    }

    world.init(scene);

    SceneIO sio{scene};

    auto this_output_path = AssetDir::output_path(UIPC_RELATIVE_SOURCE_FILE);

    write_ply(Vs, Ts, fmt::format("{}initial.ply", this_output_path));

    sio.write_surface(fmt::format("{}scene_surface{}.obj", this_output_path, 0));

    for(int i = 1; i < 200; i++)
    {
        world.advance();
        world.sync();
        world.retrieve();

        sio.write_surface(fmt::format("{}scene_surface{}.obj", this_output_path, i));
        
        //scene.geometries().find(0).vertices().find<Vector3>(builtin::position)->view();
        auto mesh = scene.geometries().find(0).geometry->geometry().as<SimplicialComplex>();
        auto currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();

        auto stress = calculate_cauchy_stress(Vs, Ts, moduli, Dm_invs, currentVs);
        write_cauchy_stress_csv(stress, fmt::format("{}stress{}.csv", this_output_path, i));
        write_ply(currentVs, Ts, fmt::format("{}tet_face_mesh{}.ply", this_output_path, i));
    }
}

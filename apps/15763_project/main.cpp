#include <app/asset_dir.h>
#include <uipc/uipc.h>
#include <uipc/constitution/affine_body_constitution.h>
#include <uipc/constitution/affine_body_fixed_joint.h>
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

void cylinder_scene() {
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
    SimplicialComplex base_mesh = cylinder(0.5, 1, 20, 40, Vs, Ts);
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

    auto mesh = scene.geometries().find(0).geometry->geometry().as<SimplicialComplex>();
    auto currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
    auto stress = calculate_cauchy_stress(Vs, Ts, moduli, Dm_invs, currentVs);
    write_ply_face_stress(currentVs, Ts, stress, fmt::format("{}high_res_cylinder_scene/cylinder_face{}.ply", this_output_path, 0));
    write_ply_vertex_stress(currentVs, Ts, stress, fmt::format("{}high_res_cylinder_scene/cylinder_vertex{}.ply", this_output_path, 0));

    sio.write_surface(fmt::format("{}high_res_cylinder_scene/scene_surface{}.obj", this_output_path, 0));

    for(int i = 1; i < 200; i++)
    {
        world.advance();
        world.sync();
        world.retrieve();

        sio.write_surface(fmt::format("{}high_res_cylinder_scene/scene_surface{}.obj", this_output_path, i));
        
        //scene.geometries().find(0).vertices().find<Vector3>(builtin::position)->view();
        auto mesh = scene.geometries().find(0).geometry->geometry().as<SimplicialComplex>();
        auto currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();

        auto stress = calculate_cauchy_stress(Vs, Ts, moduli, Dm_invs, currentVs);
        //write_cauchy_stress_csv(stress, fmt::format("{}stress{}.csv", this_output_path, i));
        write_ply_face_stress(currentVs, Ts, stress, fmt::format("{}high_res_cylinder_scene/cylinder_face{}.ply", this_output_path, i));
        write_ply_vertex_stress(currentVs, Ts, stress, fmt::format("{}high_res_cylinder_scene/cylinder_vertex{}.ply", this_output_path, i));
    }
}

void box_scene() {
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
    SimplicialComplex base_mesh = box(Vector3{0.8, 0.9, 1}, Vector3i{8, 9, 10}, Vs, Ts);
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

    auto mesh = scene.geometries().find(0).geometry->geometry().as<SimplicialComplex>();
    auto currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
    auto stress = calculate_cauchy_stress(Vs, Ts, moduli, Dm_invs, currentVs);
    write_ply_face_stress(currentVs, Ts, stress, fmt::format("{}box_scene/tet_face_mesh{}.ply", this_output_path, 0));

    sio.write_surface(fmt::format("{}box_scene/scene_surface{}.obj", this_output_path, 0));

    for(int i = 1; i < 200; i++)
    {
        world.advance();
        world.sync();
        world.retrieve();

        sio.write_surface(fmt::format("{}box_scene/scene_surface{}.obj", this_output_path, i));
        
        //scene.geometries().find(0).vertices().find<Vector3>(builtin::position)->view();
        auto mesh = scene.geometries().find(0).geometry->geometry().as<SimplicialComplex>();
        auto currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();

        auto stress = calculate_cauchy_stress(Vs, Ts, moduli, Dm_invs, currentVs);
        //write_cauchy_stress_csv(stress, fmt::format("{}stress{}.csv", this_output_path, i));
        write_ply_face_stress(currentVs, Ts, stress, fmt::format("{}box_scene/tet_face_mesh{}.ply", this_output_path, i));
    }
}

void astm_loading_scene() {
    // material parameters
    ElasticModuli moduli = ElasticModuli::youngs_poisson(82.0_GPa, 0.206);
    Float glass_density = 2510; // kg/m^3
    Float C = 2.76e-12; // m^2/N
    Float lambda = 550e-9_m;
    Float material_fringe_value = lambda / C;

    // scene dimensions
    Float scale = 4;
    Float Li = 25.0_mm * scale;
    Float L = 115.0_mm * scale;
    Float Lg = 125.0_mm * scale;
    Float b = 25.0_mm * scale;
    Float d = 10.0_mm * scale;
    Float r = 5.0_mm * scale;
    Float g = 9.8; // m/s^2
    Float initial_gap = 0.8_mm;

    // mass needed to get a certain number of fringe cycles
    Float target_num_cycles = 12.0;
    Float target_max_sigma = target_num_cycles / 2 * material_fringe_value / b;
    Float W = target_max_sigma * 2 * b * d * d / 3 / (L - Li);
    Float loading_mass = W / g;
    Float loading_mass_volume = 2 * M_PI * r * r * b;
    Float loading_mass_density = loading_mass / loading_mass_volume;

    std::cout << "target max sigma: " << target_max_sigma << std::endl;

    // scene config
    auto  config = Scene::default_config();
    config["gravity"] = Vector3{0, -g, 0};
    config["dt"] = 0.01_s;
    //config["newton"]["velocity_tol"] = 0.05_mm / 1.0_s;
    //config["linear_system"]["tol_rate"] = 1e-6;
    config["contact"]["d_hat"] = 0.5_mm;

    // scene setup
    Engine engine{"cuda"};
    World world{engine};
    Scene scene{config};

    // create constitution and contact model
    AffineBodyConstitution abd;
    scene.constitution_tabular().insert(abd);
    StableNeoHookean snh;
    scene.constitution_tabular().insert(snh);

    // friction ratio and contact resistance
    scene.contact_tabular().default_model(0.5, 1.0_GPa);
    auto default_element = scene.contact_tabular().default_element();


    // create the meshes
    // glass mesh
    vector<Vector3> glass_Vs;
    vector<Vector4i> glass_Ts;
    SimplicialComplex glass_mesh = box(Vector3{Lg, d, b}, Vector3i{75, 20, 15}, glass_Vs, glass_Ts);
    //SimplicialComplex glass_mesh = box(Vector3{Lg, d, b}, Vector3i{50, 8, 10}, glass_Vs, glass_Ts);
    //SimplicialComplex glass_mesh = box(Vector3{Lg * 3, d * 3, b * 3}, Vector3i{25 * 3, 4 * 3, 5 * 3}, glass_Vs, glass_Ts);
    default_element.apply_to(glass_mesh);
    label_surface(glass_mesh);
    label_triangle_orient(glass_mesh);

    // loading point mesh
    vector<Vector3> cylinder_Vs;
    vector<Vector4i> cylinder_Ts;
    SimplicialComplex loading_point_mesh = cylinder(r, b, 10, 1, cylinder_Vs, cylinder_Ts);
    default_element.apply_to(loading_point_mesh);
    label_surface(loading_point_mesh);
    label_triangle_orient(loading_point_mesh);

    // precompute info for stress calculation
    vector<Matrix3x3> Dm_invs;
    Dm_invs.reserve(glass_Ts.size());
    for (auto tet : glass_Ts) {
        Dm_invs.push_back(Dm_inv(glass_Vs[tet(0)], glass_Vs[tet(1)], glass_Vs[tet(2)], glass_Vs[tet(3)]));
    }

    // place the meshes in the scene
    // glass
    snh.apply_to(glass_mesh, moduli, glass_density);
    auto glass_object = scene.objects().create("glass");
    glass_object->geometries().create(glass_mesh);

    // top left loading point
    SimplicialComplex top_left_loading_point_mesh = loading_point_mesh;
    abd.apply_to(top_left_loading_point_mesh, 100.0_MPa, loading_mass_density);
    {
        Transform t = Transform::Identity();
        t.translate(Vector3{Li / 2, d / 2 + r + initial_gap, 0});
        view(top_left_loading_point_mesh.transforms())[0] = t.matrix();
    }
    auto top_left_loading_point_object = scene.objects().create("top_left_loading_point");
    auto [top_left_loading_point_geo_slot, _tl] = top_left_loading_point_object->geometries().create(top_left_loading_point_mesh);

    // top right loading point
    SimplicialComplex top_right_loading_point_mesh = loading_point_mesh;
    abd.apply_to(top_right_loading_point_mesh, 100.0_MPa, loading_mass_density);
    {
        Transform t = Transform::Identity();
        t.translate(Vector3{-Li / 2, d / 2 + r + initial_gap, 0});
        view(top_right_loading_point_mesh.transforms())[0] = t.matrix();
    }
    auto top_right_loading_point_object = scene.objects().create("top_right_loading_point");
    auto [top_right_loading_point_geo_slot, _tr] = top_right_loading_point_object->geometries().create(top_right_loading_point_mesh);

    // fix the top two loading points together so they stay a fixed width apart and don't roll
    AffineBodyFixedJoint             fixed_joint;
    vector<S<SimplicialComplexSlot>> tl_geo_slots    = {top_left_loading_point_geo_slot};
    vector<IndexT>                   tl_instance_ids = {0};
    vector<S<SimplicialComplexSlot>> tr_geo_slots    = {top_right_loading_point_geo_slot};
    vector<IndexT>                   tr_instance_ids = {0};
    vector<Float>                    strength_ratios = {100.0};
    auto joint_mesh = fixed_joint.create_geometry(
        span{tl_geo_slots}, span{tl_instance_ids},
        span{tr_geo_slots}, span{tr_instance_ids},
        span{strength_ratios});
    auto joint_object = scene.objects().create("loading_point_fixed_joint");
    joint_object->geometries().create(joint_mesh);

    // bottom left loading point
    SimplicialComplex bottom_left_loading_point_mesh = loading_point_mesh;
    abd.apply_to(bottom_left_loading_point_mesh, 100.0_MPa);
    {
        Transform t = Transform::Identity();
        t.translate(Vector3{L / 2, -d / 2 - r - initial_gap, 0});
        view(bottom_left_loading_point_mesh.transforms())[0] = t.matrix();

        auto is_fixed = bottom_left_loading_point_mesh.instances().find<IndexT>(builtin::is_fixed);
        view(*is_fixed)[0] = 1;
    }
    auto bottom_left_loading_point_object = scene.objects().create("bottom_left_loading_point");
    bottom_left_loading_point_object->geometries().create(bottom_left_loading_point_mesh);

    // bottom right loading point
    SimplicialComplex bottom_right_loading_point_mesh = loading_point_mesh;
    abd.apply_to(bottom_right_loading_point_mesh, 100.0_MPa);
    {
        Transform t = Transform::Identity();
        t.translate(Vector3{-L / 2, -d / 2 - r - initial_gap, 0});
        view(bottom_right_loading_point_mesh.transforms())[0] = t.matrix();

        auto is_fixed = bottom_right_loading_point_mesh.instances().find<IndexT>(builtin::is_fixed);
        view(*is_fixed)[0] = 1;
    }
    auto bottom_right_loading_point_object = scene.objects().create("bottom_right_loading_point");
    bottom_right_loading_point_object->geometries().create(bottom_right_loading_point_mesh);


    // initialize the scene
    world.init(scene);
    SceneIO sio{scene};
    auto this_output_path = AssetDir::output_path(UIPC_RELATIVE_SOURCE_FILE);

    auto mesh = scene.geometries().find(glass_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
    auto currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
    auto stress = calculate_cauchy_stress(glass_Vs, glass_Ts, moduli, Dm_invs, currentVs);
    write_ply_vertex_stress(glass_Vs, glass_Ts, stress, fmt::format("{}astm_loading_scene/tet_face_mesh{}.ply", this_output_path, 0));

    sio.write_surface(fmt::format("{}astm_loading_scene/scene_surface{}.obj", this_output_path, 0));

    for(int i = 1; i < 50; i++)
    {
        world.advance();
        world.sync();
        world.retrieve();

        sio.write_surface(fmt::format("{}astm_loading_scene/scene_surface{}.obj", this_output_path, i));
        
        //scene.geometries().find(0).vertices().find<Vector3>(builtin::position)->view();
        auto mesh = scene.geometries().find(0).geometry->geometry().as<SimplicialComplex>();
        auto currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();

        auto stress = calculate_cauchy_stress(glass_Vs, glass_Ts, moduli, Dm_invs, currentVs);
        //write_cauchy_stress_csv(stress, fmt::format("{}stress{}.csv", this_output_path, i));
        write_ply_face_stress(glass_Vs, glass_Ts, stress, fmt::format("{}astm_loading_scene/tet_face_mesh_face{}.ply", this_output_path, i));
        write_ply_vertex_stress(glass_Vs, glass_Ts, stress, fmt::format("{}astm_loading_scene/tet_face_mesh_vertex{}.ply", this_output_path, i));
    }
}

void table_scene() {
    // material parameters
    ElasticModuli tabletop_moduli = ElasticModuli::youngs_poisson(100.0_MPa, 0.49);
    ElasticModuli leg_moduli = ElasticModuli::youngs_poisson(150.0_MPa, 0.49);
    Float table_density = 1000; // kg/m^3
    Float cube_density = 20000; // kg/m^3

    // scene dimensions
    Float leg_height = 0.65_m;
    Float leg_radius = 0.05_m;
    Float offset_from_table_edge = 0.1_m;
    Float table_length = 1.2_m;
    Float table_width = 0.9_m;
    Float table_height = 0.04_m;
    Float cube_width = 0.1_m;
    Float cube_initial_height = 1.0_m;
    Float g = 9.8; // m/s^2
    Float initial_gap = 1.0_mm;

    Float leg_x = table_width / 2 - offset_from_table_edge;
    Float leg_y = leg_height / 2 + initial_gap;
    Float leg_z = table_length / 2 - offset_from_table_edge;

    // scene config
    auto  config = Scene::default_config();
    config["gravity"] = Vector3{0, -g, 0};
    config["dt"] = 0.01_s;
    config["contact"]["d_hat"] = 0.5_mm;

    // scene setup
    string tetmesh_dir{AssetDir::tetmesh_path()};
    Engine engine{"cuda"};
    World world{engine};
    Scene scene{config};
    bool export_everything = false;

    // create constitution and contact model
    AffineBodyConstitution abd;
    scene.constitution_tabular().insert(abd);
    StableNeoHookean snh;
    scene.constitution_tabular().insert(snh);

    // friction ratio and contact resistance
    scene.contact_tabular().default_model(0.85, 1.0_GPa);
    auto default_element = scene.contact_tabular().default_element();


    // create the meshes
    // table leg
    vector<Vector3> leg_Vs;
    vector<Vector4i> leg_Ts;
    SimplicialComplex leg_mesh = cylinder(leg_radius, leg_height, 5, 65, leg_Vs, leg_Ts);
    default_element.apply_to(leg_mesh);
    label_surface(leg_mesh);
    label_triangle_orient(leg_mesh);

    // tabletop mesh
    vector<Vector3> tabletop_Vs;
    vector<Vector4i> tabletop_Ts;
    SimplicialComplex tabletop_mesh = box(Vector3{table_width, table_height, table_length}, Vector3i{80, 4, 160}, tabletop_Vs, tabletop_Ts);
    default_element.apply_to(tabletop_mesh);
    label_surface(tabletop_mesh);
    label_triangle_orient(tabletop_mesh);

    // cube mesh
    vector<Vector3> cube_Vs;
    vector<Vector4i> cube_Ts;
    auto cube_mesh = box(Vector3{cube_width, cube_width, cube_width}, Vector3i{1, 1, 1}, cube_Vs, cube_Ts);
    label_surface(cube_mesh);
    label_triangle_orient(cube_mesh);

    // precompute info for stress calculation
    vector<Matrix3x3> leg_Dm_invs;
    leg_Dm_invs.reserve(leg_Ts.size());
    for (auto tet : leg_Ts) {
        leg_Dm_invs.push_back(Dm_inv(leg_Vs[tet(0)], leg_Vs[tet(1)], leg_Vs[tet(2)], leg_Vs[tet(3)]));
    }

    vector<Matrix3x3> tabletop_Dm_invs;
    tabletop_Dm_invs.reserve(tabletop_Ts.size());
    for (auto tet : tabletop_Ts) {
        tabletop_Dm_invs.push_back(Dm_inv(tabletop_Vs[tet(0)], tabletop_Vs[tet(1)], tabletop_Vs[tet(2)], tabletop_Vs[tet(3)]));
    }

    // place the meshes in the scene
    // back left table leg
    SimplicialComplex bl_leg_mesh = leg_mesh;
    snh.apply_to(bl_leg_mesh, leg_moduli, table_density);
    {
        Transform t = Transform::Identity();
        t.translate(Vector3{-leg_x, leg_y, -leg_z});
        t.rotate(AngleAxis(M_PI / 2, Vector3{1, 0, 0}));
        view(bl_leg_mesh.transforms())[0] = t.matrix();
    }
    auto bl_leg_object = scene.objects().create("bl_leg");
    bl_leg_object->geometries().create(bl_leg_mesh);

    // back right table leg
    SimplicialComplex br_leg_mesh = leg_mesh;
    snh.apply_to(br_leg_mesh, leg_moduli, table_density);
    {
        Transform t = Transform::Identity();
        t.translate(Vector3{leg_x, leg_y, -leg_z});
        t.rotate(AngleAxis(M_PI / 2, Vector3{1, 0, 0}));
        view(br_leg_mesh.transforms())[0] = t.matrix();
    }
    auto br_leg_object = scene.objects().create("br_leg");
    br_leg_object->geometries().create(br_leg_mesh);

    // front left table leg
    SimplicialComplex fl_leg_mesh = leg_mesh;
    snh.apply_to(fl_leg_mesh, leg_moduli, table_density);
    {
        Transform t = Transform::Identity();
        t.translate(Vector3{-leg_x, leg_y, leg_z});
        t.rotate(AngleAxis(M_PI / 2, Vector3{1, 0, 0}));
        view(fl_leg_mesh.transforms())[0] = t.matrix();
    }
    auto fl_leg_object = scene.objects().create("fl_leg");
    fl_leg_object->geometries().create(fl_leg_mesh);

    // front right table leg
    SimplicialComplex fr_leg_mesh = leg_mesh;
    snh.apply_to(fr_leg_mesh, leg_moduli, table_density);
    {
        Transform t = Transform::Identity();
        t.translate(Vector3{leg_x, leg_y, leg_z});
        t.rotate(AngleAxis(M_PI / 2, Vector3{1, 0, 0}));
        view(fr_leg_mesh.transforms())[0] = t.matrix();
    }
    auto fr_leg_object = scene.objects().create("fr_leg");
    fr_leg_object->geometries().create(fr_leg_mesh);

    // tabletop
    snh.apply_to(tabletop_mesh, tabletop_moduli, table_density);
    {
        Transform t = Transform::Identity();
        t.translate(Vector3{0, leg_height + table_height / 2 + 2 * initial_gap, 0});
        view(tabletop_mesh.transforms())[0] = t.matrix();
    }
    auto tabletop_object = scene.objects().create("tabletop");
    tabletop_object->geometries().create(tabletop_mesh);

    // cube
    abd.apply_to(cube_mesh, 100.0_MPa, cube_density);
    {
        Transform t = Transform::Identity();
        t.translate(Vector3{0, cube_initial_height, 0});
        view(cube_mesh.transforms())[0] = t.matrix();
    }
    auto cube_object = scene.objects().create("cube");
    cube_object->geometries().create(cube_mesh);

    // ground
    auto ground_object = scene.objects().create("ground");
    ground_object->geometries().create(ground(0));


    // initialize the scene
    world.init(scene);
    SceneIO sio{scene};
    auto this_output_path = AssetDir::output_path(UIPC_RELATIVE_SOURCE_FILE);

    auto mesh = scene.geometries().find(bl_leg_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
    auto currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
    auto stress = calculate_cauchy_stress(leg_Vs, leg_Ts, leg_moduli, leg_Dm_invs, currentVs);
    write_ply_face_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/bl_leg_face{}.ply", this_output_path, 0));
    write_ply_vertex_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/bl_leg_vertex{}.ply", this_output_path, 0));

    mesh = scene.geometries().find(br_leg_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
    currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
    stress = calculate_cauchy_stress(leg_Vs, leg_Ts, leg_moduli, leg_Dm_invs, currentVs);
    write_ply_face_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/br_leg_face{}.ply", this_output_path, 0));
    write_ply_vertex_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/br_leg_vertex{}.ply", this_output_path, 0));

    mesh = scene.geometries().find(fl_leg_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
    currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
    stress = calculate_cauchy_stress(leg_Vs, leg_Ts, leg_moduli, leg_Dm_invs, currentVs);
    write_ply_face_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/fl_leg_face{}.ply", this_output_path, 0));
    write_ply_vertex_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/fl_leg_vertex{}.ply", this_output_path, 0));

    mesh = scene.geometries().find(fr_leg_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
    currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
    stress = calculate_cauchy_stress(leg_Vs, leg_Ts, leg_moduli, leg_Dm_invs, currentVs);
    write_ply_face_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/fr_leg_face{}.ply", this_output_path, 0));
    write_ply_vertex_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/fr_leg_vertex{}.ply", this_output_path, 0));

    mesh = scene.geometries().find(tabletop_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
    currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
    stress = calculate_cauchy_stress(tabletop_Vs, tabletop_Ts, tabletop_moduli, tabletop_Dm_invs, currentVs);
    write_ply_face_stress(currentVs, tabletop_Ts, stress, fmt::format("{}table_scene/tabletop_face{}.ply", this_output_path, 0));
    write_ply_vertex_stress(currentVs, tabletop_Ts, stress, fmt::format("{}table_scene/tabletop_vertex{}.ply", this_output_path, 0));

    mesh = scene.geometries().find(cube_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
    currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
    stress = calculate_cauchy_stress(cube_Vs, cube_Ts, tabletop_moduli, tabletop_Dm_invs, currentVs);
    write_ply_face_stress(currentVs, cube_Ts, stress, fmt::format("{}table_scene/cube{}.ply", this_output_path, 0));

    sio.write_surface(fmt::format("{}table_scene/scene_surface{}.obj", this_output_path, 0));

    for(int i = 1; i < 50; i++)
    {
        world.advance();
        world.sync();
        world.retrieve();

        sio.write_surface(fmt::format("{}table_scene/scene_surface{}.obj", this_output_path, i));

        mesh = scene.geometries().find(bl_leg_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
        currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
        stress = calculate_cauchy_stress(leg_Vs, leg_Ts, leg_moduli, leg_Dm_invs, currentVs);
        write_ply_face_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/bl_leg_face{}.ply", this_output_path, i));
        write_ply_vertex_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/bl_leg_vertex{}.ply", this_output_path, i));

        mesh = scene.geometries().find(br_leg_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
        currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
        stress = calculate_cauchy_stress(leg_Vs, leg_Ts, leg_moduli, leg_Dm_invs, currentVs);
        write_ply_face_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/br_leg_face{}.ply", this_output_path, i));
        write_ply_vertex_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/br_leg_vertex{}.ply", this_output_path, i));

        mesh = scene.geometries().find(fl_leg_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
        currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
        stress = calculate_cauchy_stress(leg_Vs, leg_Ts, leg_moduli, leg_Dm_invs, currentVs);
        write_ply_face_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/fl_leg_face{}.ply", this_output_path, i));
        write_ply_vertex_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/fl_leg_vertex{}.ply", this_output_path, i));

        mesh = scene.geometries().find(fr_leg_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
        currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
        stress = calculate_cauchy_stress(leg_Vs, leg_Ts, leg_moduli, leg_Dm_invs, currentVs);
        write_ply_face_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/fr_leg_face{}.ply", this_output_path, i));
        write_ply_vertex_stress(currentVs, leg_Ts, stress, fmt::format("{}table_scene/fr_leg_vertex{}.ply", this_output_path, i));

        mesh = scene.geometries().find(tabletop_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
        currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
        stress = calculate_cauchy_stress(tabletop_Vs, tabletop_Ts, tabletop_moduli, tabletop_Dm_invs, currentVs);
        write_ply_face_stress(currentVs, tabletop_Ts, stress, fmt::format("{}table_scene/tabletop_face{}.ply", this_output_path, i));
        write_ply_vertex_stress(currentVs, tabletop_Ts, stress, fmt::format("{}table_scene/tabletop_vertex{}.ply", this_output_path, i));

        mesh = scene.geometries().find(cube_object->geometries().ids()[0]).geometry->geometry().as<SimplicialComplex>();
        currentVs = mesh->vertices().find<Vector3>(builtin::position)->view();
        stress = calculate_cauchy_stress(cube_Vs, cube_Ts, tabletop_moduli, tabletop_Dm_invs, currentVs);
        write_ply_face_stress(currentVs, cube_Ts, stress, fmt::format("{}table_scene/cube{}.ply", this_output_path, i));
    }
}

int main() {
    cylinder_scene();

    return 0;
}

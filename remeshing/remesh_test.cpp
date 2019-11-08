//
// Created by ziqwang on 2019-09-27.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch2/catch.hpp"
#include "viewer.h"
#include <cmath>

int evaluate_mesh_valence(surface_mesh::Surface_mesh &mesh){

    Mesh::Vertex_iterator     v_it, v_end(mesh.vertices_end());
    int variance = 0;
    for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
    {
        if(mesh.is_boundary(*v_it)){
            variance += std::abs((int)mesh.valence(*v_it) - (int)4);
        }
        else{
            variance += std::abs((int)mesh.valence(*v_it) - (int)6);
        }
    }
    return variance;
}

TEST_CASE( "Valence Checker", "max.off" )
{
    //Loading the mesh.
    nanogui::init();

    nanogui::ref<Viewer> app = new Viewer();

    mesh_processing::MeshProcessing *meshProcess = app->mesh_;
    meshProcess->load_mesh("../data/max.off");
    meshProcess->set_target_length(2.0f);

    std::cout << "Before remeshing (Valence Quality):\t" << evaluate_mesh_valence(meshProcess->mesh_) << std::endl;

    mesh_processing::REMESHING_TYPE remesh_type = mesh_processing::UNIFORM;
    meshProcess->remesh(remesh_type, 5);
    meshProcess->compute_mesh_properties();

    std::cout << "After remeshing (Valence Quality):\t" <<  evaluate_mesh_valence(meshProcess->mesh_) << std::endl;
}
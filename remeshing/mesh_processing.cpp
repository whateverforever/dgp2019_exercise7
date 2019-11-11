//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Copyright (C) 2017 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//-----------------------------------------------------------------------------
#define _USE_MATH_DEFINES
#include "mesh_processing.h"
#include <cmath>
#include <set>
namespace mesh_processing {

using std::cout;
using std::endl;
using std::max;
using std::min;
using surface_mesh::Color;
using surface_mesh::Point;
using surface_mesh::Scalar;

MeshProcessing::MeshProcessing(const string& filename)
{
    load_mesh(filename);
}

MeshProcessing::~MeshProcessing()
{
}

void MeshProcessing::remesh(const REMESHING_TYPE& remeshing_type,
    const int& num_iterations)
{
    calc_weights();
    calc_mean_curvature();
    calc_uniform_mean_curvature();
    calc_gauss_curvature();
    calc_target_length(remeshing_type);

    // main remeshing loop
    for (int i = 0; i < num_iterations; ++i) {
        split_long_edges();
        collapse_short_edges();
        equalize_valences();
        tangential_relaxation();
    }
}

void MeshProcessing::calc_target_length(const REMESHING_TYPE& remeshing_type)
{
    Mesh::Vertex_iterator v_it, v_end(mesh_.vertices_end());
    Mesh::Vertex_around_vertex_circulator vv_c, vv_end;
    Scalar length;
    Scalar mean_length;
    Scalar H;
    Scalar K;

    Mesh::Vertex_property<Scalar> curvature = mesh_.vertex_property<Scalar>("v:curvature", 0);
    Mesh::Vertex_property<Scalar> gauss_curvature = mesh_.vertex_property<Scalar>("v:gauss_curvature", 0);
    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
    Mesh::Vertex_property<Scalar> target_new_length = mesh_.vertex_property<Scalar>("v:newlength", 0);

    if (remeshing_type == UNIFORM) {
        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
            target_length[*v_it] = target_length_;

    } else if (remeshing_type == ADAPTIVE) {
        // ------------- IMPLEMENT HERE ---------
        // NOTE: this *optional* exercise requires the methods calc_[mean/gauss]_curvature() methods.
        // These methods correspond to the methods you already implemented in Exercise 5, part 2 and you should
        // copy their content, and adapt it if necessary, to fill in the required properties described in the
        // comments of these functions below.
        //
        // In this function, you'll have to implement the following:
        //
        // Compute the maximal curvature at each vertex (use the precomputed mean and gaussian curvature).
        // Calculate the desired edge length as the target_length_ divided by the maximal curvature at each vertex, and assign it to the property target_length.
        // Smooth the target_length property uniformly, use the property target_new_length to store the smoothed values intermediately.
        // Rescale the property target_new_length such that it's mean equals the user specified target_length_.
        // Finally set the target_length property to this rescaled value.
        // ------------- IMPLEMENT HERE ---------

        // calculate desired length
        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
            length = 1.0;
            if (!mesh_.is_boundary(*v_it)) {
            }
            target_length[*v_it] = length;
        }

        // smooth desired length
        for (int i = 0; i < 5; i++) {
        }

        // rescale desired length
        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
        }
    }
}

void MeshProcessing::split_long_edges()
{
    Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
    Mesh::Vertex v0, v1, v;
    bool finished;
    int i;

    Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

    for (finished = false, i = 0; !finished && i < 100; ++i) {
        finished = true;
        // ------------- IMPLEMENT HERE ---------
        // INSERT CODE:
        //  Loop through all edges:
        //  Compute the desired length as the mean between the property "target_length" of two vertices of the edge
        //  If the edge is longer than 4/3 * desired length
        //		add the midpoint to the mesh: Mesh::Vertex v_new = _mesh.add_vertex( [Position of mid point] )
        //		set the interpolated "normals" and interpolated "target_length" property of the new vertex
        //		split the edge with this vertex: mesh_.split( [Edge], [Vertex] )
        //  If any edge has been processed this way, do not leave the outer iteration (i.e. set "finished" to false)
        // ------------- IMPLEMENT HERE ---------

        for (auto e : mesh_.edges()) {
            Mesh::Vertex vertA = mesh_.vertex(e, 0);
            Mesh::Vertex vertB = mesh_.vertex(e, 1);

            float length_desired = (target_length[vertA] + target_length[vertB]) / 2;
            float length_actual = mesh_.edge_length(e);

            if (length_actual > 4 / 3 * length_desired) {
                Point midPoint = (mesh_.position(vertA) + mesh_.position(vertB)) / 2;

                Mesh::Vertex v_new = mesh_.add_vertex(midPoint);

                Point normalA = normals[vertA];
                Point normalB = normals[vertB];
                Point normal_interpolated = (normalA + normalB) / 2;

                normals[v_new] = normal_interpolated;
                target_length[v_new] = length_desired;

                mesh_.split(e, v_new);

                finished = false;
            }
        }
    }
}

void MeshProcessing::collapse_short_edges()
{
    Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
    Mesh::Vertex v0, v1;
    Mesh::Halfedge h01, h10;
    bool finished, b0, b1;
    int i;
    bool hcol01, hcol10;

    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

    for (finished = false, i = 0; !finished && i < 100; ++i) {
        finished = true;

        for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {
            if (!mesh_.is_deleted(*e_it)) // might already be deleted
            {
                // ------------- IMPLEMENT HERE ---------
                // INSERT CODE:
                // Compute the desired length as the mean between the property "target_length" of two vertices of the edge
                // If the edge is shorter than 4/5 of the desired length
                //		Check if halfedge connects a boundary vertex with a non-boundary vertex, using  mesh_.is_boundary( [Vertex] ). If so, don't collapse that halfedge.
                //		Check if halfedges are collapsible using mesh_.is_collapse_ok( [halfedge] )
                //		Select the halfedge to be collapsed if at least one halfedge can be collapsed (if both are, collapse lower valence vertex into higher valence vertex, use mesh_.valence([Vertex]) ).
                //		Collapse the halfedge, using mesh_.collapse( [Halfedge] )
                // Stay in the loop running until no collapse has been done (use the finished variable)
                // ------------- IMPLEMENT HERE ---------
                Mesh::Vertex vertA = mesh_.vertex(*e_it, 0);
                Mesh::Vertex vertB = mesh_.vertex(*e_it, 1);

                float length_desired = (target_length[vertA] + target_length[vertB]) / 2;
                float length_actual = mesh_.edge_length(*e_it);

                if (length_actual < 4 / 5 * length_desired) {
                }
            }
        }
    }

    mesh_.garbage_collection();
    if (i == 100)
        std::cerr << "collapse break\n";
}

void MeshProcessing::equalize_valences()
{
    Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
    Mesh::Vertex v0, v1, v2, v3;
    Mesh::Halfedge h;
    int val0, val1, val2, val3;
    int val_opt0, val_opt1, val_opt2, val_opt3;
    int ve0, ve1, ve2, ve3, ve_before, ve_after;
    bool finished;
    int i;

    // flip all edges
    for (finished = false, i = 0; !finished && i < 100; ++i) {
        finished = true;

        for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {
            if (!mesh_.is_boundary(*e_it)) {
                // ------------- IMPLEMENT HERE ---------
                //  Extract valences of the four vertices involved to an eventual flip.
                //  Compute the sum of the squared valence deviances before flip (see assignment for a definition of the valence deviance!)
                //  Compute the sum of the squared valence deviances after an eventual flip
                //  If valence deviance is decreased by a flip and the flip is possible (mesh_.is_flip_ok([Edge])), flip the edge: mesh_.flip([Edge])
                //  Stay in the loop until no flip has been performed (use the finished variable)
                // ------------- IMPLEMENT HERE ---------
            }
        }
    }

    if (i == 100)
        std::cerr << "flip break\n";
}

void MeshProcessing::tangential_relaxation()
{
    Mesh::Vertex_iterator v_it, v_end(mesh_.vertices_end());
    Mesh::Vertex_around_vertex_circulator vv_c, vv_end;
    int valence;
    Point u, n;
    Point laplace;

    Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
    Mesh::Vertex_property<Point> update = mesh_.vertex_property<Point>("v:update");

    // smooth
    for (int iters = 0; iters < 10; ++iters) {
        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
            if (!mesh_.is_boundary(*v_it)) {
                // ------------- IMPLEMENT HERE ---------
                //  Compute uniform Laplacian curvature approximation vector, see the slides for the simple formula
                //  Compute the tangential component of the laplacian vector.
                //  Store this tangential component of the laplacian vector in the "update" property, where it will be used in the subsequent code to move the vertex.
                // ------------- IMPLEMENT HERE ---------
            }
        }

        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
            if (!mesh_.is_boundary(*v_it))
                mesh_.position(*v_it) += update[*v_it];
    }
}

void MeshProcessing::calc_uniform_mean_curvature()
{
    Mesh::Vertex_property<Scalar> v_unicurvature = mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);

    // ------------- IMPLEMENT HERE ---------
    // If you want to be able to visualize the uniform mean curvature vector approximation,
    // please adapt your solution of Exercise 5 to fill in the property v_unicurvature
    // with the uniform laplace mean curvature vector approximation computed
    // in Exercise 5, calc_uniform_laplacian.
    // ------------- IMPLEMENT HERE ---------
}

void MeshProcessing::calc_mean_curvature()
{
    Mesh::Vertex_property<Scalar> v_curvature = mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Edge_property<Scalar> e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    // ------------- IMPLEMENT HERE ---------
    // If you want to implement the optional adaptive remeshing exercise, please
    // adapt your solution of Exercise 5 to fill in the property v_curvature
    // with the laplace-beltrami mean curvature approximation computed
    // in Exercise 5, calc_mean_curvature.
    // ------------- IMPLEMENT HERE ---------
}

void MeshProcessing::calc_gauss_curvature()
{
    Mesh::Vertex_property<Scalar> v_gauss_curvature = mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    // ------------- IMPLEMENT HERE ---------
    // If you want to implement the optional adaptive remeshing exercise, please
    // adapt your solution of Exercise 5 to fill in the property v_gauss_curvature
    // with the Gauss curvature approximation you implemented in Exercise 5,
    // calc_gauss_curvature
    // ------------- IMPLEMENT HERE ---------
}

void MeshProcessing::calc_weights()
{
    calc_edges_weights();
    calc_vertices_weights();
}

void MeshProcessing::calc_edges_weights()
{
    auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
    auto points = mesh_.vertex_property<Point>("v:point");

    Mesh::Halfedge h0, h1, h2;
    Point p0, p1, p2, d0, d1;

    for (auto e : mesh_.edges()) {
        e_weight[e] = 0.0;

        h0 = mesh_.halfedge(e, 0);
        p0 = points[mesh_.to_vertex(h0)];

        h1 = mesh_.halfedge(e, 1);
        p1 = points[mesh_.to_vertex(h1)];

        if (!mesh_.is_boundary(h0)) {
            h2 = mesh_.next_halfedge(h0);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0, d1) / norm(cross(d0, d1));
        }

        if (!mesh_.is_boundary(h1)) {
            h2 = mesh_.next_halfedge(h1);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0, d1) / norm(cross(d0, d1));
        }
    }
}

void MeshProcessing::calc_vertices_weights()
{
    Mesh::Face_around_vertex_circulator vf_c, vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar area;
    auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    for (auto v : mesh_.vertices()) {
        area = 0.0;
        vf_c = mesh_.faces(v);

        if (!vf_c) {
            continue;
        }

        vf_end = vf_c;

        do {
            fv_c = mesh_.vertices(*vf_c);

            const Point& P = mesh_.position(*fv_c);
            ++fv_c;
            const Point& Q = mesh_.position(*fv_c);
            ++fv_c;
            const Point& R = mesh_.position(*fv_c);

            area += norm(cross(Q - P, R - P)) * 0.5f * 0.3333f;

        } while (++vf_c != vf_end);

        v_weight[v] = 0.5 / area;
    }
}

void MeshProcessing::load_mesh(const string& filename)
{
    if (!mesh_.read(filename)) {
        std::cerr << "Mesh not found, exiting." << std::endl;
        exit(-1);
    }

#ifndef UNIT_TEST

    cout << "Mesh " << filename << " loaded." << endl;
    cout << "# of vertices : " << mesh_.n_vertices() << endl;
    cout << "# of faces : " << mesh_.n_faces() << endl;
    cout << "# of edges : " << mesh_.n_edges() << endl;
#endif
    // Compute the center of the mesh
    mesh_center_ = Point(0.0f, 0.0f, 0.0f);
    for (auto v : mesh_.vertices()) {
        mesh_center_ += mesh_.position(v);
    }
    mesh_center_ /= mesh_.n_vertices();

    // Compute the maximum distance from all points in the mesh and the center
    dist_max_ = 0.0f;
    for (auto v : mesh_.vertices()) {
        if (distance(mesh_center_, mesh_.position(v)) > dist_max_) {
            dist_max_ = distance(mesh_center_, mesh_.position(v));
        }
    }

    compute_mesh_properties();

    // Store the original mesh, this might be useful for some computations
    mesh_init_ = mesh_;
}

void MeshProcessing::compute_mesh_properties()
{
    Mesh::Vertex_property<Point> vertex_normal = mesh_.vertex_property<Point>("v:normal");
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
    Mesh::Vertex_property<Color> v_color_valence = mesh_.vertex_property<Color>("v:color_valence",
        Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_unicurvature = mesh_.vertex_property<Color>("v:color_unicurvature",
        Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_curvature = mesh_.vertex_property<Color>("v:color_curvature",
        Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_gaussian_curv = mesh_.vertex_property<Color>("v:color_gaussian_curv",
        Color(1.0f, 1.0f, 1.0f));

    Mesh::Vertex_property<Scalar> vertex_valence = mesh_.vertex_property<Scalar>("v:valence", 0.0f);
    for (auto v : mesh_.vertices()) {
        vertex_valence[v] = mesh_.valence(v);
    }

    Mesh::Vertex_property<Scalar> v_unicurvature = mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_curvature = mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_gauss_curvature = mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);

    calc_weights();
    calc_uniform_mean_curvature();
    calc_mean_curvature();
    calc_gauss_curvature();
    color_coding(vertex_valence, &mesh_, v_color_valence, 3 /* min */,
        8 /* max */);
    color_coding(v_unicurvature, &mesh_, v_color_unicurvature);
    color_coding(v_curvature, &mesh_, v_color_curvature);
    color_coding(v_gauss_curvature, &mesh_, v_color_gaussian_curv);

    // get the mesh attributes and upload them to the GPU
    int j = 0;
    unsigned int n_vertices(mesh_.n_vertices());

    // Create big matrices to send the data to the GPU with the required
    // format
    color_valence_ = Eigen::MatrixXf(3, n_vertices);
    color_unicurvature_ = Eigen::MatrixXf(3, n_vertices);
    color_curvature_ = Eigen::MatrixXf(3, n_vertices);
    color_gaussian_curv_ = Eigen::MatrixXf(3, n_vertices);
    normals_ = Eigen::MatrixXf(3, n_vertices);
    points_ = Eigen::MatrixXf(3, n_vertices);
    indices_ = MatrixXu(3, mesh_.n_faces());

    for (auto f : mesh_.faces()) {
        std::vector<float> vv(3);
        int k = 0;
        for (auto v : mesh_.vertices(f)) {
            vv[k] = v.idx();
            ++k;
        }
        indices_.col(j) << vv[0], vv[1], vv[2];
        ++j;
    }

    j = 0;
    for (auto v : mesh_.vertices()) {
        points_.col(j) << mesh_.position(v).x,
            mesh_.position(v).y,
            mesh_.position(v).z;

        normals_.col(j) << vertex_normal[v].x,
            vertex_normal[v].y,
            vertex_normal[v].z;

        color_valence_.col(j) << v_color_valence[v].x,
            v_color_valence[v].y,
            v_color_valence[v].z;

        color_unicurvature_.col(j) << v_color_unicurvature[v].x,
            v_color_unicurvature[v].y,
            v_color_unicurvature[v].z;

        color_curvature_.col(j) << v_color_curvature[v].x,
            v_color_curvature[v].y,
            v_color_curvature[v].z;

        color_gaussian_curv_.col(j) << v_color_gaussian_curv[v].x,
            v_color_gaussian_curv[v].y,
            v_color_gaussian_curv[v].z;
        ++j;
    }
}

void MeshProcessing::color_coding(Mesh::Vertex_property<Scalar> prop, Mesh* mesh,
    Mesh::Vertex_property<Color> color_prop, Scalar min_value,
    Scalar max_value, int bound)
{
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    if (min_value == 0.0 && max_value == 0.0) {
        // discard upper and lower bound
        unsigned int n = values.size() - 1;
        unsigned int i = n / bound;
        std::sort(values.begin(), values.end());
        min_value = values[i];
        max_value = values[n - 1 - i];
    }

    // map values to colors
    for (auto v : mesh->vertices()) {
        set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
    }
}

void MeshProcessing::set_color(Mesh::Vertex v, const Color& col,
    Mesh::Vertex_property<Color> color_prop)
{
    color_prop[v] = col;
}

Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value)
{
    Scalar v0, v1, v2, v3, v4;
    v0 = min_value + 0.0 / 4.0 * (max_value - min_value);
    v1 = min_value + 1.0 / 4.0 * (max_value - min_value);
    v2 = min_value + 2.0 / 4.0 * (max_value - min_value);
    v3 = min_value + 3.0 / 4.0 * (max_value - min_value);
    v4 = min_value + 4.0 / 4.0 * (max_value - min_value);

    Color col(1.0f, 1.0f, 1.0f);

    if (value < v0) {
        col = Color(0, 0, 1);
    } else if (value > v4) {
        col = Color(1, 0, 0);
    } else if (value <= v2) {
        if (value <= v1) { // [v0, v1]
            Scalar u = (value - v0) / (v1 - v0);
            col = Color(0, u, 1);
        } else { // ]v1, v2]
            Scalar u = (value - v1) / (v2 - v1);
            col = Color(0, 1, 1 - u);
        }
    } else {
        if (value <= v3) { // ]v2, v3]
            Scalar u = (value - v2) / (v3 - v2);
            col = Color(u, 1, 0);
        } else { // ]v3, v4]
            Scalar u = (value - v3) / (v4 - v3);
            col = Color(1, 1 - u, 0);
        }
    }
    return col;
}
}

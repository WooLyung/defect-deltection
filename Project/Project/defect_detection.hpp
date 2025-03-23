#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <iostream>
#include <vector>

#pragma comment(lib, "gmp.lib")
using namespace std;
using namespace string_literals;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;

typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Segment_3 Segment;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::IO::Color Color;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> Tree;
typedef CGAL::Side_of_triangle_mesh<Mesh, K> Inside_tester;

typedef CGAL::Search_traits_3<K> Traits;
typedef CGAL::Kd_tree<Traits> AABBTree;
typedef CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;
typedef Neighbor_search::Distance Distance;

typedef pair<Point, int> Point_with_index;
typedef CGAL::Search_traits_3<K> Base_traits;
typedef CGAL::Search_traits_adapter<Point_with_index, CGAL::First_of_pair_property_map<Point_with_index>, Base_traits> Traits2;
typedef CGAL::Kd_tree<Traits2> AABBTree2;
typedef CGAL::Orthogonal_k_neighbor_search<Traits2> Neighbor_search2;

struct MRH {
    vector<tuple<int, int, int>> H;
    vector<Point> P;
};

struct Skeleton {
    vector<Point> V;
    vector<pair<int, int>> E;
};

Mesh read_mesh(const char* path);
void write_mesh(Mesh mesh, const char* path);
bool is_segment_inside(Tree& tree, Inside_tester& inside_tester, Segment& segment);
void laplacian_smoothing(Skeleton& S, int iter, double lambda);
Skeleton medial_skeleton(Mesh& M);
MRH construct_MRH(Skeleton& S);
void adjust_skeleton(MRH& T, vector<Point>& V_o, vector<Point>& V_d, double alpha, double beta, int delta);
vector<double> calculate_defect(Mesh& M_o, Mesh& M_d, vector<Point>& V_d, vector<Point>& V_d_clone);
vector<double> defect_detection(Mesh& M_o, Mesh& M_d, double alpha, double beta, int delta);
double compute_vertex_mean_error(vector<double>& D);
double compute_weighted_mean_error(Mesh& M_d, vector<double>& D);
double compute_maximum_error(vector<double>& D);
double compute_percentages_with_error(vector<double>& D, double d);
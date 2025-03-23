#include "defect_detection.hpp"

void visualize(vector<double>& D, Mesh& M_d);
void printError(vector<double>& D, Mesh& M_d, double alpha, double beta, int delta);

int main()
{
    Mesh M_o = read_mesh("C:/_Research/mesh/A.off");
    Mesh M_d = read_mesh("C:/_Research/mesh/A_scan.off");

    double alpha = 0.1;
    double beta = 0;
    double delta = 0;

    vector<double> D = defect_detection(M_o, M_d, alpha, beta, delta);
    printError(D, M_d, alpha, beta, delta);

    visualize(D, M_d);
    write_mesh(M_d, "C:/_Research/result/A_scan_defect.off");
}

void printError(vector<double>& D, Mesh& M_d, double alpha, double beta, int delta)
{
    double e0 = compute_maximum_error(D);
    double e1 = compute_vertex_mean_error(D);
    double e2 = compute_weighted_mean_error(M_d, D);
    double e3[] = {
        compute_percentages_with_error(D, 0.3),
        compute_percentages_with_error(D, 0.4),
        compute_percentages_with_error(D, 0.5),
        compute_percentages_with_error(D, 0.6),
        compute_percentages_with_error(D, 0.7),
        compute_percentages_with_error(D, 0.8),
        compute_percentages_with_error(D, 0.9),
        compute_percentages_with_error(D, 1.0),
    };
    printf("(%4.2lf, %4.2lf, %4d) %8.6lf & %8.6lf & %8.6lf & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf\\\\\n", alpha, beta, delta, e0, e1, e2, e3[0], e3[1], e3[2], e3[3], e3[4], e3[5], e3[6], e3[7]);
}

void visualize(vector<double>& D, Mesh& M_d)
{
    auto face_color_map = M_d.add_property_map<Mesh::Face_index, Color>("f:color", Color(255, 255, 255)).first;
    for (auto face : M_d.faces())
    {
        double x = 0;
        int c = 0;

        for (auto vertex : CGAL::vertices_around_face(M_d.halfedge(face), M_d))
            x += D[vertex.idx()], c++;
        x /= c;

        if (x <= 0.3)
            face_color_map[face] = Color(255, 255, 255);
        else if (x <= 0.4)
            face_color_map[face] = Color(255, 255, 0);
        else if (x <= 0.5)
            face_color_map[face] = Color(127, 255, 255);
        else if (x <= 0.6)
            face_color_map[face] = Color(0, 255, 0);
        else if (x <= 0.7)
            face_color_map[face] = Color(0, 127, 127);
        else if (x <= 0.8)
            face_color_map[face] = Color(0, 0, 255);
        else if (x <= 0.9)
            face_color_map[face] = Color(0, 0, 127);
        else
            face_color_map[face] = Color(0, 0, 0);
    }
}
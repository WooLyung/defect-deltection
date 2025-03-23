#include "defect_detection.hpp"

Mesh read_mesh(const char* path)
{
    Mesh mesh;
    ifstream input(path);
    input >> mesh;

    return mesh;
}

void write_mesh(Mesh mesh, const char* path)
{
    ofstream output(path);
    output << mesh;
}

bool is_segment_inside(Tree& tree, Inside_tester& inside_tester, Segment& segment)
{
    return !tree.do_intersect(segment) && inside_tester(segment.source()) == CGAL::ON_BOUNDED_SIDE && inside_tester(segment.target()) == CGAL::ON_BOUNDED_SIDE;
}

void laplacian_smoothing(Skeleton& S, int iter, double lambda)
{
    cout << "laplacin_smoothing()" << endl;

    vector<set<int>> edges(S.V.size());
    for (const auto& edge : S.E)
    {
        edges[edge.first].insert(edge.second);
        edges[edge.second].insert(edge.first);
    }

    for (int i = 0; i < iter; i++)
    {
        vector<Point> vertices2 = S.V;
        for (int u = 0; u < S.V.size(); u++)
        {
            if (edges[u].size() == 0)
                continue;

            Vector dir = Vector(0, 0, 0);
            for (int v : edges[u])
                dir += S.V[v] - Point(0, 0, 0);
            dir = (dir / edges[u].size()) - (S.V[u] - Point(0, 0, 0));
            vertices2[u] += dir * lambda;
        }
        S.V = vertices2;
    }
}

Skeleton medial_skeleton(Mesh& M)
{
    cout << "medial_skeleton()" << endl;

    Skeleton S;
    Delaunay dt;
    map<Delaunay::Vertex_handle, int> t_handles;

    Tree tree(M.faces().first, M.faces().second, M);
    Inside_tester inside_tester(M);
    tree.accelerate_distance_queries();
    dt.insert(M.points().begin(), M.points().end());

    for (auto it = dt.finite_edges_begin(); it != dt.finite_edges_end(); it++)
    {
        try
        {
            Delaunay::Cell_handle cell = it->first;
            int i = it->second, j = it->third;
            Delaunay::Cell_handle neighbor = cell->neighbor(j);

            if (!dt.is_infinite(cell) && !dt.is_infinite(neighbor))
            {
                Point p1 = dt.dual(cell);
                Point p2 = dt.dual(neighbor);
                Segment segment(p1, p2);

                if (is_segment_inside(tree, inside_tester, segment))
                {
                    Delaunay::Vertex_handle v1 = cell->vertex(i);
                    Delaunay::Vertex_handle v2 = cell->vertex(j);

                    if (t_handles.find(v1) == t_handles.end())
                        t_handles.insert({ v1, t_handles.size() }), S.V.push_back(p1);
                    if (t_handles.find(v2) == t_handles.end())
                        t_handles.insert({ v2, t_handles.size() }), S.V.push_back(p2);

                    int v1i = t_handles[v1];
                    int v2i = t_handles[v2];
                    S.E.push_back({ v1i, v2i });
                }
            }
        }
        catch (exception)
        {
        }
    }

    vector<vector<int>> adj_list(S.V.size());
    vector<bool> visited(S.V.size(), false);
    vector<int> size(S.V.size(), 0);
    vector<int> component(S.V.size());
    map<int, int> compress;
    int biggest = 0;

    for (auto e : S.E)
    {
        adj_list[e.first].push_back(e.second);
        adj_list[e.second].push_back(e.first);
    }

    queue<int> q;
    for (int i = 0; i < S.V.size(); i++)
    {
        if (visited[i])
            continue;
        visited[i] = true;
        size[i]++;
        component[i] = i;
        q.push(i);

        while (!q.empty())
        {
            int x = q.front();
            q.pop();
            for (int y : adj_list[x])
            {
                if (visited[y])
                    continue;
                visited[y] = true;
                size[i]++;
                component[y] = i;
                q.push(y);

                if (size[biggest] < size[i])
                    biggest = i;
            }
        }
    }

    Skeleton filtered_S;
    for (int i = 0; i < S.V.size(); i++)
    {
        if (component[i] != biggest)
            continue;
        filtered_S.V.push_back(S.V[i]);
        compress.insert({ i, compress.size() });
    }

    for (auto e : S.E)
    {
        if (component[e.first] != biggest || component[e.second] != biggest)
            continue;
        filtered_S.E.push_back({ compress[e.first], compress[e.second] });
    }

    laplacian_smoothing(filtered_S, 100, 1.0);
    return filtered_S;
}

MRH construct_MRH(Skeleton& S)
{
    cout << "construct_MRH()" << endl;

    MRH T;

    vector<bool> isMerged;
    vector<set<int>> edges;
    priority_queue<tuple<float, int, int>> pq;
    set<int> curr;

    for (int i = 0; i < S.V.size(); i++)
    {
        T.P.push_back(S.V[i]);
        isMerged.push_back(false);
        edges.push_back(set<int>());
        curr.insert(i);
    }

    for (int i = 0; i < S.E.size(); i++)
    {
        pq.push({ -CGAL::squared_distance(T.P[S.E[i].first], T.P[S.E[i].second]), S.E[i].first, S.E[i].second });
        edges[S.E[i].first].insert(S.E[i].second);
        edges[S.E[i].second].insert(S.E[i].first);
    }

    while (!pq.empty())
    {
        auto top = pq.top();
        pq.pop();

        int i = get<1>(top);
        int j = get<2>(top);
        if (isMerged[i] || isMerged[j])
            continue;

        int u = T.P.size();
        T.H.push_back({ i, j, u });

        Point pos = Point((T.P[i].x() + T.P[j].x()) * 0.5f, (T.P[i].y() + T.P[j].y()) * 0.5f, (T.P[i].z() + T.P[j].z()) * 0.5f);
        isMerged[i] = isMerged[j] = true;
        T.P.push_back(pos);
        isMerged.push_back(false);
        edges.push_back(set<int>());

        for (int v : edges[i])
        {
            if (isMerged[v])
                continue;

            edges[u].insert(v);
            edges[v].insert(u);
            edges[v].erase(i);
            edges[v].erase(j);
            pq.push({ -CGAL::squared_distance(T.P[u], T.P[v]), u, v });
        }

        for (int v : edges[j])
        {
            if (isMerged[v])
                continue;

            edges[u].insert(v);
            edges[v].insert(u);
            edges[v].erase(i);
            edges[v].erase(j);
            pq.push({ -CGAL::squared_distance(T.P[u], T.P[v]), u, v });
        }

        edges[i].clear();
        edges[j].clear();
        edges[u].erase(i);
        edges[u].erase(j);

        curr.erase(i);
        curr.erase(j);
        curr.insert(u);
    }

    return T;
}

void adjust_skeleton(MRH& T, vector<Point>& V_o, vector<Point>& V_d, double alpha, double beta, int delta)
{
    cout << "adjust_skeleton()" << endl;

    AABBTree aabb_origin(V_o.begin(), V_o.end());
    map<int, Vector> trans{ { get<2>(*T.H.rbegin()), Vector(0, 0, 0) } };
    set<int> curV_d{ get<2>(*T.H.rbegin()) };

    for (int t = T.H.size() - 1; t >= 0; t--)
    {
        if (curV_d.size() <= delta)
        {
            vector<pair<Point, int>> V_d_with_index;
            map<int, pair<Vector, int>> tempCount;
            for (int x : curV_d)
            {
                V_d_with_index.push_back({ T.P[x] + trans[x], x });
                tempCount.insert({ x, { Vector(0, 0, 0), 0 } });
            }
            AABBTree2 aabb_points(V_d_with_index.begin(), V_d_with_index.end());

            for (Point q : V_o)
            {
                Point_with_index c = Neighbor_search2(aabb_points, q).begin()->first;
                tempCount[c.second].first += q - c.first;
                tempCount[c.second].second++;
            }

            for (int p : curV_d)
                if (tempCount[p].second != 0)
                    trans[p] += (tempCount[p].first / tempCount[p].second) * beta;
        }
        else
        {
            int w = get<2>(T.H[t]);
            Point c = Neighbor_search(aabb_origin, T.P[w] + trans[w]).begin()->first;
            trans[w] = c - T.P[w];
        }

        auto tuple = T.H[t];
        int u = get<0>(tuple);
        int v = get<1>(tuple);
        int w = get<2>(tuple);

        curV_d.erase(w);
        curV_d.insert(u);
        curV_d.insert(v);
        trans[u] = trans[v] = trans[w] * alpha;
    }

    for (int t = 0; t < V_d.size(); t++)
    {
        Point c = Neighbor_search(aabb_origin, V_d[t] + trans[t]).begin()->first;
        V_d[t] = c;
    }
}

vector<double> calculate_defect(Mesh& M_o, Mesh& M_d, vector<Point>& V_d, vector<Point>& V_d_clone)
{
    cout << "calculate_defect()" << endl;

    vector<double> D;
    vector<pair<Point, int>> vert1;
    for (int i = 0; i < V_d_clone.size(); i++)
        vert1.push_back({ V_d_clone[i], i });
    AABBTree2 aabb1(vert1.begin(), vert1.end());

    vector<pair<Point, int>> vert2;
    for (int i = 0; i < M_o.number_of_vertices(); i++)
        vert2.push_back({ M_o.point(Mesh::Vertex_index(i)) , i });
    AABBTree2 aabb2(vert2.begin(), vert2.end());

    for (auto vertex : M_d.vertices())
    {
        Neighbor_search2 search(aabb1, M_d.point(vertex));
        auto closet = search.begin()->first;
        int i = closet.second;

        Neighbor_search2 search2(aabb2, M_d.point(vertex) + (V_d[i] - V_d_clone[i]));
        auto closet2 = search2.begin()->first;
        D.push_back(sqrt((closet2.first - M_d.point(vertex)).squared_length()));
    }

    return D;
}

vector<double> defect_detection(Mesh& M_o, Mesh& M_d, double alpha, double beta, int delta)
{
    Skeleton S_o = medial_skeleton(M_o);
    Skeleton S_d = medial_skeleton(M_d);
    vector<Point> V_d_clone = S_d.V;
    MRH T = construct_MRH(S_d);
    adjust_skeleton(T, S_o.V, S_d.V, alpha, beta, delta);
    vector<double> D = calculate_defect(M_o, M_d, S_d.V, V_d_clone);
    return D;
}

double compute_vertex_mean_error(vector<double>& D)
{
    double error = 0;
    for (double e : D)
        error += e;
    return error /= D.size();
}

double compute_weighted_mean_error(Mesh& M_d, vector<double>& D)
{
    std::vector<double> vertex_areas(M_d.num_vertices(), 0.0);
    double area = 0.0;

    for (const auto& face : M_d.faces())
    {
        std::vector<Mesh::Vertex_index> vertices;
        for (const auto& vertex : CGAL::vertices_around_face(M_d.halfedge(face), M_d))
            vertices.push_back(vertex);

        Point p0 = M_d.point(vertices[0]);
        for (int i = 1; vertices.size() > i + 1; i++) {
            Point p1 = M_d.point(vertices[i]);
            Point p2 = M_d.point(vertices[i + 1]);

            double tri_area = 0.5 * CGAL::sqrt(CGAL::squared_area(p0, p1, p2));
            area += tri_area;

            vertex_areas[vertices[0]] += tri_area / 3.0;
            vertex_areas[vertices[i]] += tri_area / 3.0;
            vertex_areas[vertices[i + 1]] += tri_area / 3.0;
        }
    }

    double error = 0.0;
    for (size_t i = 0; i < D.size(); i++)
        error += D[i] * vertex_areas[i];

    return error / area;
}

double compute_maximum_error(vector<double>& D)
{
    double error = D[0];
    for (double e : D)
        error = max(e, error);
    return error;
}

double compute_percentages_with_error(vector<double>& D, double d)
{
    int cnt = 0;
    for (double e : D)
        if (e <= d)
            cnt++;
    return (double)cnt / D.size() * 100.0;
}
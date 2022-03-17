#include "plugin.hh"

#ifdef GUROBI
GRBEnv env = GRBEnv(true);
#endif

int seed = 0;
int n_constraints = 50;
double theta_triangulate = 30 * M_PI / 180;

double theta = 30 * M_PI / 180;
bool use_halfedge_collapses = true;
bool use_edge_collapses = true;
bool use_triangle_collapses = true;
bool prio = true;
bool prio_triangle_collapses = true;
int strategy = 3;

int runtime;

void AngleBounded2DSimplificationPlugin::initializePlugin() {
#ifdef GUROBI
    env.set(GRB_IntParam_LogToConsole, false);
    env.start();
#endif
    std::setlocale(LC_NUMERIC, "en_US.UTF-8");
    id = -1;
    if (OpenFlipper::Options::gui()) {
        QWidget* toolbox = new QWidget();
        QGridLayout* layout = new QGridLayout(toolbox);

        QLabel* label = new QLabel("Seed");
        layout->addWidget(label, 0, 0);
        in_seed = new QSpinBox(toolbox);
        in_seed->setValue(seed);
        layout->addWidget(in_seed, 0, 1);

        label = new QLabel("Number of constraints");
        layout->addWidget(label, 1, 0);
        in_n_constraints = new QSpinBox(toolbox);
        in_n_constraints->setValue(n_constraints);
        layout->addWidget(in_n_constraints, 1, 1);

        label = new QLabel("Triangulation angle (in degrees)");
        layout->addWidget(label, 2, 0);
        in_theta_triangulate = new QDoubleSpinBox(toolbox);
        in_theta_triangulate->setValue(theta_triangulate * 180 / M_PI);
        layout->addWidget(in_theta_triangulate, 2, 1);

        QPushButton* triangulateRandomButton = new QPushButton("Triangulate", toolbox);
        layout->addWidget(triangulateRandomButton, 3, 0);
        connect(triangulateRandomButton, SIGNAL(clicked()), this, SLOT(triangulate_random()));

        label = new QLabel("Decimation angle (in degrees)");
        layout->addWidget(label, 5, 0);
        in_theta = new QDoubleSpinBox(toolbox);
        in_theta->setValue(theta * 180 / M_PI);
        layout->addWidget(in_theta, 5, 1);

        label = new QLabel("Use halfedge collapses");
        layout->addWidget(label, 6, 0);
        in_use_halfedge_collapses = new QCheckBox(toolbox);
        in_use_halfedge_collapses->setChecked(use_halfedge_collapses);
        layout->addWidget(in_use_halfedge_collapses, 6, 1);

        label = new QLabel("Use edge collapses");
        layout->addWidget(label, 7, 0);
        in_use_edge_collapses = new QCheckBox(toolbox);
        in_use_edge_collapses->setChecked(use_edge_collapses);
        layout->addWidget(in_use_edge_collapses, 7, 1);

        label = new QLabel("Use triangle collapses");
        layout->addWidget(label, 8, 0);
        in_use_triangle_collapses = new QCheckBox(toolbox);
        in_use_triangle_collapses->setChecked(use_triangle_collapses);
        layout->addWidget(in_use_triangle_collapses, 8, 1);

        label = new QLabel("Prioritize collapses with greater minimum angle");
        layout->addWidget(label, 9, 0);
        in_prio = new QCheckBox(toolbox);
        in_prio->setChecked(prio);
        layout->addWidget(in_prio, 9, 1);

        label = new QLabel("Prioritize triangle collapses (only with general prioritization)");
        layout->addWidget(label, 10, 0);
        in_prio_triangle_collapses = new QCheckBox(toolbox);
        in_prio_triangle_collapses->setChecked(prio_triangle_collapses);
        layout->addWidget(in_prio_triangle_collapses, 10, 1);

        label = new QLabel("Collapse position strategy\n0: Manual\n1: Gurobi case 1\n2: Gurobi case 2\n3: Gradient ascent\nElse: Center");
        layout->addWidget(label, 11, 0);
        in_strategy = new QSpinBox(toolbox);
        in_strategy->setValue(strategy);
        layout->addWidget(in_strategy, 11, 1);

        QPushButton* decimateButton = new QPushButton("Decimate", toolbox);
        layout->addWidget(decimateButton, 12, 0);
        connect(decimateButton, SIGNAL(clicked()), this, SLOT(start_thread()));

        emit addToolbox("Angle-Bounded 2D Mesh Simplification", toolbox);
    }
}

void AngleBounded2DSimplificationPlugin::pluginsInitialized(const QVector<QPair<QString, QString>>& pluginOptions) {
    if (!OpenFlipper::Options::gui()) {
        std::string in_filename = "";
        std::string out_filename = "";
        for (auto it = pluginOptions.begin(); it != pluginOptions.end(); it++) {
            if (it->first.toStdString() == "in") {
                in_filename = it->second.toStdString();
            }
            if (it->first.toStdString() == "out") {
                out_filename = it->second.toStdString();
            }
        }
        if (in_filename == "") {
            triangulate_random();
        } else {
            triangulate_file(in_filename);
        }
        decimate();
        if (out_filename != "") {
            OpenMesh::IO::write_mesh(*mesh, out_filename);
        }
    }
}

void AngleBounded2DSimplificationPlugin::get_parameters() {
    if (OpenFlipper::Options::gui()) {
        seed = in_seed->value();
        n_constraints = in_n_constraints->value();
        theta_triangulate = in_theta_triangulate->value() * M_PI / 180;
        theta = in_theta->value() * M_PI / 180;
        prio = in_prio->isChecked();
        use_halfedge_collapses = in_use_halfedge_collapses->isChecked();
        use_edge_collapses = in_use_edge_collapses->isChecked();
        use_triangle_collapses = in_use_triangle_collapses->isChecked();
        prio_triangle_collapses = in_prio_triangle_collapses->isChecked();
        strategy = in_strategy->value();
    }
}

// --------------------------------------------------------------------------------------------------------------------------------

void AngleBounded2DSimplificationPlugin::intersection_test(std::list<std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d>>& data) {
    while (!data.empty()) {
        auto s1 = data.front();
        data.pop_front();

        bool is_intersection = false;
        for (unsigned int j = 0; j < m_features.size(); j++) {
            auto s2 = m_features[j];
            double t1, t2;
            if (ACG::Geometry::lineIntersection(s1.first, s1.second, s2.first, s2.second, t1, t2)) {
                if (t1 > 0.9999 && t1 < 1.0)
                    t1 = 1.0;
                if (t2 > 0.9999 && t2 < 1.0)
                    t2 = 1.0;
                if (t1 > 0.0 && t1 < 0.0001)
                    t1 = 0.0;
                if (t2 > 0.0 && t2 < 0.0001)
                    t2 = 0.0;

                if ((t1 >= 0.0 && t1 <= 1.0 && t2 > 0.0 && t2 < 1.0) ||
                    (t1 > 0.0 && t1 < 1.0 && t2 >= 0.0 && t2 <= 1.0)) {
                    auto q = s2.first + t2 * (s2.second - s2.first);
                    if (0.0 == t1)
                        s1 = std::make_pair(q, s2.second);
                    if (1.0 == t1)
                        s1 = std::make_pair(s2.first, q);
                    if (t2 > 0.0 && t2 < 1.0) {
                        m_features[j] = std::make_pair(s2.first, q);
                        m_features.push_back(std::make_pair(q, s2.second));
                    }
                    if (t1 > 0.0 && t1 < 1.0) {
                        data.push_back(std::make_pair(s1.first, q));
                        data.push_back(std::make_pair(q, s1.second));
                    } else
                        data.push_front(s1);
                    is_intersection = true;
                    break;
                }
            }
        }
        if (!is_intersection)
            m_features.push_back(s1);
    }
}

void AngleBounded2DSimplificationPlugin::get_constraints() {
    m_features.clear();

    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(-1, 1);
    auto random_point = [gen, dis]() mutable {
        double x = dis(gen);
        double y = dis(gen);
        return OpenMesh::Vec2d(x, y);
    };

    std::list<std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d>> data;
    for (int i = 0; i < n_constraints; i++) {
        OpenMesh::Vec2d first = random_point();
        OpenMesh::Vec2d second = random_point();
        std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d> s = std::make_pair(first, second);
        if (0.00001 == (s.first - s.second).sqrnorm()) {
            i--;
            continue;
        }

        bool is_overlapping = false;
        for (auto s1 : data) {
            if ((s.second[1] - s.first[1]) * (s1.second[0] - s1.first[0]) == (s1.second[1] - s1.first[1]) * (s.second[0] - s.first[0])) {
                is_overlapping = true;
                i--;
                break;
            }
        }
        if (!is_overlapping)
            m_features.push_back(s);
    }

    unsigned int old_count = 0;
    while (old_count != m_features.size()) {
        data.clear();
        data.insert(data.begin(), m_features.begin(), m_features.end());
        m_features.clear();
        old_count = data.size();
        intersection_test(data);
    }
}

int add_point_constraint(OpenMesh::Vec2d p, std::map<OpenMesh::Vec2d, int>& in_pt_map, std::vector<double>& in_points, int& id) {
    if (in_pt_map.end() == in_pt_map.find(p)) {
        in_pt_map[p] = id++;
        in_points.emplace_back(p[0]);
        in_points.emplace_back(p[1]);
        return (id - 1);
    }
    return in_pt_map[p];
}

void add_segment_constraint(int i1, int i2, int& seg_id, std::vector<int>& segments, std::vector<int>& segment_markers) {
    segments.emplace_back(i1);
    segments.emplace_back(i2);
    segment_markers.emplace_back(seg_id++);
}

std::vector<double> getBbox(std::vector<std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d>>& data) {
    std::vector<double> bbox = {0.0, 0.0, 0.0, 0.0};
    if (data.size() <= 0)
        return bbox;

    double max_x, min_x, max_y, min_y;
    min_x = max_x = data[0].first[0];
    min_y = max_y = data[0].first[1];
    for (unsigned int i = 0; i < data.size(); i++) {
        if (data[i].first[0] < min_x)
            min_x = data[i].first[0];
        else if (data[i].first[0] > max_x)
            max_x = data[i].first[0];
        if (data[i].second[0] < min_x)
            min_x = data[i].second[0];
        else if (data[i].second[0] > max_x)
            max_x = data[i].second[0];

        if (data[i].first[1] < min_y)
            min_y = data[i].first[1];
        else if (data[i].first[1] > max_y)
            max_y = data[i].first[1];
        if (data[i].second[1] < min_y)
            min_y = data[i].second[1];
        else if (data[i].second[1] > max_y)
            max_y = data[i].second[1];
    }

    double avg_x = 0.5 * (max_x + min_x);
    double d_x = 0.5 * (max_x - min_x);
    double avg_y = 0.5 * (max_y + min_y);
    double d_y = 0.5 * (max_y - min_y);

    d_x = std::max(d_x, d_y);
    d_y = d_x;

    bbox[0] = avg_x;
    bbox[1] = d_x;
    bbox[2] = avg_y;
    bbox[3] = d_y;
    return bbox;
}

void AngleBounded2DSimplificationPlugin::delaunayMeshing() {
    std::map<OpenMesh::Vec2d, int> in_pt_map;
    std::vector<double> in_points;
    int p_id = 0, seg_id = 1;
    std::vector<int> in_segments;
    std::vector<int> in_segment_markers;
    std::map<std::pair<int, int>, int> segment_set;
    for (auto h : m_features) {
        int id1 = add_point_constraint(h.first, in_pt_map, in_points, p_id);
        int id2 = add_point_constraint(h.second, in_pt_map, in_points, p_id);
        add_segment_constraint(id1, id2, seg_id, in_segments, in_segment_markers);
    }
    if (bb) {
        std::vector<double> bbox = getBbox(m_features);
        std::vector<OpenMesh::Vec2d> bbox_pt;
        double bbox_offset_factor = 1.5;
        bbox[1] *= bbox_offset_factor;
        bbox[3] *= bbox_offset_factor;
        bbox_pt.push_back(OpenMesh::Vec2d(bbox[0] + bbox[1], bbox[2] + bbox[3]));
        bbox_pt.push_back(OpenMesh::Vec2d(bbox[0] - bbox[1], bbox[2] + bbox[3]));
        bbox_pt.push_back(OpenMesh::Vec2d(bbox[0] - bbox[1], bbox[2] - bbox[3]));
        bbox_pt.push_back(OpenMesh::Vec2d(bbox[0] + bbox[1], bbox[2] - bbox[3]));
        for (int i = 0; i < 4; i++)
            add_point_constraint(bbox_pt[i], in_pt_map, in_points, p_id);
        for (int i = 0; i < 4; i++) {
            m_features.push_back(std::make_pair(bbox_pt[i], bbox_pt[(i + 1) % 4]));
            add_segment_constraint(in_pt_map[bbox_pt[i]], in_pt_map[bbox_pt[(i + 1) % 4]], seg_id, in_segments, in_segment_markers);
        }
    }
    triangulateio t_in = {}, t_out = {};
    t_in.numberofpoints = p_id;
    if (t_in.numberofpoints > 0)
        t_in.pointlist = in_points.data();
    t_in.numberofsegments = in_segments.size() / 2;
    if (t_in.numberofsegments > 0) {
        t_in.segmentlist = in_segments.data();
        t_in.segmentmarkerlist = in_segment_markers.data();
    }
    t_out = {};
    std::string flags_str = "zpQeq";
    flags_str = flags_str + std::to_string(theta_triangulate * 180 / M_PI + 0.001);
    char flags[32];
    strcpy(flags, flags_str.c_str());

    triangulate(flags, &t_in, &t_out, nullptr);

    for (int i = 0; i < t_out.numberofpoints; i++)
        mesh->add_vertex(TriMesh::Point(t_out.pointlist[2 * i], t_out.pointlist[2 * i + 1], 0));
    for (int i = 0; i < t_out.numberoftriangles; i++)
        mesh->add_face(TriMesh::VertexHandle(t_out.trianglelist[3 * i]),
                       TriMesh::VertexHandle(t_out.trianglelist[3 * i + 1]),
                       TriMesh::VertexHandle(t_out.trianglelist[3 * i + 2]));
    for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it) {
        mesh->property(is_feature_node, *v_it) = false;
        mesh->property(is_on_feature, *v_it) = false;
    }
    for (auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it) {
        mesh->property(feature_id, *e_it) = -1;
    }
    std::map<std::pair<TriMesh::VertexHandle, TriMesh::VertexHandle>, int> feature_map;
    for (int i = 0; i < t_out.numberofedges; i++) {
        if (t_out.edgemarkerlist[i] > 0) {
            auto va = TriMesh::VertexHandle(t_out.edgelist[2 * i]);
            auto vb = TriMesh::VertexHandle(t_out.edgelist[2 * i + 1]);
            mesh->property(is_on_feature, va) = true;
            mesh->property(is_on_feature, vb) = true;
            auto f1 = m_features[t_out.edgemarkerlist[i] - 1].first;
            auto f2 = m_features[t_out.edgemarkerlist[i] - 1].second;
            auto fp1 = TriMesh::Point(f1[0], f1[1], 0);
            auto fp2 = TriMesh::Point(f2[0], f2[1], 0);
            if (mesh->point(va) == fp1 || mesh->point(va) == fp2)
                mesh->property(is_feature_node, va) = true;
            if (mesh->point(vb) == fp1 || mesh->point(vb) == fp2)
                mesh->property(is_feature_node, vb) = true;
            OpenMesh::EdgeHandle e = mesh->edge_handle(mesh->find_halfedge(va, vb));
            mesh->status(e).set_feature(true);
            mesh->property(feature_id, e) = t_out.edgemarkerlist[i] - 1;
        }
    }
}

void AngleBounded2DSimplificationPlugin::triangulate_random() {
    bb = true;
    if (OpenFlipper::Options::gui()) {
        get_parameters();
    }
    std::cout << "---- Triangulate ----" << std::endl;
    std::cout << "Seed: " << seed << std::endl;
    std::cout << "Number of constraints: " << n_constraints << std::endl;
    std::cout << "Angle: " << theta_triangulate << " = " << theta_triangulate * 180 / M_PI << " degrees" << std::endl;
    get_constraints();
    BaseObjectData* base;
    emit addEmptyObject(DATA_TRIANGLE_MESH, id);
    if (!PluginFunctions::getObject(id, base)) {
        std::cerr << "Couldn't get new mesh!\n";
        return;
    }
    base->setName(QString::number(seed) + QString("-") + QString::number(n_constraints) + QString("-") + QString::number(theta_triangulate));
    base->target(false);
    base->setObjectDrawMode(ACG::SceneGraph::DrawModes::WIREFRAME, false);
    mesh = PluginFunctions::triMesh(base);
    mesh->add_property(is_feature_node, "is_feature_node");
    mesh->add_property(is_on_feature, "is_on_feature");
    mesh->add_property(feature_id, "feature_id");
    delaunayMeshing();
    mesh->garbage_collection();
    mesh->update_normals();
    if (OpenFlipper::Options::gui()) {
        emit updatedObject(id, UPDATE_ALL);
    }
    std::cout << "Number of triangles: " << mesh->n_faces() << std::endl;
}

void AngleBounded2DSimplificationPlugin::triangulate_file(std::string filename) {
    bb = false;
    if (OpenFlipper::Options::gui()) {
        get_parameters();
    }
    std::cout << "---- Triangulate ----" << std::endl;
    std::cout << filename << std::endl;
    if (id != -1) {
        emit deleteObject(id);
    }
    std::ifstream file(filename);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    std::vector<std::vector<double>> vertices;
    for (auto line : lines) {
        size_t pos = line.find(" ");
        std::string type = line.substr(0, pos);
        line.erase(0, pos + 1);
        if (type == "v") {
            std::vector<double> vertex(2);
            pos = line.find(" ");
            vertex[0] = std::stod(line.substr(0, pos));
            line.erase(0, pos + 1);
            vertex[1] = std::stod(line);
            vertices.push_back(vertex);
        }
    }
    std::vector<std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d>> features;
    for (auto line : lines) {
        size_t pos = line.find(" ");
        std::string type = line.substr(0, pos);
        line.erase(0, pos + 1);
        if (type == "l") {
            std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d> l;
            pos = line.find(" ");
            l.first = OpenMesh::Vec2d(vertices[std::stoi(line.substr(0, pos)) - 1].data());
            line.erase(0, pos + 1);
            l.second = OpenMesh::Vec2d(vertices[std::stoi(line) - 1].data());
            features.push_back(l);
        }
    }
    m_features = features;

    std::list<std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d>> data;
    unsigned int old_count = 0;
    while (old_count != m_features.size()) {
        data.clear();
        data.insert(data.begin(), m_features.begin(), m_features.end());
        m_features.clear();
        old_count = data.size();
        intersection_test(data);
    }

    BaseObjectData* base;
    emit addEmptyObject(DATA_TRIANGLE_MESH, id);
    if (!PluginFunctions::getObject(id, base)) {
        std::cerr << "Couldn't get new mesh!\n";
        return;
    }
    base->target(false);
    base->setObjectDrawMode(ACG::SceneGraph::DrawModes::WIREFRAME, false);
    mesh = PluginFunctions::triMesh(base);
    mesh->add_property(is_feature_node, "is_feature_node");
    mesh->add_property(is_on_feature, "is_on_feature");
    mesh->add_property(feature_id, "feature_id");
    delaunayMeshing();
    mesh->garbage_collection();
    mesh->update_normals();
    if (OpenFlipper::Options::gui()) {
        emit updatedObject(id, UPDATE_ALL);
    }
    std::cout << "Number of triangles: " << mesh->n_faces() << std::endl;
}

// --------------------------------------------------------------------------------------------------------------------------------

void print_vector(std::vector<double> v) {
    std::cout << "(" << v[0] << ", " << v[1] << ")" << std::endl;
}

std::vector<double> add(std::vector<double> a, std::vector<double> b) {
    return {a[0] + b[0], a[1] + b[1]};
}

std::vector<double> sub(std::vector<double> a, std::vector<double> b) {
    return {a[0] - b[0], a[1] - b[1]};
}

std::vector<double> mul(std::vector<double> v, double s) {
    return {v[0] * s, v[1] * s};
}

std::vector<double> div(std::vector<double> v, double s) {
    return {v[0] / s, v[1] / s};
}

double dot(std::vector<double> a, std::vector<double> b) {
    return a[0] * b[0] + a[1] * b[1];
}

double norm(std::vector<double> v) {
    return std::sqrt(dot(v, v));
}

std::vector<double> normalize(std::vector<double> v) {
    return div(v, norm(v));
}

std::vector<double> rotate(std::vector<double> v, double angle) {
    double sin = std::sin(angle);
    double cos = std::cos(angle);
    return {cos * v[0] - sin * v[1], sin * v[0] + cos * v[1]};
}

double angle(std::vector<double> v) {
    return std::atan2(v[1], v[0]);
}

struct Line {
    std::vector<double> p;
    std::vector<double> v;
};

struct Circle {
    std::vector<double> c;
    double r;
};

std::vector<std::vector<double>> line_line_intersection(Line a, Line b) {
    if (b.v[0] == 0) {
        double tmp = a.p[0];
        a.p[0] = a.p[1];
        a.p[1] = tmp;
        tmp = a.v[0];
        a.v[0] = a.v[1];
        a.v[1] = tmp;
        tmp = b.p[0];
        b.p[0] = b.p[1];
        b.p[1] = tmp;
        tmp = b.v[0];
        b.v[0] = b.v[1];
        b.v[1] = tmp;
    }
    if (b.v[0] == 0) {
        // b.v = (0, 0)
        return {};
    }
    double gradient = b.v[1] / b.v[0];
    double denominator = a.v[1] - gradient * a.v[0];
    if (denominator == 0) {
        // Parallel or a.v = (0, 0)
        return {};
    }
    double lambda = (b.p[1] - a.p[1] - gradient * (b.p[0] - a.p[0])) / denominator;
    return {add(a.p, mul(a.v, lambda))};
}

std::vector<std::vector<double>> line_circle_intersection(Line l, Circle c) {
    std::vector<double> a = sub(l.p, c.c);
    double av = dot(a, l.v);
    double radicant = av * av - dot(a, a) + c.r * c.r;
    if (radicant < 0) {
        return {};
    }
    if (radicant == 0) {
        return {add(l.p, mul(l.v, -av))};
    }
    return {add(l.p, mul(l.v, -av - std::sqrt(radicant))), add(l.p, mul(l.v, -av + std::sqrt(radicant)))};
}

std::vector<std::vector<double>> circle_circle_intersection(Circle a, Circle b) {
    std::vector<double> ab = mul(sub(b.c, a.c), 2);
    if (ab[0] == 0 && ab[1] == 0) {
        // Same center
        return {};
    }
    double c = dot(b.c, b.c) - dot(a.c, a.c) - b.r * b.r + a.r * a.r;
    std::vector<double> p;
    if (ab[1] == 0) {
        p = {c / ab[0], 0};
    } else {
        p = {0, c / ab[1]};
    }
    std::vector<double> v = normalize(rotate(ab, M_PI / 2));
    return line_circle_intersection({p, v}, a);
}

std::pair<TriMesh::Point, double> AngleBounded2DSimplificationPlugin::pos(std::vector<OpenMesh::HalfedgeHandle> ring, std::set<TriMesh::VertexHandle>& affected_vertices, std::vector<double> center) {
    std::vector<std::vector<double>> angles_before;
    std::vector<Line> lines;
    std::vector<Circle> circles;
    for (auto h : ring) {
        OpenMesh::Vec3d pos_a = mesh->point(mesh->from_vertex_handle(h));
        OpenMesh::Vec3d pos_b = mesh->point(mesh->to_vertex_handle(h));
        OpenMesh::Vec3d pos_c = mesh->point(mesh->to_vertex_handle(mesh->next_halfedge_handle(h)));
        std::vector<std::vector<double>> tri = {{pos_a[0], pos_a[1]}, {pos_b[0], pos_b[1]}, {pos_c[0], pos_c[1]}};
        std::vector<double> angles;
        for (size_t j = 0; j < 3; j++) {
            std::vector<double> u = normalize(sub(tri[(j + 1) % 3], tri[j]));
            std::vector<double> v = normalize(sub(tri[(j + 2) % 3], tri[j]));
            angles.push_back(std::acos(dot(u, v)));
        }
        angles_before.push_back(angles);

        std::vector<double> ab = sub(tri[1], tri[0]);
        double d = norm(ab);
        ab = div(ab, d);

        double alpha = std::min(theta, angles[0]) + 0.001;
        std::vector<double> va = rotate(ab, alpha);
        lines.push_back({tri[0], va});

        double beta = std::min(theta, angles[1]) + 0.001;
        std::vector<double> vb = rotate(ab, -beta);
        lines.push_back({tri[1], vb});

        double gamma = std::min(theta, angles[2]) + 0.001;
        std::vector<double> c = add(div(add(tri[0], tri[1]), 2), mul(rotate(ab, M_PI / 2), d / 2 / std::tan(gamma)));
        double r = d / 2 / std::sin(gamma);
        circles.push_back({c, r});
    }

    std::vector<double> p = {0, 0};

    if (strategy == 0) {

        std::vector<std::pair<std::vector<double>, std::pair<int, int>>> intersections;
        for (size_t i = 0; i < lines.size(); i++) {
            for (size_t j = i + 1; j < lines.size(); j++) {
                std::vector<std::vector<double>> ll = line_line_intersection(lines[i], lines[j]);
                for (auto s : ll) {
                    intersections.push_back({s, {i, j}});
                }
            }
            for (size_t j = 0; j < circles.size(); j++) {
                std::vector<std::vector<double>> lc = line_circle_intersection(lines[i], circles[j]);
                for (auto s : lc) {
                    intersections.push_back({s, {i, j + lines.size()}});
                }
            }
        }
        for (size_t i = 0; i < circles.size(); i++) {
            for (size_t j = i + 1; j < circles.size(); j++) {
                std::vector<std::vector<double>> cc = circle_circle_intersection(circles[i], circles[j]);
                for (auto s : cc) {
                    intersections.push_back({s, {i + lines.size(), j + lines.size()}});
                }
            }
        }

        for (auto line : lines) {
            std::vector<double> n = rotate(line.v, M_PI / 2);
            for (size_t i = 0; i < intersections.size();) {
                std::vector<double> s = intersections[i].first;
                if (dot(n, sub(s, line.p)) < -0.001) {
                    intersections.erase(intersections.begin() + i);
                } else {
                    i++;
                }
            }
        }
        for (auto circle : circles) {
            for (size_t i = 0; i < intersections.size();) {
                std::vector<double> s = intersections[i].first;
                if (norm(sub(s, circle.c)) > circle.r + 0.001) {
                    intersections.erase(intersections.begin() + i);
                } else {
                    i++;
                }
            }
        }

        if (intersections.size() == 0) {
            return {{0, 0, 0}, 0};
        }

        for (auto intersection : intersections) {
            p = add(p, intersection.first);
        }
        p = div(p, intersections.size());

    } else if (strategy == 1 || strategy == 2) {

#ifdef GUROBI
        try {
            env.set(GRB_IntParam_OutputFlag, 0);
            GRBModel model = GRBModel(env);
            GRBVar x = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x");
            GRBVar y = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y");
            GRBVar dummy = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dummy");
            GRBLinExpr obj = -dummy;
            model.set(GRB_IntParam_NonConvex, 2);
            model.setObjective(obj, GRB_MINIMIZE);
            for (auto l : lines) {
                double a = -l.v[1];
                double b = l.v[0];
                double c = l.v[1] * l.p[0] - l.v[0] * l.p[1];
                if (strategy == 1) {
                    model.addConstr(a * x + b * y + c - dummy >= 0);
                } else {
                    model.addQConstr((a * x + b * y + c) * (a * x + b * y + c) - dummy >= 0);
                    model.addConstr(a * x + b * y + c >= 0);
                }
            }
            for (auto c : circles) {
                model.addQConstr(c.r * c.r - c.c[0] * c.c[0] - c.c[1] * c.c[1] + 2 * c.c[0] * x + 2 * c.c[1] * y - x * x - y * y - dummy >= 0);
            }
            model.optimize();
            if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
                if (dummy.get(GRB_DoubleAttr_X) < 0) {
                    return {{0, 0, 0}, 0};
                }
                p = {x.get(GRB_DoubleAttr_X), y.get(GRB_DoubleAttr_X)};
            } else {
                return {{0, 0, 0}, 0};
            }
        } catch (GRBException e) {
            std::cout << "Gurobi exception " << e.getErrorCode() << " " << e.getMessage() << std::endl;
        }
#endif

    } else if (strategy == 3) {

        auto smooth_min = [this, ring](std::vector<double> c) {
            size_t n = 3 * ring.size();
            std::vector<double> phi(n);
            std::vector<std::vector<double>> J_phi(n);
            for (size_t i = 0; i < ring.size(); i++) {
                OpenMesh::Vec3d pos_a = mesh->point(mesh->from_vertex_handle(ring[i]));
                std::vector<double> a = {pos_a[0], pos_a[1]};
                OpenMesh::Vec3d pos_b = mesh->point(mesh->to_vertex_handle(ring[i]));
                std::vector<double> b = {pos_b[0], pos_b[1]};
                double ab_dot_ab = dot(sub(b, a), sub(b, a));
                double norm_ab = std::sqrt(ab_dot_ab);
                double ac_dot_ac = dot(sub(c, a), sub(c, a));
                double norm_ac = std::sqrt(ac_dot_ac);
                double bc_dot_bc = dot(sub(c, b), sub(c, b));
                double norm_bc = std::sqrt(bc_dot_bc);
                double ab_dot_ac = dot(sub(b, a), sub(c, a));
                double bc_dot_ba = dot(sub(c, b), sub(a, b));
                double ca_dot_cb = dot(sub(a, c), sub(b, c));
                phi[3 * i] = std::acos(ab_dot_ac / (norm_ab * norm_ac));
                phi[3 * i + 1] = std::acos(bc_dot_ba / (norm_bc * norm_ab));
                phi[3 * i + 2] = std::acos(ca_dot_cb / (norm_ac * norm_bc));
                J_phi[3 * i] = {
                    -((b[0] - a[0]) / (norm_ab * norm_ac) - (c[0] - a[0]) * ab_dot_ac / (norm_ab * std::pow(ac_dot_ac, 3 / 2))) / std::sqrt(1 - ab_dot_ac * ab_dot_ac / (ab_dot_ab * ac_dot_ac)),
                    -((b[1] - a[1]) / (norm_ab * norm_ac) - (c[1] - a[1]) * ab_dot_ac / (norm_ab * std::pow(ac_dot_ac, 3 / 2))) / std::sqrt(1 - ab_dot_ac * ab_dot_ac / (ab_dot_ab * ac_dot_ac))};
                J_phi[3 * i + 1] = {
                    -((a[0] - b[0]) / (norm_ab * norm_bc) - (c[0] - b[0]) * bc_dot_ba / (norm_ab * std::pow(bc_dot_bc, 3 / 2))) / std::sqrt(1 - bc_dot_ba * bc_dot_ba / (ab_dot_ab * bc_dot_bc)),
                    -((a[1] - b[1]) / (norm_ab * norm_bc) - (c[1] - b[1]) * bc_dot_ba / (norm_ab * std::pow(bc_dot_bc, 3 / 2))) / std::sqrt(1 - bc_dot_ba * bc_dot_ba / (ab_dot_ab * bc_dot_bc))};
                J_phi[3 * i + 2] = {
                    -((-a[0] - b[0] + 2 * c[0]) / (norm_ac * norm_bc) + (a[0] - c[0]) * ca_dot_cb / (std::pow(ac_dot_ac, 3 / 2) * norm_bc) + (b[0] - c[0]) * ca_dot_cb / (norm_ac * std::pow(bc_dot_bc, 3 / 2))) / std::sqrt(1 - ca_dot_cb * ca_dot_cb / (ac_dot_ac * bc_dot_bc)),
                    -((-a[1] - b[1] + 2 * c[1]) / (norm_ac * norm_bc) + (a[1] - c[1]) * ca_dot_cb / (std::pow(ac_dot_ac, 3 / 2) * norm_bc) + (b[1] - c[1]) * ca_dot_cb / (norm_ac * std::pow(bc_dot_bc, 3 / 2))) / std::sqrt(1 - ca_dot_cb * ca_dot_cb / (ac_dot_ac * bc_dot_bc))};
            }
            double ALPHA = -100;
            std::vector<double> exp(n);
            double sum = 0;
            for (size_t i = 0; i < n; i++) {
                exp[i] = std::exp(ALPHA * phi[i]);
                sum += exp[i];
            }
            double S = 1 / ALPHA * std::log(sum);
            std::vector<double> J_S(n);
            for (size_t i = 0; i < n; i++) {
                J_S[i] = exp[i] / sum;
            }
            std::vector<double> gradient(2);
            for (size_t i = 0; i < 2; i++) {
                gradient[i] = 0;
                for (size_t j = 0; j < n; j++) {
                    gradient[i] += J_S[j] * J_phi[j][i];
                }
            }
            std::pair<double, std::vector<double>> result = {S, gradient};
            return result;
        };

        double scale = std::numeric_limits<double>::max();
        for (auto h : ring) {
            OpenMesh::Vec3d pos_a = mesh->point(mesh->from_vertex_handle(h));
            OpenMesh::Vec3d pos_b = mesh->point(mesh->to_vertex_handle(h));
            double d = norm(sub({pos_b[0], pos_b[1]}, {pos_a[0], pos_a[1]}));
            if (d < scale) {
                scale = d;
            }
        }
        p = center;
        std::pair<double, std::vector<double>> f = smooth_min(p);
        double step = scale / 2;
        for (size_t i = 0; i < 100 && step > scale / 1000; i++) {
            std::vector<double> next_p = add(p, mul(normalize(f.second), step));
            std::pair<double, std::vector<double>> next_f = smooth_min(next_p);
            if (next_f.first >= f.first) {
                p = next_p;
                f = next_f;
                step *= 2;
            } else {
                step /= 2;
            }
        }

    } else {
        p = center;
    }

    double min_angle = M_PI;
    for (size_t i = 0; i < ring.size(); i++) {
        OpenMesh::Vec3d pos_a = mesh->point(mesh->from_vertex_handle(ring[i]));
        OpenMesh::Vec3d pos_b = mesh->point(mesh->to_vertex_handle(ring[i]));
        std::vector<std::vector<double>> tri = {{pos_a[0], pos_a[1]}, {pos_b[0], pos_b[1]}, p};
        std::vector<double> ab = sub(tri[1], tri[0]);
        std::vector<double> ac = sub(tri[2], tri[0]);
        double area = ab[0] * ac[1] - ac[0] * ab[1];
        if (area < 0) {
            return {{0, 0, 0}, 0};
        }
        for (size_t j = 0; j < 3; j++) {
            std::vector<double> u = normalize(sub(tri[(j + 1) % 3], tri[j]));
            std::vector<double> v = normalize(sub(tri[(j + 2) % 3], tri[j]));
            double angle = std::acos(dot(u, v));
            if (angle < std::min(theta, angles_before[i][j])) {
                return {{0, 0, 0}, 0};
            }
            if (angle < min_angle) {
                min_angle = angle;
            }
        }
    }

    for (auto h : ring) {
        affected_vertices.insert(mesh->to_vertex_handle(h));
    }
    return {{p[0], p[1], 0}, min_angle};
}

bool AngleBounded2DSimplificationPlugin::triangle_collapse_ok(OpenMesh::FaceHandle f) {
    OpenMesh::HalfedgeHandle h0 = mesh->halfedge_handle(f);
    OpenMesh::HalfedgeHandle h1 = mesh->next_halfedge_handle(h0);
    return mesh->is_collapse_ok(h0) && mesh->is_collapse_ok(h1);
}

OpenMesh::VertexHandle AngleBounded2DSimplificationPlugin::collapse_triangle(OpenMesh::FaceHandle f, OpenMesh::Vec3d p) {
    OpenMesh::HalfedgeHandle h0 = mesh->halfedge_handle(f);
    OpenMesh::HalfedgeHandle h1 = mesh->next_halfedge_handle(h0);
    OpenMesh::VertexHandle v = mesh->to_vertex_handle(h1);
    mesh->collapse(h0);
    mesh->collapse(h1);
    mesh->set_point(v, p);
    return v;
}

OpenMesh::VertexHandle AngleBounded2DSimplificationPlugin::collapse_edge(OpenMesh::EdgeHandle e, OpenMesh::Vec3d p) {
    OpenMesh::HalfedgeHandle h = mesh->halfedge_handle(e, 0);
    OpenMesh::VertexHandle v = mesh->to_vertex_handle(h);
    mesh->collapse(h);
    mesh->set_point(v, p);
    return v;
}

void AngleBounded2DSimplificationPlugin::start_thread() {
    if (algo_thread.joinable()) {
        algo_thread.join();
    }
    algo_thread = std::thread(&AngleBounded2DSimplificationPlugin::decimate, this);
}

struct Collapse {
    enum {
        TRIANGLE,
        EDGE,
        HALFEDGE
    } type;
    OpenMesh::FaceHandle face;
    OpenMesh::EdgeHandle edge;
    OpenMesh::HalfedgeHandle halfedge;
    bool collapsable;
    OpenMesh::Vec3d p;
    double min_angle;
    std::set<TriMesh::VertexHandle> affected_vertices;
};

void AngleBounded2DSimplificationPlugin::decimate() {
    if (OpenFlipper::Options::gui()) {
        get_parameters();
    }
    std::cout << "---- Decimate    ----" << std::endl;
    std::cout << "Angle: " << theta << " = " << theta * 180 / M_PI << " degrees" << std::endl;
    std::cout << "Prioritize: " << prio << std::endl;
    std::cout << "Use halfedge collapses: " << use_halfedge_collapses << std::endl;
    std::cout << "Use edge collapses: " << use_edge_collapses << std::endl;
    std::cout << "Use triangle collapses: " << use_triangle_collapses << std::endl;
    std::cout << "Prioritize triangle collapses: " << prio_triangle_collapses << std::endl;
    std::cout << "Strategy: " << strategy << std::endl;
    std::cout << "Number of triangles before: " << mesh->n_faces() << std::endl;
    std::cout << "..." << std::endl;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();

    auto get_triangle_collapse = [this](OpenMesh::FaceHandle f) {
        Collapse c;
        c.type = Collapse::TRIANGLE;
        c.face = f;
        c.collapsable = false;
        bool free = mesh->fv_iter(f).is_valid() && triangle_collapse_ok(f);
        for (auto fv_it = mesh->fv_iter(f); fv_it.is_valid() && free; fv_it++) {
            if (mesh->property(is_on_feature, *fv_it)) {
                free = false;
            }
        }
        if (!free) {
            return c;
        }
        std::vector<OpenMesh::HalfedgeHandle> ring;
        std::vector<OpenMesh::VertexHandle> vertices;
        for (auto fv_it = mesh->fv_iter(f); fv_it.is_valid(); fv_it++) {
            vertices.push_back(*fv_it);
        }
        OpenMesh::HalfedgeHandle h = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(mesh->halfedge_handle(f))))));
        OpenMesh::VertexHandle to = mesh->to_vertex_handle(h);
        if (to == vertices[0] || to == vertices[1] || to == vertices[2]) {
            h = mesh->next_halfedge_handle(mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(mesh->halfedge_handle(f))))));
        }
        OpenMesh::HalfedgeHandle s = h;
        do {
            ring.push_back(h);
            h = mesh->next_halfedge_handle(h);
            do {
                h = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(h));
                to = mesh->to_vertex_handle(h);
            } while (to == vertices[0] || to == vertices[1] || to == vertices[2]);
        } while (h != s);
        std::set<OpenMesh::VertexHandle> affected_vertices;
        std::vector<double> center = {0, 0};
        for (size_t i = 0; i < 3; i++) {
            OpenMesh::Vec3d pos_v = mesh->point(vertices[i]);
            center = add(center, {pos_v[0], pos_v[1]});
        }
        center = div(center, 3);
        std::pair<OpenMesh::Vec3d, double> pa = pos(ring, affected_vertices, center);
        if (affected_vertices.size() > 0) {
            c.collapsable = true;
            c.p = pa.first;
            c.min_angle = pa.second;
            c.affected_vertices = affected_vertices;
        }
        return c;
    };

    auto get_edge_collapse = [this](OpenMesh::EdgeHandle e) {
        Collapse c;
        c.type = Collapse::EDGE;
        c.edge = e;
        c.collapsable = false;
        OpenMesh::HalfedgeHandle h = mesh->halfedge_handle(e, 0);
        OpenMesh::VertexHandle u = mesh->from_vertex_handle(h);
        OpenMesh::VertexHandle v = mesh->to_vertex_handle(h);
        if (!mesh->is_collapse_ok(h) || mesh->property(is_on_feature, u) || mesh->property(is_on_feature, v)) {
            return c;
        }
        std::vector<OpenMesh::HalfedgeHandle> ring;
        h = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(mesh->next_halfedge_handle(h))));
        OpenMesh::HalfedgeHandle s = h;
        do {
            ring.push_back(h);
            h = mesh->next_halfedge_handle(h);
            OpenMesh::VertexHandle to;
            do {
                h = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(h));
                to = mesh->to_vertex_handle(h);
            } while (to == u || to == v);
        } while (h != s);
        std::set<OpenMesh::VertexHandle> affected_vertices;
        std::vector<double> center = {0, 0};
        OpenMesh::Vec3d pos_u = mesh->point(u);
        OpenMesh::Vec3d pos_v = mesh->point(v);
        center = div(add(add(center, {pos_u[0], pos_u[1]}), {pos_v[0], pos_v[1]}), 2);
        std::pair<OpenMesh::Vec3d, double> pa = pos(ring, affected_vertices, center);
        if (affected_vertices.size() > 0) {
            c.collapsable = true;
            c.p = pa.first;
            c.min_angle = pa.second;
            c.affected_vertices = affected_vertices;
        }
        return c;
    };

    auto get_halfedge_collapse = [this](OpenMesh::HalfedgeHandle h) {
        Collapse c;
        c.type = Collapse::HALFEDGE;
        c.halfedge = h;
        c.collapsable = false;
        OpenMesh::EdgeHandle e = mesh->edge_handle(h);
        OpenMesh::VertexHandle u = mesh->from_vertex_handle(h);
        OpenMesh::VertexHandle v = mesh->to_vertex_handle(h);
        if (!mesh->is_collapse_ok(h)) {
            return c;
        }
        if (mesh->property(feature_id, e) >= 0) {
            if (mesh->property(is_feature_node, u)) {
                return c;
            }
        } else if (mesh->property(is_on_feature, u)) {
            return c;
        }
        std::vector<OpenMesh::HalfedgeHandle> back;
        for (auto uh_it = mesh->voh_iter(u); uh_it.is_valid(); uh_it++) {
            if (!mesh->is_boundary(*uh_it) && *uh_it != h && *uh_it != mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(h))) {
                back.push_back(mesh->next_halfedge_handle(*uh_it));
            }
        }
        std::set<OpenMesh::VertexHandle> affected_vertices;
        OpenMesh::Vec3d pos_v = mesh->point(v);
        std::vector<double> p = {pos_v[0], pos_v[1]};
        double min_angle = M_PI;
        for (auto h : back) {
            OpenMesh::HalfedgeHandle hn = mesh->next_halfedge_handle(h);
            OpenMesh::HalfedgeHandle hnn = mesh->next_halfedge_handle(hn);
            OpenMesh::Vec3d pos_a = mesh->point(mesh->from_vertex_handle(h));
            OpenMesh::Vec3d pos_b = mesh->point(mesh->to_vertex_handle(h));
            OpenMesh::Vec3d pos_c = mesh->point(mesh->to_vertex_handle(hn));
            std::vector<std::vector<double>> tri_before = {{pos_a[0], pos_a[1]}, {pos_b[0], pos_b[1]}, {pos_c[0], pos_c[1]}};
            std::vector<std::vector<double>> tri_after = {{pos_a[0], pos_a[1]}, {pos_b[0], pos_b[1]}, p};
            std::vector<double> ab = sub(tri_after[1], tri_after[0]);
            std::vector<double> ac = sub(tri_after[2], tri_after[0]);
            double area = ab[0] * ac[1] - ac[0] * ab[1];
            if (area < 0) {
                return c;
            }
            std::vector<bool> skip = {
                mesh->property(feature_id, mesh->edge_handle(hnn)) >= 0 && mesh->property(feature_id, mesh->edge_handle(h)) >= 0,
                mesh->property(feature_id, mesh->edge_handle(h)) >= 0 && mesh->property(feature_id, mesh->edge_handle(hn)) >= 0,
                mesh->property(feature_id, mesh->edge_handle(hn)) >= 0 && mesh->property(feature_id, mesh->edge_handle(hnn)) >= 0};
            for (size_t j = 0; j < 3; j++) {
                if (!skip[j]) {
                    std::vector<double> u = normalize(sub(tri_before[(j + 1) % 3], tri_before[j]));
                    std::vector<double> v = normalize(sub(tri_before[(j + 2) % 3], tri_before[j]));
                    double angle_before = std::acos(dot(u, v));
                    u = normalize(sub(tri_after[(j + 1) % 3], tri_after[j]));
                    v = normalize(sub(tri_after[(j + 2) % 3], tri_after[j]));
                    double angle_after = std::acos(dot(u, v));
                    if (angle_after >= M_PI || std::isnan(angle_after) || angle_after < std::min(theta, angle_before - 0.001)) {
                        return c;
                    }
                    if (angle_after < min_angle) {
                        min_angle = angle_after;
                    }
                }
            }
            affected_vertices.insert(mesh->from_vertex_handle(h));
            affected_vertices.insert(mesh->to_vertex_handle(h));
        }
        c.collapsable = true;
        c.p = pos_v;
        c.min_angle = min_angle;
        c.affected_vertices = affected_vertices;
        return c;
    };

    auto cmp = [this](Collapse a, Collapse b) {
        double weight_a = a.min_angle;
        if (prio_triangle_collapses && a.type == Collapse::TRIANGLE) {
            weight_a += M_PI;
        }
        double weight_b = b.min_angle;
        if (prio_triangle_collapses && b.type == Collapse::TRIANGLE) {
            weight_b += M_PI;
        }
        return weight_a < weight_b;
    };
    PriorityQueue<Collapse, decltype(cmp)> oq(cmp, mesh->n_faces() + mesh->n_edges() + mesh->n_halfedges());
    std::queue<Collapse> q;
    std::vector<bool> inq(mesh->n_faces() + mesh->n_edges() + mesh->n_halfedges(), true);

    std::vector<Collapse> collapses;
    if (use_triangle_collapses) {
        for (auto f : mesh->faces()) {
            if (prio) {
                Collapse c = get_triangle_collapse(f);
                oq.push(c, f.idx());
            } else {
                Collapse c;
                c.type = Collapse::TRIANGLE;
                c.face = f;
                collapses.push_back(c);
            }
        }
    }
    if (use_edge_collapses) {
        for (auto e : mesh->edges()) {
            if (prio) {
                Collapse c = get_edge_collapse(e);
                oq.push(c, mesh->n_faces() + e.idx());
            } else {
                Collapse c;
                c.type = Collapse::EDGE;
                c.edge = e;
                collapses.push_back(c);
            }
        }
    }
    if (use_halfedge_collapses) {
        for (auto h : mesh->halfedges()) {
            if (prio) {
                Collapse c = get_halfedge_collapse(h);
                oq.push(c, mesh->n_faces() + mesh->n_edges() + h.idx());
            } else {
                Collapse c;
                c.type = Collapse::HALFEDGE;
                c.halfedge = h;
                collapses.push_back(c);
            }
        }
    }
    if (!prio) {
        std::default_random_engine rng;
        rng.seed(6364);
        std::shuffle(collapses.begin(), collapses.end(), rng);
        for (auto c : collapses) {
            q.push(c);
        }
    }

    while (!(oq.empty() && q.empty())) {
        Collapse c;
        if (prio) {
            c = oq.pop();
        } else {
            c = q.front();
            q.pop();
            switch (c.type) {
            case Collapse::TRIANGLE: {
                inq[c.face.idx()] = false;
                c = get_triangle_collapse(c.face);
            } break;
            case Collapse::EDGE: {
                inq[mesh->n_faces() + c.edge.idx()] = false;
                c = get_edge_collapse(c.edge);
            } break;
            case Collapse::HALFEDGE: {
                inq[mesh->n_faces() + mesh->n_edges() + c.halfedge.idx()] = false;
                c = get_halfedge_collapse(c.halfedge);
            } break;
            }
        }

        if (c.collapsable) {
            bool collapsed = false;
            switch (c.type) {
            case Collapse::TRIANGLE: {
                if (mesh->fv_iter(c.face).is_valid() && triangle_collapse_ok(c.face)) {
                    c.affected_vertices.insert(collapse_triangle(c.face, c.p));
                    collapsed = true;
                }
            } break;
            case Collapse::EDGE: {
                if (mesh->is_collapse_ok(mesh->halfedge_handle(c.edge, 0))) {
                    c.affected_vertices.insert(collapse_edge(c.edge, c.p));
                    collapsed = true;
                }
            } break;
            case Collapse::HALFEDGE: {
                if (mesh->is_collapse_ok(c.halfedge)) {
                    c.affected_vertices.insert(mesh->to_vertex_handle(c.halfedge));
                    mesh->collapse(c.halfedge);
                    collapsed = true;
                }
            } break;
            }

            if (collapsed) {
                if (use_triangle_collapses) {
                    std::set<OpenMesh::FaceHandle> affected_faces;
                    for (auto v : c.affected_vertices) {
                        for (auto vf_it = mesh->vf_iter(v); vf_it.is_valid(); vf_it++) {
                            affected_faces.insert(*vf_it);
                        }
                    }
                    for (auto f : affected_faces) {
                        if (prio) {
                            Collapse c = get_triangle_collapse(f);
                            oq.push(c, f.idx());
                        } else if (!inq[f.idx()]) {
                            Collapse c;
                            c.type = Collapse::TRIANGLE;
                            c.face = f;
                            q.push(c);
                            inq[f.idx()] = true;
                        }
                    }
                }
                if (use_edge_collapses) {
                    std::set<OpenMesh::EdgeHandle> affected_edges;
                    for (auto v : c.affected_vertices) {
                        for (auto ve_it = mesh->ve_iter(v); ve_it.is_valid(); ve_it++) {
                            affected_edges.insert(*ve_it);
                        }
                    }
                    for (auto e : affected_edges) {
                        if (prio) {
                            Collapse c = get_edge_collapse(e);
                            oq.push(c, mesh->n_faces() + e.idx());
                        } else if (!inq[mesh->n_faces() + e.idx()]) {
                            Collapse c;
                            c.type = Collapse::EDGE;
                            c.edge = e;
                            q.push(c);
                            inq[mesh->n_faces() + e.idx()] = true;
                        }
                    }
                }
                if (use_halfedge_collapses) {
                    std::set<OpenMesh::HalfedgeHandle> affected_halfedges;
                    for (auto v : c.affected_vertices) {
                        for (auto vh_it = mesh->voh_iter(v); vh_it.is_valid(); vh_it++) {
                            affected_halfedges.insert(*vh_it);
                        }
                    }
                    for (auto h : affected_halfedges) {
                        if (prio) {
                            Collapse c = get_halfedge_collapse(h);
                            oq.push(c, mesh->n_faces() + mesh->n_edges() + h.idx());
                        } else if (!inq[mesh->n_faces() + mesh->n_edges() + h.idx()]) {
                            Collapse c;
                            c.type = Collapse::HALFEDGE;
                            c.halfedge = h;
                            q.push(c);
                            inq[mesh->n_faces() + mesh->n_edges() + h.idx()] = true;
                        }
                    }
                }
            }
        }
    }
    mesh->garbage_collection();
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    if (OpenFlipper::Options::gui()) {
        emit updatedObject(id, UPDATE_ALL);
    }
    std::cout << "Number of triangles after: " << mesh->n_faces() << std::endl;
    std::cout << "Time: " << runtime << "ms" << std::endl;
}

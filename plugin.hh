#pragma once

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/common/Types.hh>
#include <QCheckBox>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>
#include <QSpinBox>

#include <algorithm>
#include <clocale>
#include <queue>
#include <random>
#include <thread>
#include <vector>

#include "priority_queue.hh"
extern "C" {
#include "triangle.h"
}
#include <gurobi_c++.h>

class AngleBounded2DSimplificationPlugin : public QObject, BaseInterface, LoadSaveInterface, ToolboxInterface {
    Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(LoadSaveInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_PLUGIN_METADATA(IID "org.OpenFlipper.Plugins.AngleBounded2DSimplification")

signals:
    void updatedObject(int, const UpdateType&);
    void addToolbox(QString, QWidget*);
    void addEmptyObject(DataType, int&);
    void deleteObject(int);

public:
    QString name() {
        return (QString("Angle-Bounded 2D Mesh Simplification"));
    }
    QString description() {
        return (QString("Implementation for the paper submitted to GMP 2022"));
    }

private slots:
    void noguiSupported() {}
    void initializePlugin();

    void triangulate_random();
    void start_thread();

private:
    int id;
    TriMesh* mesh;
    OpenMesh::VPropHandleT<bool> is_feature_node;
    OpenMesh::VPropHandleT<bool> is_on_feature;
    OpenMesh::EPropHandleT<int> feature_id;

    QSpinBox* in_seed;
    QSpinBox* in_n_constraints;
    QDoubleSpinBox* in_theta_triangulate;
    QDoubleSpinBox* in_theta;
    QCheckBox* in_prio;
    QCheckBox* in_use_halfedge_collapses;
    QCheckBox* in_use_edge_collapses;
    QCheckBox* in_use_triangle_collapses;
    QCheckBox* in_prio_triangle_collapses;
    QSpinBox* in_strategy;

    std::vector<std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d>> m_features;

    std::thread algo_thread;

    void get_parameters();

    void intersection_test(std::list<std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d>>& data);
    void get_constraints();
    void delaunayMeshing();

    std::pair<TriMesh::Point, double> pos(std::vector<OpenMesh::HalfedgeHandle> ring, std::set<TriMesh::VertexHandle>& neighbours, std::vector<double> center);
    bool triangle_collapse_ok(OpenMesh::FaceHandle f);
    OpenMesh::VertexHandle collapse_triangle(OpenMesh::FaceHandle f, OpenMesh::Vec3d p);
    OpenMesh::VertexHandle collapse_edge(TriMesh::EdgeHandle e, OpenMesh::Vec3d p);
    void decimate();
};

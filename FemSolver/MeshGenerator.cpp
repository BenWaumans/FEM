#include "MeshGenerator.hpp"
#include <GC_MakeArcOfCircle.hxx>
#include <GC_MakeSegment.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Pnt.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Trsf.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <TopoDS_Face.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <gp_Vec.hxx>
#include <cmath>
#include <TopoDS.hxx>
#include <gp.hxx>

static double CalcPointOnCircle(const gp_Pnt &point, const double &x, const double &radius)
{
    return sqrt(pow(radius, 2) - pow(x - point.X(), 2)) + point.Y();
}

static TopoDS_Wire MirrorWire(const TopoDS_Wire &wire, const gp_Ax1 &axis)
{
    gp_Trsf trsf;
    trsf.SetMirror(axis);
    BRepBuilderAPI_Transform bRepTrsf(wire, trsf);
    TopoDS_Shape mirroredShape = bRepTrsf.Shape();
    TopoDS_Wire mirroredWire = TopoDS::Wire(mirroredShape);
    BRepBuilderAPI_MakeWire mkWire;
    mkWire.Add(wire);
    mkWire.Add(mirroredWire);
    return mkWire.Wire();
}

MeshGenerator::MeshGenerator(const double &hight,
                             const double &width,
                             const double &tickness_s,
                             const double &tickness_g,
                             const double &radius,
                             const double &length)
{
    gp_Pnt pnt1(hight / 2, 0, 0);
    gp_Pnt pnt2(hight / 2, width / 2, 0);
    gp_Pnt pnt3(hight / 2 - tickness_g, width / 2, 0);
    gp_Pnt pnt4(hight / 2 - tickness_g, tickness_s + radius, 0);
    gp_Pnt pnt5(
        hight / 2 - tickness_g - radius / 2,
        CalcPointOnCircle(pnt4, hight / 2 - tickness_g - radius / 2, radius),
        0);
    gp_Pnt pnt6(hight / 2 - tickness_g - radius, tickness_s, 0);
    gp_Pnt pnt7(0, tickness_s, 0);

    Handle(Geom_TrimmedCurve) segments[5];
    segments[0] = GC_MakeSegment(pnt1, pnt2);
    segments[1] = GC_MakeSegment(pnt2, pnt3);
    segments[2] = GC_MakeSegment(pnt3, pnt4);
    segments[3] = GC_MakeArcOfCircle(pnt4, pnt5, pnt6);
    segments[4] = GC_MakeSegment(pnt6, pnt7);

    BRepBuilderAPI_MakeWire wire = BRepBuilderAPI_MakeWire();
    for(Handle(Geom_TrimmedCurve) segment : segments) {
        wire.Add(BRepBuilderAPI_MakeEdge(segment));
    }

    TopoDS_Wire wireProfile = MirrorWire(MirrorWire(wire, gp::OX()), gp::OY());

    TopoDS_Face faceProfile = BRepBuilderAPI_MakeFace(wireProfile);
    gp_Vec prismVec(0, 0, length);
    TopoDS_Shape body = BRepPrimAPI_MakePrism(faceProfile, prismVec);
}


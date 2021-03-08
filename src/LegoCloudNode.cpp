#include "LegoCloudNode.h"

#include <qmath.h>
#include <QGraphicsScene>
#include <iostream>
#include <fstream>
#include <utility>

#include "LegoDimensions.h"
#include "LegoCloud.h"
#include "LegoGraph.h"

#include <boost/optional.hpp>

#ifdef WIN32
#include <windows.h>
#include <gl/GLU.h>
#elif __APPLE__
#include <glu.h>
#else
#include <GL/glu.h>
#endif

#define KNOB_RESOLUTION_DISPLAY 15
#define KNOB_RESOLUTION_OBJ_EXPORT 8

template<typename Id>
struct VertexCache {
    QMap<Vector3, int> known;
    Id id;

    VertexCache(Id idI) : known(), id(idI) {}

    template<typename NewF>
    int getIndex(const Vector3& v, NewF newF) {
        const auto lb = known.lowerBound(v);
        if(lb == known.end() || !(lb.key() == v)) {
            int newIdx = 1 + known.size();
            known.insert(lb, v, newIdx);
            newF(id);
            return newIdx;
        } else {
            return *lb;
        }
    }

};

struct OBJVertStream {
    std::ostream& out;
    const char* spec;
    int degree;
};

struct LazyVertex {
    Vector3 loc;
    VertexCache<OBJVertStream>* cache;
    operator int() {
        return cache->getIndex(loc, [this](const OBJVertStream& id){
            id.out << "v" << id.spec;
            for(int i = 0; i < id.degree; ++i){
                id.out << " " << loc.data()[i];
            }
            id.out << std::endl;
        });
    }
};

template<size_t ...I> struct integer_sequence {};
template<size_t N, size_t ...I>
struct partial_index_sequence {
    using type = typename partial_index_sequence<N - 1, N, I...>::type;
};
template<size_t ...I>
struct partial_index_sequence<0, I...> {
    using type = integer_sequence<0, I...>;
};

template<size_t N> using index_sequence = typename partial_index_sequence<N - 1>::type;

template<size_t N>
std::array<int, N> intCollection(const std::array<int, N>& i) { return i; }

template<typename Collection, size_t ...I>
std::array<int, sizeof...(I)> intCollection(QPair<Collection, integer_sequence<I...>>& pair) {
    return {pair.first[I] ...};
}

template<typename Elem, typename std::enable_if<!std::is_same<Elem, int>::value, size_t>::type N>
std::array<int, N> intCollection(std::array<Elem, N>& e) {
    auto conversionPair = qMakePair(e, index_sequence<N>());
    return intCollection(conversionPair);
}

QVector<int> intCollection(const QVector<int>& v) { return v; }

QVector<int> intCollection(boost::none_t) { return QVector<int>(); }

template<typename F, size_t ... I>
auto unrollLoop1D(F f, integer_sequence<I...>) -> std::array<decltype(f((size_t)0)), sizeof...(I)> {
    return { f(I) ... };
}

template<typename F, size_t ...I, size_t ...J>
auto unrollLoop2D(F f, integer_sequence<I...> iS, integer_sequence<J...>) -> std::array<std::array<decltype(f((size_t)0, (size_t)0)), sizeof...(J)>, sizeof...(I)> {
    return unrollLoop1D([&f](size_t i) -> std::array<decltype(f((size_t)0, (size_t)0)), sizeof...(J)> {
        return {
          f(i, J) ...
        };
    }, iS);
}



template<typename Vp, typename Vt = boost::none_t, typename Vn = boost::none_t>
std::ostream& writeFace(std::ostream& out, Vp vp, Vt vt = Vt(), Vn vn = Vn()) {
    auto vpIndices = intCollection(vp);
    auto vtIndices = intCollection(vt);
    auto vnIndices = intCollection(vn);
    out << "f";
    for(size_t i = 0; i < vpIndices.size(); ++i){
        out << " " << vpIndices[i];
        if(vtIndices.size() || vnIndices.size()) {
            out << "/";
            if(vtIndices.size()) {
                out << vtIndices[i];
            }
            if(vnIndices.size()){
                out << "/";
                out << vnIndices[i];
            }
        }
    }
    return (out << std::endl);
}

template<typename Vp, typename Vt = boost::none_t, typename Vn = boost::none_t>
struct LazyFace {
    Vp vp;
    Vt vt;
    Vn vn;
private:
    bool converted;
public:
    LazyFace(Vp vpI, Vt vtI = Vt(), Vn vnI = Vn())
    : vp(vpI), vt(vtI), vn(vnI), converted(false) {}
    template<typename Vp_, typename Vt_, typename Vn_>
    friend std::ostream& operator<<(std::ostream& os, LazyFace<Vp_, Vt_, Vn_>& lf);
};

template<typename Vp, typename Vt, typename Vn>
std::ostream& operator<<(std::ostream& os, LazyFace<Vp, Vt, Vn>& lf) {
    return (lf.converted
            ? os
            : ((lf.converted = true), writeFace(os, lf.vp, lf.vt, lf.vn)));
}

std::array<LazyVertex, 4> getUVs(int width, int heightN, int heightD, int colorId, bool rotate90, VertexCache<OBJVertStream>& uvCache) {
    constexpr static const int UNIT = 100;
    int colorU = colorId % 16;
    int colorV = colorId / 16;

    int baseU = 16;
    int baseV = 16;

    int brickU;
    int brickV;
    switch(width) {
        case 1: {
            brickU = 0;
            brickV = (UNIT*heightN)/heightD;
        } break;
        case 2: {
            brickU = 0;
            brickV = (2*UNIT*heightN)/heightD;
        } break;
        case 3: {
            brickU = 0;
            brickV = (3*UNIT*heightN)/heightD;
        } break;
        case 4:{
            brickU = (10-width)*UNIT;
            brickV = (3*UNIT*heightN)/heightD;
        } break;
        case 6: {
            brickU = (10-width)*UNIT;
            brickV = (2*UNIT*heightN)/heightD;
        } break;
        case 8: {
            brickU = (10-width)*UNIT;
            brickV = (UNIT*heightN)/heightD;
        } break;
        case 10:{
            brickU = (10-width)*UNIT;
            brickV = 0;
        } break;
        default: throw std::out_of_range("Unknown width");
    }
    if(heightN == 1 && heightD == 1) {
        ; // no-op
    } else if(heightN == 2 && heightD == 1) {
        baseU += 1024;
    } else if(heightN == 6 && heightD == 5) {
        baseV += 1024;
    } else {
        throw std::out_of_range("Unknown height");
    }

    Vector3 color(colorU, colorV, 0);
    Vector3 base(baseU, baseV, 0);
    // CCW around, starting from origin in +v
    std::array<Vector3, 4> rawUVs{
        Vector3(brickU, brickV, 0),
        Vector3(brickU, brickV + (UNIT * heightN) / heightD, 0),
        Vector3(brickU + (UNIT * width), brickV + (UNIT * heightN) / heightD, 0),
        Vector3(brickU + (UNIT * width), brickV, 0)};
    std::cout << width << "x" << (heightN/float(heightD)) << " face has CCW coordinates: {";
    std::cout << "[" << (base+rawUVs[0]) << "], [" << (base+rawUVs[1]) << "], ["  << (base+rawUVs[2]) << "], ["  << (base+rawUVs[3]) << "]" ;
    std::cout << std::endl;
    if(rotate90) {
        /* now start along a +u edge */
        std::rotate(rawUVs.begin(), rawUVs.begin() + 1, rawUVs.end());
        //std::swap(rawUVs[1], rawUVs[3]);
    }
    return unrollLoop1D([&](size_t i) -> LazyVertex {
        return {color + (base + rawUVs[i]) / 2048, &uvCache};
    }, index_sequence<4>());
}

LegoCloudNode::LegoCloudNode()
  : legoCloud_(new LegoCloud()), renderLayerByLayer_(false), renderLayer_(0), knobList_(glGenLists(1)),
    renderBricks_(true), renderGraph_(false), colorRendering_(RealColor), drawDirty_(true)
{
    baseScale_ = 1.6;
    basePoint_ = Vector3(0,0,0);//(2653477.260585, -3695084.806065, 0);
}

LegoCloudNode::~LegoCloudNode()
{
  delete legoCloud_;
  //std::cout << "Node destroyed" << std::endl;
}

void LegoCloudNode::recomputeAABB()
{
  int minX = INT_MAX;
  int minLevel = 0;
  int minY = INT_MAX;

  int maxX = INT_MIN;
  int maxLevel = legoCloud_->getLevelNumber();
  int maxY = INT_MIN;

  for(int level=0; level < legoCloud_->getLevelNumber(); level++)
  {
    for(QList<LegoBrick>::const_iterator brick = legoCloud_->getBricks(level).begin(); brick != legoCloud_->getBricks(level).constEnd(); brick++)//QTL
    {
      if(brick->getPosX() < minX)
        minX = brick->getPosX();
      if(brick->getPosY() < minY)
        minY = brick->getPosY();

      if(brick->getPosX()+1 > maxX)
        maxX = brick->getPosX()+1;
      if(brick->getPosY()+1 > maxY)
        maxY = brick->getPosY()+1;
    }
  }

  boundsMin_ = Vector3(minX*LEGO_KNOB_DISTANCE, minLevel*LEGO_HEIGHT, minY*LEGO_KNOB_DISTANCE);
  boundsMax_ = Vector3(maxX*LEGO_KNOB_DISTANCE, maxLevel*LEGO_HEIGHT, maxY*LEGO_KNOB_DISTANCE);
}

void LegoCloudNode::render()
{
  /*float specularity_[4] = {0.7, 0.7, 0.7, 1.0};
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularity_);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 32);
  glLineWidth(getDefaultMaterial()->getEdgeWidth());
  glPointSize(getDefaultMaterial()->getVertexRadius());*/

//  glPushAttrib(GL_CULL_FACE | GL_POLYGON_MODE | GL_LINE_SMOOTH | GL_BLEND);
//  glEnable(GL_CULL_FACE);
//  glCullFace(GL_BACK);

//  glDisable(GL_LIGHTING);

//  glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
  //glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

  /*
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);                 // Set Line Antialiasing
  glEnable(GL_BLEND);                         // Enable Blending
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  */


  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(GL_SMOOTH);

  const LegoGraph& graph = legoCloud_->getLegoGraph();

  if(renderBricks_)
  {
//    glPushAttrib(GL_POLYGON_OFFSET_FILL);
//    glEnable(GL_POLYGON_OFFSET_FILL);
//    glPolygonOffset(1.0, 1.0);

    LegoGraph::vertex_iterator vertexIt, vertexItEnd;
    for (boost::tie(vertexIt, vertexItEnd) = boost::vertices(graph); vertexIt != vertexItEnd; ++vertexIt)
    {
      LegoBrick* brick = graph[*vertexIt].brick;
      if((!renderLayerByLayer_ && brick->isOuter())|| (renderLayerByLayer_ && brick->getLevel() == renderLayer_))
      {
        setColor(*vertexIt);
        drawLegoBrick(*brick);
        //drawNeighbourhood(*brick, legoCloud_->getNeighbours(brick));
      }
    }

    glColor3i(0,0,0);
    for (boost::tie(vertexIt, vertexItEnd) = boost::vertices(graph); vertexIt != vertexItEnd; ++vertexIt)
    {
      LegoBrick* brick = graph[*vertexIt].brick;
      if((!renderLayerByLayer_ && brick->isOuter())|| (renderLayerByLayer_ && brick->getLevel() == renderLayer_))
      {

        drawBrickOutline(*brick);
        //drawNeighbourhood(*brick, legoCloud_->getNeighbours(brick));
      }
    }
//    glPopAttrib();
  }

  drawDirty_ = false;

  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHTING);

  glDisable(GL_DEPTH_TEST);

  if(renderGraph_)
    drawLegoGraph(graph);

//  glDisable(GL_LINE_SMOOTH);
//  glDisable(GL_BLEND);
//  glDisable(GL_CULL_FACE);
  glPopAttrib();
}




void LegoCloudNode::drawLegoBrick(const LegoBrick &brick) const
{
  Vector3 p1;//Back corner down left
  p1[0] = brick.getPosX()*LEGO_KNOB_DISTANCE + LEGO_HORIZONTAL_TOLERANCE;
  p1[1] = brick.getLevel()*LEGO_HEIGHT;
  p1[2] = brick.getPosY()*LEGO_KNOB_DISTANCE + LEGO_HORIZONTAL_TOLERANCE;

  Vector3 p2;//Front corner up right
  p2[0] = p1[0] + brick.getSizeX()*LEGO_KNOB_DISTANCE - LEGO_HORIZONTAL_TOLERANCE;
  p2[1] = p1[1] + LEGO_HEIGHT;
  p2[2] = p1[2] + brick.getSizeY()*LEGO_KNOB_DISTANCE - LEGO_HORIZONTAL_TOLERANCE;

  drawBox(p1, p2);

  drawKnobs(brick, p1);
}

void LegoCloudNode::drawBox(const Vector3 &p1, const Vector3 &p2) const
{
  glBegin(GL_QUADS);
  //Bottom
  glNormal3f(0.0, -1.0, 0.0);
  glVertex3f(p1[0], p1[1], p1[2]);
  glVertex3f(p2[0], p1[1], p1[2]);
  glVertex3f(p2[0], p1[1], p2[2]);
  glVertex3f(p1[0], p1[1], p2[2]);

  //Back
  glNormal3f(0.0, 0.0, -1.0);
  glVertex3f(p1[0], p1[1], p1[2]);
  glVertex3f(p1[0], p2[1], p1[2]);
  glVertex3f(p2[0], p2[1], p1[2]);
  glVertex3f(p2[0], p1[1], p1[2]);

  //Right
  glNormal3f(1.0, 0.0, 0.0);
  glVertex3f(p2[0], p1[1], p1[2]);
  glVertex3f(p2[0], p2[1], p1[2]);
  glVertex3f(p2[0], p2[1], p2[2]);
  glVertex3f(p2[0], p1[1], p2[2]);

  //Front
  glNormal3f(0.0, 0.0, 1.0);
  glVertex3f(p2[0], p2[1], p2[2]);
  glVertex3f(p1[0], p2[1], p2[2]);
  glVertex3f(p1[0], p1[1], p2[2]);
  glVertex3f(p2[0], p1[1], p2[2]);

  //Left
  glNormal3f(-1.0, 0.0, 0.0);
  glVertex3f(p1[0], p1[1], p1[2]);
  glVertex3f(p1[0], p1[1], p2[2]);
  glVertex3f(p1[0], p2[1], p2[2]);
  glVertex3f(p1[0], p2[1], p1[2]);

  //Top
  glNormal3f(0.0, 1.0, 0.0);
  glVertex3f(p2[0], p2[1], p2[2]);
  glVertex3f(p2[0], p2[1], p1[2]);
  glVertex3f(p1[0], p2[1], p1[2]);
  glVertex3f(p1[0], p2[1], p2[2]);
  glEnd();

  //____________
}

void LegoCloudNode::drawBrickOutline(const LegoBrick &brick) const
{
  Vector3 p1;//Back corner down left
  p1[0] = brick.getPosX()*LEGO_KNOB_DISTANCE + LEGO_HORIZONTAL_TOLERANCE;
  p1[1] = brick.getLevel()*LEGO_HEIGHT;
  p1[2] = brick.getPosY()*LEGO_KNOB_DISTANCE + LEGO_HORIZONTAL_TOLERANCE;

  Vector3 p2;//Front corner up right
  p2[0] = p1[0] + brick.getSizeX()*LEGO_KNOB_DISTANCE - LEGO_HORIZONTAL_TOLERANCE;
  p2[1] = p1[1] + LEGO_HEIGHT;
  p2[2] = p1[2] + brick.getSizeY()*LEGO_KNOB_DISTANCE - LEGO_HORIZONTAL_TOLERANCE;

  glBegin(GL_LINE_LOOP);
  //Right
  glVertex3f(p2[0], p1[1], p1[2]);
  glVertex3f(p2[0], p2[1], p1[2]);
  glVertex3f(p2[0], p2[1], p2[2]);
  glVertex3f(p2[0], p1[1], p2[2]);
  glEnd();

  glBegin(GL_LINE_LOOP);
  //Left
  glVertex3f(p1[0], p1[1], p1[2]);
  glVertex3f(p1[0], p1[1], p2[2]);
  glVertex3f(p1[0], p2[1], p2[2]);
  glVertex3f(p1[0], p2[1], p1[2]);
  glEnd();

  glBegin(GL_LINES);
  glVertex3f(p1[0], p1[1], p1[2]);
  glVertex3f(p2[0], p1[1], p1[2]);

  glVertex3f(p1[0], p1[1], p2[2]);
  glVertex3f(p2[0], p1[1], p2[2]);

  glVertex3f(p1[0], p2[1], p2[2]);
  glVertex3f(p2[0], p2[1], p2[2]);

  glVertex3f(p1[0], p2[1], p1[2]);
  glVertex3f(p2[0], p2[1], p1[2]);
  glEnd();


  //Draw knob outline:
  Vector3 p;//Center of back left knob (top)
  p[0] = p1[0] + LEGO_KNOB_DISTANCE/2.0;
  p[1] = p1[1] + LEGO_HEIGHT + LEGO_KNOB_HEIGHT;
  p[2] = p1[2] + LEGO_KNOB_DISTANCE/2.0;

  for(int x = 0; x < brick.getSizeX(); ++x)
  {
    for(int y = 0; y < brick.getSizeY(); ++y)
    {

      //Draw loop on top
      glBegin(GL_LINE_LOOP);
      for(int i = 0; i <= KNOB_RESOLUTION_DISPLAY; ++i)
      {
        double angle = -i*(2*M_PI/double(KNOB_RESOLUTION_DISPLAY));
        glVertex3d(p[0] + x*LEGO_KNOB_DISTANCE + cos(angle)*LEGO_KNOB_RADIUS,
                   p[1],
                   p[2] + y*LEGO_KNOB_DISTANCE + sin(angle)*LEGO_KNOB_RADIUS);
      }
      glEnd();

      glBegin(GL_LINE_LOOP);
      for(int i = 0; i <= KNOB_RESOLUTION_DISPLAY; ++i)
      {
        double angle = -i*(2*M_PI/double(KNOB_RESOLUTION_DISPLAY));
        glVertex3d(p[0] + x*LEGO_KNOB_DISTANCE + cos(angle)*LEGO_KNOB_RADIUS,
                   p[1] - LEGO_KNOB_HEIGHT,
                   p[2] + y*LEGO_KNOB_DISTANCE + sin(angle)*LEGO_KNOB_RADIUS);
      }
      glEnd();
    }
  }

}

void LegoCloudNode::drawKnobs(const LegoBrick &brick, const Vector3 &p1) const
{
  Vector3 p;//Center of back left knob (top)
  p[0] = p1[0] + LEGO_KNOB_DISTANCE/2.0;
  p[1] = p1[1] + LEGO_HEIGHT + LEGO_KNOB_HEIGHT;
  p[2] = p1[2] + LEGO_KNOB_DISTANCE/2.0;

  //glColor3fv(brick.getColor());

  for(int x = 0; x < brick.getSizeX(); ++x)
  {
    for(int y = 0; y < brick.getSizeY(); ++y)
    {

      //Draw cylinder
      glBegin(GL_QUAD_STRIP);
      for(int i = 0; i <= KNOB_RESOLUTION_DISPLAY; ++i)
      {
        double angle = -i*(2*M_PI/double(KNOB_RESOLUTION_DISPLAY));
        glNormal3d(cos(angle), 0.0, sin(angle));
        glVertex3d(p[0] + x*LEGO_KNOB_DISTANCE + cos(angle)*LEGO_KNOB_RADIUS,
                   p[1],
                   p[2] + y*LEGO_KNOB_DISTANCE + sin(angle)*LEGO_KNOB_RADIUS);

        glVertex3d(p[0] + x*LEGO_KNOB_DISTANCE + cos(angle)*LEGO_KNOB_RADIUS,
                   p[1] - LEGO_KNOB_HEIGHT,
                   p[2] + y*LEGO_KNOB_DISTANCE + sin(angle)*LEGO_KNOB_RADIUS);
      }
      glEnd();

      //Draw disk on top
      glBegin(GL_TRIANGLE_FAN);
      glNormal3d(0.0, 1.0, 0.0);
      glVertex3d(p[0]+x*LEGO_KNOB_DISTANCE, p[1], p[2]+y*LEGO_KNOB_DISTANCE);
      for(int i = 0; i <= KNOB_RESOLUTION_DISPLAY; ++i)
      {
        double angle = -i*(2*M_PI/double(KNOB_RESOLUTION_DISPLAY));
        glVertex3d(p[0] + x*LEGO_KNOB_DISTANCE + cos(angle)*LEGO_KNOB_RADIUS,
                   p[1],
                   p[2] + y*LEGO_KNOB_DISTANCE + sin(angle)*LEGO_KNOB_RADIUS);
      }
      glEnd();
    }
  }
}

void LegoCloudNode::drawNeighbourhood(const LegoBrick &brick, const QSet<LegoBrick *> &neighbours) const
{

  const double DRAW_HEIGHT = 0.01;
  Vector3 brickCenter;
  brickCenter[0] = brick.getPosX()*LEGO_KNOB_DISTANCE + (brick.getSizeX()*LEGO_KNOB_DISTANCE)/2.0;
  brickCenter[1] = brick.getLevel()*LEGO_HEIGHT + LEGO_KNOB_HEIGHT + DRAW_HEIGHT;
  brickCenter[2] = brick.getPosY()*LEGO_KNOB_DISTANCE + (brick.getSizeY()*LEGO_KNOB_DISTANCE)/2.0;

  //glLineWidth(2);
  glColor3fv(brick.getRandColor().data());

  Vector3 neighbourCenter;

  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);

  foreach(const LegoBrick* neighbour, neighbours)
  {
    neighbourCenter[0] = neighbour->getPosX()*LEGO_KNOB_DISTANCE + (neighbour->getSizeX()*LEGO_KNOB_DISTANCE)/2.0;
    neighbourCenter[1] = neighbour->getLevel()*LEGO_HEIGHT + LEGO_KNOB_HEIGHT + DRAW_HEIGHT;
    neighbourCenter[2] = neighbour->getPosY()*LEGO_KNOB_DISTANCE + (neighbour->getSizeY()*LEGO_KNOB_DISTANCE)/2.0;
    glVertex3fv(brickCenter.data());
    glVertex3fv(neighbourCenter.data());
  }

  glEnd();
  glEnable(GL_LIGHTING);
}

void LegoCloudNode::drawLegoGraph(const LegoGraph & graph) const
{
  //Draw vertices

  glPushAttrib(GL_LIGHTING);
  glDisable(GL_LIGHTING);
  glBegin(GL_POINTS);
  LegoGraph::vertex_iterator vertexIt, vertexItEnd;
  for (boost::tie(vertexIt, vertexItEnd) = boost::vertices(graph); vertexIt != vertexItEnd; ++vertexIt)
  {
    const LegoBrick* brick = graph[*vertexIt].brick;
    if(!renderLayerByLayer_ || (renderLayerByLayer_ && (brick->getLevel() == renderLayer_ || brick->getLevel() == renderLayer_+1)))
    {
      Vector3 brickCenter;
      brickCenter[0] = brick->getPosX()*LEGO_KNOB_DISTANCE + (brick->getSizeX()*LEGO_KNOB_DISTANCE)/2.0;
      brickCenter[1] = brick->getLevel()*LEGO_HEIGHT + LEGO_HEIGHT;
      brickCenter[2] = brick->getPosY()*LEGO_KNOB_DISTANCE + (brick->getSizeY()*LEGO_KNOB_DISTANCE)/2.0;

      //glColor3d(1,0,0);
      setColor(*vertexIt);
      glVertex3fv(brickCenter.data());
    }
  }
  glEnd();

  //Draw edges
  LegoGraph::edge_iterator edgeIt, edgeEnd;
  //glColor3d(0,0,1);

  glBegin(GL_LINES);
  for (boost::tie(edgeIt, edgeEnd) = boost::edges(graph); edgeIt != edgeEnd; ++edgeIt)
  {
    const LegoBrick* source = graph[boost::source(*edgeIt, graph)].brick;
    const LegoBrick* target = graph[boost::target(*edgeIt, graph)].brick;
    if(!renderLayerByLayer_ || (renderLayerByLayer_ &&
                                (source->getLevel() == renderLayer_ || source->getLevel() == renderLayer_+1) &&
                                (target->getLevel() == renderLayer_ || target->getLevel() == renderLayer_+1)))
    {
      Vector3 sourceCenter;
      sourceCenter[0] = source->getPosX()*LEGO_KNOB_DISTANCE + (source->getSizeX()*LEGO_KNOB_DISTANCE)/2.0;
      sourceCenter[1] = source->getLevel()*LEGO_HEIGHT + LEGO_HEIGHT;
      sourceCenter[2] = source->getPosY()*LEGO_KNOB_DISTANCE + (source->getSizeY()*LEGO_KNOB_DISTANCE)/2.0;

      Vector3 targetCenter;
      targetCenter[0] = target->getPosX()*LEGO_KNOB_DISTANCE + (target->getSizeX()*LEGO_KNOB_DISTANCE)/2.0;
      targetCenter[1] = target->getLevel()*LEGO_HEIGHT + LEGO_HEIGHT;
      targetCenter[2] = target->getPosY()*LEGO_KNOB_DISTANCE + (target->getSizeY()*LEGO_KNOB_DISTANCE)/2.0;

      glColor3d(0,0,1);

      glVertex3fv(sourceCenter.data());
      glVertex3fv(targetCenter.data());
    }

  }
  glEnd();

  glPopAttrib();
}

void LegoCloudNode::setColor(const LegoGraph::vertex_descriptor& vertex) const
{
  const LegoGraph& graph = legoCloud_->getLegoGraph();

  switch(colorRendering_)
  {
    case RealColor:
      glColor3fv(legoCloud_->getLegalColor()[graph[vertex].brick->getColorId()].data());
      break;

    case Random:
      //glColor3d(boost::out_degree(vertex, graph)/10.0, 0.0, 0.0);

      //glColor3fv(graph[vertex].brick->getRandColor());
      glColor3fv(legoCloud_->getLegalColor()[graph[vertex].brick->getHash() % legoCloud_->getLegalColor().size()].data());
      /*if(graph[vertex].brick->isOuter())
        glColor3fv(graph[vertex].brick->getRandColor());
      else
        glColor3fv(graph[vertex].brick->getRandColor() + Color3(0.5, 0,0));*/
      break;

    case ConnectedComp:
      {
        int red = 31 + graph[vertex].connected_comp;
        int green = 31*red + graph[vertex].connected_comp;
        int blue = 31*green + graph[vertex].connected_comp;
        int mod = 50;

        glColor3d((red%mod)/double(mod),
                  (green%mod)/double(mod),
                  (blue%mod)/double(mod));
      }
      break;

    case BiconnectedComp:
      if(graph[vertex].badArticulationPoint)
      {
        glColor3d(1.0,0.0,0.0);
      }
      else
      {
        glColor3fv(graph[vertex].brick->getRandColor().data());
      }
      break;

    default:
      glColor3d(0.0,0.0,0.0);
  }
}

void LegoCloudNode::drawInstructions(QGraphicsScene *scene, bool hintLayerBelow)
{
  const int BRICK_PIXEL_SIZE = 20;

  scene->clear();
  scene->setSceneRect(0, 0, legoCloud_->getWidth()*BRICK_PIXEL_SIZE, legoCloud_->getDepth()*BRICK_PIXEL_SIZE);

  for(QList<LegoBrick>::const_iterator brick = legoCloud_->getBricks(renderLayer_).begin(); brick != legoCloud_->getBricks(renderLayer_).constEnd(); brick++)//QTL
  {
    Color3 color = legoCloud_->getLegalColor()[brick->getColorId()];
    scene->addRect(brick->getPosX()*BRICK_PIXEL_SIZE, brick->getPosY()*BRICK_PIXEL_SIZE, brick->getSizeX()*BRICK_PIXEL_SIZE, brick->getSizeY()*BRICK_PIXEL_SIZE, QPen(),
                   QBrush(QColor(color[0]*255, color[1]*255, color[2]*255), Qt::SolidPattern));

    for(int x = 0; x < brick->getSizeX(); ++x)
    {
      for(int y = 0; y < brick->getSizeY(); ++y)
      {
        double knobRadius = LEGO_KNOB_RADIUS / LEGO_KNOB_DISTANCE;
        scene->addEllipse((brick->getPosX() + x + (0.5 - knobRadius))*BRICK_PIXEL_SIZE, (brick->getPosY() + y + (0.5 - knobRadius))*BRICK_PIXEL_SIZE, 2.0*knobRadius*BRICK_PIXEL_SIZE, 2.0*knobRadius*BRICK_PIXEL_SIZE);
      }
    }
  }

  //Draw the underneath layer (if there is one)
  if(renderLayer_ >= 1 && hintLayerBelow)
  {
    //foreach(const LegoBrick& brick, legoCloud_->getBricks(renderLayer_-1))
    for(QList<LegoBrick>::const_iterator brick = legoCloud_->getBricks(renderLayer_-1).begin(); brick != legoCloud_->getBricks(renderLayer_-1).constEnd(); brick++)
    {
      QColor color(0, 0, 0, 200);
      //scene->addRect(brick.getPosX(), brick.getPosY(), brick.getSizeX(), brick.getSizeY(), QPen(Qt::NoPen), QBrush(color, Qt::SolidPattern));
      scene->addRect(brick->getPosX()*BRICK_PIXEL_SIZE, brick->getPosY()*BRICK_PIXEL_SIZE, brick->getSizeX()*BRICK_PIXEL_SIZE,
                     brick->getSizeY()*BRICK_PIXEL_SIZE, QPen(Qt::NoPen), QBrush(color, Qt::Dense5Pattern));
    }
  }
}

void LegoCloudNode::exportToLdr(QString filename)
{
    std::ofstream ldrFile(filename.toStdString().c_str());
    if(!ldrFile.is_open())
    {
        std::cerr << "LegoCloudNode: unable to create or open the file: " << filename.toStdString().c_str() << std::endl;
    }
    else
    {
        ldrFile << legoCloud_->toLDraw();
    }
    ldrFile.close();
}

void LegoCloudNode::exportToObj(QString filename)
{
  if(KNOB_RESOLUTION_OBJ_EXPORT % 4){
      std::cerr << "This code path doesn't work if the knob resolution is not a multiple of 4" << std::endl;
      std::cerr << "But it should be easy to fix!" << std::endl;
  } else {
      const float LEGO_VERTICAL_TOLERANCE = 0.0001f;

      Vector3 trueBasePoint = (basePoint_ / baseScale_)* LEGO_KNOB_DISTANCE;
      float oldY = trueBasePoint.y();
      //trueBasePoint.x() = -trueBasePoint.x();
      trueBasePoint.y() = trueBasePoint.z();
      trueBasePoint.z() = -oldY;

      std::cout << "True base point: " << trueBasePoint << std::endl;

      //const QList<LegoBrick*>& outterBricks = legoCloud_->getOuterBricks();
      std::cout << "Flood filling open space markers...";
      LegoBrick *outside = new LegoBrick(-1, 0, 0, legoCloud_->getWidth(), legoCloud_->getDepth());
      QMap<VoxelCoord, QVector<VoxelCoord>> queue;
      auto visit = [&queue,&outside,this](const VoxelCoord& vc) {
          auto& voxelGrid = legoCloud_->getVoxelGrid();
          LegoBrick*& voxel = voxelGrid[vc.z][vc.x][vc.y];
          if(nullptr == voxel) {
              voxel = outside;
              QVector<VoxelCoord> neighbors = legoCloud_->getL1Neighbors(vc);
              brickr::erase_if(neighbors, [&voxelGrid](const VoxelCoord& vc){
                  return nullptr != voxelGrid[vc.z][vc.x][vc.y];
              });
              queue.insert(vc,neighbors);
          }
      };
      for(int x = 0; x < legoCloud_->getWidth(); ++x) {
          for(int y = 0; y < legoCloud_->getDepth(); ++y) {
              visit({x, y, 0});
              visit({x, y, legoCloud_->getHeight() - 1});
          }
      }
      for(int z = 0; z < legoCloud_->getHeight(); ++z) {
          for(int y = 0; y < legoCloud_->getDepth(); ++y) {
              visit({0, y, z});
              visit({legoCloud_->getWidth() - 1, y, z});
          }
      }
      for(int z = 0; z < legoCloud_->getHeight(); ++z) {
          for(int x = 0; x < legoCloud_->getWidth(); ++x) {
              visit({x, 0, z});
              visit({x, legoCloud_->getDepth() - 1, z});
          }
      }
      while(queue.size()) {
          QVector<VoxelCoord> neighbors = queue.take(queue.firstKey());
          for(auto it = neighbors.begin(); it != neighbors.end(); ++it){
              visit(*it);
          }
      }
      std::cout << "Done" << std::endl;
      //Create the file
      std::ofstream objFile (filename.toStdString().c_str());
      if (!objFile.is_open())
      {
        std::cerr << "LegoCloudNode: unable to create or open the file: " << filename.toStdString().c_str() << std::endl;
      }
      objFile << "mtllib ldraw_palette.mtl" << std::endl;

      for(int i = 1; i <= KNOB_RESOLUTION_OBJ_EXPORT; ++i){
          // must be multiple of 4 to get cardinal directions encoded as well!
          double angle = -(2*M_PI)*(i/double(KNOB_RESOLUTION_OBJ_EXPORT));
          Vector3 normal(cos(angle),0,sin(angle)); // iterate clockwise: right, back, left, front
          objFile << "vn " << normal << std::endl;
      }
      objFile << "vn " << Vector3(0, -1, 0) << std::endl; // down
      objFile << "vn " << Vector3(0, 1, 0) << std::endl; // up

      /*
      for(int i = 0; i < legoCloud_->getLegalColor().size(); ++i) {
          float u = (0.5 + (i % 16)) / 16.0;
          float v = (0.5 + (i / 16)) / 16.0;
          objFile << "vt " << u << " " << v << std::endl;
      }
      */

      VertexCache<OBJVertStream> uvCache({objFile, "t", 2});

      VertexCache<OBJVertStream> xyzCache({objFile, "", 3});
      using LazyQuadIdx = std::array<LazyVertex,4>;
      using QuadIdx = std::array<int,4>;
      using TriIdx = std::array<int, 3>;
      using LazyTriIdx = std::array<LazyVertex, 3>;

      // left/right, back/front, bottom/top
      auto faceIdxToNormIdx = [](size_t i){
          static constexpr size_t QUARTER_TURN = KNOB_RESOLUTION_OBJ_EXPORT / 4;
          static constexpr size_t HALF_TURN = KNOB_RESOLUTION_OBJ_EXPORT / 2;
          static constexpr size_t FULL_TURN = KNOB_RESOLUTION_OBJ_EXPORT;
          static constexpr size_t BASES[3] = {HALF_TURN, QUARTER_TURN, FULL_TURN + 1};
          static constexpr size_t OFFSETS[3] = {HALF_TURN, HALF_TURN, 1};
          size_t axis = i / 2;
          bool offset = i % 2;
          return BASES[axis] + offset*OFFSETS[axis];
      };

      auto faceIdxToFirstVertexIdxs = [](size_t i) {
          size_t axis = i / 2;
          bool offset = (i % 2);
          size_t fixedBit = offset << ((axis + 2) % 3);
          size_t firstBit = 1 << ((axis + offset) % 3);
          size_t secondBit = 1 << ((axis + !offset) % 3);
          return std::array<size_t, 4>{fixedBit, fixedBit | firstBit, fixedBit | firstBit | secondBit, fixedBit | secondBit};
      };

      auto faceIdxToSecondVertexIdx = [](size_t i) -> size_t {
          static constexpr size_t idxMap[3] = {1, 0, 2}; // execute y<->z swap going from voxels to obj
          return 0;//idxMap[i / 2];
      };

      auto faceIdxToFirstEdgeIdxs = [](size_t face) -> std::array<size_t, 4> {
          size_t axis = face / 2;
          bool offset = face % 2;
          size_t r0 = 2 - offset;
          size_t r1 = 1 + offset;

          return {(axis+r0)%3,(axis+r1)%3,(axis+r0)%3,(axis+r1)%3};
      };

      auto faceIdxToSecondEdgeIdxs = [](size_t face) -> std::array<size_t, 4> {
          size_t offset = face % 2;

          return {2 * offset, 2 * offset, 3, 1};
      };

      objFile << "g bricks" << std::endl;
      objFile << "usemtl ldraw_palette" << std::endl;
      legoCloud_->forEachBrick([&](const LegoBrick* brick, const Vector3 (&pLocal)[2], int brickColor){
          Vector3 p[2] = {trueBasePoint + pLocal[0], trueBasePoint + pLocal[1]};
          Vector3 o[2] = {
              // back/down/left corner gets offset front/up/right
              Vector3(LEGO_HORIZONTAL_TOLERANCE, LEGO_VERTICAL_TOLERANCE, LEGO_HORIZONTAL_TOLERANCE),
              // front/up/right corner gets offset back/down/left
              Vector3(-LEGO_HORIZONTAL_TOLERANCE, -LEGO_VERTICAL_TOLERANCE, -LEGO_HORIZONTAL_TOLERANCE)
          };

          TriIdx triColors = {brickColor, brickColor, brickColor};
          QuadIdx quadColors = {brickColor, brickColor, brickColor, brickColor};

          BrickSize size = brick->getSize();
          int sizeX = brick->getSizeX();
          int sizeY = brick->getSizeY();
          auto faceIdxToUVs = [&](size_t face) -> LazyQuadIdx {
              size_t axis = face / 2;
              bool offset = face % 2;

              /*
              float transform[3][4] = {
                  {1, 0, 0, 10.0f*size.second}, // All the basic LDraw bricks are longer in X than Z
                  {0, 1, 0, -24.0f*brick.getLevel()}, // -Y up
                  {0, 0, 1, 10.0f*size.first}, // +Z into screen
              };
              if(size.first == brick.getSizeX()) {
                  std::swap(transform[0],transform[2]);
              }
              */
              LazyQuadIdx uvQuad;
              switch(axis) {
              case 0: {
                  uvQuad = getUVs(sizeY,6,5,brickColor,!offset,uvCache);
              } break;
              case 1: {
                  uvQuad = getUVs(sizeX,6,5,brickColor,offset,uvCache);
              } break;
              case 2: {
                  uvQuad = getUVs(size.second, size.first,1,brickColor,sizeX == size.first, uvCache);
                  //std::swap(uvQuad[1],uvQuad[3]);
              } break;
              default: throw std::out_of_range("3d geometry has no axis greater than 2");
              }
              return uvQuad;
          };
          auto lazyBoxUVs = unrollLoop1D(faceIdxToUVs, index_sequence<6>());


          //Write 8 vertices of box * 3 way split for bevel (x/y off, y/z off, z/x off)
          auto lazyBoxVerts = [&xyzCache,&p, &o](size_t corner, size_t offDir) -> LazyVertex {
              bool pM[3] = {(bool)(corner & 4), (bool)(corner & 2), (bool)(corner & 1)};

              bool oM[3];
              oM[(offDir + 0) % 3] = true;
              oM[(offDir + 1) % 3] = true;
              oM[(offDir + 2) % 3] = false;

              Vector3 retVec;
              for(int i = 0; i < 3; ++i) {
                  retVec.data()[i] = p[pM[i]][i];// + oM[i]*o[pM[i]][i];
              }

              return {retVec, &xyzCache};
          };

          std::array<std::array<LazyVertex, /*3*/1>, 8> boxVerts = unrollLoop2D(lazyBoxVerts, index_sequence<8>(), index_sequence</*3*/1>());

          using BrickSizeF = int(LegoBrick::*)() const;
          auto faceIdxToNeighbor = [&brick](size_t i ) -> VoxelCoord {
              static constexpr BrickSizeF axisSizeFs[3] = {&LegoBrick::getSizeX, &LegoBrick::getSizeY, &LegoBrick::getSizeLevels};
              VoxelCoord vc{0, 0, 0};
              size_t axis = i / 2;
              bool offset = i % 2;
              vc[axis] = 1;
              return vc * (offset ? (brick->*axisSizeFs[axis])() : -1);
          };

          /*
          std::array<LazyFace<LazyTriIdx,TriIdx,TriIdx>,8> bevelCorners = unrollLoop1D([&boxVerts, &triColors, &faceIdxToSecondVertexIdx, &faceIdxToNormIdx](size_t i) -> LazyFace<LazyTriIdx,TriIdx,TriIdx> {
              std::array<size_t, 3> faceJs = {0, 1, 2};
              for(size_t j = 0; j < 3; ++j){
                  if((bool)(i & 1) ^ (bool)(i & 2) ^ (bool)(i & 4)){
                      faceJs[j] = 2 - faceJs[j];
                  }

                  faceJs[j] = 2*faceJs[j] + (bool)(i & (1 << ((faceJs[j] + 2) % 3)));
              }
              return LazyFace<LazyTriIdx,TriIdx,TriIdx>(
                          unrollLoop1D([&i,&boxVerts,&faceIdxToSecondVertexIdx,&faceJs](size_t j) -> LazyVertex {
                            return boxVerts[i][faceIdxToSecondVertexIdx(faceJs[j])];
                          }, index_sequence<3>()),
                          triColors,
                          unrollLoop1D([&faceIdxToNormIdx,&faceJs](size_t i) -> int{
                            return faceIdxToNormIdx(faceJs[i]);
                          }, index_sequence<3>()));
          }, index_sequence<8>());

          std::array<std::array<LazyFace<LazyQuadIdx,QuadIdx,QuadIdx>, 4>, 3> bevelEdges = unrollLoop2D([&](size_t axis, size_t corner) -> LazyFace<LazyQuadIdx,QuadIdx,QuadIdx>{
              size_t faceIs[2] = {2*axis,2*axis+1};
              std::array<std::array<size_t, 4>,2> vertexIs{
                  faceIdxToFirstVertexIdxs(faceIs[0]),
                  faceIdxToFirstVertexIdxs(faceIs[1])
              };
              std::swap(vertexIs[1][1],vertexIs[1][3]);
              size_t faceJs[4] = {2*((axis+1)%3),2*((axis+2)%3),2*((axis+1)%3)+1,2*((axis+2)%3)+1};
              return LazyFace<LazyQuadIdx,QuadIdx,QuadIdx>(
                          unrollLoop1D([&](size_t i) -> LazyVertex {
                              bool offset = ((i+1) % 4)/2;
                              return boxVerts[vertexIs[offset][corner]][faceIdxToSecondVertexIdx(faceJs[(corner+i/2)%4])];
                          }, index_sequence<4>()),
                          quadColors,
                          unrollLoop1D([&](size_t i) -> int{
                            return faceIdxToNormIdx(faceJs[(corner + i / 2) % 4]);
                          }, index_sequence<4>()));
           }, index_sequence<3>(), index_sequence<4>());
           //*/

          for(size_t face = 0; face < 6; ++face) {
              if(face == 4) {
                  continue; // skip the bottom
              }
              VoxelCoord neighbor = faceIdxToNeighbor(face);
              if(legoCloud_->visible(brick, outside, neighbor)) {
                  int normIdx = faceIdxToNormIdx(face);
                  std::array<size_t, 4> is = faceIdxToFirstVertexIdxs(face);
                  size_t j = faceIdxToSecondVertexIdx(face);
                  writeFace(objFile,
                            LazyQuadIdx{boxVerts[is[0]][j],boxVerts[is[1]][j],boxVerts[is[2]][j],boxVerts[is[3]][j]},
                            lazyBoxUVs[face],//quadColors,
                            QuadIdx{normIdx, normIdx, normIdx, normIdx});
                  std::array<size_t, 4> ks = faceIdxToFirstEdgeIdxs(face);
                  std::array<size_t, 4> ls = faceIdxToSecondEdgeIdxs(face);
                  /*
                  for(size_t i = 0; i < 4; ++i){
                      objFile << bevelCorners[is[i]];
                      objFile << bevelEdges[ks[i]][ls[i]];
                  }
                  //*/
              }
          }
      });

      objFile << "g studs" << std::endl;
      objFile << "usemtl ldraw_palette" << std::endl;
      legoCloud_->forEachBrick([&](const LegoBrick* brick, const Vector3 (&pLocal)[2], int brickColor){
          Vector3 p[2] = {trueBasePoint + pLocal[0], trueBasePoint + pLocal[1]};
          Vector3 knobCenter;//Center of back left knob (top)
          knobCenter[0] = p[0][0] + LEGO_KNOB_DISTANCE/2.0;
          knobCenter[1] = p[0][1] + LEGO_HEIGHT + LEGO_KNOB_HEIGHT;
          knobCenter[2] = p[0][2] + LEGO_KNOB_DISTANCE/2.0;

          Vector3 colorUV(brickColor % 16, brickColor / 16, 0);
          LazyVertex studUV{colorUV + Vector3(0.75, 0.75,0), &uvCache};

          LazyQuadIdx quadColors = {studUV, studUV, studUV, studUV};

          QVector<QVector<bool>> visibleStud(brick->getSizeX());
          QVector<QVector<QVector<LazyVertex>>> studIndices(brick->getSizeX());

          int nextLevel = 1 + brick->getLevel();
          //Write knobs vertices
          for(int x = 0; x < brick->getSizeX(); ++x)
          {
              visibleStud[x] = QVector<bool>(brick->getSizeY());
              studIndices[x] = QVector<QVector<LazyVertex>>(brick->getSizeY());
            for(int y = 0; y < brick->getSizeY(); ++y)
            {
                visibleStud[x][y] = (nextLevel == legoCloud_->getLevelNumber()) || (outside == legoCloud_->getVoxelGrid()[nextLevel][x+brick->getPosX()][y+brick->getPosY()]);
                if(visibleStud[x][y])
                {
                    studIndices[x][y] = QVector<LazyVertex>(2 * KNOB_RESOLUTION_OBJ_EXPORT);

                    for(int i = 0; i < KNOB_RESOLUTION_OBJ_EXPORT; ++i)
                    {
                      double angle = -i*(2*M_PI/double(KNOB_RESOLUTION_OBJ_EXPORT));
                      studIndices[x][y][i] = LazyVertex{
                          Vector3(knobCenter[0] + x*LEGO_KNOB_DISTANCE + cos(angle)*LEGO_KNOB_RADIUS,
                                  knobCenter[1],
                                  knobCenter[2] + y*LEGO_KNOB_DISTANCE + sin(angle)*LEGO_KNOB_RADIUS), &xyzCache};
                      studIndices[x][y][KNOB_RESOLUTION_OBJ_EXPORT+i] = LazyVertex{
                          Vector3(knobCenter[0] + x*LEGO_KNOB_DISTANCE + cos(angle)*LEGO_KNOB_RADIUS,
                                  knobCenter[1] - LEGO_KNOB_HEIGHT,
                                  knobCenter[2] + y*LEGO_KNOB_DISTANCE + sin(angle)*LEGO_KNOB_RADIUS), &xyzCache};
                    }

                } else {
                    studIndices[x][y] = QVector<LazyVertex>(0);
                }
            }
          }

          using CapIdx = std::array<int, KNOB_RESOLUTION_OBJ_EXPORT>;
          using LazyCapIdx = std::array<LazyVertex, KNOB_RESOLUTION_OBJ_EXPORT>;
          CapIdx capNorm;
          capNorm.fill(KNOB_RESOLUTION_OBJ_EXPORT + 2);
          LazyCapIdx capColor;
          capColor.fill(studUV);
          //Write top caps face
          for(int x = 0; x < brick->getSizeX(); ++x)
          {
            for(int y = 0; y < brick->getSizeY(); ++y)
            {
                if(visibleStud[x][y]){
                    auto& knobIndex = studIndices[x][y];
                    writeFace(objFile,
                              qMakePair(knobIndex, index_sequence<KNOB_RESOLUTION_OBJ_EXPORT>()),
                              capColor,
                              capNorm);
                }
            }
          }

          //Write cylinder faces
          for(int x = 0; x < brick->getSizeX(); ++x)
          {
            for(int y = 0; y < brick->getSizeY(); ++y)
            {
                if(visibleStud[x][y]){
                    auto& knobIndex = studIndices[x][y];

                    for(int i = 0; i < KNOB_RESOLUTION_OBJ_EXPORT; ++i)
                    {
                        int j = (i + 1) % KNOB_RESOLUTION_OBJ_EXPORT;
                        int iNorm = i ? i : KNOB_RESOLUTION_OBJ_EXPORT; // un-mod
                        int jNorm = j ? j : KNOB_RESOLUTION_OBJ_EXPORT; // un-mod
                        writeFace(objFile,
                                  LazyQuadIdx{knobIndex[i], knobIndex[KNOB_RESOLUTION_OBJ_EXPORT + i], knobIndex[KNOB_RESOLUTION_OBJ_EXPORT + j], knobIndex[j]},
                                  quadColors,
                                  QuadIdx{iNorm, iNorm, jNorm, jNorm});

                    }
                }
            }
          }
      });

      objFile.close();

      std::cout << "Clearing open space markers...";
      for(int z = 0; z < legoCloud_->getHeight(); ++z) {
          for(int x = 0; x < legoCloud_->getWidth(); ++x) {
              for(int y = 0; y < legoCloud_->getDepth(); ++y) {
                  LegoBrick*& voxel = legoCloud_->getVoxelGrid()[z][x][y];
                  if(outside == voxel){
                      voxel = nullptr;
                  }
              }
          }
      }
      delete outside;
      std::cout << "Done" << std::endl;
  }
}

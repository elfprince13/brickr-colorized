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

template<size_t N> using index_sequence = typename partial_index_sequence<N>::type;

template<size_t N>
std::array<int, N> intCollection(const std::array<int, N>& i) { return i; }

template<typename Collection, size_t ...I>
std::array<int, sizeof...(I)> intCollection(QPair<Collection, integer_sequence<I...>>& pair) {
    return {pair.first[I] ...};
}

QVector<int> intCollection(const QVector<int>& v) { return v; }

QVector<int> intCollection(boost::none_t) { return QVector<int>(); }



template<typename Vp, typename Vt = boost::none_t, typename Vn = boost::none_t>
void writeFace(std::ostream& out, Vp vp, Vt vt = Vt(), Vn vn = Vn()) {
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
            } else {
                out << "/";
            }
            if(vnIndices.size()){
                out << vnIndices[i];
            }
        }
    }
    out << std::endl;
}

LegoCloudNode::LegoCloudNode()
  : legoCloud_(new LegoCloud()), renderLayerByLayer_(false), renderLayer_(0), knobList_(glGenLists(1)),
    renderBricks_(true), renderGraph_(false), colorRendering_(RealColor), drawDirty_(true)
{

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
      int brickIndex = 0;
      int vertexIndex = 1;
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

      for(int i = 1; i <= KNOB_RESOLUTION_OBJ_EXPORT; ++i){
          // must be multiple of 4 to get cardinal directions encoded as well!
          double angle = -(2*M_PI)*(i/double(KNOB_RESOLUTION_OBJ_EXPORT));
          Vector3 normal(cos(angle),0,sin(angle)); // iterate clockwise: right, back, left, front
          objFile << "vn " << normal << std::endl;
      }
      objFile << "vn " << Vector3(0, 1, 0) << std::endl; // up
      objFile << "vn " << Vector3(0, -1, 0) << std::endl; // down

      VertexCache<OBJVertStream> xyzCache({objFile, "", 3});
      for(int level = 0; level < legoCloud_->getLevelNumber(); level++)
      {
        for(QList<LegoBrick>::const_iterator brickIt = legoCloud_->getBricks(level).begin(); brickIt != legoCloud_->getBricks(level).constEnd(); brickIt++)
        {
          const LegoBrick* brick = &(*brickIt);

          Vector3 p1;//Back corner down left
          p1[0] = brick->getPosX()*LEGO_KNOB_DISTANCE;// + LEGO_HORIZONTAL_TOLERANCE;
          p1[1] = brick->getLevel()*LEGO_HEIGHT;// + LEGO_VERTICAL_TOLERANCE;
          p1[2] = brick->getPosY()*LEGO_KNOB_DISTANCE;// + LEGO_HORIZONTAL_TOLERANCE;
          //p1 *= 10.0;

          Vector3 p2;//Front corner up right
          p2[0] = p1[0] + brick->getSizeX()*LEGO_KNOB_DISTANCE;// - LEGO_HORIZONTAL_TOLERANCE;
          p2[1] = p1[1] + LEGO_HEIGHT;// - LEGO_VERTICAL_TOLERANCE;
          p2[2] = p1[2] + brick->getSizeY()*LEGO_KNOB_DISTANCE;// - LEGO_HORIZONTAL_TOLERANCE;
          //p2 *= 10.0;

          Vector3 knobCenter;//Center of back left knob (top)
          knobCenter[0] = p1[0] + LEGO_KNOB_DISTANCE/2.0;
          knobCenter[1] = p1[1] + LEGO_HEIGHT + LEGO_KNOB_HEIGHT;
          knobCenter[2] = p1[2] + LEGO_KNOB_DISTANCE/2.0;


          //objFile << "g default" << std::endl;
          //Write 8 vertices of box
          LazyVertex boxVerts[8] = {
              {Vector3(p1[0], p1[1], p1[2]), &xyzCache},
              {Vector3(p1[0], p1[1], p2[2]), &xyzCache},
              {Vector3(p1[0], p2[1], p1[2]), &xyzCache},
              {Vector3(p1[0], p2[1], p2[2]), &xyzCache},
              {Vector3(p2[0], p1[1], p1[2]), &xyzCache},
              {Vector3(p2[0], p1[1], p2[2]), &xyzCache},
              {Vector3(p2[0], p2[1], p1[2]), &xyzCache},
              {Vector3(p2[0], p2[1], p2[2]), &xyzCache}
          };

          QVector<QVector<bool>> visibleStud(brick->getSizeX());
          QVector<QVector<QVector<LazyVertex>>> studIndices(brick->getSizeX());
          int studCount = 0;

          //Write knobs vertices
          for(int x = 0; x < brick->getSizeX(); ++x)
          {
              visibleStud[x] = QVector<bool>(brick->getSizeY());
              studIndices[x] = QVector<QVector<LazyVertex>>(brick->getSizeY());
            for(int y = 0; y < brick->getSizeY(); ++y)
            {
                visibleStud[x][y] = ((1 + level) == legoCloud_->getLevelNumber()) || (outside == legoCloud_->getVoxelGrid()[level+1][x+brickIt->getPosX()][y+brickIt->getPosY()]);
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
          //objFile << "s off" << std::endl;
          //Write box faces
          using QuadIdx = std::array<int,4>;
          if(legoCloud_->visible(brick, outside, {-1, 0, 0})) {
              //left
              int normIdx = KNOB_RESOLUTION_OBJ_EXPORT / 2;
              writeFace(objFile,
                        QuadIdx{boxVerts[0], boxVerts[1], boxVerts[3], boxVerts[2]},
                        boost::none,
                        QuadIdx{normIdx, normIdx, normIdx});
          }
          if(legoCloud_->visible(brick, outside, {brick->getSizeX(), 0, 0})) {
              //right
              int normIdx = KNOB_RESOLUTION_OBJ_EXPORT;
              writeFace(objFile,
                        QuadIdx{boxVerts[4], boxVerts[6], boxVerts[7], boxVerts[5]},
                        boost::none,
                        QuadIdx{normIdx, normIdx, normIdx});
          }
          if(legoCloud_->visible(brick, outside, {0, 0, -1})) {
              //bottom
              int normIdx = KNOB_RESOLUTION_OBJ_EXPORT + 2;
              writeFace(objFile,
                        QuadIdx{boxVerts[0], boxVerts[4], boxVerts[5], boxVerts[1]},
                        boost::none,
                        QuadIdx{normIdx, normIdx, normIdx});
          }
          if(legoCloud_->visible(brick, outside, {0, 0, 1})) {
              //top
              int normIdx = KNOB_RESOLUTION_OBJ_EXPORT + 1;
              writeFace(objFile,
                        QuadIdx{boxVerts[2], boxVerts[3], boxVerts[7], boxVerts[6]},
                        boost::none,
                        QuadIdx{normIdx, normIdx, normIdx});
          }
          if(legoCloud_->visible(brick, outside, {0, -1, 0})) {
              //back
              int normIdx = KNOB_RESOLUTION_OBJ_EXPORT / 4;
              writeFace(objFile,
                        QuadIdx{boxVerts[0], boxVerts[2], boxVerts[6], boxVerts[4]},
                        boost::none,
                        QuadIdx{normIdx, normIdx, normIdx});
          }
          if(legoCloud_->visible(brick, outside, {0, brick->getSizeY(), 0})) {
              //front
              int normIdx = 3 * (KNOB_RESOLUTION_OBJ_EXPORT / 4);
              writeFace(objFile,
                        QuadIdx{boxVerts[1], boxVerts[5], boxVerts[7], boxVerts[3]},
                        boost::none,
                        QuadIdx{normIdx, normIdx, normIdx});
          }

          using CapIdx = std::array<int, KNOB_RESOLUTION_OBJ_EXPORT>;
          CapIdx capNorm;
          capNorm.fill(KNOB_RESOLUTION_OBJ_EXPORT + 1);
          //Write top caps face
          for(int x = 0; x < brick->getSizeX(); ++x)
          {
            for(int y = 0; y < brick->getSizeY(); ++y)
            {
                if(visibleStud[x][y]){
                    auto& knobIndex = studIndices[x][y];
                    writeFace(objFile,
                              qMakePair(knobIndex, index_sequence<KNOB_RESOLUTION_OBJ_EXPORT>()),
                              boost::none,
                              capNorm);
                }
            }
          }

          //Write cylinder faces
          //objFile << "s 1" << std::endl;
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
                                  QuadIdx{knobIndex[i], knobIndex[KNOB_RESOLUTION_OBJ_EXPORT + i], knobIndex[KNOB_RESOLUTION_OBJ_EXPORT + j], knobIndex[j]},
                                  boost::none,
                                  QuadIdx{iNorm, iNorm, jNorm, jNorm});

                    }
                }
            }
          }

          brickIndex++;
          vertexIndex += 8 + studCount*(2*KNOB_RESOLUTION_OBJ_EXPORT);
        }
      }


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

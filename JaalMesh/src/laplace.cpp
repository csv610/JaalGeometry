#include "Mesh.hpp"

using namespace Jaal;


///////////////////////////////////////////////////////////////////////////////
Point3D
LaplaceSmoothing::displace( const Point3D &head, const Point3D &tail, double r )
{
    double d = JMath::length( head, tail);
    if( d < 1.0E-10) return head;

    Vec3D  vec = JMath::unit_vector(head, tail);

    Point3D newpos;
    newpos[0] = tail[0] + d*r*vec[0];
    newpos[1] = tail[1] + d*r*vec[1];
    newpos[2] = tail[2] + d*r*vec[2];
    return newpos;
}

///////////////////////////////////////////////////////////////////////////////

double
LaplaceSmoothing::update_vertex_position(Vertex *vertex, const Point3D &newpos)
{
    double area0, area1, localerror = 0.0;

    JNodeSequence vneighs;
    vertex->getRelations( vneighs );

    Point3D oldpos = vertex->getXYZCoords();

    JFaceSequence vfaces;
    vertex->getRelations( vfaces );

    area0 = 0.0;
    for (size_t i = 0; i < vfaces.size(); i++)
    {
        MeshEntity::idtype id = vfaces[i]->getID();
        area0 += facearea[id];
    }

    // Change the coordinate ( temporary )
    vertex->setXYZCoords(newpos);

    // Check for validity of the sounding patch...
    int some_face_inverted = 0;
    area1 = 0.0;
    for (size_t i = 0; i < vfaces.size(); i++)
    {
        MeshEntity::idtype id = vfaces[i]->getID();
        facearea[id] = vfaces[i]->getArea(); // Recalculate new face area ...
        area1 += facearea[id];
        if (vfaces[i]->concaveAt() >= 0)
        {
            some_face_inverted = 1;
            break;
        }
    }

    // Both the area must be invariant and no face must inverte in order to update the
    // postion.
    if (fabs(area1 - area0) > 1.0E-10 || some_face_inverted == 1)
    {
        vertex->setXYZCoords(oldpos); // Revert back
        localerror = 0.0;
        return localerror;
    }

    // Update the area of the surrounding faces ...
    localerror = max(localerror, fabs(newpos[0] - oldpos[0]));
    localerror = max(localerror, fabs(newpos[1] - oldpos[1]));
    localerror = max(localerror, fabs(newpos[2] - oldpos[2]));
    maxerror = max(maxerror, localerror);

    return localerror;
}

///////////////////////////////////////////////////////////////////////////////

double
LaplaceSmoothing::update_vertex_position(Vertex *vertex)
{
    if (vertex->isConstrained()) return 0.0;

    const Point3D &oldpos = vertex->getXYZCoords();

    JNodeSequence neighs;
    vertex->getRelations( neighs );
    assert(neighs.size());

    Point3D newpos, xyz;
    newpos[0] = 0.0;
    newpos[1] = 0.0;
    newpos[2] = 0.0;

    double weight, sum_weight = 0.0;
    for (size_t i = 0; i < neighs.size(); i++)
    {
        weight = lapweight->get(vertex, neighs[i]);
        xyz = neighs[i]->getXYZCoords();
        newpos[0] += weight * xyz[0];
        newpos[1] += weight * xyz[1];
        newpos[2] += weight * xyz[2];
        sum_weight += weight;
    }

    newpos[0] /= sum_weight;
    newpos[1] /= sum_weight;
    newpos[2] /= sum_weight;

    newpos[0] = (1.0 - lambda) * oldpos[0] + lambda * newpos[0];
    newpos[1] = (1.0 - lambda) * oldpos[1] + lambda * newpos[1];
    newpos[2] = (1.0 - lambda) * oldpos[2] + lambda * newpos[2];

    return update_vertex_position(vertex, newpos);
}

///////////////////////////////////////////////////////////////////////////////
int
LaplaceSmoothing::execute()
{
    return global_smoothing();
}

///////////////////////////////////////////////////////////////////////////////

int
LaplaceSmoothing::global_smoothing()
{
    assert(lapweight);
    mesh->getTopology()->search_boundary();

    int relexist0 = mesh->buildRelations(0, 0);
    int relexist2 = mesh->buildRelations(0, 2);

    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    facearea.resize(numfaces);
    size_t num_inverted_before_smoothing = 0;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setID(i); // Needed because of facearea is an  external array...
        facearea[i] = face->getArea();
        if (face->concaveAt() >= 0) num_inverted_before_smoothing++;
    }

    for (int iter = 0; iter < numIters; iter++)
    {
        maxerror = 0.0;
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *v = mesh->getNodeAt(i);
            if( !v->isRemoved() )
                maxerror = max(maxerror, update_vertex_position(v));
        }
        cout << " iter " << maxerror << endl;

        if (maxerror < 1.0E-06) break;
    }

    if (!relexist0) mesh->clearRelations(0, 0);
    if (!relexist2) mesh->clearRelations(0, 2);

    /*
        size_t num_inverted_after_smoothing = 0;
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = mesh->getFaceAt(i);
            if (face->invertedAt() >= 0) num_inverted_after_smoothing++;
        }
        assert(num_inverted_after_smoothing <= num_inverted_before_smoothing);
    */

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
LaplaceSmoothing::localized_at(const JNodeSequence &vertexQ)
{
    assert(lapweight);

    assert( mesh->getAdjTable(0,0)  );
    assert( mesh->getAdjTable(0,2)  );

    /*
        int relexist0 = mesh->buildRelations(0, 0);
        int relexist2 = mesh->buildRelations(0, 2);
    */

    //How many faces will be affected by changing the coordinates of vertexQ

    JFaceSequence vfaces;
    size_t numfaces = mesh->getSize(2);

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i); // assert( face );
        face->setVisitBit(0);
    }

    for (size_t i = 0; i < vertexQ.size(); i++)
    {
        Vertex *v = vertexQ[i];
        v->getRelations( vfaces );
        for (size_t j = 0; j < vfaces.size(); j++)
            vfaces[j]->setVisitBit(1);
    }

    facearea.resize(numfaces);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setID(i); // Needed because of facearea is an  external array...
        facearea[i] = 0.0;
        if (face->isVisited())
            facearea[i] = face->getArea();
    }

    for (int iter = 0; iter < numIters; iter++)
    {
        maxerror = 0.0;
        for (size_t i = 0; i < vertexQ.size(); i++)
            maxerror = max(maxerror, update_vertex_position(vertexQ[i]));
        if (maxerror < 1.0E-06) break;
    }

    int status = 0;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if (face->isVisited() && face->concaveAt() >= 0)
        {
            status = 1;
            break;
        }
    }

    /*
        if (!relexist0) mesh->clearRelations(0, 0);
        if (!relexist2) mesh->clearRelations(0, 2);
    */

    return status;
}

///////////////////////////////////////////////////////////////////////////////

int
LaplaceSmoothing::convexify()
{
    if (mesh->getTopology()->getElementsType(2) != JFace::QUADRILATERAL) return 1;
    assert( mesh->isPruned() );

    int relexist0 = mesh->buildRelations(0, 0);
    int relexist2 = mesh->buildRelations(0, 2);

    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    facearea.resize(numfaces);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setID(i); // Needed because of facearea is an  external array...
        facearea[i] = face->getArea();
    }

    JNodeSequence fixnodes;
    double r = 0.5;

    size_t ncount = 0;
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = mesh->getFaceAt(iface);
        int pos = face->concaveAt();
        if (pos >= 0)
        {
            Vertex *v0 = face->getNodeAt((pos + 0) % 4);
            if (!v0->isBoundary())
            {
                ncount++;
                fixnodes.clear();
                // Constrained the face nodes ..
                for (int j = 0; j < face->getSize(0); j++)
                {
                    Vertex *vertex = mesh->getNodeAt(j);
                    if (!vertex->isConstrained())
                    {
                        vertex->setConstrainedMark(1);
                        fixnodes.push_back(vertex);
                    }
                }

                Vertex *v1 = face->getNodeAt((pos + 1) % 4);
                Vertex *v2 = face->getNodeAt((pos + 2) % 4);
                Vertex *v3 = face->getNodeAt((pos + 3) % 4);

                Point3D pm ;
                Vertex::mid_point(v1, v3, pm, r);
                Point3D newpos = displace( pm, v2->getXYZCoords(), 1.0001);

                v0->setXYZCoords(newpos);

                // Upate the remaining nodes...
                for (size_t i = 0; i < numnodes; i++)
                    update_vertex_position(mesh->getNodeAt(i));
                // unconstrained the face nodes ..
                for (size_t j = 0; j < fixnodes.size(); j++)
                    fixnodes[j]->setConstrainedMark(0);
            }
        }
    }

    if (ncount) {
        cout << "#of Concave faces " << ncount << endl;
        global_smoothing();
    }

    if (!relexist0) mesh->clearRelations(0, 0);
    if (!relexist2) mesh->clearRelations(0, 2);

    /*
        ncount = mesh->count_concave_faces();

        if (ncount)
        {
            cout << "Warning: Laplacian could not eliminate concave faces " << ncount << endl;
            return 1;
        }
    */
    return 0;
}

////////////////////////////////////////////////////////////////////////////////



/*
double
LaplaceSmoothing::angle_based_update (const LVertex &lnode)
{
  //
  // Presently this works only for 2D case where z = 0.0;
  //
  Vertex *vertex = lnode.apex;
  if (vertex->isConstrained ()) return 0.0;

  map<Vertex*, vector<Vertex*> >::const_iterator it;
  Point3D p0 = vertex->getXYZCoords ();
  double x = p0[0];
  double y = p0[1];
  double xsum = 0.0;
  double ysum = 0.0;
  for (it = lnode.neighs.begin (); it != lnode.neighs.end (); ++it)
    {
      Vertex *v1 = it->first;
      vector<Vertex*> vneighs = it->second;
      assert (vneighs.size () == 2);
      Point3D pj = v1->getXYZCoords ();
      Point3D pp = vneighs[0]->getXYZCoords ();   // Right side of pj
      Point3D pm = vneighs[1]->getXYZCoords ();   // Left  side of  pj
      Vec3D vec0 = JMath::create_vector (p0, pj);
      Vec3D vec1 = JMath::create_vector (pp, pj);
      Vec3D vec2 = JMath::create_vector (pm, pj);
      double a1 = JMath::getVectorAngle (vec0, vec1, ANGLE_IN_RADIANS);
      double b1 = JMath::getVectorAngle (vec0, vec2, ANGLE_IN_RADIANS);
      double t  = 0.5 * (b1 - a1);
      double x0 = pj[0];
      double y0 = pj[1];
      double xn = x0 + (x - x0) * cos (t) - (y - y0) * sin (t);
      double yn = y0 + (x - x0) * sin (t) + (y - y0) * cos (t);
      xsum += xn;
      ysum += xn;
    }
  int nsize = lnode.neighs.size ();
  Point3D newpos;
  newpos[0] = xsum / (double) nsize;
  newpos[1] = ysum / (double) nsize;
  newpos[2] = 0.0;
  return update_vertex_position (vertex, newpos);
}

///////////////////////////////////////////////////////////////////////////////

void
LaplaceSmoothing::build_angle_based_neighbors ()
{
  size_t numnodes = mesh->getSize (0);
  lnodes.clear ();
  lnodes.reserve(numnodes);
  LVertex lnode;
  for (size_t inode = 0; inode < numnodes; inode++)
    {
      Vertex *vertex = mesh->getNodeAt (inode);
      if (!vertex->isBoundary ())
        {
          lnode.apex = vertex;
          lnode.neighs.clear();
          vector<Face*> vfaces = vertex->getRelations2 ();
          for (size_t j = 0; j < vfaces.size (); j++)
            {
              Face *face = vfaces[j];
              int pos = face->queryPosOf (vertex);
              Vertex *v1 = face->getNodeAt ((pos + 1) % 4);
              Vertex *v2 = face->getNodeAt ((pos + 2) % 4);
              Vertex *v3 = face->getNodeAt ((pos + 3) % 4);
              lnode.neighs[v1].push_back (v2);
              lnode.neighs[v2].push_back (v1);
              lnode.neighs[v3].push_back (v2);
              lnode.neighs[v2].push_back (v3);
            }
            lnodes.push_back (lnode);
        }
    }
}
 */


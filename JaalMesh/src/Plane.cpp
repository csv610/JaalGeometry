#include "Plane.hpp"

//////////////////////////////////////////////////////////////////////////////////

JPlane::JPlane()
{
}
//////////////////////////////////////////////////////////////////////////////////

JPlane::JPlane(const JPlane &P)
{
    a = P.a;
    b = P.b;
    c = P.c;
    d = P.d;
}
//////////////////////////////////////////////////////////////////////////////////

JPlane::JPlane(double _a, double _b, double _c, double _d)
{
    a = _a;
    b = _b;
    c = _c;
    d = _d;
}
//////////////////////////////////////////////////////////////////////////////////

JPlane::JPlane(const Vec3D &normal, double _d)
{
    a = normal[0];
    b = normal[1];
    c = normal[2];
    d = _d;
}

//////////////////////////////////////////////////////////////////////////////////

JPlane JPlane::fromPointNormal(const Vec3D &point, const Vec3D &normal)
{
    JPlane Result;
    Vec3D pnormal = JMath::unit_vector(normal);
    double d = -JMath::dot_product(point, pnormal);
    return JPlane(pnormal, d);
}
//////////////////////////////////////////////////////////////////////////////////

JPlane JPlane::fromPointVectors(const Point3D &point, const Vec3D &V1, const Vec3D &V2)
{
    Vec3D normal;
    JMath::cross_product(V1, V2, normal);
    return fromPointNormal(point, normal);
}

//////////////////////////////////////////////////////////////////////////////////

JPlane JPlane::normalize()
{
    JPlane Result;
    double distance = sqrt(a * a + b * b + c * c);
    Result.a = a / distance;
    Result.b = b / distance;
    Result.c = c / distance;
    Result.d = d / distance;
    return Result;
}

JPlane JPlane::fromPoints(const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
    Vec3D vec0, vec1;
    JMath::make_vector(p1,p0, vec0);
    JMath::make_vector(p2,p0, vec1);
    return fromPointVectors(p0, vec0, vec1);
}

Vec3D JPlane::intersectLine(const JLine3D &line) const
{
//   return intersectLine(Line.P0, Line.P0 + Line.D);
}

Vec3D JPlane::intersectLine(const Vec3D &V1, const Vec3D &V2) const
{
    /*
        Vec3D Diff = V1 - V2;
        double Denominator = a * Diff.x + b * Diff.y + c * Diff.z;
        if(Denominator == 0.0f)
        {
            return (V1 + V2) * 0.5f;
        }
        double u = (a * V1.x + b * V1.y + c * V1.z + d) / Denominator;
        return (V1 + u * (V2 - V1));
    */
}

Vec3D JPlane::intersectLine(const Vec3D &V1, const Vec3D &V2, bool &Hit) const
{
    /*
        Hit = true;
        Vec3D Diff = V2 - V1;
        double denominator = a * Diff.x + b * Diff.y + c * Diff.z;
        if(denominator == 0) {Hit = false; return V1;}
        double u = (a * V1.x + b * V1.y + c * V1.z + d) / denominator;

        return (V1 + u * (V2 - V1));
    */
}

double JPlane::intersectLineRatio(const Vec3D &V1, const Vec3D &V2)
{
    /*
        Vec3D Diff = V2 - V1;
        double Denominator = a * Diff.x + b * Diff.y + c * Diff.z;
        if(Denominator == 0.0f)
        {
            return 0.0f;
        }
        return (a * V1.x + b * V1.y + c * V1.z + d) / -Denominator;
    */
}

double JPlane::getSignedDistance(const Vec3D &Pt) const
{
//    return (a * Pt.x + b * Pt.y + c * Pt.z + d);
}

double JPlane::getUnsignedDistance(const Vec3D &Pt) const
{
//    return JMath::Abs(a * Pt.x + b * Pt.y + c * Pt.z + d);
}

Vec3D JPlane::getClosestPoint(const Vec3D &Point)
{
//   return (Point - Normal() * SignedDistance(Point));
}

/*
bool JPlane:: PlanePlaneIntersection(const JPlane &P1, const JPlane &P2, JLine3D &L)
{
    double Denominator = P1.a * P2.b - P1.b * P2.a;
    if(Denominator == 0.0f)
    {
        // this case should be handled by switching axes...
        return false;
    }
    L.P0 = Vec3D((P2.d * P1.b - P1.d * P2.b) / Denominator, (P1.d * P2.a - P2.d * P1.a) / Denominator, 0.0f);

    L.D = Vec3D::Cross(P1.Normal(), P2.Normal());
    if(L.D.Length() == 0.0f)
    {
        return false;
    }
    L.D = Vec3D::Normalize(L.D);
    return true;
}
*/


double JPlane::dot(const JPlane &P, const Vec4D &V)
{
//    return P.a * V.x + P.b * V.y + P.c * V.z + P.d * V.w;
}

double JPlane::dotCoord(const JPlane &P, const Vec3D &V)
{
//    return P.a * V.x + P.b * V.y + P.c * V.z + P.d;
}

double JPlane::dotNormal(const JPlane &P, const Vec3D &V)
{
//    return P.a * V.x + P.b * V.y + P.c * V.z;
}

/*
The remainder of the code in this file implements the function JPlane::FitToPoints.
Part of this is taken from:
http://www.google.com/codesearch?hl=en&q=+best+fit+plane+show:HZpynAfUBdA:AwtEUprnjf8:-KEKP5FuK5I&sa=N&cd=1&ct=rc&cs_p=ftp://figment.csee.usf.edu/pub/segmentation-comparison/USF-segger-code.tar.Z&cs_f=calcplane.c#first
*/

/*    This routine finds the scatter matrix of a number of points equal
**    to "Points_Total".  The x,y and z coordinates of the points are
**    stored in the "X_Coord", "Y_Coord" and "Z_Coord" arrays.  "Centroid"
**    is a 3-element array containing the center of gravity of the set
**    of points.  The scatter matrix will be returned in the 3x3 array
**    called "ScatterMatrix".  The xyz placement in the Scatter Matrix
**    will be returned by the 3 element array "Order", where the index
**    of Order indicates the column (and row) in "ScatterMatrix", and
**    a value of 0 means x, 1 means y, and 2 means z.
*/

void getScatterMatrix( const vector<Vec4D> &Points, const Vec3D &centroid,
                       double scatterMatrix[3][3], int order[3])
{

#ifdef CSV

    int    i,TempI;
    double    TempD;

    /*    To compute the correct scatter matrix, the centroid must be
    **    subtracted from all points.  If the plane is to be forced to pass
    **    through the origin (0,0,0), then the centroid was earlier set
    **    equal to (0,0,0).  This, of course, is NOT the true centroid of
    **    the set of points!  Since the matrix is symmetrical about its
    **    diagonal, one-third of it is redundant and is simply found at
    **    the end.
    */
    for (i=0; i<3; i++) {
        scatterMatrix[i][0]=scatterMatrix[i][1]=scatterMatrix[i][2]=0;
    }

    for(size_t i = 0; i < Points.Length(); i++)
    {
        const Vec4f &P = Points[i];
        Vec3D d = Vec3D(P.x, P.y, P.z) - centroid;
        double weight = P.w;
        scatterMatrix[0][0] += d.x*d.x*weight;
        scatterMatrix[0][1] += d.x*d.y*weight;
        scatterMatrix[0][2] += d.x*d.z*weight;
        scatterMatrix[1][1] += d.y*d.y*weight;
        scatterMatrix[1][2] += d.y*d.z*weight;
        scatterMatrix[2][2] += d.z*d.z*weight;
    }
    scatterMatrix[1][0]=scatterMatrix[0][1];
    scatterMatrix[2][0]=scatterMatrix[0][2];
    scatterMatrix[2][1]=scatterMatrix[1][2];

    /*    Now, perform a sort of "Matrix-sort", whereby all the larger elements
    **    in the matrix are relocated towards the lower-right portion of the
    **    matrix.  This is done as a requisite of the tred2 and tqli algorithms,
    **    for which the scatter matrix is being computed as an input.
    **    "Order" is a 3 element array that will keep track of the xyz order
    **    in the scatterMatrix.
    */
    order[0]=0;        /* Beginning order is x-y-z, as found above */
    order[1]=1;
    order[2]=2;
    if (scatterMatrix[0][0] > scatterMatrix[1][1]) {
        TempD=scatterMatrix[0][0];
        scatterMatrix[0][0]=scatterMatrix[1][1];
        scatterMatrix[1][1]=TempD;
        TempD=scatterMatrix[0][2];
        scatterMatrix[0][2]=scatterMatrix[2][0]=scatterMatrix[1][2];
        scatterMatrix[1][2]=scatterMatrix[2][1]=TempD;
        TempI=order[0];
        order[0]=order[1];
        order[1]=TempI;
    }

    if (scatterMatrix[1][1] > scatterMatrix[2][2]) {
        TempD=scatterMatrix[1][1];
        scatterMatrix[1][1]=scatterMatrix[2][2];
        scatterMatrix[2][2]=TempD;
        TempD=scatterMatrix[0][1];
        scatterMatrix[0][1]=scatterMatrix[1][0]=scatterMatrix[0][2];
        scatterMatrix[0][2]=scatterMatrix[2][0]=TempD;
        TempI=order[1];
        order[1]=order[2];
        order[2]=TempI;
    }
    if (scatterMatrix[0][0] > scatterMatrix[1][1]) {
        TempD=scatterMatrix[0][0];
        scatterMatrix[0][0]=scatterMatrix[1][1];
        scatterMatrix[1][1]=TempD;
        TempD=scatterMatrix[0][2];
        scatterMatrix[0][2]=scatterMatrix[2][0]=scatterMatrix[1][2];
        scatterMatrix[1][2]=scatterMatrix[2][1]=TempD;
        TempI=order[0];
        order[0]=order[1];
        order[1]=TempI;
    }
#endif
}

#ifdef CSV


/*
**    This code is taken from ``Numerical Recipes in C'', 2nd
**    and 3rd editions, by Press, Teukolsky, Vetterling and
**    Flannery, Cambridge University Press, 1992, 1994.
**
*/

/*
**    tred2 Householder reduction of a double, symmetric matrix a[1..n][1..n].
**    On output, a is replaced by the orthogonal matrix q effecting the
**    transformation. d[1..n] returns the diagonal elements of the
**    tridiagonal matrix, and e[1..n] the off-diagonal elements, with
**    e[1]=0.
**
**    For my problem, I only need to handle a 3x3 symmetric matrix,
**    so it can be simplified.
**    Therefore n=3.
**
**    Attention: in the book, the index for array starts from 1,
**    but in C, index should start from zero. so I need to modify it.
**    I think it is very simple to modify, just substract 1 from all the
**    index.
*/

#define    SIGN(a,b)    ((b)<0? -fabs(a):fabs(a))

void tred2(double a[3][3], double d[3], double e[3]) {
    int        l,k,i,j;
    double    scale,hh,h,g,f;

    for(i=3; i>=2; i--)
    {
        l=i-1;
        h=scale=0.0;
        if(l>1)
        {
            for(k=1; k<=l; k++)
                scale+=fabs(a[i-1][k-1]);
            if(scale==0.0)        /* skip transformation */
                e[i-1]=a[i-1][l-1];
            else
            {
                for(k=1; k<=l; k++)
                {
                    a[i-1][k-1]/=scale;    /* use scaled a's for transformation. */
                    h+=a[i-1][k-1]*a[i-1][k-1];    /* form sigma in h. */
                }
                f=a[i-1][l-1];
                g=f>0? -sqrt(h):sqrt(h);
                e[i-1]=scale*g;
                h-=f*g;    /* now h is equation (11.2.4) */
                a[i-1][l-1]=f-g;    /* store u in the ith row of a. */
                f=0.0;
                for(j=1; j<=l; j++)
                {
                    a[j-1][i-1]=a[i-1][j-1]/h; /* store u/H in ith column of a. */
                    g=0.0;    /* form an element of A.u in g */
                    for(k=1; k<=j; k++)
                        g+=a[j-1][k-1]*a[i-1][k-1];
                    for(k=j+1; k<=l; k++)
                        g+=a[k-1][j-1]*a[i-1][k-1];
                    e[j-1]=g/h; /* form element of p in temorarliy unused element of e. */
                    f+=e[j-1]*a[i-1][j-1];
                }
                hh=f/(h+h);    /* form K, equation (11.2.11) */
                for(j=1; j<=l; j++) /* form q and store in e overwriting p. */
                {
                    f=a[i-1][j-1]; /* Note that e[l]=e[i-1] survives */
                    e[j-1]=g=e[j-1]-hh*f;
                    for(k=1; k<=j; k++) /* reduce a, equation (11.2.13) */
                        a[j-1][k-1]-=(f*e[k-1]+g*a[i-1][k-1]);
                }
            }
        }
        else
            e[i-1]=a[i-1][l-1];
        d[i-1]=h;
    }


    /*
    **    For computing eigenvector.
    */
    d[0]=0.0;
    e[0]=0.0;

    for(i=1; i<=3; i++) /* begin accumualting of transfomation matrices */
    {
        l=i-1;
        if(d[i-1]) /* this block skipped when i=1 */
        {
            for(j=1; j<=l; j++)
            {
                g=0.0;
                for(k=1; k<=l; k++) /* use u and u/H stored in a to form P.Q */
                    g+=a[i-1][k-1]*a[k-1][j-1];
                for(k=1; k<=l; k++)
                    a[k-1][j-1]-=g*a[k-1][i-1];
            }
        }
        d[i-1]=a[i-1][i-1];
        a[i-1][i-1]=1.0; /* reset row and column of a to identity matrix for next iteration */
        for(j=1; j<=l; j++)
            a[j-1][i-1]=a[i-1][j-1]=0.0;
    }
}



/*
**    QL algo with implicit shift, to determine the eigenvalues and
**    eigenvectors of a double,symmetric  tridiagonal matrix, or of a double,
**    symmetric matrix previously reduced by algo tred2.
**    On input , d[1..n] contains the diagonal elements of the tridiagonal
**    matrix. On output, it returns the eigenvalues. The vector e[1..n]
**    inputs the subdiagonal elements of the tridiagonal matrix, with e[1]
**    arbitrary. On output e is destroyed. If the eigenvectors of a
**    tridiagonal matrix are desired, the matrix z[1..n][1..n] is input
**    as the identity matrix. If the eigenvectors of a matrix that has
**    been reduced by tred2 are required, then z is input as the matrix
**    output by tred2. In either case, the kth column of z returns the
**    normalized eigenvector corresponding to d[k].
**
*/
void tqli(double d[3], double e[3], double z[3][3]) {
    int        m,l,iter,i,k;
    double    s,r,p,g,f,dd,c,b;

    for(i=2; i<=3; i++)
        e[i-2]=e[i-1];    /* convenient to renumber the elements of e */
    e[2]=0.0;
    for(l=1; l<=3; l++)
    {
        iter=0;
        do
        {
            for(m=l; m<=2; m++)
            {
                /*
                **    Look for a single small subdiagonal element
                **    to split the matrix.
                */
                dd=fabs(d[m-1])+fabs(d[m]);
                if(fabs(e[m-1])+dd == dd)
                    break;
            }
            if(m!=l)
            {
                if(iter++ == 30)
                {
                    printf("\nToo many iterations in TQLI");
                }
                g=(d[l]-d[l-1])/(2.0f*e[l-1]); /* form shift */
                r=sqrt((g*g)+1.0f);
                g=d[m-1]-d[l-1]+e[l-1]/(g+SIGN(r,g)); /* this is dm-ks */
                s=c=1.0;
                p=0.0;
                for(i=m-1; i>=l; i--)
                {
                    /*
                    **    A plane rotation as in the original
                    **    QL, followed by Givens rotations to
                    **    restore tridiagonal form.
                    */
                    f=s*e[i-1];
                    b=c*e[i-1];
                    if(fabs(f) >= fabs(g))
                    {
                        c=g/f;
                        r=sqrt((c*c)+1.0f);
                        e[i]=f*r;
                        c*=(s=1.0f/r);
                    }
                    else
                    {
                        s=f/g;
                        r=sqrt((s*s)+1.0f);
                        e[i]=g*r;
                        s*=(c=1.0f/r);
                    }
                    g=d[i]-p;
                    r=(d[i-1]-g)*s+2.0f*c*b;
                    p=s*r;
                    d[i]=g+p;
                    g=c*r-b;
                    for(k=1; k<=3; k++)
                    {
                        /*
                        **    Form eigenvectors
                        */
                        f=z[k-1][i];
                        z[k-1][i]=s*z[k-1][i-1]+c*f;
                        z[k-1][i-1]=c*z[k-1][i-1]-s*f;
                    }
                }
                d[l-1]=d[l-1]-p;
                e[l-1]=g;
                e[m-1]=0.0f;
            }
        } while(m != l);
    }
}

bool JPlane::fitToPoints(const Vector<Vec3D> &Points, double &ResidualError)
{
    Vec3D Basis1, Basis2;
    Vector<Vec4f> weightedPoints(Points.Length());
    for(size_t i = 0; i < Points.Length(); i++)
    {
        weightedPoints[i] = Vec4f(Points[i], 1.0f);
    }
    double NormalEigenvalue;
    return FitToPoints(weightedPoints, Basis1, Basis2, NormalEigenvalue, ResidualError);
}

bool JPlane:: fitToPoints(const Vector<Vec4f> &Points, Vec3D &Basis1, Vec3D &Basis2, double &NormalEigenvalue, double &ResidualError)
{
    Vec3D centroid, Normal;

    double scatterMatrix[3][3];
    int  order[3];
    double diagonalMatrix[3];
    double offDiagonalMatrix[3];

    // Find centroid
    centroid = Vec3D::Origin;
    double Totalweight = 0.0f;
    for(size_t i = 0; i < Points.Length(); i++)
    {
        Totalweight += Points[i].w;
        centroid += Vec3D(Points[i].x, Points[i].y, Points[i].z) * Points[i].w;
    }
    centroid /= Totalweight;

    // Compute scatter matrix
    getScatterMatrix(Points, centroid, scatterMatrix, order);

    tred2(scatterMatrix,diagonalMatrix,offDiagonalMatrix);
    tqli(diagonalMatrix,offDiagonalMatrix,scatterMatrix);

    /*
    **    Find the smallest eigenvalue first.
    */
    double Min = diagonalMatrix[0];
    double Max = diagonalMatrix[0];

    size_t MinIndex = 0;
    size_t MiddleIndex = 0;
    size_t MaxIndex = 0;
    for(size_t i = 1; i < 3; i++)
    {
        if(diagonalMatrix[i] < Min)
        {
            Min = diagonalMatrix[i];
            MinIndex = i;
        }
        if(diagonalMatrix[i] > Max)
        {
            Max = diagonalMatrix[i];
            MaxIndex = i;
        }
    }
    for(size_t i = 0; i < 3; i++)
    {
        if(MinIndex != i && MaxIndex != i)
        {
            MiddleIndex = i;
        }
    }
    /*
    **    The normal of the plane is the smallest eigenvector.
    */
    for(size_t i = 0; i < 3; i++)
    {
        Normal[order[i]] = scatterMatrix[i][MinIndex];
        Basis1[order[i]] = scatterMatrix[i][MiddleIndex];
        Basis2[order[i]] = scatterMatrix[i][MaxIndex];
    }
    NormalEigenvalue = JMath::Abs(diagonalMatrix[MinIndex]);
    Basis1.SetLength(diagonalMatrix[MiddleIndex]);
    Basis2.SetLength(diagonalMatrix[MaxIndex]);

    if(!Basis1.Valid() || !Basis2.Valid() || !Normal.Valid())
    {
        *this = ConstructFromPointNormal(centroid, Vec3D::eX);
        Basis1 = Vec3D::eY;
        Basis2 = Vec3D::eZ;
    }
    else
    {
        *this = ConstructFromPointNormal(centroid, Normal);
    }

    ResidualError = 0.0f;
    for(size_t i = 0; i < Points.Length(); i++)
    {
        ResidualError += UnsignedDistance(Vec3D(Points[i].x, Points[i].y, Points[i].z));
    }
    ResidualError /= Points.Length();

    return true;
}
#endif

#include <bits/stdc++.h>
#include <windows.h>
#include "include/glut.h"
using namespace std;
#define pi (2*acos(0.0))
#include "bitmap_image.h"
//#include "lib.h"
#define WINDOW_WIDTH 700
#define WINDOW_HEIGHT 700
#define VIEW_ANGLE 90
//.....................mritunjoy...................////
#define MAX_T 100000000.00
#define EPS 1e-6
int recursions;
int pixels;
double angle,near_dist,far_dist,aspect_ratio,w,h ;
double source_factor=1.0;

void clip(double &value, double min_val, double max_val)
{
    if(value < min_val)
        value = min_val;
    if(value > max_val)
        value = max_val;
}
class point
{
public:
    double x,y,z;
    point()
    {

    }
    point(double x1,double y1,double z1)
    {
        x=x1;
        y=y1;
        z=z1;
    }

    point operator * (const point &p)const
    {
        return point(x*p.x,y*p.y,z*p.z) ;
    }
    point operator * (const double &t)const
    {
        return point(x*t,y*t,z*t) ;
    }
    point operator + (const point &p)const
    {
        return point(x+p.x,y+p.y,z+p.z) ;
    }
    point operator - (const point &p)const
    {
        return point(x-p.x,y-p.y,z-p.z) ;
    }
    double value()
    {
        return sqrt(x*x+y*y+z*z) ;
    }


} pos,u,l,r,nu,nl,nr;

double dotProduct( point p1,point p2 )
{
    return p1.x * p2.x + p1.y*p2.y + p1.z * p2.z;
}
point crossProduct(point p1,point p2)
{

    double x = p1.y*p2.z-p1.z*p2.y ;
    double y = p2.x*p1.z-p1.x*p2.z ;
    double z = p1.x*p2.y-p2.x*p1.y ;
    return point(x,y,z) ;
}

point normalizePoint(point p)
{
    double scal = 1.00/sqrt(p.x*p.x+p.y*p.y+p.z*p.z) ;
    return p*scal ;
}

point Rotate(point rotat,point fixed,double rotateAngle)
{
    point perpen;
    perpen.x = fixed.y*rotat.z - fixed.z*rotat.y;
    perpen.y = fixed.z*rotat.x - fixed.x*rotat.z;
    perpen.z = fixed.x*rotat.y - fixed.y*rotat.x;
    point ret  =  rotat*cos(rotateAngle)+perpen*sin(rotateAngle) ;
    return  ret;
}


void drawSquare(double a)
{
    glBegin(GL_QUADS);
    {
        glVertex3f( a, a,0);
        glVertex3f( a,-a,0);
        glVertex3f(-a,-a,0);
        glVertex3f(-a, a,0);
    }
    glEnd();
}

class color
{
public:
    double red,green,blue;
public:
    color()
    {
        red=green=blue=0 ;
    }
    color(double r,double g,double b)
    {
        red=r ;
        green=g;
        blue=b ;
    }
    color operator+(color c)
    {
        double r,g,b;
        r=red+c.red;
        b=blue+c.blue;
        g=green+c.green;

        return color(r,g,b) ;
    }

    color operator*(double c)
    {
        double r,g,b;
        r=red*c;
        b=blue*c;
        g=green*c;

        return color(r,g,b) ;
    }
    color operator*(color c)
    {
        double r,g,b;
        r=red*c.red;
        b=blue*c.blue;
        g=green*c.green;

        return color(r,g,b) ;
    }

};

class Ray
{
public:
    point direction;
    point p;

    Ray()
    {
        p = point( 0, 0, 0 );
        direction = point( 0, 0, 0 );
    }
    Ray(point P0,point dir )
    {
        p = P0;
        direction = dir;

    }
};

class Sphere
{
public:
    color col;
    point cen;
    double rad,amb,diff,spec,reflec,exponent;
public:
    Sphere()
    {
    }
    Sphere(point c,double r,color _col,double _amb,double _diff,double _spec,double _reflec,double _exponent)
    {
        diff=_diff ;
        spec=_spec ;
        reflec = _reflec;
        exponent=_exponent ;
        cen=c ;
        rad = r ;
        col = _col ;
        amb = _amb ;
    }

    void draw()
    {
        glColor3f(col.red, col.green, col.blue);
        glPushMatrix();
        glTranslatef(cen.x, cen.y, cen.z);
        glutSolidSphere(rad, 100, 100);
        glPopMatrix();
    }
    point getNormal(point& P)
    {
        return (P - cen) *rad;
    }

    double intersect(Ray ray)
    {
        point L = ray.p - cen;
        double b = -dotProduct(L,ray.direction );
        double d = (b * b) - dotProduct(L,L) + rad*rad;
        double t= MAX_T;
        if (d > 0)
        {
            d = sqrt(d);
            double  t1 = b - d;
            double t2 = b + d;
            if (t2 > 0)
            {
                if (t1 < 0)
                {
                    t = t2 ;
                }
                else
                {
                    t = t1;
                }
            }
        }
        return t;
    }
};

class Triangle
{
public:
    double amb,diff,spec,reflec,exponent;
    color col ;
    point p1,p2,p3 ;


    Triangle()
    {

    }

    Triangle(point p_1,point p_2,point p_3,color colr,double ambient,double diffuse,double specular,double reflection,double expt)
    {
        p1=p_1 ;
        p2=p_2 ;
        p3=p_3 ;
        col=colr ;
        amb=ambient;
        diff=diffuse;
        spec=specular;
        reflec=reflection;
        exponent=expt;
    }
    void draw()
    {
        glColor3f(col.red,col.green,col.blue) ;
        glBegin(GL_TRIANGLES);
        glVertex3f(p1.x,p1.y,p1.z);
        glVertex3f(p2.x,p2.y,p2.z);
        glVertex3f(p3.x,p3.y,p3.z);
        glEnd();
    }
    point getNormal(point p)
    {
        point aa = p2 - p1;
        point bb = p3 - p1;
        point N = crossProduct(aa,bb);
        return normalizePoint(N) ;
    }
    int insideTrian(point P)
    {
        point a,b,c;
        a=p1;
        b=p2;
        c=p3;
        point W = P - a;
        point U = b - a;
        point V = c - a;
        point vCrosU = crossProduct(V,U);
        point vCrosW = crossProduct(V,W);
        if (dotProduct(vCrosW,vCrosU) < 0) //r
            return false;

        point uCrosV = crossProduct(U,V);
        point uCrosW = crossProduct(U,W);

        // t
        if (dotProduct(uCrosW,uCrosV) < 0)
            return false;

        // r> 0, t>0.
        // r+t <= 1, r,t <= 1
        double deno = uCrosV.value();
        double t = uCrosW.value() / deno;
        double r = vCrosW.value()/ deno;
        return(r+t<=1.00) ;
    }

    double intersect(Ray &ray)
    {/*
        point N = getNormal(p1) ;
        double d = dotProduct(N,ray.direction) ;
        if(fabs(d)<EPS)
            return MAX_T ;
        double t = (-(dotProduct( N,ray.p)) + dotProduct(N,p1)) / d;
        point intp = ray.p+ray.direction*t;
        if(insideTrian(intp))
        {
            return t ;
        }
        return MAX_T ;
        */

        point normal = getNormal(point(0, 0, 0));

        double t;
        bool flag = false;
        double denom = dotProduct(normal, ray.direction);

        if (denom < 0.0)
        {
            normal = normal * -1.0;
            denom = dotProduct(normal, ray.direction);
        }

        // checking whether the ray intersects with the triangle's plane
        if (abs(denom) < EPS)
            return MAX_T;

        t = (dotProduct(normal, p1 - ray.p)) / denom;
        if (t >= 0)
            flag = true;

        if (!flag)
            return MAX_T;

        // checking whether the point of intersection is inside the triangle
        bool b1, b2, b3;
        point intersectionPoint = ray.p + ray.direction * t;

        point N = crossProduct(p2 - p1, p3 - p1);
        double area2 = N.value();

        point C, edge1, edge2;

        edge1 = p2 - p1;
        edge2 = intersectionPoint - p1;
        C = crossProduct(edge1, edge2);
        if (dotProduct(N, C) < 0)
            return MAX_T;

        edge1 = p3 - p2;
        edge2 = intersectionPoint - p2;
        C = crossProduct(edge1, edge2);
        if (dotProduct(N, C) < 0)
            return MAX_T;

        edge1 = p1 - p3;
        edge2 = intersectionPoint - p3;
        C = crossProduct(edge1, edge2);
        if (dotProduct(N, C) < 0)
            return MAX_T;

        return t;
    }
};

class Square
{
public:
    int i,j,n,m ;
    double reflec,diff,spec,amb,exponent ;
    double side ;
    color col = color(0,0,0) ;
    point Dir,N ;

    Square(int ii,int jj,int nn,int mm,point normal,point dir,double len,double expt)
    {
        reflec = 0.3 ;
        diff=1.00 ;
        spec =0 ;
        exponent=expt;
        i = ii ;
        j=jj ;
        n = nn ;
        m=mm ;
        N = normal ;
        Dir = dir ;
        side=len ;
        if (i%2==0 && j%2==0)
            col =  color( 0.0, 0.0, 0.0)  ;
        else if(i%2==0 && j%2!=0)
            col =  color( 1.0, 1.0, 1.0) ;
        else if(i%2!=0 &&  j%2==0)
            col =  color( 1.0, 1.0, 1.0) ;
        else if (i%2!=0 && j%2!=0)
            col =  color( 0.0, 0.0, 0.0) ;
    }
    Square(point mid,color col,double Len,double a,double d,double s,double r,double expt)
    {
        Dir = mid ;
        N = point(0,0,mid.z) ;
        side=Len ;
        diff = d;
        reflec=r ;
        spec=s ;
        exponent=expt;
        col = col ;
        amb=a ;
    }
    point getNormal(point &P)
    {
        return normalizePoint(point(0,0,1)) ;
    }
    void draw()
    {
        glPushMatrix();
        {
            glTranslatef(Dir.x,Dir.y,Dir.z);
            glColor3f(col.red,col.green,col.blue) ;
            drawSquare(side/2);
        }
        glPopMatrix() ;
    }

    double intersect( Ray& ray)
    {
        double d = dotProduct(N,ray.direction) ;
        point intp;
        if (d != 0)
        {
            double dist = (-( dotProduct(N,ray.p)) + dotProduct(N,Dir)) / d;
            if (dist > 0)
            {
                intp =ray.p + ray.direction * dist;

                int ff=intp.x>=Dir.x-1 && intp.x<=Dir.x+side && intp.y>=Dir.y-1 && intp.y<=Dir.y+side &&
                       intp.z>=Dir.z-1 && intp.z<=Dir.z+side;
                if(ff)
                {
                    return dist;

                }
            }
        }
        return MAX_T;
    }

};

class Light
{
public:
    point cent ;
    color col ;
    double diff ;
    double rad,spec ;
    Light(double x,double y,double z)
    {
        cent = point(x,y,z) ;
        col  = color(1.0,1.0,1.0) ;
        rad = 3.00 ;
        diff=1 ;
        spec=1 ;
    }
    Light()
    {

    }

    void draw()
    {
        glPushMatrix();
        {
            glTranslatef(cent.x,cent.y,cent.z );
            glColor3f(col.red,col.green,col.blue );
            glutSolidSphere(rad,5,5);
        }
        glPopMatrix();
    }
    point getNormal(point& p)
    {
        return (p - cent) *rad;
    }
    double intersect(Ray ray)
    {
        point L = ray.p - cent;
        double b = -dotProduct(L,ray.direction );
        double d = (b * b) - dotProduct(L,L) + rad*rad;
        double t= MAX_T;
        if (d > 0)
        {
            d = sqrt(d);
            double  t1 = b - d;
            double t2 = b + d;
            if (t2 > 0)
            {
                if (t1 < 0)
                {
                    t = t2 ;
                }
                else
                {
                    t = t1;
                }
            }
        }
        return t;
    }
};

vector<Sphere>sphs ;
vector<Triangle>trians ;
vector<Square>sqs ;
vector<Light>Lights ;

point get_reflected_ray_direction(Ray &ray, point normal)
{
    point out_dir = ray.direction - normal * 2.0 * dotProduct(ray.direction, normal);
    return normalizePoint(out_dir);
}
double getReflection(int i,int j)
{
    if(i==1)
    {
        return sqs[j].reflec ;
    }
    else if(i==2)
    {
        return sphs[j].reflec ;
    }

    else if(i==3)
    {
        return 0 ;
    }

    else if(i==4)
    {
        return trians[j].reflec ;
    }

    return 0 ;
}
double getspec(int i,int j)
{
    if(i==1)
    {
        return sqs[j].spec ;
    }
    else if(i==2)
    {
        return sphs[j].spec ;
    }

    else if(i==3)
    {
        return Lights[j].spec ;
    }

    else if(i==4)
    {
        return trians[j].spec ;
    }
    return 0 ;
}
color getColor(int i,int j)
{
    if(i==1)
    {
        return sqs[j].col;
    }
    else if(i==2)
    {
        return sphs[j].col ;
    }
    else if(i==3)
    {
        return Lights[j].col ;
    }
    else if(i==4)
    {
        return trians[j].col ;
    }
    return color(0,0,0) ;
}

double getambient(int i,int j)
{
    if(i==1)
    {
        return sqs[j].amb ;
    }
    else if(i==2)
    {
        return sphs[j].amb ;
    }
    else if(i==3)
    {
        return 0 ;
    }
    else if(i==4)
    {
        return trians[j].amb ;
    }
    return 0 ;
}

double getdiffuse(int i,int j)
{
    if(i==1)
    {
        return sqs[j].diff ;
    }
    else if(i==2)
    {
        return sphs[j].diff ;
    }
    else if(i==3)
    {
        return Lights[j].diff ;
    }
    else if(i==4)
    {
        return trians[j].diff ;
    }
    return 0 ;
}

double getexponent(int i,int j)
{
    if(i==1)
    {
        return sqs[j].exponent ;
    }
    else if(i==2)
    {
        return sphs[j].exponent ;
    }
    else if(i==3)
    {
        return 1 ;
    }
    else if(i==4)
    {
        return trians[j].exponent ;
    }
    return 0 ;
}

point getNormal(int i,int j,point P)
{
    if(i==1)
    {
        return sqs[j].getNormal(P) ;
    }
    else if(i==2)
    {
        return sphs[j].getNormal(P) ;
    }
    else if(i==3)
    {
        return Lights[j].getNormal(P) ;
    }
    else if(i==4)
    {
        return trians[j].getNormal(P) ;
    }
    return point(0,0,0) ;
}
class nearest
{
public:
    int I;
    int J;
    double min_t;
    nearest()
    {
    }
    nearest(int ii,int jj,double tt)
    {
        I=ii;
        J=jj;
        min_t=tt;
    }
};
nearest getNearest(Ray& ray,int Level,color& col)
{
    if(Level>recursions)
    {
        col = color(0,0,0) ;
    }
    double min_t = MAX_T ;
    int I=-1,J=-1 ;
    //intersect square
    for(int i=0; i<(int)sqs.size(); i++)
    {
        double t = sqs[i].intersect(ray) ;
        if(t<min_t)
        {
            I = 1 ;
            J = i ;
            min_t=t ;
        }
    }
    ///         Now with the spheres
    for(int i=0; i<(int)sphs.size(); i++)
    {
        double  t = sphs[i].intersect(ray) ;
        if(t<min_t)
        {
            I = 2 ;
            J = i ;
            min_t=t ;
        }
    }

    for(int i=0; i<(int)trians.size(); i++)
    {
        double  t = trians[i].intersect(ray) ;
        if(t<min_t)
        {
            I = 4 ;
            J = i ;
            min_t=t ;
        }
    }
    //
    /*
    for(int i=0; i<(int)Lights.size(); i++)
    {
        double  t = Lights[i].intersect(ray) ;
        if(t<min_t)
        {
            //col = color(1,1,1) ;
            //return  ;
            I = 3 ;
            J = i ;
            min_t=t ;
        }
    }
    */
    return nearest(I,J,min_t);
}
void fill_color(Ray& ray,int ii,int jj,double t,int Level,color& col)
{
    /// determine the color at the point
    double min_t=t;
    int I=ii;
    int J=jj;
    ray.direction=normalizePoint(ray.direction) ;
    point intp = ray.p+ray.direction*min_t ;
    point normal = getNormal(I,J,intp);
    point reflection = get_reflected_ray_direction(ray, normal);
    /// Trace light
    //col=col+getColor(I,J)*0.4 ;
    for(int i=0; i<(int)Lights.size(); i++)
    {
        double ambient =getambient(I,J), lambert = 0.0, phong = 0.0;
        point dir =normalizePoint(Lights[i].cent - intp);
        point start = intp + dir * EPS;

        Ray L(start, dir);
        point R = normalizePoint(normal * (dotProduct(L.direction, normal) * 2.0) - L.direction);
        point V = normalizePoint(intp * -1.0);

        bool flag = false;
        for(int k=0; k<(int)sqs.size(); k++)
        {
            if(sqs[k].intersect(L)>0 && sqs[k].intersect(L)!=MAX_T)
            {
                flag=true ;
                break ;
            }
        }
        for(int k=0; k<(int)sphs.size(); k++)
        {
            if(sphs[k].intersect(L)>0 && sphs[k].intersect(L)!=MAX_T)
            {
                flag=true ;
                break ;
            }
        }

        for(int k=0; k<(int)trians.size(); k++)
        {
            if(trians[k].intersect(L)>0 && trians[k].intersect(L)!=MAX_T)
            {
                flag=true ;
                break ;
            }
        }

        if (!flag)
        {
            lambert = source_factor * getdiffuse(I,J) * dotProduct(L.direction, normal);
            phong = getspec(I,J) * pow(dotProduct(R, V), getexponent(I,J));

            clip(lambert, 0.0, 1.0);
            clip(phong, 0.0, 1.0);
        }
        col.red += ((ambient + lambert + phong) * getColor(I,J).red);
        col.green += ((ambient + lambert + phong) * getColor(I,J).green);
        col.blue += ((ambient + lambert + phong) * getColor(I,J).blue);

        if (Level < recursions)
        {
            point start = intp + reflection * EPS;

            Ray reflectionRay(start, reflection);
            color reflected_color;

            nearest neer = getNearest(reflectionRay,Level,reflected_color);
            if (neer.I != -1 || neer.J!=-1)
            {
                fill_color(reflectionRay,neer.I,neer.J,neer.min_t, Level + 1,reflected_color);
                col.red += reflected_color.red * getReflection(neer.I,neer.J);
                col.green += reflected_color.green * getReflection(neer.I,neer.J);
                col.blue+= reflected_color.blue * getReflection(neer.I,neer.J);
            }

            clip(col.red, 0.0, 1.0);
            clip(col.green, 0.0, 1.0);
            clip(col.blue, 0.0, 1.0);

        }
    }
}



void image_capture()
{

    double plane_distance = (WINDOW_HEIGHT/2)/tan(VIEW_ANGLE*pi/360);
    point top_left = pos + (l * plane_distance - r * (WINDOW_WIDTH/2) + u * (WINDOW_HEIGHT/2));

    cout << "Plane distance : " << plane_distance << endl;
    cout << "top_left : "<< top_left.x<<"  "<<top_left.y<<" "<<top_left.z<<endl;
    cout << "Saving...";

    double du = (WINDOW_WIDTH*1.0) / pixels;
    double dv = (WINDOW_HEIGHT*1.0) / pixels;
    bitmap_image image(pixels,pixels);

    for (int i = 0; i < pixels; i++)
    {
        for (int j = 0; j < pixels; j++)
        {

            point direction_to_top_left = top_left + r*i*du - u*j*dv;

            Ray ray(pos, direction_to_top_left - pos);
            ray.direction=normalizePoint(ray.direction);

            color draw_color;
            nearest nearr=getNearest(ray,1,draw_color);
            fill_color(ray,nearr.I,nearr.J,nearr.min_t,1,draw_color);
            image.set_pixel(i,j,draw_color.red*255.00,draw_color.green*255.00,draw_color.blue*255.00);
        }
        cout<<"Processed Row::"<<i<<endl ;
    }
    cout<<"Saving the image"<<endl ;
    image.save_image("out.bmp");
    cout<<"finish"<<endl;


}

void keyboardListener(unsigned char key, int x,int y)
{
    switch(key)
    {

    case '1'://rotate/look left
        l = Rotate(l,u,(M_PI/180)*3) ;
        r = Rotate(r,u,(M_PI/180)*3) ;
        break;
    case '2'://rotate/look right
        l = Rotate(l,u,-(M_PI/180)*3) ;
        r = Rotate(r,u,-(M_PI/180)*3) ;
        break;
    case '3': //look up
        l = Rotate(l,r,(M_PI/180)*3) ;
        u = Rotate(u,r,(M_PI/180)*3) ;
        break ;
    case '4': // look down
        l = Rotate(l,r,-(M_PI/180)*3) ;
        u = Rotate(u,r,-(M_PI/180)*3) ;
        break ;
    case '5'://tilt clockwise
        u = Rotate(u,l,(-M_PI/180)*3) ;
        r = Rotate(r,l,(-M_PI/180)*3) ;
        break;
    case '6'://tilt anticlockwise
        u = Rotate(u,l,(M_PI/180)*3) ;
        r = Rotate(r,l,(M_PI/180)*3) ;
        break;
    case '0':
        image_capture();
        break;

    default:
        break;
    }
}


void specialKeyListener(int key, int x,int y)
{
    switch(key)
    {
    case GLUT_KEY_UP:		// up arrow key
        //cameraHeight += 3.0;
        pos = pos+l*3 ;
        break;
    case GLUT_KEY_DOWN:		//down arrow key
        //cameraHeight -= 3.0;
        pos = pos-l*3 ;
        break;

    case GLUT_KEY_RIGHT:
        //cameraAngle += 0.03;
        pos = pos+r*3 ;
        break;
    case GLUT_KEY_LEFT:
        //cameraAngle -= 0.03;
        pos = pos-r*3 ;
        break;

    case GLUT_KEY_PAGE_UP:
        pos = pos+u*3 ;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos = pos-u*3 ;
        break;
    default:
        break;
    }
}


void display()
{

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    gluLookAt(pos.x,pos.y,pos.z, pos.x+l.x,pos.y+l.y,pos.z+l.z,u.x,u.y,u.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    //drawAxes();
    int n=sqs.size();
    for(int i=0; i<n; i++)
    {
        sqs[i].draw() ;
    }
    n=trians.size();
    for(int i=0; i<n; i++)
    {
        trians[i].draw() ;
    }
    n=Lights.size();
    for(int i=0; i<n; i++)
    {
        Lights[i].draw() ;
    }
    n=sphs.size();
    for(int i=0; i<n; i++)
    {
        sphs[i].draw() ;
    }

    glutSwapBuffers();
}


void animate()
{
    //angle+=0.05;
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void InitData()
{
    Sphere sph ;
    int n_objects ;
    double rad,amb,diff,spec,reflec,exponent,x,y,z,red,green,blue,height,base;
    string ob;
    freopen("scene.txt", "r", stdin);
    cin >> recursions>>pixels>>n_objects ;
    cout<<"recursions: "<<recursions<<" pixels: "<<pixels<<" number of objects: "<<n_objects<<endl;
    while(n_objects>0)
    {
        cin>>ob;
        if(ob=="sphere")
        {
            cin>>x>>y>>z>>rad>>red>>green>>blue>>amb>>diff>>spec>>reflec>>exponent;
            //cout<<exponent;
            sph = Sphere(point(x,y,z),rad,color(red,green,blue),amb,diff,spec,reflec,exponent) ;
            sphs.push_back(sph) ;
        }
        else if(ob=="pyramid")
        {
            cin>>x>>y>>z>>base>>height>>red>>green>>blue>>amb>>diff>>spec>>reflec>>exponent;
            Triangle trian  ;
            point m = point(x+base/2.0,y+base/2.0,z) ;
            point top = m+point(0,0,1)*height ;
            point a = point(x,y,z) ;
            point b = point(x+base,y,z) ;
            point c = point(x+base,y+base,z) ;
            point d = point(x,y+base,z) ;

            trian = Triangle(top,a,b,color(red,green,blue),amb,diff,spec,reflec,exponent) ;
            trians.push_back(trian) ;
            trian = Triangle(top,b,c,color(red,green,blue),amb,diff,spec,reflec,exponent) ;
            trians.push_back(trian) ;
            trian = Triangle(top,c,d,color(red,green,blue),amb,diff,spec,reflec,exponent) ;
            trians.push_back(trian) ;
            trian = Triangle(top,d,a,color(red,green,blue),amb,diff,spec,reflec,exponent) ;
            trians.push_back(trian) ;

            Square sq = Square((a+c),color(red,green,blue),base,amb,diff,spec,reflec,exponent) ;
            sqs.push_back(sq) ;
        }
        n_objects--;
    }
    int n_lights;
    cin>>n_lights;
    cout<<n_lights<<endl ;
    for(int i=0; i<n_lights; i++)
    {
        cin>>x>>y>>z ;
        Lights.push_back(Light(x,y,z)) ;
    }

    for(int i=0,m=-10; i<20; i++,m++)
    {
        for(int j=0,n=-10; j<20; j++,n++)
        {
            Square sq = Square(i,j,n,m,point( 0, 0, 1 ), point(m*30,n*30,0 ),30,1) ;
            sqs.push_back(sq) ;
        }
    }
}
void init()
{
    //codes for initialization
    u.x=0;
    u.y=0;
    u.z=1;
    r.x=1 ;
    r.y=0;
    r.z =0 ;
    l.x=0 ;
    l.y=1 ;
    l.z=0 ;
    pos.x=0;
    pos.y=-200;
    pos.z=30;
    InitData();

    //clear the screen
    glClearColor(0, 0, 0, 0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(90, 1, 1, 1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color
    glutCreateWindow("Ray Tracing");
    init();
    glEnable(GL_DEPTH_TEST);	//enable Depth Testing
    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMainLoop();		//The main loop of OpenGL
    return 0;
}

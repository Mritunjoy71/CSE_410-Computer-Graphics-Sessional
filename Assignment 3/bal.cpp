#include<bits/stdc++.h>
using namespace std;
#include <windows.h>
#include <GL/glut.h>
#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
int recursions;
int pixels;

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




void drawAxes()
{

    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    {
        glVertex3f( 100,0,0);
        glVertex3f(-100,0,0);

        glVertex3f(0,-100,0);
        glVertex3f(0, 100,0);

        glVertex3f(0,0, 100);
        glVertex3f(0,0,-100);
    }
    glEnd();

}


void drawGrid()
{
    int i;
    if(drawgrid==1)
    {
        glColor3f(0.6, 0.6, 0.6);	//grey
        glBegin(GL_LINES);
        {
            for(i=-8; i<=8; i++)
            {

                if(i==0)
                    continue;	//SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i*10, -90, 0);
                glVertex3f(i*10,  90, 0);

                //lines parallel to X-axis
                glVertex3f(-90, i*10, 0);
                glVertex3f( 90, i*10, 0);
            }
        }
        glEnd();
    }
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
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

    color operator*(color c)
    {
        double r,g,b;
        r=red*c.red;
        b=blue*c.blue;
        g=green*c.green;

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



};




void draw_1_8thSphere(double radius,int slices,int stacks)
{
    struct point points[100][100];
    int i,j;
    double h,r;
    //generate points
    for(i=0; i<=stacks; i++)
    {
        h=radius*sin(((double)i/(double)stacks)*(pi/2));
        r=radius*cos(((double)i/(double)stacks)*(pi/2));
        for(j=0; j<=slices; j++)
        {
            points[i][j].x=r*cos(((double)j/(double)slices)*pi/2);
            points[i][j].y=r*sin(((double)j/(double)slices)*pi/2);
            points[i][j].z=h;
        }
    }
    //draw quads using generated points
    for(i=0; i<stacks; i++)
    {
        glColor3f(1, 0, 0);
        for(j=0; j<slices; j++)
        {
            glBegin(GL_QUADS);
            {
                //upper hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
            }
            glEnd();
        }
    }
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
        pos = pos+l ;
        break;
    case GLUT_KEY_DOWN:		//down arrow key
        //cameraHeight -= 3.0;
        pos = pos-l ;
        break;

    case GLUT_KEY_RIGHT:
        //cameraAngle += 0.03;
        pos = pos+r ;
        break;
    case GLUT_KEY_LEFT:
        //cameraAngle -= 0.03;
        pos = pos-r ;
        break;

    case GLUT_KEY_PAGE_UP:
        pos = pos+u ;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos = pos-u ;
        break;

    case GLUT_KEY_INSERT:
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

    drawAxes();
    //drawGrid();

    //glColor3f(1,0,0);

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate()
{
    //angle+=0.05;
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init()
{


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



    //clear the screen
    glClearColor(0,0,0,0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(80,	1,	1,	1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv)
{
    freopen("description.txt","r",stdin);
    Sphere sph ;
    int n_objects ;
    double rad,amb,diff,spec,reflec,exponent,x,y,z,red,green,blue,high,wid ;
    cin>>rad;
    cout<<rad;


    init();
    glutInit(&argc,argv);
    glutInitWindowSize(500,500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("Ray Tracing");


    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    // glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}

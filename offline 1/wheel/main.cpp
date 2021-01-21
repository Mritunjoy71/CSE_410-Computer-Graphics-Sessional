#include<stdio.h>
#include<stdlib.h>
#include<math.h>
using namespace std;


#include <windows.h>
//#include <GL/glut.h>
#define pi (2*acos(0.0))
#include "glut.h"

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double W_angle;

struct point
{
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


} pos,u,l,r,nu,nl,nr,Centre,V_dir,V_fixed;

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
    if(drawaxes==1)
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



void drawWheel(double x, double y, double radius)
{
    //glTranslated(Centre.x,Centre.y,Centre.z);
    int i;
    int segments=360;
    struct point points[360][2];
    //generate points
    for(i=0; i<segments; i++)
    {
        points[i][0].x=x+radius*cos(((double)i/(double)segments)*2*pi);
        points[i][0].y=y+radius*sin(((double)i/(double)segments)*2*pi);
        points[i][0].z=0;

        points[i][1].x=x+radius*cos(((double)i/(double)segments)*2*pi);
        points[i][1].y=y+radius*sin(((double)i/(double)segments)*2*pi);
        points[i][1].z=5;
    }


    for(int i=0; i<segments; i++)
    {
        glBegin(GL_QUADS);
        {
            glVertex3f(points[i][0].x,points[i][0].y,points[i][0].z);
            glVertex3f(points[i][1].x,points[i][1].y,points[i][1].z);
            glVertex3f(points[i+1][1].x,points[i+1][1].y,points[i+1][1].z);
            glVertex3f(points[i+1][0].x,points[i+1][0].y,points[i+1][0].z);
        }
        glEnd();
    }
    glBegin(GL_QUADS);
    {
        glVertex3f(points[0][0].x,points[0][0].y,points[0][0].z);
        glVertex3f(points[segments/2][1].x,points[segments/2][1].y,points[segments/2][1].z);
        glVertex3f(points[segments/2][0].x,points[segments/2][0].y,points[segments/2][0].z);
        glVertex3f(points[0][1].x,points[0][1].y,points[0][1].z);
    }
    glEnd();

    glBegin(GL_QUADS);
    {
        glVertex3f(points[segments/4][0].x,points[segments/4][0].y,points[segments/4][0].z);
        glVertex3f(points[3*segments/4][1].x,points[3*segments/4][1].y,points[3*segments/4][1].z);
        glVertex3f(points[3*segments/4][0].x,points[3*segments/4][0].y,points[3*segments/4][0].z);
        glVertex3f(points[segments/4][1].x,points[segments/4][1].y,points[segments/4][1].z);
    }
    glEnd();

     //glTranslated(Centre.x,Centre.y,Centre.z);



}

void keyboardListener(unsigned char key, int x,int y)
{
    switch(key)
    {
    case 'w':
        Centre.x =Centre.x+V_dir.x*2 ;
        break ;
    case 's':
        Centre.x =Centre.x-V_dir.x*2 ;
        break ;
    case 'a':
        W_angle++ ;
        V_dir =Rotate(V_dir,V_fixed,W_angle*(M_PI/180.0)) ;
        break ;
    case 'd':
        W_angle-- ;
        V_dir =Rotate(V_dir,V_fixed,W_angle*(M_PI/180.0)) ;
        break ;

    default:
        break;
    }
}


void specialKeyListener(int key, int x,int y)
{
    switch(key)
    {
    case GLUT_KEY_UP:		// up arrow key
        cameraHeight += 3.0;
        //pos = pos+l ;
        break;
    case GLUT_KEY_DOWN:		//down arrow key
        cameraHeight -= 3.0;
        //pos = pos-l ;
        break;

    case GLUT_KEY_RIGHT:
        cameraAngle += 0.03;
        //pos = pos+r ;
        break;
    case GLUT_KEY_LEFT:
        cameraAngle -= 0.03;
        //pos = pos-r ;
        break;

    case GLUT_KEY_PAGE_UP:
        //pos = pos+u ;
        break;
    case GLUT_KEY_PAGE_DOWN:
        //pos = pos-u ;
        break;

    case GLUT_KEY_INSERT:
        break;

    case GLUT_KEY_HOME:
        break;
    case GLUT_KEY_END:
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
    gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    //gluLookAt(pos.x,pos.y,pos.z+50, pos.x+l.x,pos.y+l.y,pos.z+l.z+50,u.x,u.y,u.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

    //drawSphere(30,24,20);




    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    //glTranslated(Centre.x,Centre.y,Centre.z);
    //glTranslated(Centre.x,Centre.y,0) ;
    glRotated(W_angle,0,0,1) ;
    glTranslated(0,0,15) ;
    glRotated(90,1,0,0) ;
    drawWheel(Centre.x,Centre.y,15) ;
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
    //codes for initialization
    drawgrid=1;
    drawaxes=1;
    cameraHeight=150.0;
    cameraAngle=1.0;
    angle=0;

    u.x=0;
    u.y=0;
    u.z=1;

    r.x=(-1.00)/sqrt(2) ;
    r.y=(1.00)/sqrt(2) ;
    r.z =0 ;

    l.x=(-1.00)/sqrt(2) ;
    l.y=(-1.00)/sqrt(2) ;
    l.z=0 ;

    pos.x=100;
    pos.y=100;
    pos.z=0;

    nl = normalizePoint(l) ;
    nu = normalizePoint(u) ;
    nr = normalizePoint(r) ;
    //clear the screen
    glClearColor(0,0,0,0);

    Centre.x=0;
    Centre.y=0;
    Centre.z =0;
    W_angle= 0 ;
    V_dir.x=1 ;
    V_dir.y=0;
    V_dir.z=0 ;
    V_fixed.x=0;
    V_fixed.y=1 ;
    V_fixed.z=0 ;

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
    glutInit(&argc,argv);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("Wheel");

    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
  //  glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}

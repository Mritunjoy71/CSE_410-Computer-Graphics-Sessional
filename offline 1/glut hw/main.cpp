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
double an;
double side;
double max_side;

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





void drawSquares()
{
    glColor3f(1, 1,0);

    //glTranslated(0, 0, max_side);
    drawSquare(side);

    glPushMatrix();
    {
        glColor3f(0, 0, 1);
        glTranslated(max_side,0,0);
        drawSquare(side);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glColor3f(0, 1, 0);
        glTranslated(-max_side,0,0);
        drawSquare(side);
    }
    glPopMatrix();
    glPushMatrix();
    {
        glColor3f(0, 1, 1);
        glTranslated(0,max_side,0);
        drawSquare(side);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glColor3f(1, 0, 0);
        glTranslated(0,-max_side,0);
        drawSquare(side);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glColor3f(1, 0, 1);
        glTranslated(2*max_side,0,0);
        drawSquare(side);
    }
    glPopMatrix();



}



void keyboardListener(unsigned char key, int x,int y)
{
    switch(key)
    {

    case '1':
        an=an+M_PI;
        break;

    case '2':

        break;

    case '3':

        break;

    case '4':

        break;

    case '5':

        break;

    case '6':

        break;

    case '7':

        break;

    case '8':

        break;

    case '9':

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
        pos = pos+u ;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos = pos-u ;
        break;

    case GLUT_KEY_INSERT:
        break;

    case GLUT_KEY_HOME:
        side -= 0.2;
        side = max(side, 0.0);
        break;
    case GLUT_KEY_END:
        side += 0.2;
        side = min(side,max_side);
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
    drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

    //drawOneEighthSphere(30,24,20);
    //drawSpheres();//8 ta
    glRotated(an,0,1,0);
    drawSquares();//6 ta
    //drawCylinders();//12 ta
    //draw_1_4thCylinder(10,10,60,60);




    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate()
{
    angle+=0.05;
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init()
{
    //codes for initialization
    max_side=20;
    side=10;
    drawgrid=0;
    drawaxes=1;
    cameraHeight=150.0;
    cameraAngle=1.0;
    angle=0;

    an=0;

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
    glutInit(&argc,argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("fully controllable camera sphere to/from cube");

    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    // glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}

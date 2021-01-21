#include <bits/stdc++.h>
using namespace std;
#define PI acos(-1.0)
stack<pair<double **, bool> > Stack;
double eyex, eyey, eyez, lookx, looky, lookz, upx, upy, upz, fovy, aspectratio, near, far,ang,ax,ay,az;
double p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, sx, sy, sz, tx, ty, tz, lx, ly, lz, rx, ry, rz, ux, uy, uz, fovx, t, r;
ofstream stage1;
ofstream stage2;
ofstream stage3;
string com_line;

double **T, **R, **V, **P;

double **identity()
{
    //cout<<"in identity function\n";
    double **I = new double *[4];
    for (int i = 0; i < 4; i++)
    {
        I[i] = new double[4];
    }

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
                I[i][j] = 1.0;
            else
                I[i][j] = 0.0;
        }
    }
    return I;
}

double **mult(double **mat1, double **mat2)
{
    double **mul = identity();
    for (int m = 0; m < 4; m++)
    {
        mul[m][m] = 0.0;
    }
    int i, j, k;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            for (k = 0; k < 4; k++)
                mul[i][j] += mat1[i][k] * mat2[k][j];
        }
    }
    return mul;
}

void push()
{
    double **push_mat = Stack.top().first;
    Stack.push(make_pair(push_mat,true));
}
void pop()
{
    while(Stack.top().second == false)
    {
        Stack.pop();
    }
    Stack.pop();
}

void transformation(double** mat)
{
    double** top_mat = Stack.top().first;
    double** new_mat = mult(top_mat, mat);

    Stack.push(make_pair(new_mat, false));
}

int main()
{

    freopen("scene.txt", "r", stdin);
    cout << fixed << showpoint << setprecision(7);

    stage1 << fixed << showpoint << setprecision(7);
    stage2 << fixed << showpoint << setprecision(7);
    stage3 << fixed << showpoint << setprecision(7);

    Stack.push(make_pair(identity(), false));
    stage1.open("stage1.txt");
    stage2.open("stage2.txt");
    stage3.open("stage3.txt");
    cin >> eyex >> eyey >> eyez >> lookx >> looky >> lookz >> upx >> upy >> upz >> fovy >> aspectratio >> near >> far;
    //cout<<"far:"<<far<<endl;

    lx = lookx - eyex;
    ly = looky - eyey;
    lz = lookz - eyez;
    double val;
    val = sqrt(lx * lx + ly * ly + lz * lz);
    lx = lx / val;
    ly = ly / val;
    lz = lz / val;

    rx = ly * upz - lz * upy;
    ry = lz * upx - lx * upz;
    rz = lx * upy - ly * upx;
    val = sqrt(rx * rx + ry * ry + rz * rz);
    rx = rx / val;
    ry = ry / val;
    rz = rz / val;

    ux = ry * lz - rz * ly;
    uy = rz * lx - rx * lz;
    uz = rx * ly - ry * lx;
    val = sqrt(ux * ux + uy * uy + uz * uz);
    ux = ux / val;
    uy = uy / val;
    uz = uz / val;

    T = identity();
    T[0][3] = -eyex;
    T[1][3] = -eyey;
    T[2][3] = -eyez;

    R = identity();
    R[0][0] = rx;
    R[0][1] = ry;
    R[0][2] = rz;
    R[1][0] = ux;
    R[1][1] = uy;
    R[1][2] = uz;
    R[2][0] = -lx;
    R[2][1] = -ly;
    R[2][2] = -lz;
    V = mult(R, T);
    fovx = fovy * aspectratio;
    t = near * tan(fovy / 2.0 * PI / 180.0);
    r = near * tan(fovy / 2.0 * PI / 180.0);
    P = identity();
    P[0][0] = near / r;
    P[1][1] = near / t;
    P[2][2] = -(far + near) / (far - near);
    P[2][3] = -(2 * far * near) / (far - near);
    P[3][2] = -1.0;
    P[3][3]=0.0;

    while (1)
    {
        cin >> com_line;
        if (com_line == "triangle")
        {
            cin >> p1x >> p1y >> p1z >> p2x >> p2y >> p2z >> p3x >> p3y >> p3z;
            double **trian = identity();
            trian[0][0] = p1x;
            trian[1][0] = p1y;
            trian[2][0] = p1z;
            trian[3][0] = 1.0;

            trian[0][1] = p2x;
            trian[1][1] = p2y;
            trian[2][1] = p2z;
            trian[3][1] = 1.0;

            trian[0][2] = p3x;
            trian[1][2] = p3y;
            trian[2][2] = p3z;
            trian[3][2] = 1.0;

            trian[0][3] = 1.0;
            trian[1][3] = 1.0;
            trian[2][3] = 1.0;
            trian[3][3] = 1.0;

            double **mod_trian = mult(Stack.top().first, trian);
            int i, j;
            cout << "modelling\n";
            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    if (j != 2)
                    {
                        stage1 << mod_trian[j][i] << " ";
                        cout << mod_trian[j][i] << " ";
                        // cout<<"hellooooo";
                    }
                    else
                    {
                        stage1 << mod_trian[j][i];
                        cout << mod_trian[j][i];
                    }
                }
                stage1 << endl;
                cout << endl;
            }
            stage1 << endl;
            cout << endl;

            double **view_trian = mult(V, mod_trian);
            cout << "view\n";
            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    if (j != 2)
                    {
                        stage2 << view_trian[j][i] << " ";
                        cout << view_trian[j][i] << " ";
                        // cout<<"hellooooo";
                    }
                    else
                    {
                        stage2 << view_trian[j][i];
                        cout << view_trian[j][i];
                    }
                }
                stage2 << endl;
                cout << endl;
            }
            stage2 << endl;
            cout << endl;


            double **pro_trian = mult(P, view_trian);
            cout<<"Projection\n";
            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    if (j != 2)
                    {
                        stage3 << pro_trian[j][i]/pro_trian[3][i] << " ";
                        cout  << pro_trian[j][i]/pro_trian[3][i] << " ";
                        // cout<<"hellooooo";
                    }
                    else
                    {
                        stage3 << pro_trian[j][i]/pro_trian[3][i];
                        cout << pro_trian[j][i]/pro_trian[3][i];
                    }
                }
                stage3 << endl;
                cout << endl;
            }
            stage3 << endl;
            cout<<endl;

        }

        else if (com_line == "translate")
        {
            cin >> tx >> ty >> tz;
            double **mat=identity();
            mat[0][3]=tx;
            mat[1][3]=ty;
            mat[2][3]=tz;

            transformation(mat);
        }

        else if (com_line == "scale")
        {
            cin >> sx >> sy >> sz;
            double **mat=identity();
            mat[0][0]=sx;
            mat[1][1]=sy;
            mat[2][2]=sz;

            transformation(mat);
        }
        else if (com_line == "rotate")
        {
            cin >> ang >> ax >> ay >> az;
            double sin_t = sin(PI * ang / 180.0);
            double cos_t = cos(PI * ang / 180.0);
            double val = sqrt(ax * ax + ay * ay + az * az);
            double aax = ax /val;
            double aay = ay / val;
            double aaz = az / val;
            double** mat = identity();
            //rodrigues formula
            mat[0][0] = cos_t + (1 - cos_t)*(aax * aax);
            mat[0][1] = (1 - cos_t)*(aax * aay)  - sin_t*aaz ;
            mat[0][2] = (1 - cos_t)*(aax * aaz) + sin_t*aay;
            mat[1][0] = (1 - cos_t)*(aax * aay) + sin_t*aaz;
            mat[1][1] = cos_t + (1 - cos_t)*(aay * aay);
            mat[1][2] = (1 - cos_t)*(aay * aaz) - sin_t*aax;
            mat[2][0] = (1 - cos_t)*(aax * aaz) - sin_t*aay ;
            mat[2][1] = (1 - cos_t)*(aay * aaz)  + sin_t*aax ;
            mat[2][2] = cos_t + (1 - cos_t)*(aaz * aaz) ;

            transformation(mat);
        }

         else if (com_line == "push")
        {
            push();
        }
        else if (com_line == "pop")
        {
            pop();
        }
        else if (com_line == "end")
        {
            break;
        }
    }

    return 0;
}

#include<bits/stdc++.h>
using namespace std;
int main()
{
    float eyex,eyey,eyez,lookx,looky,lookz,upx,upy,upz,foy,aspectratio,near,far;
    float p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z,sx,sy,sz,tx,ty,tz;
    ifstream fin;
    string line;
    fin.open("scene.txt");

    int wc=0;
    int lc=0;
    while (fin)
    {
        getline(fin, line);
        istringstream ss(line);
        if(lc<4)
        {
            do
            {
                string word;
                ss >> word;
                stringstream geek(word);
                {
                    switch (wc)
                    {
                    case 0:
                        geek >> eyex;
                        cout << "ex:"<<eyex << endl;
                        break;
                    case 1:
                        geek >> eyey;
                        cout << "ey:"<<eyey << endl;
                        break;
                    case 2:
                        geek >> eyez;
                        cout << "ez:"<<eyez << endl;
                        break;
                    case 4:
                        geek >> lookx;
                        cout << "lx:"<<lookx << endl;
                        break;
                    case 5:
                        geek >> looky;
                        cout << "ly:"<<looky << endl;
                        break;
                    case 6:
                        geek >> lookz;
                        cout << "lz:"<<lookz << endl;
                        break;
                    case 8:
                        geek >> upx;
                        cout << "ux:"<<upx << endl;
                        break;
                    case 9:
                        geek >> upy;
                        cout << "uy:"<<upy << endl;
                        break;
                    case 10:
                        geek >> upz;
                        cout << "uz:"<<upz << endl;
                        break;
                    case 12:
                        geek >> foy;
                        cout << "foy:"<<foy << endl;
                        break;
                    case 13:
                        geek >> aspectratio;
                        cout << "aspect ratio:"<<aspectratio << endl;
                        break;
                    case 14:
                        geek >> near;
                        cout << "near:"<<near << endl;
                        break;
                    case 15:
                        geek >> far;
                        cout << "far:"<<far << endl;
                        break;
                    }
                }

                wc++;
            }
            while (ss);
        }
        else if(line=="triangle")
        {

        }

        lc++;
    }
    fin.close();
    return 0;
}

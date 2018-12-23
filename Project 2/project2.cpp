//
//  project2.cpp
//  Graphics
//
//  Created by Neha Pullabhotla on 10/21/18.
//  Copyright © 2018 Neha Pullabhotla. All rights reserved.
//

#include <iostream>
//#include <OpenGL/gl.h>
//#include <OpenGL/glu.h>
#include <GL/glut.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

//bounding box struct
typedef struct
{
    float x_max, x_min, y_max, y_min, z_max, z_min;
    
}BoundingBox;

BoundingBox _boundingBox;

//initialize bounding box;
BoundingBox makeBoundingBox(float x_max, float x_min, float y_max, float y_min, float z_max, float z_min)
{
    _boundingBox.x_max = x_max;
    _boundingBox.x_min = x_min;
    _boundingBox.y_max = y_max;
    _boundingBox.y_min = y_min;
    _boundingBox.z_max = z_max;
    _boundingBox.z_min = z_min;
    
    
    return _boundingBox;
}

//global functions
void drawScene();
void BB_to_NDC();
void displayRotationAxis(float x1, float y1, float z1, float x2, float y2, float z2);
float windowWidth = 500;
float windowHeight = 500;
GLvoid *font_style = GLUT_BITMAP_TIMES_ROMAN_24;
static long font = (long)GLUT_BITMAP_8_BY_13; // Font selection.

//function to intialize draw string
void drawstr(int x, int y, const char* format, size_t length)
{
    glRasterPos3f(x,y,0);
    for(int i=0; i<length; ++i)
        glutBitmapCharacter(font_style, format[i]);
    glEnd();
}

// Routine to draw a bitmap character string.
void writeBitmapString(void *font, char *string)
{
    char *c;
    
    for (c = string; *c != '\0'; c++) glutBitmapCharacter(font, *c);
}

// Routine to convert floating point to char string.
void floatToString(char * destStr, int precision, float val)
{
    sprintf(destStr, "%f", val);
    destStr[precision] = '\0';
}

//class Point
class Point
{
public:
    float x, y, z, num_points;
    
    Point() : x(0), y(0), z(0) {};
    Point(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {};
    void translate3D(float tx, float ty, float tz);
    BoundingBox boundingBox;
    void setBB(BoundingBox boundingBox);
    void computeBoundingBox(vector<vector<Point> >&coords);
    void rotate3D(float x1, float y1, float z1, float x2, float y2, float z2, float angle);
    float xNDC = 0, yNDC = 0, zNDC = 0;
    float xPix = 0, yPix = 0, zPix = 0;
    float x_rot = 0, y_rot = 0, z_rot = 0;
    void scale3D(float sx, float sy, float sz);
    float scale_x = 0, scale_y = 0, scale_z = 0;
    float transl_x = 0, transl_y = 0, transl_z = 0;
};

//struct initialize edges
typedef struct{
    int edge1, edge2;
} Edges;

void Point::setBB(BoundingBox boundingBox)
{
    _boundingBox = boundingBox;
}

//global vectors
vector<vector<Point> > threeDpolygon;
vector<vector<Edges> > all_edges;
int polygon3D_ID;
void outPutFile(vector<vector<Point> > &coords);

//read input from file
void readFromFile(string file, vector<vector<Point> > &list, vector<vector<Edges> > &list2)
{
    ifstream inputFile;
    inputFile.open(file.c_str());
    
    int num_polygons;
    
    inputFile >> num_polygons;
    int num_points;
    int num_edges;
    float x, y, z;
    int edge1, edge2;
    Point points;
    Edges edges;
    //vector<Vertex> coordinates;
    
    for(int i = 0; i < num_polygons; i++)
    {
        inputFile >> num_points;
        
        vector<Point> polygon;
        
        for(int j = 0; j < num_points; j++)
        {
            inputFile >> x >> y >> z;
            
            points.num_points = num_points;
            points.x = x;
            points.y = y;
            points.z = z;
            
            polygon.push_back(points);
        }
        
        list.push_back(polygon);
        
        inputFile >> num_edges;
        //cout << "num_edges = " << num_edges << endl;
        
        vector<Edges> read_edges;
        for(int k = 0; k < num_edges; k++)
        {
            inputFile >> edge1 >> edge2;
            edges.edge1 = edge1;
            edges.edge2 = edge2;
            read_edges.push_back(edges);
        }
        
        list2.push_back(read_edges);
    }
    
    
    
    //cout << "edges[0][0] = " << read_edges[0].edge1 << " " << all_edges[0][0].edge2 << endl;
    inputFile.close();
}

//find min and max of bounding box
void Point::computeBoundingBox(vector<vector<Point> >&coords)
{
    //cout << "in here woohoo" << endl;
    float min_x = 0, min_y = 0, min_z = 0, max_x = 0, max_y = 0, max_z = 0;
    //cout << "coords[0][0].x = " << coords[0][0].x << endl;
    for(int i = 1; i < threeDpolygon.size(); i++)
    {
        min_x = (coords[0][0].x < coords[i][0].x) ? coords[0][0].x : coords[i][0].x;
        min_y = (coords[0][0].y < coords[i][0].y) ? coords[0][0].y : coords[i][0].y;
        min_z = (coords[0][0].z < coords[i][0].z) ? coords[0][0].z : coords[i][0].z;
        max_x = (coords[0][0].x > coords[i][0].x) ? coords[0][0].x : coords[i][0].x;
        max_y = (coords[0][0].y > coords[i][0].y) ? coords[0][0].y : coords[i][0].y;
        max_z = (coords[0][0].z > coords[i][0].z) ? coords[0][0].z : coords[i][0].z;
    }
    for(int i = 0; i < coords.size(); i++)
    {
        for(int j = 1; j < coords[i].size(); j++)
        {
            if(coords[i][j].x < min_x)
            {
                min_x = coords[i][j].x;
            }
            
            if(coords[i][j].y < min_y)
            {
                min_y = coords[i][j].y;
            }
            if(coords[i][j].z < min_z)
            {
                min_z = coords[i][j].z;
            }
            
            max_x = (coords[i][j].x > max_x) ? coords[i][j].x : max_x;
            max_y = (coords[i][j].y > max_y) ? coords[i][j].y : max_y;
            max_z = (coords[i][j].z > max_z) ? coords[i][j].z : max_z;
        }
    }
    
    //cout << "max_x = " << max_x << endl;
    makeBoundingBox(max_x+0.5, min_x-0.5, max_y+0.5, min_y-0.5, max_z+0.5, min_z-0.5);
    //BB_to_NDC();
}


//3d translation
//create identity matrix
void matrix4x4SetIdentity(float matIdent4x4[4][4])
{
    int row, col;
    
    for(row = 0; row < 4; row++)
    {
        for(col = 0; col < 4; col++)
        {
            matIdent4x4[row][col] = (row == col);
        }
    }
}

//function for translation
void Point::translate3D(float tx, float ty, float tz)
{
    vector<vector<Point> > all_translate;
    float matTransl3D[4][4];
    Point translate_point;
    
    matrix4x4SetIdentity(matTransl3D);
    matTransl3D[0][3] = tx;
    matTransl3D[1][3] = ty;
    matTransl3D[2][3] = tz;
    
        vector<Point> newTransl;
    
        for(int j = 0; j < threeDpolygon[polygon3D_ID-1].size(); j++)
        {
            translate_point.transl_x = threeDpolygon[polygon3D_ID-1][j].x + tx;
            translate_point.transl_y = threeDpolygon[polygon3D_ID-1][j].y + ty;
            translate_point.transl_z = threeDpolygon[polygon3D_ID-1][j].z + tz;
            
            newTransl.push_back(translate_point);
            
        }
        
        all_translate.push_back(newTransl);
    
    for(int i = 0; i < all_translate.size(); i++)
    {
        for(int j = 0; j < all_translate[i].size(); j++)
        {
            threeDpolygon[polygon3D_ID-1][j].x = all_translate[i][j].transl_x;
            threeDpolygon[polygon3D_ID-1][j].y = all_translate[i][j].transl_y;
            threeDpolygon[polygon3D_ID-1][j].z = all_translate[i][j].transl_z;
        }
    }
    
    outPutFile(threeDpolygon);
    
    computeBoundingBox(threeDpolygon);
    BB_to_NDC();
}

//intialize matrices for rotation
float matRot[4][4];

void matrixPreMultiplyMatrix(float matrix1[4][4], float matrix2[4][4])
{
    int row, col;
    
    float matTemp[4][4];
    
    for(row = 0; row < 4; row++)
    {
        for(col = 0; col < 4; col++)
        {
            matTemp[row][col] = matrix1[row][0] * matrix2[0][col] + matrix1[row][1] * matrix2[1][col] + matrix1[row][2] * matrix2[2][col] + matrix1[row][3] * matrix2[3][col];
        }
    }
    
    for(row = 0; row < 4; row++)
    {
        for(col = 0; col < 4; col++)
        {
            matrix2[row][col] = matTemp[row][col];
            
            //cout << "matrix = " << matrix2[row][col] << " ";
        }
        //cout << endl;
    }
}

void TranslRot3D(float tx, float ty, float tz)
{
    float matTranslRot3D[4][4];
    matrix4x4SetIdentity(matTranslRot3D);
    
    matTranslRot3D[0][3] = tx;
    matTranslRot3D[1][3] = ty;
    matTranslRot3D[2][3] = tz;
    
    matrixPreMultiplyMatrix(matTranslRot3D, matRot);
}

//3d rotation
void Point::rotate3D(float x1, float y1, float z1, float x2, float y2, float z2, float angle)
{
    vector<vector<Point> > all_rotate;
    float matQuaternionRot[4][4];
    
    float axisVecLength = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
    float cosA = cos(angle);
    float oneC = 1 - cosA;
    float sinA = sin(angle);
    float ux = (x2-x1)/(axisVecLength);
    float uy = (y2-y1)/(axisVecLength);
    float uz = (z2-z1)/(axisVecLength);
    
    TranslRot3D(-x1, -y1, -z1);
    
    matrix4x4SetIdentity(matQuaternionRot);
    
    matQuaternionRot[0][0] = ux * ux * oneC + cosA;
    matQuaternionRot[0][1] = ux * uy * oneC - uz * sinA;
    matQuaternionRot[0][2] = ux * uz * oneC + uy * sinA;
    matQuaternionRot[1][0] = uy * ux * oneC + uz * sinA;
    matQuaternionRot[1][1] = uy * uy * oneC + cosA;
    matQuaternionRot[1][2] = uy * uz * oneC - ux * sinA;
    matQuaternionRot[2][0] = uz * ux * oneC - uy * sinA;
    matQuaternionRot[2][1] = uz * uy * oneC + ux * sinA;
    matQuaternionRot[2][2] = uz * uz * oneC + cosA;
    
    matrixPreMultiplyMatrix(matQuaternionRot, matRot);
    TranslRot3D(x1, y1, z1);
    Point rot_points;
    
        vector<Point> rot_coords;
        //cout << "vertices = " << threeDpolygon[polygon3D_ID-1].size() << endl;
        for(int j = 0; j < threeDpolygon[polygon3D_ID-1].size(); j++)
        {
            //cout << "3d polygon x = " << threeDpolygon[i][j].x << endl;
            rot_points.x_rot = (matRot[0][0] * threeDpolygon[polygon3D_ID-1][j].x) + (matRot[0][1] * threeDpolygon[polygon3D_ID-1][j].y) + (matRot[0][2] * threeDpolygon[polygon3D_ID-1][j].z);
            rot_points.y_rot = (matRot[1][0] * threeDpolygon[polygon3D_ID-1][j].x) + (matRot[1][1] * threeDpolygon[polygon3D_ID-1][j].y) + (matRot[1][2] * threeDpolygon[polygon3D_ID-1][j].z);
           rot_points.z_rot = (matRot[2][0] * threeDpolygon[polygon3D_ID-1][j].x) + (matRot[2][1] * threeDpolygon[polygon3D_ID-1][j].y) + (matRot[2][2] * threeDpolygon[polygon3D_ID-1][j].z);
            
            rot_coords.push_back(rot_points);
        }
        
        all_rotate.push_back(rot_coords);

    //cout << "all rotate x = " << all_rotate[0][0].x_rot << endl;
    for(int i = 0; i < all_rotate.size(); i++)
    {
        for(int j = 0; j < all_rotate[i].size(); j++)
        {
            threeDpolygon[polygon3D_ID-1][j].x = all_rotate[i][j].x_rot;
            threeDpolygon[polygon3D_ID-1][j].y = all_rotate[i][j].y_rot;
            threeDpolygon[polygon3D_ID-1][j].z = all_rotate[i][j].z_rot;
        }
    }
    
    outPutFile(threeDpolygon);

    computeBoundingBox(threeDpolygon);
    BB_to_NDC();
    displayRotationAxis(x1, y1, z1, x2, y2, z2);
    
}

//function to display rotation axis
void displayRotationAxis(float x1, float y1, float z1, float x2, float y2, float z2)
{
    float xdelta = _boundingBox.x_max-_boundingBox.x_min;
    float ydelta = _boundingBox.y_max-_boundingBox.y_min;
    float zdelta = _boundingBox.z_max-_boundingBox.z_min;
    float delta = max(xdelta, ydelta);
    delta = max(delta, zdelta);
    
    float x1NDC = ((x1 - _boundingBox.x_min)/(delta))*2 + (-1);
    float x2NDC = ((x2 - _boundingBox.x_min)/(delta))*2 + (-1);
    float y1NDC = ((y1 - _boundingBox.x_min)/(delta))*2 + (-1);
    float y2NDC = ((y2 - _boundingBox.x_min)/(delta))*2 + (-1);
    float z1NDC = ((z1 - _boundingBox.x_min)/(delta))*2 + (-1);
    float z2NDC = ((z2 - _boundingBox.x_min)/(delta))*2 + (-1);
    
    
    // Xw = (xndc+1)*(width/2)+x
    float pixX1 = (x1NDC - (-1))/2 * 250;
    float pixY1 = (y1NDC - (-1))/2 * 250;
    float pixX2 = (x2NDC - (-1))/2 * 250;
    float pixY2 = (y2NDC - (-1))/2 * 250;
    float pixZ1 = (z1NDC - (-1))/2 * 250;
    float pixZ2 = (z2NDC - (-1))/2 * 250;
    
    //XY Projection
    glViewport(0, 250, 250, 250);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 250, 0, 250, -1, 1);
    
    glColor3f(1.0, 0.0, 0.0);
    
    glBegin(GL_LINES);
    glVertex2f(pixX1, pixY1);
    glVertex2f(pixX2, pixY2);
    glEnd();

    //XZ Projection
    glViewport(250, 250, 250, 250);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 250, 0, 250, -1, 1);
    
    glColor3f(0.0, 0.0, 1.0);
    
    glBegin(GL_LINES);
    glVertex2f(pixX1, pixZ1);
    glVertex2f(pixX2, pixZ2);
    glEnd();
    
    //YZ Projection
    glViewport(0, 0, 250, 250);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 250, 0, 250, -1, 1);
    
    glColor3f(1.0, 0.11, 0.68);
    
    glBegin(GL_LINES);
    glVertex2f(pixY1, pixZ1);
    glVertex2f(pixY2, pixZ2);
    glEnd();
    
    glFlush();
}

//3d scaling
void Point::scale3D(float sx, float sy, float sz)
{
    vector<vector<Point> > all_scale;
    float centroid_x = 0, centroid_y = 0, centroid_z = 0;
    //find the centroid -- average of all vertices
    //for(int i = 0; i < threeDpolygon.size(); i++)
    //{
        for(int j = 0; j < threeDpolygon[polygon3D_ID-1].size(); j++)
        {
            centroid_x = centroid_x + threeDpolygon[polygon3D_ID-1][j].x;
            centroid_y += threeDpolygon[polygon3D_ID-1][j].y;
            centroid_z += threeDpolygon[polygon3D_ID-1][j].z;
        }
        
        centroid_x = centroid_x/threeDpolygon[polygon3D_ID-1].size();
        //cout << "centroid_x = " << centroid_x << endl;
        centroid_y = centroid_y/threeDpolygon[polygon3D_ID-1].size();
        centroid_z = centroid_z/threeDpolygon[polygon3D_ID-1].size();
    //}
    
    float matScale3D[4][4];
    
    matrix4x4SetIdentity(matScale3D);
    
    matScale3D[0][0] = sx;
    matScale3D[0][3] = (1 - sx) * centroid_x;
    matScale3D[1][1] = sy;
    matScale3D[1][3] = (1 - sy) * centroid_y;
    matScale3D[2][2] = sz;
    matScale3D[2][3] = (1 - sz) * centroid_z;
    
    Point scale_point;
   
        vector<Point> scale_coords;
        for(int j = 0; j < threeDpolygon[polygon3D_ID-1].size(); j++)
        {
            scale_point.scale_x = (matScale3D[0][0] * threeDpolygon[polygon3D_ID-1][j].x) + (matScale3D[0][1] * threeDpolygon[polygon3D_ID-1][j].y) + (matScale3D[0][2] * threeDpolygon[polygon3D_ID-1][j].z) + (matScale3D[0][3] * 1);
            //cout << "x = " << scale_point.scale_x << endl;
            scale_point.scale_y = (matScale3D[1][0] * threeDpolygon[polygon3D_ID-1][j].x) + (matScale3D[1][1] * threeDpolygon[polygon3D_ID-1][j].y) + (matScale3D[1][2] * threeDpolygon[polygon3D_ID-1][j].z) + (matScale3D[1][3] * 1);
            //cout << "y = " << scale_point.scale_y << endl;
            scale_point.scale_z = (matScale3D[2][0] * threeDpolygon[polygon3D_ID-1][j].x) + (matScale3D[2][1] * threeDpolygon[polygon3D_ID-1][j].y) + (matScale3D[2][2] * threeDpolygon[polygon3D_ID-1][j].z) + (matScale3D[2][3] * 1);
            //cout << "z = " << scale_point.scale_z << endl;
            scale_coords.push_back(scale_point);
            //cout << "scale x = " << scale_coords[0].scale_x << endl;
        }
        
        all_scale.push_back(scale_coords);
   
    for(int k = 0; k < all_scale.size(); k++)
    {
        for(int l = 0; l < all_scale[k].size(); l++)
        {
            threeDpolygon[polygon3D_ID-1][l].x = all_scale[k][l].scale_x;
            threeDpolygon[polygon3D_ID-1][l].y = all_scale[k][l].scale_y;
            threeDpolygon[polygon3D_ID-1][l].z = all_scale[k][l].scale_z;
        }
    }
    outPutFile(threeDpolygon);
    
    //cout << " threeDpolygon[k][l].x = " <<  threeDpolygon[1][0].x << endl;
    computeBoundingBox(threeDpolygon);
    BB_to_NDC();
}


//function to convert world to NDC coords
void BB_to_NDC()
{
    vector<vector<Point> > NDC_vec;
    vector<vector<Point> > pixVector;
//    cout << "bounding box x_max = " << _boundingBox.x_max << endl;
//    cout << "bounding box x_min = " << _boundingBox.x_min << endl;
//    cout << "bounding box y_min = " << _boundingBox.y_min << endl;
//    cout << "bounding box y_max = " << _boundingBox.y_max << endl;
    float xdelta = _boundingBox.x_max-_boundingBox.x_min;
    float ydelta = _boundingBox.y_max-_boundingBox.y_min;
    float zdelta = _boundingBox.z_max-_boundingBox.z_min;
    float delta = max(xdelta, ydelta);
    delta = max(delta, zdelta);
   // cout << "Delta = " << delta << endl;
    Point NDC_val;

    for(int i = 0; i < threeDpolygon.size(); i++)
    {
        vector<Point> NDC;
       for(int j = 0; j < threeDpolygon[i].size(); j++)
       {
           NDC_val.xNDC = ((threeDpolygon[i][j].x - _boundingBox.x_min)/delta) * 2 + (-1) * 0.9;
           NDC_val.yNDC = ((threeDpolygon[i][j].y - _boundingBox.x_min)/delta) * 2 + (-1) * 0.9;
           NDC_val.zNDC = ((threeDpolygon[i][j].z - _boundingBox.x_min)/delta) * 2 + (-1) * 0.9;
           //cout << "x = " << NDC_val.xNDC << endl;
           NDC.push_back(NDC_val);
       }
        
        NDC_vec.push_back(NDC);
    }
    
    //cout << "ndc = " << NDC_vec[0][0].xNDC << endl;
    
    glViewport(0, 0, 500, 500);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 500, 0.0, 500, -1, 1);
    
    //draw grid
    //glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    glVertex2f(250,0 );
    glVertex2f(250, 500);
    glVertex2f(0,250 );
    glVertex2f(500, 250);
    glEnd();
    
    Point pixel_point;
   // Xw = (xndc+1)*(width/2)+x
    for(int i = 0; i < NDC_vec.size(); i++)
    {
        vector<Point> pix_coords;
        for(int j = 0; j < NDC_vec[i].size(); j++)
        {
            pixel_point.xPix = (NDC_vec[i][j].xNDC - (-1))/2 * 250;
            pixel_point.yPix = (NDC_vec[i][j].yNDC - (-1))/2 * 250;
            pixel_point.zPix = (NDC_vec[i][j].zNDC - (-1))/2 * 250;
            
//            cout << "pixX = " << pixel_point.xPix << endl;
//            cout << "pixY = " << pixel_point.yPix << endl;
//            cout << "pixZ = " << pixel_point.zPix << endl;
            
            pix_coords.push_back(pixel_point);
        }
        
        pixVector.push_back(pix_coords);
    }
    
    //XY Projection
    glViewport(0, 250, 250, 250);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 250, 0, 250, -1, 1);
    
    string str = "xy plane";
    drawstr(5,5, str.c_str(), str.length());
    //XY Projection
    glBegin(GL_LINES);
    
    for(int i = 0; i < pixVector.size(); i++)
    {
        for(int j = 0; j < all_edges[i].size(); j++)
        {
            glVertex2f(pixVector[i][all_edges[i][j].edge1-1].xPix, pixVector[i][all_edges[i][j].edge1-1].yPix);
            glVertex2f(pixVector[i][all_edges[i][j].edge2-1].xPix, pixVector[i][all_edges[i][j].edge2-1].yPix);
        }
    }
    
    glEnd();
    
    //XZ Projection
    glViewport(250, 250, 250, 250);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 250, 0, 250, -1, 1);
    
    str = "xz plane";
    drawstr(5, 5, str.c_str(), str.length());
    //XZ Projection
    glBegin(GL_LINES);
    
    for(int i = 0; i < pixVector.size(); i++)
    {
        for(int j = 0; j < all_edges[i].size(); j++)
        {
            glVertex2f(pixVector[i][all_edges[i][j].edge1-1].xPix, pixVector[i][all_edges[i][j].edge1-1].zPix);
            glVertex2f(pixVector[i][all_edges[i][j].edge2-1].xPix, pixVector[i][all_edges[i][j].edge2-1].zPix);
        }
    }
    
    glEnd();

    //YZ Projection
    glViewport(0, 0, 250, 250);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 250, 0, 250, -1, 1);
    
    str = "yz plane";
    drawstr(5, 5, str.c_str(), str.length());
    glBegin(GL_LINES);
    
    for(int i = 0; i < pixVector.size(); i++)
    {
        for(int j = 0; j < all_edges[i].size(); j++)
        {
            glVertex2f(pixVector[i][all_edges[i][j].edge1-1].zPix, pixVector[i][all_edges[i][j].edge1-1].yPix);
            glVertex2f(pixVector[i][all_edges[i][j].edge2-1].zPix, pixVector[i][all_edges[i][j].edge2-1].yPix);
        }
    }
    
    glEnd();
    
    glFlush();
}


// Routine to output interaction instructions to the C++ window.
void printInteraction(void)
{
    std::cout << "Hello:" << std::endl;
    std::cout << "Please enter 0 for default/original display of objects" << endl;
    std::cout << "Please select any of the following transformations: " << std::endl;
    std::cout << "1. Translation" << std::endl;
    std::cout << "2. Rotation" << std::endl;
    std::cout << "3. Scaling" << std::endl;
    std::cout << "or if you would like to exit:" << endl;
    std::cout << "4. Quit" << std::endl;
}

//function to write to output file
void outPutFile(vector<vector<Point> > &coords)
{
    ofstream outputFile;
    outputFile.open("output_transform.txt");
    
    outputFile << coords.size() << endl;
    outputFile << endl;
    
    for(int i = 0; i < coords.size(); i ++)
    {
        outputFile << coords[i].size() << endl;
        
        for(int j = 0; j < coords[i].size(); j++)
        {
            outputFile << coords[i][j].x << " " << coords[i][j].y << " " << coords[i][j].z << endl;
        }
        
        outputFile << all_edges[i].size() << endl;
        
        for(int k = 0; k < all_edges[i].size(); k++)
        {
            outputFile << all_edges[i][k].edge1 << " " << all_edges[i][k].edge2 << endl;
        }
        
        outputFile << endl;
    }
    
    outputFile.close();
}

//interaction
void read_input(int menu_num)
{
    float tx, ty, tz;
    float x1, y1, z1, x2, y2, z2, angle, sx, sy, sz;
    string display;
    
    if(menu_num == 0)
    {
        threeDpolygon[0][0].computeBoundingBox(threeDpolygon);
        BB_to_NDC();
    }
    
    if(menu_num == 1)
    {
        cout << "Please enter a translation vector of the form tx ty tz: " << endl;
        cin >> tx >> ty >> tz;
        cout << "Please select ID of object to be translated (in range from 1 to 3): " << endl;
        cin >> polygon3D_ID;
//        cout << "Do you want to display coordinates of the object?" << endl;
//        cin >> display;
        
        if(polygon3D_ID == 1)
        {
            threeDpolygon[0][0].translate3D(tx, ty, tz);            
        }
        
        else if(polygon3D_ID == 2)
        {
            threeDpolygon[1][0].translate3D(tx, ty, tz);
        }
        
        else
        {
            threeDpolygon[2][0].translate3D(tx, ty, tz);
        }       
    }
    
    if(menu_num == 2)
    {
        
        cout << "Please enter two points to define a rotation axis: " << endl;
        cout << "Point 1 of the form x1 y1 z1: ";
        cin >> x1 >> y1 >> z1;
        cout << "Point 2 of the form x2 y2 z2: ";
        cin >> x2 >> y2 >> z2;
        cout << "Please enter any angle in radians: " << endl;
        cin >> angle;
        cout << "Please select ID of object to be rotated (in range from 1 to 3): " << endl;
        cin >> polygon3D_ID;
        
        if(polygon3D_ID == 1)
        {
            matrix4x4SetIdentity(matRot);
         threeDpolygon[0][0].rotate3D(x1, y1, z1, x2, y2, z2, angle);
        }
        
        else if(polygon3D_ID == 2)
        {
            matrix4x4SetIdentity(matRot);
          threeDpolygon[1][0].rotate3D(x1, y1, z1, x2, y2, z2, angle);
        }
        
        else{
            matrix4x4SetIdentity(matRot);
            threeDpolygon[2][0].rotate3D(x1, y1, z1, x2, y2, z2, angle);
        }
    }
    
    if(menu_num == 3)
    {
        cout << "Please enter scaling factor of the form sx sy sz: " << endl;
        cin >> sx >> sy >> sz;
        cout << "Please select ID of object to be scaled (in range from 1 to 3): " << endl;
        cin >> polygon3D_ID;
        
        if(polygon3D_ID == 1)
        {
            threeDpolygon[0][0].scale3D(sx, sy, sz);
        }
        
        else if(polygon3D_ID == 2)
        {
            threeDpolygon[1][0].scale3D(sx, sy, sz);
        }
        
        else{
            //cout << "in here " << endl;
            threeDpolygon[2][0].scale3D(sx, sy, sz);
        }
    }
    
    if(menu_num == 4)
    {
        exit(0);
    }
   // glutPostRedisplay();
}

int main(int argc, char *argv[])
{
    string inputFile;
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);

    //set window size to 200*200
    glutInitWindowSize(windowWidth, windowHeight);
    //set window position
    glutInitWindowPosition(100, 100);
    
    //create and set main window title
    glutCreateWindow("Project 2!!");
    //glClearColor(0.0, 0.0, 0.0, 0.0); //clears the buffer of OpenGL
    
    cout << "Enter input file name in .txt form: " << endl;
    getline(cin, inputFile);

    readFromFile(inputFile, threeDpolygon, all_edges);

    //sets display function
    glutDisplayFunc(drawScene);
    //glutReshapeFunc(resize);
    //glutSpecialFunc(specialKeyFunction);
    
    glutMainLoop();//main display loop, will display until terminate
    return 0;
}

//main display loop, this function will be called again and again by OpenGL
void drawScene()
{

    //Misc.
    int menu_num;
    
    glClear(GL_COLOR_BUFFER_BIT);

    printInteraction();
    cin >> menu_num;
    
    glColor3f(1.0f, 1.0f, 1.0f); //white color
    
    read_input(menu_num);
    
    
    glutPostRedisplay();
    
    glFlush();
}


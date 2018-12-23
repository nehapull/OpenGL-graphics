#include <iostream>
// #include <GLUT/glut.h>
// #include <OpenGL/glu.h>
// #include <OpenGL/gl.h>
#include <GL/glut.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <list>
#include <string>
#include <cstring>
#include <algorithm>

using namespace std;

//initialize color
struct Color{
    float r, g, b;
    Color() {};
    Color(float red, float green, float blue) : r(red), g(green), b(blue) {};
    Color(float x) : r(x), g(x), b(x) {};
    int r_int = 0, g_int = 0, b_int = 0;
};

//point class
class Point
{
public:
    float x, y, z, num_points;
    Color color;
    float xNDC = 0, yNDC = 0, zNDC = 0;
    float xPix = 0, yPix = 0, zPix = 0;
    Point() : x(0), y(0), z(0) {};
    Point(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {};
    Point(float x, float y, float z, Color color = Color()) : x(x), y(y), color(color) {};
    float triangle_face;
    void computeBoundingBox(vector<vector<Point> >&coords);

    //function to calculate normals
    vector<Point> computeNormal();
    Point normalize(const Point &p);
    
    Point operator-(const Point &p)
    {
        Point point;
        
        point.x = this->x - p.x;
        point.y = this->y - p.y;
        point.z = this->z - p.z;
        
        return point;
    }
};

//initialize global pix buff
float *PixelBuffer;

void DDA(int x0, int y0, int xEnd, int yEnd);
void MakePixel(int x, int y);


vector<vector<Point> > threeDpolygon;
vector<vector<int> > all_faces;

//read input from file
void readFromFile(string file, vector<vector<Point> > &list, vector<vector<int> > &list2)
{
    ifstream inputFile;
    inputFile.open(file.c_str());
    
    int num_polygons;
    
    inputFile >> num_polygons;
    //cout << "num_polygons = " << num_polygons << endl;
    int num_points;
    int num_faces;
    float x, y, z;
    //int p1, p2, p3;
    Point points;
    int indices;
    string line;
    //vector<Vertex> coordinates;
    for(int i = 0; i < num_polygons; i++)
    {
        inputFile >> num_points;
        //cout << "num points = " << num_points << endl;
        
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
        
        //read number of faces of polygon
        inputFile >> num_faces;
        
        //read each face's edges connecting triangles
        vector<int> index;
        for(int k = 0; k < num_faces; k++)
        {
            for(int i = 0; i < 3; i++)
            {
                inputFile >> indices;
                index.push_back(indices);
            }
        }
        //cout << "all faces = " << list2.size() << endl;
        list2.push_back(index);
    }
    
    //cout << "list[0] size = " << list2[0].size() << endl;
    inputFile.close();
}


//bounding box
typedef struct
{
    float x_max, x_min, y_max, y_min, z_max, z_min;
    
}BoundingBox;

BoundingBox _boundingBox;

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

void Point::computeBoundingBox(vector<vector<Point> >&coords)
{
    float min_x = 0, min_y = 0, min_z = 0, max_x = 0, max_y = 0, max_z = 0;
 
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
    
    makeBoundingBox(max_x, min_x, max_y, min_y, max_z, min_z);
}

//MakePixel for DDA & wireframe ortho proj
void MakePixel(int x, int y)
{
    PixelBuffer[(400 * y + x) * 3 + 0] = 1; //Stores the pixel value into Array
    PixelBuffer[(400 * y + x) * 3 + 1] = 1;
    PixelBuffer[(400 * y + x) * 3 + 2] = 1;
    return;    
}

//DDA algorithm
void DDA(int x0, int y0, int xEnd, int yEnd)
{
    int dx = xEnd - x0;
    int dy = yEnd - y0;
    int numPixels, k;
    
    float xIncrement, yIncrement, x = x0, y = y0;
    
    //m <= 1 case
    if(abs(dx) > abs(dy))
    {
        numPixels = abs(dx);
    }
    
    //m > 1 case
    else
    {
        numPixels = abs(dy);
    }
    
    xIncrement = float(dx)/float(numPixels);
    yIncrement = float(dy)/float(numPixels);
    
    MakePixel(round(x), round(y));
    
    for(k = 0; k < numPixels; k++)
    {
        x += xIncrement;
        y += yIncrement;
        
        MakePixel(round(x), round(y));
    }
}

//convert world to NDC coords
vector<vector<Point> > BB_to_NDC()
{
    vector<vector<Point> > NDC_vec;
    vector<vector<Point> > pixVector;

    float xdelta = _boundingBox.x_max-_boundingBox.x_min;
    float ydelta = _boundingBox.y_max-_boundingBox.y_min;
    float zdelta = _boundingBox.z_max-_boundingBox.z_min;
    float delta = max(xdelta, ydelta);
    delta = max(delta, zdelta);
    float all_min = min(_boundingBox.x_min, min(_boundingBox.y_min, _boundingBox.z_min));
  
    Point NDC_val;
    
    for(int i = 0; i < threeDpolygon.size(); i++)
    {
        vector<Point> NDC;
        for(int j = 0; j < threeDpolygon[i].size(); j++)
        {
            //NDC_val.xNDC = (((threeDpolygon[i][j].x - _boundingBox.x_min)/delta) * 2) + (-1);
            NDC_val.xNDC = (((threeDpolygon[i][j].x - all_min)/delta) * 2) + (-1);
           //cout << "x NDC = " << NDC_val.xNDC << endl;
            //NDC_val.yNDC = (((threeDpolygon[i][j].y - _boundingBox.x_min)/delta) * 2) + (-1);
           NDC_val.yNDC = (((threeDpolygon[i][j].y - all_min)/delta) * 2) + (-1);
           //cout << "y NDC = " << NDC_val.yNDC << endl;
           // NDC_val.zNDC = (((threeDpolygon[i][j].z - _boundingBox.x_min)/delta) * 2) + (-1);
           NDC_val.zNDC = (((threeDpolygon[i][j].z - all_min)/delta) * 2) + (-1);
            //cout << "z NDC = " << NDC_val.zNDC << endl;
          
            NDC.push_back(NDC_val);
        }
        
        NDC_vec.push_back(NDC);
    }
  
    Point pixel_point;
    // Xw = (xndc+1)*(width/2)+x
    for(int i = 0; i < NDC_vec.size(); i++)
    {
        vector<Point> pix_coords;
        for(int j = 0; j < NDC_vec[i].size(); j++)
        {
            pixel_point.xPix = round((NDC_vec[i][j].xNDC - (-1))/2 * 350);
            pixel_point.yPix = round((NDC_vec[i][j].yNDC - (-1))/2 * 350);
            pixel_point.zPix = round((NDC_vec[i][j].zNDC - (-1))/2 * 350);
            
                        // cout << "pixX = " << pixel_point.xPix << endl;
                        // cout << "pixY = " << pixel_point.yPix << endl;
                        // cout << "pixZ = " << pixel_point.zPix << endl;
            
            pix_coords.push_back(pixel_point);
        }
        
        pixVector.push_back(pix_coords);
    }

    return pixVector;
}

//draw wireframe ortho projections
void drawWireFrame()
{
	PixelBuffer = new float[400 * 400 * 3];
  
  vector<vector<Point> > draw = BB_to_NDC();

        //XY Projection
   memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));
   glViewport(0, 400, 400, 400);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   //glOrtho(0, 400, 0, 400, -1, 1);
   
   for(int i = 0; i < draw.size(); i++)
   {	
       for(int j = 0; j < all_faces[i].size(); j += 3)
       {
           DDA(draw[i][all_faces[i][j]-1].xPix, draw[i][all_faces[i][j]-1].yPix, draw[i][all_faces[i][j+1]-1].xPix, draw[i][all_faces[i][j+1]-1].yPix);
           DDA(draw[i][all_faces[i][j]-1].xPix, draw[i][all_faces[i][j]-1].yPix, draw[i][all_faces[i][j+2]-1].xPix, draw[i][all_faces[i][j+2]-1].yPix);
           DDA(draw[i][all_faces[i][j+2]-1].xPix, draw[i][all_faces[i][j+2]-1].yPix, draw[i][all_faces[i][j+1]-1].xPix, draw[i][all_faces[i][j+1]-1].yPix);
       }
   }

   glRasterPos2f(-1, -1);
   glDrawPixels(400, 400, GL_RGB, GL_FLOAT, PixelBuffer);
   memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));

   glFlush();

  
   //XZ Projection
   memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));
   glViewport(400, 400, 400, 400);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   //glOrtho(0, 400, 0, 400, -1, 1);
   //PixelBuffer = new float[PixelBufferSize/2];
   //cout << "here" << endl;
   for(int i = 0; i < draw.size(); i++)
   {	
       for(int j = 0; j < all_faces[i].size(); j += 3)
       {
           DDA(draw[i][all_faces[i][j]-1].xPix, draw[i][all_faces[i][j]-1].zPix, draw[i][all_faces[i][j+1]-1].xPix, draw[i][all_faces[i][j+1]-1].zPix);
           DDA(draw[i][all_faces[i][j]-1].xPix, draw[i][all_faces[i][j]-1].zPix, draw[i][all_faces[i][j+2]-1].xPix, draw[i][all_faces[i][j+2]-1].zPix);
           DDA(draw[i][all_faces[i][j+2]-1].xPix, draw[i][all_faces[i][j+2]-1].zPix, draw[i][all_faces[i][j+1]-1].xPix, draw[i][all_faces[i][j+1]-1].zPix);
       }
   }
   
   glRasterPos2f(-1, -1);
   glDrawPixels(400, 400, GL_RGB, GL_FLOAT, PixelBuffer);
   memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));
   
   glFlush();

   //YZ Projection
   memset(PixelBuffer, 0, 400 * 400 * 3);
   glViewport(0, 0, 400, 400);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   //glOrtho(0, 400, 0, 400, -1, 1);
   //PixelBuffer = new float[PixelBufferSize/2];
   for(int i = 0; i < draw.size(); i++)
   {	
       for(int j = 0; j < all_faces[i].size(); j += 3)
       {
           DDA(draw[i][all_faces[i][j]-1].zPix, draw[i][all_faces[i][j]-1].yPix, draw[i][all_faces[i][j+1]-1].zPix, draw[i][all_faces[i][j+1]-1].yPix);
           DDA(draw[i][all_faces[i][j]-1].zPix, draw[i][all_faces[i][j]-1].yPix, draw[i][all_faces[i][j+2]-1].zPix, draw[i][all_faces[i][j+2]-1].yPix);
           DDA(draw[i][all_faces[i][j+2]-1].zPix, draw[i][all_faces[i][j+2]-1].yPix, draw[i][all_faces[i][j+1]-1].zPix, draw[i][all_faces[i][j+1]-1].yPix);
       }
   }
   
   glRasterPos2f(-1, -1);
   glDrawPixels(400, 400, GL_RGB, GL_FLOAT, PixelBuffer);
   memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));

   glFlush();
}


int polygon_ID = 0;
bool halfToneOn;

//compute normals for each face
vector<Point> Point::computeNormal()
{
	//Point obj to store vector x y z
	Point point;
	Point normal;

	vector<Point> normal_vector;
	for(int j = 0; j < threeDpolygon.size(); j++)
	{
		//vector to store normal vectors
		vector<Point> normal_face;
		for(int i = 0; i < all_faces[j].size(); i += 3)
		{
			//vector to store normal for each face

			//vector1: P2-P1
			float vector1_x = (threeDpolygon[j][all_faces[j][i+1]-1].x) - (threeDpolygon[j][all_faces[j][i]-1].x);
			float vector1_y = (threeDpolygon[j][all_faces[j][i+1]-1].y) - (threeDpolygon[j][all_faces[j][i]-1].y);
			float vector1_z = (threeDpolygon[j][all_faces[j][i+1]-1].z) - (threeDpolygon[j][all_faces[j][i]-1].z);
      // cout << "vector1_x = " << vector1_x << " " << "vector1_y = " << vector1_y << " " << "vector1_z = " << vector1_z << endl; 

			//vector2: P3-P1
			float vector2_x = (threeDpolygon[j][all_faces[j][i+2]-1].x) - (threeDpolygon[j][all_faces[j][i]-1].x);
			float vector2_y = (threeDpolygon[j][all_faces[j][i+2]-1].y) - (threeDpolygon[j][all_faces[j][i]-1].y);
			float vector2_z = (threeDpolygon[j][all_faces[j][i+2]-1].z) - (threeDpolygon[j][all_faces[j][i]-1].z); 
      // cout << "vector2_x = " << vector2_x << " " << "vector2_y = " << vector2_y << " " << "vector2_z = " << vector2_z << endl; 

			//cross product = v.x = aybz − azby, v.y = azbx − axbz, v.z = axby − aybx
			float vector_x = (vector1_y * vector2_z) - (vector1_z * vector2_y);
			float vector_y = (vector1_z * vector2_x) - (vector1_x * vector2_z);
			float vector_z = (vector1_x * vector2_y) - (vector1_y * vector2_x);

			//||v1 x v2|| = sqrt(pow(v, 2))
			float magnitude = sqrt(pow(vector_x, 2) + pow(vector_y, 2) + pow(vector_z, 2));

			//normalize the vector
			vector_x = vector_x/magnitude;
			vector_y = vector_y/magnitude;
			vector_z = vector_z/magnitude; 

			point.x = vector_x;
			point.y = vector_y;
			point.z = vector_z;

			//store each vectors x y z values in a Point vector
			normal_face.push_back(point);

			/*print normal for each face*/
        // cout << "normal x = " << normal_face[count].x << " " << "normal y = " << normal_face[count].y << " " << "normal z = "
        // << normal_face[count].z << endl; 

        // count++;

		}

		float normalx_sum = 0, normaly_sum = 0, normalz_sum = 0;

		//take sum of each normal for each triangle
		for(int i = 0; i < normal_face.size(); i++)
		{
			//sum of normal x, normal y, normal z
			normalx_sum += normal_face[i].x;
			//cout << "normalx_sum = " << normalx_sum << endl;
			normaly_sum += normal_face[i].y;
			//cout << "normaly_sum = " << normaly_sum << endl;
			normalz_sum += normal_face[i].z;
			//cout << "normalz_sum = " << normalz_sum << endl;
		}

		//divide the sum of normals by the number of faces
		normalx_sum = normalx_sum/all_faces.size();
		normaly_sum = normaly_sum/all_faces.size();
		normalz_sum = normalz_sum/all_faces.size();

		//find magnitude of above values
		float normal_magnitude = sqrt(pow(normalx_sum, 2) + pow(normaly_sum, 2) + pow(normalz_sum, 2));

		//divide by normal vector magnitude
		float normalx_vector = normalx_sum/normal_magnitude;
		float normaly_vector = normaly_sum/normal_magnitude;
		float normalz_vector = normalz_sum/normal_magnitude;
		
		//store each normal in point and store in point vector
		normal.x = normalx_vector;
		normal.y = normaly_vector;
		normal.z = normalz_vector;

		normal_vector.push_back(normal); 

	}

	return normal_vector;
}

//implement phong model
void computeIP(Color ka, Color kd, Color ks, float IA, float IL, Point fromPoint, float K, Point lightSource, int n, vector<vector<Color> > &allIP)
{
	Color KaVec, KdVec, KsVec, Ip;
	float maxI;
	Point fromPointVec, lightVec, reflVec, viewVec;
	float fromPoint_mag = 0, coefficient = 0, lightpos_mag = 0, calpha = 0, cbeta = 0;
	vector<Point> normal = threeDpolygon[0][0].computeNormal();
	
	//vector to store Ip values

   for(int j = 0; j < threeDpolygon.size(); j++)
   {	

    vector<Color> Ip_vals;
	for(int i = 0; i < threeDpolygon[j].size(); i++)
	{
		//obtain the normal vector
		//vector(ka) * IA
		KaVec.r = ka.r * IA;
		KaVec.g = ka.g * IA;
		KaVec.b = ka.b * IA;

		//||f-p||
		fromPointVec.x = fromPoint.x - threeDpolygon[j][i].x;
		fromPointVec.y = fromPoint.y - threeDpolygon[j][i].y;
		fromPointVec.z = fromPoint.z - threeDpolygon[j][i].z;

		//magnitude
		fromPoint_mag = sqrt(pow(fromPointVec.x, 2) + pow(fromPointVec.y, 2) + pow(fromPointVec.z, 2)); 

		//IL/||f-p||+K
		coefficient = IL/(fromPoint_mag + K);

		//light vector = x-p/||x-p||
		lightVec.x = lightSource.x - threeDpolygon[j][i].x;
		lightVec.y = lightSource.y - threeDpolygon[j][i].y;
		lightVec.z = lightSource.z - threeDpolygon[j][i].z;

		//light-pos magnitude
		lightpos_mag = sqrt(pow(lightVec.x, 2) + pow(lightVec.y, 2) + pow(lightVec.z, 2)); 

		//calpha ndotl = ax*bx + ay*by + az*bz
		calpha = (normal[j].x * lightVec.x) + (normal[j].y * lightVec.y) + (normal[j].z * lightVec.z);

		//vector(kd) * calpha
		KdVec.r = kd.r * max((float)0.0, calpha);
		KdVec.g = kd.g * max((float)0.0, calpha);
		KdVec.b = kd.g * max((float)0.0, calpha);

		//reflection vector = -l + 2(ndotl)n
		reflVec.x = -lightVec.x + 2 * (calpha) * normal[j].x;
		reflVec.y = -lightVec.y + 2 * (calpha) * normal[j].y;
		reflVec.z = -lightVec.z + 2 * (calpha) * normal[j].z;

		//view vector = f-p/||f-p||
		viewVec.x = fromPointVec.x/fromPoint_mag;
		viewVec.y = fromPointVec.y/fromPoint_mag;
		viewVec.z = fromPointVec.z/fromPoint_mag;

		//cbeta = rdotv = ax*bx + ay*by + az*bz
		cbeta = (reflVec.x * viewVec.x) + (reflVec.y * viewVec.y) + (reflVec.z * viewVec.z);
		cbeta = pow(max((float)0.0, cbeta), n);

		//vector(ks) * pow(cbeta, n);
		KsVec.r = ks.r * cbeta;
		KsVec.g = ks.g * cbeta;
		KsVec.b = ks.b * cbeta;

		//tie everything together
		Ip.r = KaVec.r + coefficient * (KdVec.r + KsVec.r);
		Ip.g = KaVec.g + coefficient * (KdVec.g + KsVec.g);
		Ip.b = KaVec.b + coefficient * (KdVec.b + KsVec.b);

		//normalize Ip values
		//find max Ip
		maxI = max(Ip.r, Ip.g);
		maxI = max(maxI, Ip.b);
		Ip.r = Ip.r/maxI;
		Ip.g = Ip.g/maxI;
		Ip.b = Ip.b/maxI;
		//cout << "Ip.r = " << Ip.r << " " << "Ip.g = " << Ip.g << " " << "Ip.b = " << Ip.b << endl;
		Ip_vals.push_back(Ip);
	}

    allIP.push_back(Ip_vals);

    // cout << "IP vals of r = " << allIP[j][i].r << endl;
    // cout << "IP vals of g = " << allIP[j][i].g << endl;
    // cout << "IP vals of b = " << allIP[j][i].b << endl;
  }
}


float magnitude(Point p)
{
    return sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2));
}

void externalInput();

struct Vertex
{
	float a = 0;
    float b = 0;
    Color color;
    
    Vertex(float a, float b, Color color = Color()) : a(a), b(b), color(color) {};
    //Vertex(float a, float b) : a(a), b(b);

    Vertex(){};
};

struct edge_bucket {
    int yMax;       
    int yMin;       
    int xIntercept;          
    int sign;       
    int dX;         
    int dY;         
    double sum;     
    Vertex p1, p2;
};

//setPixel for gouraud shading
void setPixel(int x, int y, Color color){

    PixelBuffer[(y * 400 + x) * 3 + 0] = color.r;
    PixelBuffer[(y * 400 + x) * 3 + 1] = color.g;
    PixelBuffer[(y * 400 + x) * 3 + 2] = color.b;

    return;   
}

bool findYMin(edge_bucket* edge1, edge_bucket* edge2) {
    return edge1->yMin < edge2->yMin;
}

bool findXIntercept (edge_bucket* edge1, edge_bucket* edge2) {
    if (edge1->xIntercept < edge2->xIntercept) {
        return true;
    } 
    else if (edge1->xIntercept > edge2->xIntercept) {
        return false;
    } 
    else {
        return ((edge1->dX/edge1->dY) < (edge2->dX/edge2->dY));
    }
}

int vertices;

list<edge_bucket*> createEdge(string proj_plane, vector<vector<Color> > IP_value) {

        //convert world to NDC
        vector<vector<Point> > pixel = BB_to_NDC();

        //new bucket values to initialize a bucket object
        int index = (vertices+3) - 1;
        int yMax;
        int yMin;
        int xIntercept;
        int sign;
        int dX;
        int dY;

        Vertex vertex1;
        Vertex vertex2;

        list<edge_bucket*> edgeTable;  

        	for(int j = vertices; j < vertices + 3; j++)
        	{
        		if(proj_plane == "xy")
        		{
            	vertex1.a = pixel[polygon_ID][all_faces[polygon_ID][index]-1].xPix;
            	//cout << "v1.x = " << vertex1.a << endl;
            	vertex1.b = pixel[polygon_ID][all_faces[polygon_ID][index]-1].yPix;
            	//cout << "v1.y = " << vertex1.b << endl;
            	vertex1.color = IP_value[polygon_ID][all_faces[polygon_ID][index]-1];
            	vertex2.a = pixel[polygon_ID][all_faces[polygon_ID][j]-1].xPix;
            	//cout << "v2.x = " << vertex2.a << endl;
            	vertex2.b = pixel[polygon_ID][all_faces[polygon_ID][j]-1].yPix;
            	//cout << "v2.y = " << vertex2.b << endl;
            	vertex2.color = IP_value[polygon_ID][all_faces[polygon_ID][j]-1];
            	}

            	if(proj_plane == "xz")
            	{

            	vertex1.a = pixel[polygon_ID][all_faces[polygon_ID][index]-1].xPix;
            	//cout << "v1.x = " << v1.a << endl;
            	vertex1.b = pixel[polygon_ID][all_faces[polygon_ID][index]-1].zPix;
            	//cout << "v1.y = " << v1.b << endl;
            	vertex1.color = IP_value[polygon_ID][all_faces[polygon_ID][index]-1];
              //cout << "v1.color.r = " << v1.color.r << endl;
            	vertex2.a = pixel[polygon_ID][all_faces[polygon_ID][j]-1].xPix;
            	//cout << "v2.x = " << pixVector[polygon_ID][all_faces[polygon_ID][j]-1].xPix << endl;
            	vertex2.b = pixel[polygon_ID][all_faces[polygon_ID][j]-1].zPix;
            	//cout << "v2.y = " << pixVector[polygon_ID][all_faces[polygon_ID][j]-1].zPix << endl;
            	vertex2.color = IP_value[polygon_ID][all_faces[polygon_ID][j]-1];
              //cout << "v2.color.r = " << v2.color.r << endl;
            	}

            	if(proj_plane == "yz")
            	{
            	vertex1.a = pixel[polygon_ID][all_faces[polygon_ID][index]-1].zPix;
            	//cout << "v1.x = " << vertex1.a << endl;
            	vertex1.b = pixel[polygon_ID][all_faces[polygon_ID][index]-1].yPix;
            	//cout << "v1.y = " << vertex1.b << endl;
            	vertex1.color = IP_value[polygon_ID][all_faces[polygon_ID][index]-1];
            	vertex2.a = pixel[polygon_ID][all_faces[polygon_ID][j]-1].zPix;
            	//cout << "v2.x = " << vertex2.a << endl;
            	vertex2.b = pixel[polygon_ID][all_faces[polygon_ID][j]-1].yPix;
            	//cout << "v2.y = " << vertex2.b << endl;
            	vertex2.color = IP_value[polygon_ID][all_faces[polygon_ID][j]-1];
            	}
            	
            	if (vertex2.a < vertex1.a) {
                swap(vertex1, vertex2);
            	}

            	yMax = max(vertex1.b, vertex2.b);
            	yMin = min(vertex1.b, vertex2.b);

              //x value associated with min y value for each edge
            	xIntercept = (vertex1.b < vertex2.b) ? vertex1.a : vertex2.a;

            	sign = (vertex2.b - vertex1.b) < 0 ? -1 : 1;
            	dX = abs(vertex2.a - vertex1.a);
            	dY = abs(vertex2.b - vertex1.b);

              //add edge to edge table if slope is not equal to 0 (not horizontal)
            	if (dY != 0) 
              {

              //populate new edge table  
              edge_bucket *newBucket = new edge_bucket;
                
              newBucket->yMax = yMax;
            	newBucket->yMin = yMin;
            	newBucket->xIntercept = xIntercept;
            	newBucket->sign = sign;
            	newBucket->dX = dX;
            	newBucket->dY = dY;
            	newBucket->sum = 0;                
              newBucket->p1 = vertex1;
              newBucket->p2 = vertex2; 

               edgeTable.push_back(newBucket);

            	}

            	index = j;
        	}
 
 	return edgeTable;
}

//gouraud shading 
void fillPolyhedra (list<edge_bucket*> edgeTable) {

		Color edgePoint_left, edgePoint_right, pixels_between;
		
        int scanline = edgeTable.front()->yMin;
        edge_bucket bucket1;
        edge_bucket bucket2;
        list<edge_bucket*> activeList;
        
        while (!edgeTable.empty()) {
            
          if (!activeList.empty()) {
            for (auto i = activeList.begin(); i != activeList.end();) {
                edge_bucket* curr = *i;
                    
              if (curr->yMax == scanline) {
              i = activeList.erase(i);
              edgeTable.remove (curr);
              delete (curr);
             } 

             else {
                        i++;
                    }
                }
            }
            

            for (auto i = edgeTable.begin (); i != edgeTable.end(); i++) {
                edge_bucket* edge = *i;
                
                if (edge->yMin == scanline) {
                    activeList.push_back(edge);
                }
            }

            //sort the active edge table based on min x value
            activeList.sort(findXIntercept);

            
            //fill polyhedra
            for (auto i = activeList.begin(); i != activeList.end(); i++) {
                bucket1 = **i;
                advance(i, 1);
                bucket2 = **i;
                	
                  //linear interpolation
                  edgePoint_left.r = ((float)(scanline - bucket1.yMax)/(bucket1.yMin - bucket1.yMax)) * bucket1.p1.color.r + ((float)(bucket1.yMin - scanline)/(bucket1.yMin - bucket1.yMax)) * bucket1.p2.color.r;
                  edgePoint_left.g = ((float)(scanline - bucket1.yMax)/(bucket1.yMin - bucket1.yMax)) * bucket1.p1.color.g + ((float)(bucket1.yMin - scanline)/(bucket1.yMin - bucket1.yMax)) * bucket1.p2.color.g;
                  edgePoint_left.b = ((float)(scanline - bucket1.yMax)/(bucket1.yMin - bucket1.yMax)) * bucket1.p1.color.b + ((float)(bucket1.yMin - scanline)/(bucket1.yMin - bucket1.yMax)) * bucket1.p2.color.b;

          
                  edgePoint_right.r = ((float)(scanline - bucket2.yMin)/(bucket2.yMax - bucket2.yMin)) * bucket2.p1.color.r + ((float)(bucket2.yMax - scanline)/(bucket2.yMax - bucket2.yMin)) * bucket2.p2.color.r;
                  edgePoint_right.g = ((float)(scanline - bucket2.yMin)/(bucket2.yMax - bucket2.yMin)) * bucket2.p1.color.g + ((float)(bucket2.yMax - scanline)/(bucket2.yMax - bucket2.yMin)) * bucket2.p2.color.g;
                  edgePoint_right.b = ((float)(scanline - bucket2.yMin)/(bucket2.yMax - bucket2.yMin)) * bucket2.p1.color.b + ((float)(bucket2.yMax - scanline)/(bucket2.yMax - bucket2.yMin)) * bucket2.p2.color.b;
          

                for (int x = bucket1.xIntercept; x < bucket2.xIntercept; x++) {

                	pixels_between.r = edgePoint_left.r * ((float)(bucket2.xIntercept-x)/(bucket2.xIntercept-bucket1.xIntercept)) + edgePoint_right.r * ((float)(x - bucket1.xIntercept)/(bucket2.xIntercept-bucket1.xIntercept));
                	pixels_between.g = edgePoint_left.g * ((float)(bucket2.xIntercept-x)/(bucket2.xIntercept-bucket1.xIntercept)) + edgePoint_right.g * ((float)(x - bucket1.xIntercept)/(bucket2.xIntercept-bucket1.xIntercept));
                	pixels_between.b = edgePoint_left.b * ((float)(bucket2.xIntercept-x)/(bucket2.xIntercept-bucket1.xIntercept)) + edgePoint_right.b * ((float)(x - bucket1.xIntercept)/(bucket2.xIntercept-bucket1.xIntercept));

                      //draw the pixels
                      setPixel(x, scanline, pixels_between);
                }
            }
            
            //increment scanline
            scanline++;
            
            for (auto i = activeList.begin(); i != activeList.end(); i++) {

                edge_bucket* increment = *i;
                
                if (increment->dX != 0) {
                    increment->sum += increment->dX;
                    
                    while (increment->sum >= increment->dY) {
                        increment->xIntercept += (increment->sign);
                        increment->sum -= increment->dY;
                    }
                }
            }
        }
  }

void HalfTone(float* PixelBuff);
void HalfTonesetPixel(int x, Color color);

//draw polygon using gouraud shading
void drawPolygon(string plane, vector<vector<Color> > illum)
{
  PixelBuffer = new float[400 * 400 * 3];

	if(plane == "xy")
	{
    memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));

    polygon_ID = 0;

      //XY Projection
    glViewport(0, 400, 400, 400);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho(0, 400, 0, 400, -1, 1);

   for(int j = 0; j < threeDpolygon.size(); j++)
   { 

	   for(int i = 0; i < all_faces[j].size(); i += 3)
	   {
		    vertices = i;

		  list<edge_bucket*> finalEdgeTable = createEdge("xy", illum);
		  finalEdgeTable.sort(findYMin);
		  if(!finalEdgeTable.empty())
      { 

      fillPolyhedra(finalEdgeTable);

      }
	   }

     polygon_ID++;
   }

     if(halfToneOn)
     {
       HalfTone(PixelBuffer);
     }
	   glRasterPos2f(-1, -1);
	   glDrawPixels(400, 400, GL_RGB, GL_FLOAT, PixelBuffer);

     memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));

     glFlush();
  }

    //XZ Projection
  	if(plane == "xz")
	{
    memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));
    polygon_ID = 0;

      //XZ Projection
    glViewport(400, 400, 400, 400);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho(0, 100, 0, 100, -1, 1);
   
    for(int j = 0; j < threeDpolygon.size(); j++)
   { 
    for(int i = 0; i < all_faces[j].size(); i += 3)
  {
    
    vertices = i;
    //cout << "in vertice = " << vertices << endl;

    list<edge_bucket*> finalEdgeTable = createEdge("xz", illum);
    finalEdgeTable.sort(findYMin);
    if(!finalEdgeTable.empty())
    {
      // cout << "finaledgetable not empty" << endl;
      fillPolyhedra(finalEdgeTable);
    }
  }

    polygon_ID++;
  }
     if(halfToneOn)
     {
       HalfTone(PixelBuffer);
     }

     glRasterPos2f(-1, -1);
     glDrawPixels(400, 400, GL_RGB, GL_FLOAT, PixelBuffer);

     memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));

     glFlush();

  }

  //YZ Projection
 if(plane == "yz")
  {
    memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));
      polygon_ID = 0;

    glViewport(0, 0, 400, 400);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho(0, 400, 0, 400, -1, 1);
   
   for(int j = 0; j < threeDpolygon.size(); j++)
   { 

    for(int i = 0; i < all_faces[j].size(); i += 3)
    {
      vertices = i;

      list<edge_bucket*> finalEdgeTable = createEdge("yz", illum);
      finalEdgeTable.sort(findYMin);
      if(!finalEdgeTable.empty())
      {
      fillPolyhedra(finalEdgeTable);
      }
    }

     polygon_ID++;
   }

     if(halfToneOn)
     {
       HalfTone(PixelBuffer);
     }
    glRasterPos2f(-1, -1);
    glDrawPixels(400, 400, GL_RGB, GL_FLOAT, PixelBuffer);

     memset(PixelBuffer, 0, 400 * 400 * 3 * sizeof(float));

     glFlush();
  }
}   

//set pixel for half toning
void HalfTonesetPixel(int x, Color color){

    PixelBuffer[x + 0] = color.r;
    PixelBuffer[x + 1] = color.g;
    PixelBuffer[x + 2] = color.b;

    return;   
}

//get pixel for half toning 
Vertex getPixel(int i, float* PixelBuff)
{
	Vertex vertex;
	vertex.color.r = PixelBuff[i-2];
	vertex.color.g = PixelBuff[i-1];
	vertex.color.b = PixelBuff[i];

	return vertex;
}

//Step1: divide the screen into 3X3 pixel blocks
//Step2: get the pixels from the original 3X3 pixel buffer
//Step3: implement half toning on each megapixel
void HalfTone(float* PixelBuff)
{
	for(int i = 400 * 400 * 3; i >= 0; i-=27)
	{
		vector<Vertex> pixels;
		vector<float> maxRGB;
		float sum = 0;
		//vector<int> random_vec;

		//each pixel in megapixel
		for(int j = i; j > i-27; j-=3)
		{
			int count = 0;
			float max_rgb; 
			Vertex pixel;
			pixel = getPixel(j, PixelBuff);
			pixel.a = j;

			max_rgb = max(pixel.color.r, max(pixel.color.g, pixel.color.b));

			pixels.push_back(pixel);
			maxRGB.push_back(max_rgb);
		}

		//cout << "maxrgb [0] = " << maxRGB[0] << endl;
		for(int k = 0; k < maxRGB.size(); k++)
		{
			sum = sum + maxRGB[k];
		}

		float average_rgb = (sum/maxRGB.size());
		int numPixelsOn = round(9 * average_rgb);

		for(int i = 0; i < numPixelsOn; i++)
		{
			int random_num = rand() % 9;

			pixels[random_num].color.r = 0.1;
			pixels[random_num].color.g = 0.1;
			pixels[random_num].color.b = 0.1;

			HalfTonesetPixel(pixels[random_num].a, pixels[random_num].color);
		}

	}

}

//painter's algo: find min z value for each face and sort from back to front
void painter_algo(vector<vector<Color> > phong)
{
	vector<Point> normal = threeDpolygon[0][0].computeNormal();
  Point min_value;
  vector<vector<Point> > min_vals;
  vector<vector<Point> > pixel = BB_to_NDC();

	//for each face
	
  for(int j = 0; j < threeDpolygon.size(); j++)
  {	
    vector<Point> localMin;

	for(int i = 0; i < all_faces[j].size(); i += 3)
	{
		float min_z = min(pixel[j][all_faces[j][i]-1].zPix, pixel[j][all_faces[j][i+1]-1].zPix);
		min_z = min(min_z, pixel[j][all_faces[j][i+2]-1].zPix);
		float min_x = min(pixel[j][all_faces[j][i]-1].xPix, pixel[j][all_faces[j][i+1]-1].xPix);
		min_x = min(min_x, pixel[j][all_faces[j][i+2]-1].xPix);
		float min_y = min(pixel[j][all_faces[j][i]-1].yPix, pixel[j][all_faces[j][i+1]-1].yPix);
		min_y = min(min_y, pixel[j][all_faces[j][i+2]-1].yPix);

		//cout << "min_z = " << min_z << endl;
    min_value.z = min_z;
    min_value.x = min_x;
    min_value.y = min_y;

    localMin.push_back(min_value);
		
    }
    
    min_vals.push_back(localMin);
  }


  //XY Projection
	//view pos: (0, 0, 1) for xy plane
   //for each polygon 
	for(int i = 0; i < min_vals.size(); i++)
	{
    int point = 0;
    //for each face
		for(int j = 0; j < min_vals[i].size(); j++)
	  {

		if(min_vals[i][j].z < min_vals[i][(j+1)%min_vals[i].size()].z)
		{	

      swap(all_faces[i][point], all_faces[i][(point+3)%all_faces[i].size()]);
			swap(all_faces[i][(point+1)%all_faces[i].size()], all_faces[i][(point+1+3)%all_faces[i].size()]);
      swap(all_faces[i][(point+2)%all_faces[i].size()], all_faces[i][(point+2+3)%all_faces[i].size()]);

		}

    // cout << "all faces = " << all_faces[i][point] << " " << all_faces[i][(point+1)%all_faces[i].size()] << " " <<
    // all_faces[i][(point+2)%all_faces[i].size()] << endl;

		point += 3;

	  }
	}

	  drawPolygon("xy", phong);


  //XZ Projection

	//view pos: (0, -1, 0) for xz plane

for(int i = 0; i < min_vals.size(); i++)
  {
    int point = 0;
    //for each face
    for(int j = 0; j < min_vals[i].size(); j++)
    {

    if(min_vals[i][(j+1)%min_vals[i].size()].y < min_vals[i][j].y)
    { 
      swap(all_faces[i][(point+3)%all_faces[i].size()], all_faces[i][point]);
      swap(all_faces[i][(point+1+3)%all_faces[i].size()], all_faces[i][(point+1)%all_faces[i].size()]);
      swap(all_faces[i][(point+2+3)%all_faces[i].size()], all_faces[i][(point+2)%all_faces[i].size()]);
    }

    point += 3;

    }
  }

    drawPolygon("xz", phong);

  //YZ Projection

  for(int i = 0; i < min_vals.size(); i++)
  {
    int point = 0;
    //for each face
    for(int j = 0; j < min_vals[i].size(); j++)
    {

    if(min_vals[i][j].x < min_vals[i][(j+1)%min_vals[i].size()].x)
    { 
      swap(all_faces[i][point], all_faces[i][(point+3)%all_faces[i].size()]);
      swap(all_faces[i][(point+1)%all_faces[i].size()], all_faces[i][(point+1+3)%all_faces[i].size()]);
      swap(all_faces[i][(point+2)%all_faces[i].size()], all_faces[i][(point+2+3)%all_faces[i].size()]);
    }

    //   cout << "all faces = " << all_faces[i][point] << " " << all_faces[i][(point+1)%all_faces[i].size()] << " " <<
    // all_faces[i][(point+2)%all_faces[i].size()] << endl;
    point += 3;

    }
  }

  //memset(PixelBufferYZ, 0, 400 * 400 * 3 * sizeof(float));
  drawPolygon("yz", phong);

}

//user interface
void externalInput()
{
	  
    int choice;

    cout << "Hello!" << std::endl;
    cout << "Please enter 0 for default/wireframe display of objects" << endl;
    cout << "Please enter 1 for Phong model and Gouraud shading" << endl;
    cout << "or if you would like to exit, please enter 2:" << endl;
    cin >> choice;

    if(choice == 0)
    {
    	threeDpolygon[0][0].computeBoundingBox(threeDpolygon);
        drawWireFrame();	
    }

    if(choice == 1)
    {
      vector<vector<Color> > phong;
      Color ka, kd, ks;
      float K = 0, IA, IL;
      int n = 0;
      Point lightSource, fromPoint;
      char half_tone;

    cout << "Please indicate if you want half toning on (Y or N): " << endl;
    cin >> half_tone;
    if(half_tone == 'Y')
    {
    	halfToneOn = true;
    }
    else{
    	halfToneOn = false;
    }
    
    cout << "Enter ka values :" << endl;
    cin >> ka.r >> ka.g >> ka.b;
    cout << "Enter kd values :" << endl;
    cin >> kd.r >> kd.g >> kd.b;
    cout << "Enter ks values :" << endl;
    cin >> ks.r >> ks.g >> ks.b;

    int numVertices = 0;
    
   for(int i = 0; i < threeDpolygon.size(); i++)
   { 
    for(int j = 0; j < threeDpolygon[i].size(); j++)
    {
        numVertices++;
        K += magnitude(lightSource - threeDpolygon[i][j]);
    }
   }
   
    K = K/numVertices;
    
    cout << "Enter n value :" << endl;
    cin >> n;
    cout << "Enter light source values :" << endl;
    cin >> lightSource.x >> lightSource.y >> lightSource.z;
    cout << "Enter from point values :" << endl;
    cin >> fromPoint.x >> fromPoint.y >> fromPoint.z;
    cout << "Enter IA value :" << endl;
    cin >> IA;
    cout << "Enter IL value :" << endl;
    cin >> IL;

    cout << "K value = " << K << endl;

    threeDpolygon[0][0].computeBoundingBox(threeDpolygon);
    //BB_to_NDC();
    computeIP(ka, kd, ks, IA, IL, fromPoint, K, lightSource, n, phong);
    painter_algo(phong);

  }



  if(choice == 2)
    {
      exit(0);
    }
}

//main display loop, this function will be called again and again by OpenGL
void display()
{
    externalInput();

    glutPostRedisplay();
        
    glFlush();
}


int main(int argc, char *argv[])
{
    string inputFile;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);

    //set window size to 200*200
    glutInitWindowSize(800, 800);
    //set window position
    glutInitWindowPosition(100, 100);
    
    //create and set main window title
    glutCreateWindow("Project 3!!");
    glClearColor(0.0, 0.0, 0.0, 0.0); //clears the buffer of OpenGL
    
    readFromFile("proj3_test.txt", threeDpolygon, all_faces);

    //sets display function
    glutDisplayFunc(display);
    
    glutMainLoop();//main display loop, will display until terminate
    return 0;
}

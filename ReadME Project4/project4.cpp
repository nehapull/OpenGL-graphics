// #include <OpenGL/gl.h>
// #include <OpenGL/glu.h>
// #include <GLUT/glut.h>
#include <GL/glut.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>

#include <iostream>

using namespace std;

float *PixelBuffer;
void display();

struct Color{
    float r, g, b;
    Color() {};
    Color(float red, float green, float blue) : r(red), g(green), b(blue) {};
    Color(float x) : r(x), g(x), b(x) {};
    int r_int = 0, g_int = 0, b_int = 0;
};

class Point
{

public:

	float _x, _y;
	Point(float x, float y) : _x(x), _y(y) {}; 
	Point() {};
	//Point find_bezier(float t);
};

int numControlPoints; 
//float resolution;

//MakePixel for DDA
void MakePixel(int x, int y, Color color)
{
    PixelBuffer[(500 * y + x) * 3 + 0] = color.r; //Stores the pixel value into Array
    PixelBuffer[(500 * y + x) * 3 + 1] = color.g;
    PixelBuffer[(500 * y + x) * 3 + 2] = color.b;
    return;    
}

void DDA(int x0, int y0, int xEnd, int yEnd, Color color)
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
    
    MakePixel(round(x), round(y), color);

    for(k = 0; k < numPixels; k++)
    {
        x += xIncrement;
        y += yIncrement;

        MakePixel(round(x), round(y), color);
    }
}

//Find bezier points
/*algo: for(j = 1 to n)
		for(i = 0 to n - j) do
		bji = (1-t)b(j-1)i + t*b(i+1)(j-1)
*/		
Point find_bezier(float t, vector<Point>& control_points)
{
	int j = control_points.size()-1;

	while(j > 0)
	{
		for(int i = 0 ; i < j; i++)
		{
			control_points[i]._x = (1-t) * control_points[i]._x + t * control_points[i+1]._x;
			control_points[i]._y = (1-t) * control_points[i]._y + t * control_points[i+1]._y;
		}

		j--;
	}

	return control_points[0];
}

vector<Point> bezierPoints;
//line segments for bezier
void computeBezierPoints(vector<Point>& control_points, int& resolution)
{
	//reset bezierPoints vector
	bezierPoints.clear();

	//Point pointOnCurve;
	//resolution
	float dt = ((float)1)/((float)resolution);
	//cout << "dt = " << dt << endl;

	for(float t = 0; t <= 1; t += dt)
	{
		Point pointOnCurve = find_bezier(t, control_points);
		bezierPoints.push_back(pointOnCurve);
	}
}

void drawBezier(vector<Point>& control_points, int& resolution)
{
	//cout << "control points size = " << control_points.size() << endl;
	computeBezierPoints(control_points, resolution);
	Color color;
    color.r = 0.7;
    color.g = 0.56;
    color.b = 0.56;

	//cout << "bezier point size = " << bezierPoints.size() << endl;
	for(int i = 0; i < bezierPoints.size() - 1; i++)
	{
		DDA(bezierPoints[i]._x, bezierPoints[i]._y, bezierPoints[i+1]._x, bezierPoints[i+1]._y, color);
	}
}

vector<Point> restorePoints;
int glcontrolPoints;
int res_bezier;

void get_pointsBezier()
{
	Point user_controlPoints, getNDControlPoints, getPixControlPoints;
	vector<Point> holdNDControlPoints;
	vector<Point> holdPixControlPoints;
	int choice;
	float min_x = 0, min_y = 0, max_x = 0, max_y = 0, minBox = 0, maxBox = 0;

	Color color;
    color.r = 0;
    color.g = 1;
    color.b = 0;
	int newPoint;
	int i;
	int deletePoint;
	vector<Point> hold_ctrlPoints;
	//cout << "hold points size = " << hold_ctrlPoints.size() << endl;

	cout << "Please enter 0 to draw the bezier curve" << endl;
	cout << "Press 1 if you would like to add a new control point" << endl;
	cout << "Press 2 if you would like to delete a control point" << endl;
	cout << "Press 3 if you would like to modify a control point" << endl;
	cin >> choice;

	if(choice == 0)
	{
	cout << "Please enter the number of control points: " << endl;
	cin >> numControlPoints;
	cout << "Please enter the coordinates for control points: " << endl;

	for(i = 0; i < numControlPoints; i++)
	{
		cin >> user_controlPoints._x >> user_controlPoints._y;

		hold_ctrlPoints.push_back(user_controlPoints);
	}	


	cout << "Please enter the resolution" << endl;
	cin >> res_bezier;

	min_x = hold_ctrlPoints[0]._x;
	min_y = hold_ctrlPoints[0]._y;
	//find min
	for(int i = 1; i < hold_ctrlPoints.size(); i++)
	{
		min_x = min(min_x, hold_ctrlPoints[i]._x);
		min_y = min(min_y, hold_ctrlPoints[i]._y);
	}

	minBox = min(min_x, min_y);
	max_x = hold_ctrlPoints[0]._x;
	max_y = hold_ctrlPoints[0]._y;

	//find max
	for(int i = 1; i < hold_ctrlPoints.size(); i++)
	{
		max_x = max(max_x, hold_ctrlPoints[i]._x);
		max_y = max(max_y, hold_ctrlPoints[i]._y);
	}

	maxBox = max(max_x, max_y);
	float delta_x = (max_x - min_x);
	float delta_y = (max_y - min_y);
	float delta = max(delta_x, delta_y);

	//do ndc mapping to pixel coords
	for(int l = 0; l < hold_ctrlPoints.size(); l++)
	{
		getNDControlPoints._x = (float)((hold_ctrlPoints[l]._x - min_x)/(delta));
		getNDControlPoints._y = (float)((hold_ctrlPoints[l]._y - min_y)/(delta));

		holdNDControlPoints.push_back(getNDControlPoints);

		//cout << "ndc coords = " << holdNDControlPoints[l]._x << " " << holdNDControlPoints[l]._y << endl;
	}


	//map ndc to pixel coords
	for(int k = 0; k < holdNDControlPoints.size(); k++)
	{
		getPixControlPoints._x = (float)round(holdNDControlPoints[k]._x * 350); 
		getPixControlPoints._y = (float)round(holdNDControlPoints[k]._y * 350);

		holdPixControlPoints.push_back(getPixControlPoints);
		//cout << "pix coords = " << holdPixControlPoints[k]._x << " " << holdPixControlPoints[k]._y << endl;
	}	

	for(int j = 0; j < numControlPoints-1; j++)
	{
		DDA(holdPixControlPoints[j]._x, holdPixControlPoints[j]._y, holdPixControlPoints[j+1]._x, holdPixControlPoints[j+1]._y, color);
		//cout << "size = " << holdControlPoints.size() << endl;
	}

	restorePoints = holdPixControlPoints;
	glcontrolPoints = i;

	drawBezier(holdPixControlPoints, res_bezier);
   }

	//glDrawPixels(500, 500, GL_RGB, GL_FLOAT, PixelBuffer);

   if(choice == 1)
   {
	//add new control point
		Point newCntrlPoint;
		memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
		holdPixControlPoints = restorePoints;
		cout << "Please enter coordinates for new control point: " << endl;
		cin >> newCntrlPoint._x >> newCntrlPoint._y;
		//cout << "i = " << i << endl;
		holdPixControlPoints.insert(holdPixControlPoints.begin()+glcontrolPoints, newCntrlPoint);
		restorePoints = holdPixControlPoints;
		for(i = 0; i < restorePoints.size()-1; i++)
	{
		DDA(restorePoints[i]._x, restorePoints[i]._y, restorePoints[i+1]._x, restorePoints[i+1]._y, color);
	}
		drawBezier(holdPixControlPoints, res_bezier);
  }

   if(choice == 2)
   {
	//delete control point
		memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
		holdPixControlPoints = restorePoints;
		cout << "Please enter the control point number you want to delete: " << endl;
		cin >> deletePoint;
		holdPixControlPoints.erase(holdPixControlPoints.begin()+(deletePoint-1));
		// cout << "elements in vector = " << hold_ctrlPoints[0]._x << " " << hold_ctrlPoints[0]._y << endl;
		// cout << "elements in vector = " << hold_ctrlPoints[1]._x << " " << hold_ctrlPoints[1]._y << endl; 
		// cout << "elements in vector = " << hold_ctrlPoints[2]._x << " " << hold_ctrlPoints[2]._y << endl;
		restorePoints = holdPixControlPoints;
		glcontrolPoints--;

		for(i = 0; i < restorePoints.size()-1; i++)
	{
		DDA(restorePoints[i]._x, restorePoints[i]._y, restorePoints[i+1]._x, restorePoints[i+1]._y, color);
	}
		drawBezier(holdPixControlPoints, res_bezier);
  }

  	if(choice == 3)
  	{
		int modifyPoint;
		Point modifycoords;
		memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
		holdPixControlPoints = restorePoints;
		cout << "Please enter the control point number you want to modify: " << endl;
		cin >> modifyPoint;
		cout << "Please enter new coordinates in order to modify control point: " << endl;
		cin >> modifycoords._x >> modifycoords._y;
		holdPixControlPoints[modifyPoint-1]._x = modifycoords._x;
		holdPixControlPoints[modifyPoint-1]._y = modifycoords._y;
		restorePoints = holdPixControlPoints;
		for(i = 0; i < restorePoints.size()-1; i++)
	{
		DDA(restorePoints[i]._x, restorePoints[i]._y, restorePoints[i+1]._x, restorePoints[i+1]._y, color);
	}
		drawBezier(holdPixControlPoints, res_bezier);	
	}

	glDrawPixels(500, 500, GL_RGB, GL_FLOAT, PixelBuffer);
	//glutPostRedisplay();
}

//initialize a vector to store knot values on u
vector<float> knotVals;

//calculate the distribution of knots
void calcKnots(int numControlPoints, int orderNumber)
{
	//re-initialize knot vector for every curve
	knotVals.clear();

	//number of knots = order number + num of control points
	int numKnots = orderNumber + numControlPoints;
	//cout << "num knots = " << numKnots << endl;

	//calculate the parameter line "u"
	for(int i = 0; i < numKnots; i++)
	{
		knotVals.push_back(i);
		//cout << "knot " << i << " = " << knotVals[i] << endl;
	}

	//cout << "knot size = " << knotVals.size() << endl;
}

//calculate points on the curve for b-spline using de Boor alg
Point deBoorAlgo(float uBar, int I, int& orderNumber, vector<Point>& control_points)
{
	vector<vector<Point> > storeCtrlPoints;
	storeCtrlPoints.push_back(control_points); 
	//cout << "store = " << storeCtrlPoints[0][0]._x << endl;

	for(int j = 1; j <= orderNumber-1; j++)
	{
		vector<Point> findControlPoint;
		for(int i = I-(orderNumber-1); i <= I-j; i++)
		{
			//initialize vector
			while(findControlPoint.size() != i)
			{
				findControlPoint.push_back(Point());
			}

			Point newCtrlPoint;
			newCtrlPoint._x = ((knotVals[i+orderNumber] - uBar)/(knotVals[i+orderNumber] - knotVals[i+j])) 
			* storeCtrlPoints[j-1][i]._x + ((uBar - knotVals[i+j])/(knotVals[i+orderNumber] - knotVals[i+j])) 
			* storeCtrlPoints[j-1][i+1]._x; 
			newCtrlPoint._y = ((knotVals[i+orderNumber] - uBar)/(knotVals[i+orderNumber] - knotVals[i+j])) 
			* storeCtrlPoints[j-1][i]._y + ((uBar - knotVals[i+j])/(knotVals[i+orderNumber] - knotVals[i+j])) 
			* storeCtrlPoints[j-1][i+1]._y; 

			findControlPoint.push_back(newCtrlPoint);
		}

		storeCtrlPoints.push_back(findControlPoint);
	}

	return storeCtrlPoints[orderNumber-1][I-(orderNumber-1)];
}

//function to find interval uI <= uBar <= uI+1
int findInterval(int uI, int domainInitial)
{
	for(int I = domainInitial; I < knotVals.size()-1; I++)
	{
		//check if i within interval uI <= uBar < uI+1
		if(uI >= knotVals[I] && uI < knotVals[I+1])
		{
			// cout << "i = " << i << endl;
			return I; 
		} 
	}

	return -1;
} 

//implement de boor algorithm
void drawSpline(int& orderNumber, vector<Point>& control_points, float& resolution)
{	

	Color color;
    color.r = 0.7;
    color.g = 0.56;
    color.b = 0.56;

	//initialize vector deBoorPoints
     vector<Point> deBoorPoints;

	//intialize Point deBoor
	Point deBoor;

	//intitialize number of control points
	numControlPoints = control_points.size()-1;

	//find I such that uI <= uBar <= uI+1
	float domainInitial = knotVals[orderNumber - 1];
	//cout << "domain init = " << domainInitial << endl;
	float domainEnd = knotVals[numControlPoints + 1];
	//cout << "domain end = " << domainEnd << endl;
	
	//resolution
	float dt = ((float)domainEnd-domainInitial)/((float)resolution);
	//cout << "dt = " << dt << endl;

	for(float uI = domainInitial; uI <= domainEnd; uI += dt)
	{
		int I = findInterval(uI, floor(domainInitial));
		deBoor = deBoorAlgo(uI, I, orderNumber, control_points);
		deBoorPoints.push_back(deBoor);
	}

	for(int i = 0; i < deBoorPoints.size()-1; i++)
	{
		DDA(deBoorPoints[i]._x, deBoorPoints[i]._y, deBoorPoints[i+1]._x, deBoorPoints[i+1]._y, color);
	}

}

//draw originial curve
void drawOrignialSpline(int& numControlPoints, int& orderNumber, vector<Point>& controlPoints, float& resolution)
{
	calcKnots(numControlPoints, orderNumber);
	drawSpline(orderNumber, controlPoints, resolution);

	//memset(PixelBuffer, 0, 500 * 500 * sizeof(float));
}

//draw curve with modified knot values
void drawKnotSpline(int& orderNumber, vector<Point>& controlPoints, float& resolution)
{
	drawSpline(orderNumber, controlPoints, resolution);
}

//add new control point for b-spline curve
void addPointSpline(vector<Point>& newSetCntrlPoints, int& orderNumber, float& resolution)
{
	//if user wants to add new control point, first ask for new coordinates
	//and pass in new vector of control points
	//recalculate knot vector
	calcKnots(newSetCntrlPoints.size(), orderNumber);
	drawSpline(orderNumber, newSetCntrlPoints, resolution);
}

//modify existing control point
void modifySpline(vector<Point>& newSetCntrlPoints, int& orderNumber, float& resolution)
{
	calcKnots(newSetCntrlPoints.size(), orderNumber);
	drawSpline(orderNumber, newSetCntrlPoints, resolution);
}

//delete control point
void deleteSpline(vector<Point>& newSetCntrlPoints, int& orderNumber, float& resolution)
{
	calcKnots(newSetCntrlPoints.size(), orderNumber);
	drawSpline(orderNumber, newSetCntrlPoints, resolution);
}

vector<Point> resetControlPoints;
int globalControlPoints;
int orderNumber;
bool changeOrder;
float res_spline;

void getPointsSpline()
{
	int numControlPoints;
	bool clearBuffer = false;
	float min_x = 0, min_y = 0, minBox = 0, max_x = 0, max_y = 0, maxBox = 0;
	// int globalControlPoints;
	Point getControlPoints;
	Point getNDControlPoints;
	Point getPixControlPoints;
	vector<Point> holdControlPoints;
	vector<Point> holdNDControlPoints;
	vector<Point> holdPixControlPoints;
	Color color;
	color.r = 0;
	color.g = 1;
	color.b = 0;
	// vector<Point> resetControlPoints;
	int choice = 0;
	int i;

	//initialize input order number to 2 iff user does not want to change order
	if(!changeOrder)
	{
	 orderNumber = 2;
	}

	cout << "Please enter 0 to draw the b-spline curve" << endl;
	cout << "Press 1 if you would like to add a new control point" << endl;
	cout << "Press 2 if you would like to modify the control point" << endl;
	cout << "Press 3 if you would like to delete a control point" << endl;
	cout << "Press 4 to change the order of the curve and/or modify the knot values" << endl;
	cin >> choice;

	
	//cout << "Please enter 0 to draw the b-spline curve" << endl;
	//cin >> choice;
	if(choice == 0)
	{
	//memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
	//input num control points
	cout << "Please enter the number of control points" << endl;
	cin >> numControlPoints;

	//coordinates of control points
	//cout << "control points = " << numControlPoints << endl;
	cout << "Please enter the coordinates of control points" << endl;
	
	for(i = 0; i < numControlPoints; i++)
	{
		cin >> getControlPoints._x >> getControlPoints._y;

		holdControlPoints.push_back(getControlPoints);
	}

	cout << "Please enter the resolution" << endl;
	cin >> res_spline;

	min_x = holdControlPoints[0]._x;
	min_y = holdControlPoints[0]._y;
	//find min
	for(int i = 1; i < holdControlPoints.size(); i++)
	{
		min_x = min(min_x, holdControlPoints[i]._x);
		min_y = min(min_y, holdControlPoints[i]._y);
	}

	minBox = min(min_x, min_y);
	max_x = holdControlPoints[0]._x;
	max_y = holdControlPoints[0]._y;

	//find max
	for(int i = 1; i < holdControlPoints.size(); i++)
	{
		max_x = max(max_x, holdControlPoints[i]._x);
		max_y = max(max_y, holdControlPoints[i]._y);
	}

	maxBox = max(max_x, max_y);
	float delta_x = (max_x - min_x);
	float delta_y = (max_y - min_y);
	float delta = max(delta_x, delta_y);

	//do ndc mapping to pixel coords
	for(int l = 0; l < holdControlPoints.size(); l++)
	{
		getNDControlPoints._x = (float)((holdControlPoints[l]._x - min_x)/(delta));
		getNDControlPoints._y = (float)((holdControlPoints[l]._y - min_y)/(delta));

		holdNDControlPoints.push_back(getNDControlPoints);

		//cout << "ndc coords = " << holdNDControlPoints[l]._x << " " << holdNDControlPoints[l]._y << endl;
	}


	//map ndc to pixel coords
	for(int k = 0; k < holdNDControlPoints.size(); k++)
	{
		getPixControlPoints._x = (float)round(holdNDControlPoints[k]._x * 350); 
		getPixControlPoints._y = (float)round(holdNDControlPoints[k]._y * 350);

		holdPixControlPoints.push_back(getPixControlPoints);
		//cout << "pix coords = " << holdPixControlPoints[k]._x << " " << holdPixControlPoints[k]._y << endl;
	}	

	for(int j = 0; j < numControlPoints-1; j++)
	{
		DDA(holdPixControlPoints[j]._x, holdPixControlPoints[j]._y, holdPixControlPoints[j+1]._x, holdPixControlPoints[j+1]._y, color);
		//cout << "size = " << holdControlPoints.size() << endl;
	}

	drawOrignialSpline(numControlPoints, orderNumber, holdPixControlPoints, res_spline);
	//glDrawPixels(500, 500, GL_RGB, GL_FLOAT, PixelBuffer);

	//store original control points in case of reset
	resetControlPoints = holdPixControlPoints;
	globalControlPoints = i;

	}

	//user input to add new control point
	if(choice == 1)
	{
		memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
		Point newcoords;
		cout << "Please enter coordinates for new control point" << endl;
		cin >> newcoords._x >> newcoords._y;
		resetControlPoints.insert(resetControlPoints.begin() + globalControlPoints, newcoords);
		for(int j = 0; j < resetControlPoints.size()-1; j++)
		{
			DDA(resetControlPoints[j]._x, resetControlPoints[j]._y, resetControlPoints[j+1]._x, resetControlPoints[j+1]._y, color);
		//cout << "size = " << holdControlPoints.size() << endl;
		}
		addPointSpline(resetControlPoints, orderNumber, res_spline);
		//glDrawPixels(500, 500, GL_RGB, GL_FLOAT, PixelBuffer); 
	}

	//user input to modify control point
	if(choice == 2)
	{
		memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
		int modifyPoint;
		Point modifycoords;
		cout << "Please enter the control point number you want to modify" << endl;
		cin >> modifyPoint;
		cout << "Please enter the new coordinates for the control point to be modified" << endl;
		cin >> modifycoords._x >> modifycoords._y;
		resetControlPoints[modifyPoint-1]._x = modifycoords._x;
		resetControlPoints[modifyPoint-1]._y = modifycoords._y;

		for(int j = 0; j < resetControlPoints.size()-1; j++)
		{
			DDA(resetControlPoints[j]._x, resetControlPoints[j]._y, resetControlPoints[j+1]._x, resetControlPoints[j+1]._y, color);
		//cout << "size = " << holdControlPoints.size() << endl;
		}

		modifySpline(resetControlPoints, orderNumber, res_spline);
		//glDrawPixels(500, 500, GL_RGB, GL_FLOAT, PixelBuffer);
		//memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
	}

	//user input to delete control point
	if(choice == 3)
	{
		memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
		int deletePoint;
		cout << "Please enter the control point number you want to delete" << endl;
		cin >> deletePoint;
		resetControlPoints.erase(resetControlPoints.begin() + (deletePoint-1));

		for(int j = 0; j < resetControlPoints.size()-1; j++)
		{
			DDA(resetControlPoints[j]._x, resetControlPoints[j]._y, resetControlPoints[j+1]._x, resetControlPoints[j+1]._y, color);
		//cout << "size = " << holdControlPoints.size() << endl;
		}

		deleteSpline(resetControlPoints, orderNumber, res_spline);
		//glDrawPixels(500, 500, GL_RGB, GL_FLOAT, PixelBuffer);
		//memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
		globalControlPoints--;
	}

	if(choice == 4)
	{
		memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
		char yesKnot;
		int newKnotVals;
		int drawChoice;
		changeOrder = true;
		cout << "Please enter the order value" << endl;
		cin >> orderNumber;

		cout << "Please enter 'Y' if you want to modify knot values or 'N' if you want to use the default knot values" << endl;
		cin >> yesKnot;
		if(yesKnot == 'Y')
		{
			cout << "Please enter new knot values" << endl;
			for(int i = 0; i < (orderNumber + resetControlPoints.size()); i++)
			{
				cout << "Knot value " << i << " = " << endl;
				cin >> newKnotVals; 

				knotVals.push_back(newKnotVals);
			}
		}  
		
		cout << "Please enter 1 to draw the new b-spline curve" << endl;
		cin >> drawChoice;
		if(drawChoice == 1)
		{
			for(int j = 0; j < resetControlPoints.size()-1; j++)
		{
			DDA(resetControlPoints[j]._x, resetControlPoints[j]._y, resetControlPoints[j+1]._x, resetControlPoints[j+1]._y, color);
		//cout << "size = " << holdControlPoints.size() << endl;
		}
			drawKnotSpline(orderNumber, resetControlPoints, res_spline);
			// glDrawPixels(500, 500, GL_RGB, GL_FLOAT, PixelBuffer);
			//memset(PixelBuffer, 0, 500 * 500 * 3 * sizeof(float));
		}

	}
	   glDrawPixels(500, 500, GL_RGB, GL_FLOAT, PixelBuffer);
   		//glutPostRedisplay();
}

bool curveType;

void writeToOutput()
{
	ofstream outputFile;
    outputFile.open("outputCurves.txt");

    // if(curveType == true)
    // {
    	outputFile << "Bezier curve" << endl; 
    
    outputFile << restorePoints.size() << endl;
    for(int i = 0; i < restorePoints.size(); i++)
    {
    	outputFile << restorePoints[i]._x << " " << restorePoints[i]._y << endl;
    }

    //}
    cout << endl;
    // if(curveType == false)
    // {
    	outputFile << "B-Spline curve" << endl; 
    
    outputFile << resetControlPoints.size() << endl;
    for(int i = 0; i < resetControlPoints.size(); i++)
    {
    	outputFile << resetControlPoints[i]._x << " " << resetControlPoints[i]._y << endl;
    }

    outputFile << orderNumber << endl;
    outputFile << knotVals.size() << endl;

    for(int j = 0; j < knotVals.size(); j++)
    {
    	outputFile << knotVals[j] << endl;
    }

    //}

    outputFile.close();
}

//main display loop, this function will be called again and again by OpenGL
void display()
{
	//Misc.
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();
	int selectCurve;

	cout << "Please enter 1 to draw/modify Bezier curves" << endl;
	cout << "Please enter 2 to draw/modify B-Spline curves" << endl;
	cout << "Please enter 3 to exit" << endl;
	cin >> selectCurve;

	if(selectCurve == 1)
	{
		get_pointsBezier();
	}

	if(selectCurve == 2)
	{
		getPointsSpline();
	}

	writeToOutput();

	if(selectCurve == 3)
	{
		exit(0);
	}

	//draws pixel on screen, width and height must match pixel buffer dimension
	//glDrawPixels(500, 500, GL_RGB, GL_FLOAT, PixelBuffer);

	//window refresh
	glutPostRedisplay();
	glFlush();
}

int main(int argc, char *argv[])
{
	//allocate new pixel buffer, need initialization!!
	PixelBuffer = new float[500 * 500 * 3];

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE);
	//set window size to 200*200
	glutInitWindowSize(500, 500);
	//set window position
	glutInitWindowPosition(100, 100);

	//create and set main window title
	glutCreateWindow("Project 4!!");
	glClearColor(0, 0, 0, 0); //clears the buffer of OpenGL

	//sets display function
	glutDisplayFunc(display);

	glutMainLoop();//main display loop, will display until terminate
	return 0;
}

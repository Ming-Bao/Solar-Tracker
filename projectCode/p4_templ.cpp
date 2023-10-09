#include "image_pr4.h" 
#include <thread>
#include <chrono>
using namespace std;

//kernals for findng vertical and horizontal edges
double kernal_x[3][3] = {
	{-1, 0, 1},
	{-2, 0, 2},
	{-1, 0, 1},
};
double kernal_y[3][3] = {
	{1, 2, 1},
	{0, 0, 0},
	{-1, -2, -1},
};

//struct for passing vector of x and y positions
struct Edge {
	vector<int> x;
	vector<int> y;
};

//struct for passing x and y positions of a single pixel
struct Pixel {
	int x;
	int y;
};

// you can use that
struct Orbit {
	// logged position and time
	vector<int> x;
	vector<int> y;
	vector<int> t;
	int xc,yc,r;  // center and radius
	int x_sunrise,y_sunrise = 900;
	double omega = 0.1;
} orbit;

//returns if the position contains a red pixel
bool is_red(int row, int col){
	return (int)get_pixel(image, row, col,0) > (int)get_pixel(image, row, col,1)*1.5 &&
		   (int)get_pixel(image, row, col,0) > (int)get_pixel(image, row, col,2)*1.5;
}

//filters the x and y positions of any red edges and returns them in an Edge struct
Edge filter_edge(){
	Edge edge;
	 
	for (int col = 1; col < image.width-1; col++){
		for (int row = 1; row < image.height-1; row++){
			//set_pixel(image1, row, col, (int)get_pixel(image, row, col,0), (int)get_pixel(image, row, col,1), (int)get_pixel(image, row, col,2));
			if (is_red(row, col)){
				double rx = 0;
				double ry = 0;

				for(int yk = -1; yk <= 1; yk++){
					for(int xk = -1; xk <= 1; xk++){
						rx += (int)get_pixel(image, row+yk, col+xk, 0) * kernal_x[yk+1][xk+1];
						ry += (int)get_pixel(image, row+yk, col+xk, 0) * kernal_y[yk+1][xk+1];
					}
				}
				
				double r_final = sqrt((rx*rx) + (ry*ry));

				if(r_final > 75){
					//cout<<"r_final: "<<r_final<<endl;
					edge.x.push_back(col);
					edge.y.push_back(row);
				}
			}
		}
	}
	return edge;
}

//stores the x, y position and time of the sun
void get_sun_position_with_formula(int t){
	orbit.x.push_back(orbit.xc + orbit.r*cos(orbit.omega*t));
	orbit.y.push_back(orbit.yc + orbit.r*sin(orbit.omega*t));
	orbit.t.push_back(t);
}

//calculates the determinate of the 3x3 matrix given
double det(double matrix[3][3]){
	double a = matrix[0][0];
	double b = matrix[0][1];
	double c = matrix[0][2];

	double d = matrix[1][0];
	double e = matrix[1][1];
	double f = matrix[1][2];

	double g = matrix[2][0];
	double h = matrix[2][1];
	double i = matrix[2][2];

	return (a*((e*i)-(f*h))) - (b*((d*i)-(f*g))) + (c*((d*h)-(e*g)));
}

//calculates the orbit of the sun
void calculate_orbit(){
	double a11=0, a12=0, a13=0, a21=0, a22=0, a23=0, a31=0, a32=0, a33=0, b1=0, b2=0, b3=0, b3_1=0, b3_2=0;
	int n = orbit.x.size();
	double omega = 0.1;

	//calculating vairable values
	a11 = n;
	a12 = 0;

	a21 = 0;
	a22 = n;
	
	a33 = n;

	for (int i = 0; i < n; i++){
		double xi = orbit.x.at(i);
		double yi = orbit.y.at(i);
		double ti = orbit.t.at(i);

		a13 += cos(omega * ti);
		b1  += xi;

		a23 += sin(omega * ti);
		b2  += yi;

		a31 += cos(omega * ti);
		a32 += sin(omega * ti);

		b3_1 += xi * cos(omega * ti);
		b3_2 += yi * sin(omega * ti);
	}
	b3 = b3_1 + b3_2;

	//initilising matricies
	double xc_top[3][3] ={
		{b1, a12, a13},
		{b2, a22, a23},
		{b3, a32, a33},
	};

	double xc_bottom[3][3] ={
		{a11, a12, a13},
		{a21, a22, a23},
		{a31, a32, a33},
	};

	double yc_top[3][3] ={
		{a11, b1, a13},
		{a21, b2, a23},
		{a31, b3, a33},
	};

	double yc_bottom[3][3] ={
		{a11, a12, a13},
		{a12, a22, a23},
		{a31, a32, a33},
	};

	double r_top[3][3] ={
		{a11, a12, b1},
		{a21, a22, b2},
		{a31, a32, b3},
	};

	double r_bottom[3][3] ={
		{a11, a12, a13},
		{a21, a22, a23},
		{a31, a32, a33},
	};

	//calculating orbit x, y, and r
	orbit.xc = det(xc_top)/det(xc_bottom);
	orbit.yc = det(yc_top)/det(yc_bottom);
	orbit.r  = det(r_top)/det(r_bottom);
	cout<<"xc: "<<orbit.xc<<"yc: "<<orbit.yc<<"r: "<<orbit.r<<endl;

	orbit.x_sunrise = orbit.xc - sqrt(pow(orbit.r, 2) - pow(900-orbit.yc, 2));
}

//checks if the sun is visiable
bool sun_visiable(){
	int red = 0;
	for (int row = 0; row < image.height; row ++){
		for (int col = 0; col < image.width; col ++){
			if (is_red(row, col)){
				red++;
			}
		}
	}
	return red > 300;
}

//go towards the sun
void go_to_sun(){
	int x_sun = orbit.x.at(orbit.x.size()-1);
	int y_sun = orbit.y.at(orbit.y.size()-1);
	int x_square;
	int y_square;

	if(y_sun > 900){
		x_sun = orbit.x_sunrise;
		y_sun = orbit.y_sunrise;
	}


	get_aim(x_square, y_square);

	double top = y_sun-y_square;
	double bottom = x_sun-x_square;

	move_aim(atan2(top, bottom));
}

//checks if 2 points are close enough to be part of the circle
bool close_enough(int x1, int y1, int x2, int y2, int threshold) {
    int distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    return distance < threshold;
}

//draws the edges and center of the sun
void draw_edge_and_center(Pixel pixel, Edge edge, int edge_count){
	ImagePPM image1;
	image1.width = 900;
	image1.height = 900;
    image1.n_bytes =  image1.width*image1.height*3;
    image1.data = new char[image1.n_bytes];
    
    for (int i = 0; i < edge_count; i++){
		int x = edge.x.at(i);
		int y = edge.y.at(i);
		
		set_pixel(image1, y, x, 0, 255, 0);
	}
	
	for(int yk = -1; yk <= 1; yk++){
		for(int xk = -1; xk <= 1; xk++){
			set_pixel(image1, pixel.y + xk , pixel.x + yk, 255, 0, 0);
		}
	}
	
    save_bmp_file("test.bmp", image1);
}

//vote for the center of the circle
Pixel vote_center(){
	Edge edge = filter_edge();
	int vote[900][900];
	int edge_count = edge.x.size();
	
	for (int x = 0; x < 900; x++){
		for (int y = 0; y < 900; y++){
			vote[x][y] = 0;
		}
	} 

	for (int pix_1 = 0; pix_1 < edge_count; pix_1++){
		int x1 = edge.x.at(pix_1);
		int y1 = edge.y.at(pix_1);
		
		for (int pix_2 = pix_1; pix_2 < edge_count; pix_2++){
			//check if the pixel 2 is within 33 pixels from pixel 1
			if(!close_enough(x1, y1, edge.x.at(pix_2), edge.y.at(pix_2), 33)){continue;}
			int x2 = edge.x.at(pix_2);
			int y2 = edge.y.at(pix_2);

			for (int pix_3 = pix_2; pix_3 < edge_count; pix_3++){
				//check if the pixel 3 is within 33 pixels from pixel 1
				if(!close_enough(x1, y1, edge.x.at(pix_3), edge.y.at(pix_3), 33)){continue;}
				int x3 = edge.x.at(pix_3);
				int y3 = edge.y.at(pix_3);

				//calculating circle center positions
				double a = (2*x1) - (2*x2);   
				double b = (2*y1) - (2*y2);
				double c = pow(x2, 2) - pow(x1, 2) + pow(y2, 2) - pow(y1, 2);

				double d = (2*x1) - (2*x3);
				double e = (2*y1) - (2*y3);
				double f = pow(x3, 2) - pow(x1, 2) + pow(y3, 2) - pow(y1, 2);

				int xc = ((-c*e) + (f*b))/((a*e) - (b*d));
				int yc = ((-a*f) + (c*d))/((a*e) - (b*d));

				if(xc < 900 && xc > 0 && yc < 900 && yc > 0){
					vote[xc][yc]++;
				}
			}
		}
	}

	Pixel pixel;
	double score = 0;

	for (int x = 0; x < 900; x++){
		for (int y = 0; y < 900; y++){
			if(vote[x][y] > score){
				score = vote[x][y];
				
				pixel.x = x;
				pixel.y = y;
			}
		}
	} 

	draw_edge_and_center(pixel, edge, edge_count);
    
    cout<<score<<endl;
	if(score < 25000){
		pixel.x = -1000;
		pixel.y = -1000;
	}
	return pixel;
}

//store position for orbit
void store_orbit_position(int t){
	Pixel pixel = vote_center();
	int x = pixel.x;
	int y = pixel.y;
	
	if(x > 0 && y > 0){
		orbit.t.push_back(t);
		orbit.x.push_back(x);
		orbit.y.push_back(y);
	}
}

int main(){        
	cout<<"start..."<<endl;
	init(2);
	//int x_sun, y_sun; // current position of the sun
    for ( int time = 0 ; time < 950; time++){
		cout<<endl;

		//vote_center();
		//filter_edge();
		if(time < 41) 		{store_orbit_position(time);}
		else if(time == 41) {calculate_orbit();}
		else{
			get_sun_position_with_formula(time);
		}
		go_to_sun();

      	draw_all(time);  // image is ready, 
      	cout<<" time="<<time<<endl;
      	this_thread::sleep_for(chrono::milliseconds(300));
    }
	return 0;
}

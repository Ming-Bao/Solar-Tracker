#include "image_pr4.h" 
#include <thread>
#include <chrono>
using namespace std;

// you can use that
struct Orbit {
	// logged position and time
	vector<int> x;
	vector<int> y;
	vector<int> t;
	int xc,yc,r;  // center and radius
	int x_sunrise,y_sunrise;
	double omega = 0.1;
} orbit;

//stores the x, y position and time of the sun
void store_sun_position(int t){
	// vector<int> x;
	// vector<int> y;

	// for (int row = 0; row < image.height; row ++){
	// 	for (int col = 0; col < image.width; col ++){
	// 		if ((int)get_pixel(image, row,col,0) > ((int)get_pixel(image, row,col,2) + (int)get_pixel(image, row,col,2))){
	// 			x.push_back(col);
	// 			y.push_back(row);
	// 		}
	// 	}
	// }

	// if(x.size() != 0 && y.size() != 0 && t <= 41){
	// 	orbit.t.push_back(t);
	// 	orbit.x.push_back(x.at((int)x.size()/2));
	//  	orbit.y.push_back(y.at((int)y.size()/2));	
	// }

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
	double a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3, b3_1, b3_2;
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
}

//checks if the sun is visiable
bool sun_visiable(){
	int red = 0;
	for (int row = 0; row < image.height; row ++){
		for (int col = 0; col < image.width; col ++){
			if ((int)get_pixel(image, row,col,0) > ((int)get_pixel(image, row,col,2) + (int)get_pixel(image, row,col,2))){
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

	get_aim(x_square, y_square);

	double top = y_sun-y_square;
	double bottom = x_sun-x_square;

	move_aim(atan2(top, bottom));
}

//store position for orbit
void store_orbit_position(int t){
	vector<int> x;
	vector<int> y;

	for (int row = 0; row < image.height; row ++){
		for (int col = 0; col < image.width; col ++){
			if ((int)get_pixel(image, row,col,0) > ((int)get_pixel(image, row,col,2) + (int)get_pixel(image, row,col,2))){
				x.push_back(col);
				y.push_back(row);
			}
		}
	}

	if(x.size() != 0 && y.size() != 0 && t <= 41){
		orbit.t.push_back(t);
		orbit.x.push_back(x.at((int)x.size()/2));
	 	orbit.y.push_back(y.at((int)y.size()/2));	
	}
}

int main(){        
	cout<<"start..."<<endl;
	init(3);
	//int x_sun, y_sun; // current position of the sun
    for ( int time = 0 ; time < 950; time++){
		cout<<endl;
		
		if(time < 41) 		{store_orbit_position(time);}
		else if(time == 41) {calculate_orbit();}
		else  				{store_sun_position(time);}

		go_to_sun();

      	draw_all(time);  // image is ready, 
      	cout<<" time="<<time<<endl;
      	this_thread::sleep_for(chrono::milliseconds(300));
    }
	return 0;
}
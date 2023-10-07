/*
 * 
 * */

#include "image_pr4.h" 
#include <thread>
#include <chrono>

// you can use that
struct Orbit {
	// logged position and time
	std::vector<int> x;
	std::vector<int> y;
	std::vector<int> t;
	int xc,yc,r;  // center and radius
	int x_sunrise,y_sunrise;
	double omega = 0.1;
} orbit;





int main()
{        
	std::cout<<"start..."<<std::endl;
	init(2);
	//int x_sun, y_sun; // current position of the sun
    for ( int time = 0 ; time < 950; time++){
       draw_all(time);  // image is ready, 
       std::cout<<" time="<<time<<std::endl;
       std::this_thread::sleep_for(std::chrono::milliseconds(500));
   }
    return 0;
}


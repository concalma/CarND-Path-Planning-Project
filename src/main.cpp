#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include <algorithm>
                    

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}


class SensorFusion {
    void estimateCarsPositions(int i) {
        // car position at step i. s and d
    }
};

struct Point {
    double x;
    double y;
};

class Path {
    public:
        Path() {
            idx = 0;
        }

        vector<Point> path;
        void add(Point p) {
            path.push_back(p);
        }

        void add(double x, double y) {
            Point p = {x,y};
            path.push_back(p);
        }

        void removeLastN( int n ) {
            if( n<=path.size() ) {
                path.resize( path.size() - n );
            }
        }

        Point getCurrent() {
            return path[idx];
        }

        int size() { return path.size(); }

        Point& operator[] (int x) {
          return path[x];
        }

        /* advanced d distance over the norm of the path starting at idx */
        void advance(double d) {
            if(path.size()<2) {
                cout << "Path::advance. path is less than 2\n";
            }

            double accum = 0, norm;
            for( int i=idx; accum<d && idx<(path.size()-1); i++ ) {
                double xx = path[i+1].x - path[i].x;
                double yy = path[i+1].y - path[i].y;
                norm = sqrt( xx*xx + yy*yy );
                accum += norm;
                idx++;
            }

            // we choose the point before the distance overflow.
            idx--;
            idx = std::max(idx,  0); 
        }

    private:
        int idx;
};

class PathSimulator {
    public:
        PathSimulator( 
                vector<double> map_waypoints_x,
                vector<double> map_waypoints_y,
                vector<double> map_waypoints_s,
                vector<double> map_waypoints_dx,
                vector<double> map_waypoints_dy
                ) {
            steps = 50;
            this->map_waypoints_x = map_waypoints_x;
            this->map_waypoints_y = map_waypoints_y;
            this->map_waypoints_s = map_waypoints_s;
            this->map_waypoints_dx = map_waypoints_dx;
            this->map_waypoints_dx = map_waypoints_dx;
        }


        // calculates path in x,y car coordinates
        //
        Path calculatePath(double car_x, double car_y, double car_yaw, int dest_lane ) {
            int spline_points = 5000;
            Path path;
            tk::spline s = computeSpline( car_x, car_y, car_yaw, dest_lane );

            // Oversample spline
            for(int i=0; i<spline_points; i++) {
                double x = spline_validity[0] + (double)i*(spline_validity[1]-spline_validity[0])/(double)spline_points;
                path.add( x, s(x) );
            }
            return path;
        }


        tk::spline computeSpline(double car_x, double car_y, double car_yaw, int dest_lane) {
            int lane = dest_lane;
            // Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
            // later we will interpolate these wps with a spline and fill it with more points that control speed
            vector<double> ptsx;
            vector<double> ptsy;

            // ref x,y,yaw states
            // either we will reference the starting point as where the car is or at the previous paths end point
            ref_x = car_x;
            ref_y = car_y;
            ref_yaw = deg2rad(car_yaw);

            int prev_size = historian.size();
            // if previous size is almost empty, use car as starting ref
            if( prev_size < 2 ) {
                //use two points that make the path tangent to the car
                double prev_car_x = car_x - cos(car_yaw);
                double prev_car_y = car_y - sin(car_yaw);

                ptsx.push_back(prev_car_x);
                ptsx.push_back(car_x);

                ptsy.push_back(prev_car_y);
                ptsy.push_back(car_y);

                spline_validity[0] = car_x; 
            }
            else {
                ref_x = historian[prev_size-1].x;
                ref_y = historian[prev_size-1].y;

                double ref_x_prev = historian[prev_size-2].x;
                double ref_y_prev = historian[prev_size-2].y;
                ref_yaw = atan2( ref_y - ref_y_prev, ref_x - ref_x_prev );

                //use two points that make the path tangent to the previous path's endpoint
                ptsx.push_back(ref_x_prev);
                ptsx.push_back(ref_x);
                spline_validity[0] = ref_x; 

                ptsy.push_back(ref_y_prev);
                ptsy.push_back(ref_y);
            }

            //In Frenet add evenly 30m spaced points ahead of the starting reference
            vector<double> next_wp0 = getXY( car_s +30/*meters*/, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y );
            vector<double> next_wp1 = getXY( car_s +60/*meters*/, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y );
            vector<double> next_wp2 = getXY( car_s +90/*meters*/, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y );

            ptsx.push_back(next_wp0[0]);
            ptsx.push_back(next_wp1[0]);
            ptsx.push_back(next_wp2[0]);
            spline_validity[1] = next_wp2[0]; 

            ptsy.push_back(next_wp0[1]);
            ptsy.push_back(next_wp1[1]);
            ptsy.push_back(next_wp2[1]);


            // transformation to local car coordinates
            for( int i=0; i<ptsx.size(); i++ )
            {
                double shift_x = ptsx[i]-ref_x;
                double shift_y = ptsy[i]-ref_y;

                ptsx[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
                ptsy[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));
            }


            // create a spline
            tk::spline s;
            // set (x,y) points to the spline
            s.set_points(ptsx, ptsy);

            return s;
        }

        double spline_validity[2];
            

        void estimateCarsPositions(int i) {
        }

        double mph_to_metersps( double d ) {
            return d*1609.34/3600.0;
        }

        double metersps_to_mph( double d ) {
            return (d*3600.0)/1609.34;
        }

        // if there is a car a car ahead, returns the speed of it
        bool isCarAhead(int lane_of_interest, int time_step, double &speed ) {
            int lane = lane_of_interest;
            // find ref_v to use
            for(int i=0; i<sensor_fusion.size(); i++) {
                // if car is in my lane ?
                float d = sensor_fusion[i][6];
                if(d< (2+4*lane+2) && d>(2+4*lane-2))
                {
                    // velocity and 'distance' of the car ahead
                    double vx = sensor_fusion[i][3];
                    double vy = sensor_fusion[i][4];
                    double check_speed = sqrt( vx*vx + vy*vy );
                    double check_car_s = sensor_fusion[i][5]; // s value of the car

                    check_car_s +=  ((double)time_step * 0.02*check_speed ); // if using previous points project s value
                    // check s values greater than ours and s gap
                    if((check_car_s > car_s) && ((check_car_s-car_s) <30))
                    {
                        // do some logic here, lower reference velocity so we don't crash into the car in front fo us
                        // we can also flag lange changing
                        //ref_vel = 29.5;
                        //
                        speed = mph_to_metersps( check_speed );
                        return true;
                    }
                }
            }
            return false;
        }

        // returns a path in car coordinates
        Path simulate(int destlane) {
            Path final_path;
            
            // Calculate tentative fine path
            Path p = calculatePath(car_x, car_y, car_yaw, destlane);
            for( int i=0; i<steps; i++ ) {
                // sf.estimateCarsPositions(i);
                 // checkColisions(); // false -> invalid path
               
                // if car ahead decelerate to match speed

                double carahead_speed = 0;
                double v0 = mph_to_metersps( car_speed );
                if ( isCarAhead( destlane, i, carahead_speed ) && carahead_speed < speedmax  ) {
                    double vt1 = v0 - accmax*step_s; 
                    vt1 = std::max( vt1, carahead_speed ); // in case car is stopped
                    double d = vt1*step_s; // distance traveled at vt1
                    p.advance(d);
                } else {
                    // no car ahead accelerate and clip
                    double vt1 = v0 + accmax*step_s;
                    vt1 = std::max( vt1, speedmax );
                    double d = vt1*step_s; // distance traveled
                    p.advance(d);
                }

                final_path.add(p.getCurrent());
            }

            return final_path;
        }

        void initSim( int previously_unused, double car_x, double car_y, double car_s, double car_d, double car_yaw, double car_speed, vector<vector<double>> sf ) {
            this->car_x = car_x;
            this->car_y = car_y;
            this->car_s = car_s;
            this->car_d = car_d;
            this->car_yaw = car_yaw;
            this->car_speed = car_speed;
            this->sensor_fusion = sf;

            // removed previously points that were not used from the historian
            if( previously_unused > 0 ) {
                historian.removeLastN( previously_unused );
            }
        }

        // calculates the next 50 x,y points in map coordinates
        void getNextPoints( vector<double> &next_x_vals, vector<double> &next_y_vals ) {
            int lane = 0;
            Path path = simulate(lane);

            for( int i=0; i<path.size(); i++ ) {
                double x_ref = path[i].x;
                double y_ref = path[i].y;
                


                // convert to map coordinates
                double x_point = (x_ref * cos(ref_yaw) - y_ref*sin(ref_yaw));
                double y_point = (x_ref * sin(ref_yaw) + y_ref*cos(ref_yaw));

                x_point += ref_x;
                y_point += ref_y;

                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);

                historian.add(x_point, y_point);
            }
        }


        Path historian;

    private:
        double accmax = 6.95; // max acceleration is 10 m/s^2
        double step_s = 0.02; 
        double speedmax = mph_to_metersps( 40.9 ); // m/s 

        int steps;

        // sim pars
        double car_x, car_y, car_s, car_d, car_yaw, car_speed;

        double ref_x, ref_y, ref_yaw; // updated in computeSpline()
        // Sensor Fusion Data, a list of all other cars on the same side of the road.
        //  vector<vector<double>>.   [N][0,1,2,vx,vy,car_s,car_d]
        vector<vector<double>> sensor_fusion;

        vector<double> map_waypoints_x;
        vector<double> map_waypoints_y;
        vector<double> map_waypoints_s;
        vector<double> map_waypoints_dx;
        vector<double> map_waypoints_dy;
};


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  // string map_file_ = "../data/highway_map.csv";
  string map_file_ = "../data/highway_map_bosch1.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  PathSimulator ps(map_waypoints_x, map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy );


  double ref_vel = 0; //mph
  h.onMessage([&ps, &ref_vel, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    int lane = 1;


    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner. x,y map coordinates
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
            //  vector<vector<double>>.   [N][0,1,2,vx,vy,car_s,car_d]
          	auto sensor_fusion = j[1]["sensor_fusion"];


            ps.initSim( previous_path_x.size(),  car_x, car_y, car_s, car_d, car_yaw, car_speed, sensor_fusion );

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

            ps.getNextPoints( next_x_vals, next_y_vals );
            
            // int prev_size = previous_path_x.size();


            /*
            // SENSOR FUSION
            if(prev_size>0) {
                car_s = end_path_s;
            }

            bool too_close = false;

            // find ref_v to use
            for(int i=0; i<sensor_fusion.size(); i++) {
                // if car is in my lane ?
                float d = sensor_fusion[i][6];
                if(d< (2+4*lane+2) && d>(2+4*lane-2))
                {
                    // velocity and 'distance' of the car ahead
                    double vx = sensor_fusion[i][3];
                    double vy = sensor_fusion[i][4];
                    double check_speed = sqrt( vx*vx + vy*vy );
                    double check_car_s = sensor_fusion[i][5]; // s value of the car

                    check_car_s +=  ((double)prev_size * 0.02*check_speed ); // if using previous points project s value
                    // check s values greater than ours and s gap
                    if((check_car_s > car_s) && ((check_car_s-car_s) <30))
                    {
                        // do some logic here, lower reference velocity so we don't crash into the car in front fo us
                        // we can also flag lange changing
                        //ref_vel = 29.5;
                        //
                        too_close = true;
                    }
                }
            }

            if(too_close)
            {
                ref_vel -= .224; // 5m/s2 
            }
            else if( ref_vel <49.5)
            {
                ref_vel += .224;
            }

            */

            /*

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

            // Start with all of the previous path points from last time
            for( int i =0; i<previous_path_x.size(); i++ ) {
                next_x_vals.push_back( previous_path_x[i] );
                next_y_vals.push_back( previous_path_y[i] );
            }

            // Calculate how to break up spline points so that we travel at our desired reference velocity
            double target_x = 30.0;
            double target_y = s(target_x);
            double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));

            double x_add_on = 0;

            //fill up the rest of our path planner after filling it with previous points, here we will always output 50 points
            for( int i=1; i<=50-previous_path_x.size(); i++ ) {
                double N = (target_dist/(.02*ref_vel/2.24));
                double x_point = x_add_on + (target_x)/N;
                double y_point = s(x_point); 

                x_add_on = x_point;

                double x_ref = x_point;
                double y_ref = y_point;

                // car -> map xform
                x_point = (x_ref * cos(ref_yaw) - y_ref*sin(ref_yaw));
                y_point = (x_ref * sin(ref_yaw) + y_ref*cos(ref_yaw));

                x_point += ref_x;
                y_point += ref_y;

                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
            }
*/

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	json msgJson;
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

















































































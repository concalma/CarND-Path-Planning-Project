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
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
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

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
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
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
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
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
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


class PathPlanner {
    enum State {
        FREE_DRIVING, EXECUTING_PATH
    };

    enum LaneCandidate {
        LEFT, CURRENT, RIGHT
    };



    State state;


    int prev_size = 0;
    double car_x, car_y, car_s, car_d, car_yaw, car_speed;
    vector<double> previous_path_x, previous_path_y;
    vector<vector<double>> sensor_fusion;

    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;
    double seconds_lookahead = 3;
    double max_speed = 49;
    double max_acc = .160;


    public:
       PathPlanner(
               vector<double> &map_waypoints_x,
               vector<double> &map_waypoints_y,
               vector<double> &map_waypoints_s,
               vector<double> &map_waypoints_dx,
               vector<double> &map_waypoints_dy) {

           this->map_waypoints_x = map_waypoints_x;
           this->map_waypoints_y = map_waypoints_y;
           this->map_waypoints_s = map_waypoints_s;
           this->map_waypoints_dx = map_waypoints_dx;
           this->map_waypoints_dy = map_waypoints_dy;

           this->state = FREE_DRIVING;
       }

       struct Path {
           vector<double> x_vals;
           vector<double> y_vals;
           vector<double> s_vals;
           vector<double> d_vals;
           double target_car_s; // s value from which lane changes has been accomplished
           int locked_limit; // when locked path, this is the index up to when this needs to be executed
           int locked_idx; // current playback point for locked path
       };
       struct PathCost {
           Path path;
           double cost;
       };
       Path lockedPath;
       // returns a path plan 
       Path plan( double car_x, double car_y, double car_s, double car_d,  double car_yaw, double car_speed, vector<double> previous_path_x, vector<double> previous_path_y, vector<vector<double>> sensor_fusion ) {
           this->car_x = car_x;
           this->car_y = car_y;
           this->car_s = car_s;
           this->car_d = car_d;
           this->car_yaw = car_yaw;
           this->car_speed = car_speed;
           prev_size = previous_path_x.size();
           this->previous_path_x = previous_path_x;
           this->previous_path_y = previous_path_y;
           this->sensor_fusion = sensor_fusion;
           lockedPath.locked_idx += 50- prev_size;

           Path ret;
           PathCost pc, pcl, pcr;

           switch(state) {
               case FREE_DRIVING:
                   pcl = check_candidate_path(LEFT);
                   pc  = check_candidate_path(CURRENT);
                   pcr = check_candidate_path(RIGHT);


                   printf("%f(%d), %f, %f(%d), car_d: %f\n", pcl.cost, pcl.path.locked_idx, pc.cost, pcr.cost, pcr.path.locked_limit, car_d );

                   // if current path is still lowest cost. keep going
                   if(pc.cost < pcl.cost && pc.cost < pcr.cost ) {
                       // keep going
                   } else if( pcr.cost < pc.cost && pcr.cost < pcl.cost ) {
                    // lock in path and execute change right
                       printf("change right\n");
                       lockedPath = pcr.path;
                       state = EXECUTING_PATH;
                       playbackLocked(ret);
                       return ret;
                   } else if( pcl.cost < pc.cost && pcl.cost < pcr.cost ) {
                       printf("change left\n");
                       lockedPath = pcl.path;
                       state = EXECUTING_PATH;
                       playbackLocked(ret);
                       return ret;
                   } else {
                       // keep going  
                   }

                   // copy just 50 points
                   for(int i=0;i<50;i++) {
                    ret.x_vals.push_back( pc.path.x_vals[i] );
                    ret.y_vals.push_back( pc.path.y_vals[i] );
                   }

                   break;
               case EXECUTING_PATH:
                   playbackLocked(ret);
                   break;
           }

           return ret;
       }

       void playbackLocked(Path &ret) {
           printf("playing back locked path: %d-%d, %d\n", lockedPath.locked_idx, lockedPath.locked_idx+50, lockedPath.locked_limit );
           // copy just 50 points
           for(int i=0;i<50;i++) {
               ret.x_vals.push_back( lockedPath.x_vals[lockedPath.locked_idx + i ] );
               ret.y_vals.push_back( lockedPath.y_vals[lockedPath.locked_idx + i ] );
           }

           if( lockedPath.locked_idx > lockedPath.locked_limit ) {
               //disengage lock
               state = FREE_DRIVING;
           }
       }

  
       Path get_path( int dest_lane, double target_speed ) {
           printf("%d, %f, %d\n",dest_lane, target_speed, prev_size);
          	vector<double> ptsx;
            vector<double> ptsy;
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);
            double ref_vel = 0;
            //
            // Do I have have previous points
            if ( prev_size < 2 ) {
                // Initial state
                double prev_car_x = car_x - cos(car_yaw);
                double prev_car_y = car_y - sin(car_yaw);

                ptsx.push_back(prev_car_x);
                ptsx.push_back(car_x);

                ptsy.push_back(prev_car_y);
                ptsy.push_back(car_y);

                printf("initial state\n");

            } else {
                // Use the last two points.
                ref_x = previous_path_x[prev_size - 1];
                ref_y = previous_path_y[prev_size - 1];

                double ref_x_prev = previous_path_x[prev_size - 2];
                double ref_y_prev = previous_path_y[prev_size - 2];
                ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

                ptsx.push_back(ref_x_prev);
                ptsx.push_back(ref_x);

                ptsy.push_back(ref_y_prev);
                ptsy.push_back(ref_y);

                // get last velocity from these two points
#define DIST(a,b) (sqrt((a)*(a)+(b)*(b)))
                ref_vel = 2.24*DIST(ref_y-ref_y_prev, ref_x-ref_x_prev) / 0.02;  // distance last step / time_step
                // printf("%f\n", ref_vel);
            }

            // Setting up target points in the future.
            vector<double> next_wp0 = getXY(car_s + 30, 2 + 4*dest_lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp1 = getXY(car_s + 60, 2 + 4*dest_lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp2 = getXY(car_s + 90, 2 + 4*dest_lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);

            ptsx.push_back(next_wp0[0]);
            ptsx.push_back(next_wp1[0]);
            ptsx.push_back(next_wp2[0]);

            ptsy.push_back(next_wp0[1]);
            ptsy.push_back(next_wp1[1]);
            ptsy.push_back(next_wp2[1]);

            // Making coordinates to local car coordinates.
            for ( int i = 0; i < ptsx.size(); i++ ) {
              double shift_x = ptsx[i] - ref_x;
              double shift_y = ptsy[i] - ref_y;

              ptsx[i] = shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw);
              ptsy[i] = shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw);
            }

            // Create the spline.
            tk::spline s;
            s.set_points(ptsx, ptsy);



            // Output path points from previous path for continuity.
            Path ret; //TODO
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
            for ( int i = 0; i < prev_size; i++ ) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            // Calculate distance y position on 30 m ahead.
            double target_x = 30.0;
            double target_y = s(target_x);
            double target_dist = sqrt(target_x*target_x + target_y*target_y);

            double x_add_on = 0;

            // TODO. make sure car_speed is proper

            for( int i = 1; i < 1250 - prev_size; i++ ) {
              if( target_speed > ref_vel ) {
                  ref_vel += max_acc;
                  if(ref_vel > max_speed)
                      ref_vel = max_speed;
              } else { 
                  ref_vel -= max_acc;
                  if(ref_vel < target_speed)
                      ref_vel = target_speed;
              }

              double N = target_dist/(0.02*ref_vel/2.24);
              double x_point = x_add_on + target_x/N;
              double y_point = s(x_point);

              x_add_on = x_point;

              double x_ref = x_point;
              double y_ref = y_point;

              x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
              y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);

              x_point += ref_x;
              y_point += ref_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);
            }

            ret.x_vals = next_x_vals;
            ret.y_vals = next_y_vals;

            // calculate freent
            ret.target_car_s = car_s + 60;
            ret.locked_limit = 0;
            ret.locked_idx = 0;

            for( int i=0; i<ret.x_vals.size(); i++ ) {
                vector<double> fr = getFrenet( ret.x_vals[i], ret.y_vals[i], ref_yaw /*TODO*/, map_waypoints_x, map_waypoints_y );

                ret.s_vals.push_back(fr[0]);
                //printf("%f,", ret.s_vals[i] );
                ret.d_vals.push_back(fr[1]);

                if( ret.locked_limit==0 && ret.s_vals[i]> ret.target_car_s )  {
                    printf("set locked_limit=%d, %f,%f \n", i, ret.s_vals[i], ret.target_car_s );
                    ret.locked_limit = i;
                }
            }

            //printf("\n");
            return ret;
       }



       PathCost check_candidate_path( LaneCandidate target_lane_candidate ) {
           PathCost ret;
           ret.cost = 0;

           // TODO. manage target lanes here
           int current_lane = (int)(car_d/4);   // 0..4 -> 0, 4..8->1, 8..12->2
           int target_lane = 0; //current_lane;

           switch(target_lane_candidate) {
               case CURRENT:
                   target_lane = current_lane;
                   break;
               case LEFT:
                   switch(current_lane) {
                    case 0:
                        ret.cost = 2;
                        return ret;
                    case 1:
                        target_lane = 0;
                        break;
                    case 2:
                        target_lane = 1;
                        break;
                   }
                   break;
               case RIGHT:
                   switch(current_lane) {
                    case 0:
                        target_lane = 1;
                        break;
                    case 1:
                        target_lane = 2;
                        break;
                    case 2:
                        ret.cost = 2;
                        return ret;
                   }
                   break;
           }

           double target_speed = 49;
           // check to see if there's a car ahead of us in target lane 
           // and update target speed accordingly
           double min_car_s = std::numeric_limits<double>::max();
            for ( int i = 0; i < sensor_fusion.size(); i++ ) {
                double d = sensor_fusion[i][6];
                if( d<0 || d>12 ) 
                    continue;  // d out of bounds

                int sensed_lane = (int)(d/4);   // 0..4 -> 0, 4..8->1, 8..12->2


                if ( sensed_lane == target_lane ) {
                    //  speed.
                    double vx = sensor_fusion[i][3];
                    double vy = sensor_fusion[i][4];
                    double check_speed = sqrt(vx*vx + vy*vy);
                    double check_car_s = sensor_fusion[i][5];
                    //check_car_s += ((double)prev_size*0.02*check_speed);

                    // figure out car speed of car ahead of us
                    if( check_car_s > car_s &&  check_car_s-car_s < 70 ) {
                        if( check_car_s < min_car_s ) {
                            min_car_s = check_car_s;
                            target_speed = check_speed;
                            // printf("car ahead speed: %f, %d, %f, %f\n", target_speed, sensed_lane, check_car_s, car_s );
                        }
                    }


                } 
            }

            // TODO. get target_speed from vehicle ahead of us in target lane

            ret.path = get_path( target_lane, target_speed );

            double cost = 0;
            cost += target_lane_candidate!=CURRENT ? 0.8 : 1 ; // penalize lane change a bit

            double sdiff = ret.path.s_vals[1000]-ret.path.s_vals[0];
            ret.cost = 1/ (sdiff*cost);  // the higher the car_s, the lower the cost


            // printf("car_s: %f, car_d: %f... %f,%f\n", car_s, car_d, ret.path.s_vals[0], ret.path.d_vals[0]);


            // simulate path in search of colisions
            // path is valid
            // path is invalid
            // for( step in path ) {
            //    for( every car sensed ) 
            //        check for collision
            // }

            return ret;
          }
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
  string map_file_ = "../data/highway_map.csv";
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

  // Car's lane. Stating at middle lane.
  int lane = 1;

  // Reference velocity.
  double ref_vel = 0.0; // mph

  PathPlanner pp(map_waypoints_x, map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy );

  h.onMessage([&ref_vel, &lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &pp]
    (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {



    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
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

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            PathPlanner::Path path = pp.plan( car_x, car_y, car_s, car_d, car_yaw, car_speed, previous_path_x, previous_path_y, sensor_fusion );

            json msgJson;

          	msgJson["next_x"] = path.x_vals;
          	msgJson["next_y"] = path.y_vals;

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

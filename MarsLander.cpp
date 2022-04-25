#include <iostream>
#include <string>
#include <vector>
#include <cmath>


#define DEBUG
//#define DEEP_DEBUG
using namespace std;

// Global variable - gravity
double const g = 3.72;
// Global variable - max vertical speed for landing
double const landing_v_speed_limit = -35; 
// Global variable - max horizontal speed for landing
double const landing_h_speed_limit = 18;
// Global variable - rotation angle limit
int angle_limit = 90;
// Global variable epsilon - math error
double eps = 0.03;
// Global variable ve - speed epsilon
double ve = 3;

int accuracy = 330;

void PP(int const & x, int const & y, bool last = true) {
    cerr<<"["<<x<<";"<<y<<"]";
    if (!last) cerr<<"<-";
}

template<class Z>
ostream & operator<<(ostream & os, vector<Z> const & v) {
    for (auto const & item : v) os<<"-->"<<item;
    return os;
}

struct Point {
    Point() {}
    Point(int x, int y) : x(x), y(y) {}

    int x;
    int y;
};

struct RouteParameters {
    RouteParameters(int x0, int x1, int y0, int y1,
    int vx0=0, int vx1=0, int vy0=0, int vy1=0, 
    int rotate0=0, int power1=0
    ) : x0(x0), x1(x1), y0(y0), y1(y1), 
    vx0(vx0),vx1(vx1),vy0(vy0),vy1(vy1),
    rotate0(rotate0), power1(power1)
    {}


    int x0,x1;
    int y0,y1;
    int vx0,vx1;
    int vy0,vy1;
    int power0,power1;
    int rotate0, rotate1;
    int rotation0,rotation1;
};

void p_swap(Point & p1, Point & p2) {
    Point t = p1;
    p1.x = p2.x;
    p1.y = p2.y;
    p2.x = t.x;
    p2.y = t.y;
}

ostream & operator<<(ostream & os, Point const & p) {
    os << "[" << p.x << ";" << p.y<<"]"<<endl;
    return os;
}

class Surface {
    public:
    explicit Surface(int n=0) : _s_size(n) {
        _max_height = -3000;
        _plato_length = 0;
        _plato_height = 3000;
    }

    Surface(Surface const & sf) {
        Surface * n = new Surface(sf.get_size());
        this->points = sf.points;
        this->_s_size = sf._s_size;
        this->_max_height = sf._max_height ;
        this->_plato_start = sf._plato_start;
        this->_plato_length = sf._plato_length; 
        this->_plato_height = sf._plato_height;

    }

    Surface & operator=(Surface const & sf) {
        Surface * n = new Surface(sf.get_size());
        return *n;
    }

    void add_point(Point const & p) {
        points.push_back(p);
    }
    int get_size() const {
        return _s_size;
    }

    int get_max_height() const {
        return _max_height;
    }

    int get_max_height_x() const {
        return _max_height_x;
    }

    int get_plato_start() const {
        return _plato_start;
    }

    int get_plato_end() const {
        return _plato_start+_plato_length;
    }

    int get_plato_height() const {
        return _plato_height;
    }

    int get_height(int x) const {
        int _prev_x = 0;
        int _prev_y = 0;
        double angle = 0;
        int height = 0;
        for (auto const & item : points) {
            if (item.x > x) {
                if (item.y == _prev_y) {
                    height = item.y;
                    break;
                }
                else
                    angle = (item.x-_prev_x)/((double)(item.y-_prev_y));
                double dx = x - _prev_x;
                height = (_prev_y + dx/angle);
                break;
            }
            _prev_x = item.x;
            _prev_y = item.y;
        }
        return height;
    }

    vector<Point> const & get_points() const {
        return points;
    }
    void build_surface();

    private:
    void print_surface_result() {
        cerr << "Max: " << _max_height
        << " [" << _plato_start
        << "->" << _plato_start+_plato_length 
        << "] at " << _plato_height
        << endl;
    }

    private:
    vector<Point> points;
    int _s_size;
    int _max_height;
    int _max_height_x;
    int _plato_start;
    int _plato_length; // 1000 - minimum
    int _plato_height;
    

};

void Surface::build_surface() {
    int _prev_y = -1;
    int _prev_x = -1;
    bool isFound = false;
    for (auto & p : points) {
        if (p.y > _max_height) {
            _max_height = p.y;
            _max_height_x = p.x;
        }
        if (!isFound && _prev_y == p.y) { // start plato
            int len = p.x - _prev_x;
            if (len >= 1000) { // we expect for an UNIQUE plato
                _plato_start = _prev_x;
                _plato_length = len;
                _plato_height = p.y;
                isFound = true;
            }
        }
        _prev_x = p.x;
        _prev_y = p.y;
    }
    print_surface_result();
}

class Ship {
    public:
    Ship(Surface const & sf) : surface(new Surface(sf)) {
        surface->build_surface();
        isDestinationReached = true;
        max_angle = acos(g/4.) * 180.0 / M_PI; // arccos(max_acc/g)
        cerr<<max_angle<<endl;
        isObjectBuilded = false;
        calls_1=calls_2=0;
    }

    void Systems(int x, int y, int h_speed,
    int v_speed, int fuel, int rotate, int power) {
        this->x = x;
        this->y = y;
        this->h_speed = h_speed; 
        this->v_speed = v_speed; 
        this->fuel = fuel; 
        this->rotate = rotate; 
        this->power = power; 
        Processing();
    }

    int GetPower() const;
    int GetRotate() const;

    ~Ship() {
        delete surface;
    }
    private:
    
    void Processing();
    void CurrentState();
    void Landing() {
        cerr <<"Landing"<<endl;
        _prev_power = 4;
        
        if (v_speed > 0.3*landing_v_speed_limit)
            _prev_power = 3;
        if (abs(h_speed) > landing_h_speed_limit-1) {
            angle_limit = max_angle;
            // from the right to left
            if (h_speed < 0)
                _prev_rotate = ((-1) * angle_limit );
            // from the left to right
            else
                _prev_rotate = (angle_limit );
        }
        else {
            _prev_rotate = 0;
        }
        if ((y - route[route_point].y) < accuracy * 0.5)
            _prev_rotate = 0;

    }
    void RouteCorrection(int vx, int vy);
    void CalculateDestination(int & sx, int & sy, int & angle, int & a, int & t) {
        
        sx = x + h_speed * t + a * t * t * sin(angle*M_PI/180 - 90);
        sx = sx > 0 ? sx : 0;
        sy = (y + v_speed * t + (a * cos(angle*M_PI/180) - g) * t * t/ 2);
        sy = sy > 0 ? sy : 0;

    }
    bool RouteModule() ;
    void CalculateRoute();
    void CalculateStrategy();
    void SqrtCalculating(double defPower = 2) {

        if (defPower > 4) return;
        double a,b,c,x1=1,x2=1;
        double pwr = defPower;
        double tga= (route[route_point].y - y)/(route[route_point].x - x);
        double g = -3.71;
        double vv = 0;
        double vh = 0;
        
        
        double p = (vh * tga - vv + g) / pwr;
        //cerr<< "parameter p="<<p<<endl;
        a = tga*tga +1; 
        //cerr<<"a = "<<a<<endl;
        b = 2 * p * tga;
        //cerr<<"b = "<<b<<endl;
        c = p*p - 1;
        //cerr<<"c = "<<c<<endl;
        
        double sd = (b*b - 4*a*c);
        if (sd < 0 ) {
            cerr<< "Cant be reached : D="<<sd<<"\n###\n";
            SqrtCalculating(++pwr);
        }
        else {
            sd = sqrt(sd);
            cerr<<"sqrt(D)="<<sd<<endl;
            x1 = (0 - b - sd) / (2 * a);
            cerr << "x1=" << x1 << endl;
            cerr << "compare "<<tga<<" && "<<(vv + pwr * cos(asin(x1))-g)/(vh + pwr * x1 )<<endl;
            cerr << "The angle is "<<asin(x1) * 180 / M_PI<<endl;
            if (sd != 0) {
                x = (0 - b + sd) / (a * 2);
                cerr << "x2=" << x2 <<endl;
                cerr << "compare "<<tga<<" && "<<(vv + pwr * cos(asin(x2))-g)/(vh + pwr * x2 )<<endl;
                cerr << "The angle is "<<asin(x2) * 180 / M_PI<<endl;
                
            }
        _prev_rotate = min(x1*180/M_PI,x2*180/M_PI);
        }
        _prev_power = pwr;
    }
    void CorrectRoute() ;
    bool BuildRoute() {
        cerr<<"BuildRoute()\n";
        route.push_back(Point{x,y}); 
        for (int i = 0; i < route.size()-1; ++i) {
            // cerr<<"==========\n"<<" Distance == "
            //     <<DistanceToPoint(route[i+1], route[i])<<endl;
            if (DistanceToPoint(route[i+1], route[i]) > 5*accuracy) {
                Point np = {
                    (route[i].x + route[i+1].x)/2, GetLineHeight(route[i], route[i+1],
                    (route[i].x + route[i+1].x)/2) + accuracy/3
                };
            //    cerr<<route[i]<<route[i+1]<<"added:"<<np<<"##^^##"<<endl;
                vector<Point>::iterator it = route.begin() + i;
                route.insert(it+1,np);
            }
        }
        cerr<<"@@@@@@@@@@@@@@@@@@@@\n"<<route<<"@@@@@@@@@@@@@@\n";
        
        return true;
    }
    void BuildCurve(Point const & lhs, Point const & rhs) { 
        // set others points by BuildCurve
        // principals: 
        // 1) we get two point: ship position and landing position
        // bild line between them
        // 2) Going throught obstacles and check: do we need correct route?
        // 3) If yes, then call <build curve> with
        // obstacle and ship position recursive.
        // 3.1) Obstacle point added to route
        // 4) If not, then return (no need to add ship point to route)
        //
        cerr<<"BuildCurve()\n";
        // cycle direction
        int i = (isNeedRight ?  obstacles.size()-1 : 0 );
        int n = (isNeedRight ?  0 : obstacles.size() );
        for (; (i!=n) ; (isNeedRight ? --i : ++i)) {
            if ((rhs.x > obstacles[i].x) && (obstacles[i].x > lhs.x) ) {
                if (obstacles[i].y > GetLineHeight(lhs,rhs,obstacles[i].x)) {
                    route.push_back({obstacles[i].x, obstacles[i].y+accuracy});
                    if (isNeedRight)
                        BuildCurve(lhs, obstacles[i]);
                    else
                        BuildCurve(obstacles[i], rhs);
                    break;
                }
            }
        }
        
    }
    int GetLineHeight(Point lhs, Point rhs, int const x) {
        // also depends on bypass direction
        // just swap lhs<->rhs
        double xtg = 0;
        //if (!isNeedRight) p_swap(lhs,rhs);
        xtg = lhs.y + (x-lhs.x)*(rhs.y - lhs.y)/static_cast<double>(rhs.x - lhs.x);
        #ifdef DEEP_DEBUG
            cerr <<"{ xtg = "<<xtg <<"}";
        #endif
        return (xtg);
    }
    void SetLandDestination() {
        Point landing;
        landing.x = (surface->get_plato_end() - surface->get_plato_start())/3;
        if (!isNeedRight)
            landing.x = 2*landing.x;
        landing.x = landing.x + surface->get_plato_start();
        landing.y = surface->get_height(landing.x) + accuracy/2;
        route.push_back(landing);
    
    }
    void CheckDestination() {

        if ((DistanceToPoint(route[route_point]) < accuracy))
        {
            cerr << "#"<<route_point<<" Route point achieved!"<<endl;
            if (route.size()>1) {
                --route_point;
                route.pop_back();
            }
            if (route.size()==1) {
                cerr<<"set isLowHeight=true"<<endl;
                isLowHeight = true;
            }

            cerr<<"Next point: "<<route[route_point];
            distance = DistanceToPoint(route[route_point]);
        }
    }
    void CalculateMovementDirection() {
        // this initial point for analysis
        if (x < surface->get_plato_start())
            isNeedRight = true;
        else if (x > surface->get_plato_end())
            isNeedRight = false;
        cerr<<"IsNeedRight = "<<isNeedRight<<endl;
    }
    void CalculateObstacles() {
            // if under the max AND max is obstacle
            // compare x and sureface->_max_height_x
            // ship righter: plato->obstacle->ship
            // obstacle have more height than platoe
            cerr<<"CalculateObstacles()\n";
            int obstacle_start = 0;
            int obstacle_end = 0;
            int _prev_y = 0, _prev_x = 0;
            // check all points
            for (auto const & p : surface->get_points()) {
                // if found first upstage
                if (p.y > _prev_y && !obstacle_start)
                        obstacle_start = p.x;
                else {
                    // if found second upstage:
                    // point height less than previous x height
                    if (p.y <= _prev_y && obstacle_start) {
                        obstacle_end = _prev_x; 
                        obstacle_start = surface->get_height(_prev_x);
                        obstacles.push_back({obstacle_end,obstacle_start});
                        obstacle_start = 0;
                        obstacle_end = 0;
                    }
                }
                _prev_x = p.x;
                _prev_y = p.y;
            }
            isObjectBuilded = true;
    }
    int DistanceToPoint(Point const & p) {
        int d = sqrt((x - p.x)*(x - p.x) + (y - p.y)*(y - p.y));
        return d;
    }
    int DistanceToPoint(Point const & p1, Point const & p2) {
        int d = sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
        return d;
    }
    int FuelResources(int _time, int power = 4) {
        return _time*power;
    }
    void CorrectLastPoint();


    private:
    Surface * surface;
    Point destination;
    int route_point;
    vector<Point> route;
    vector<RouteParameters> RP;
    int x;
    int y;
    int max_angle;
    int h_speed; 
    int v_speed; 
    int fuel; 
    int rotate; 
    int power; 
    int distance;
    int _prev_rotate, _prev_power, _prev_time;
    vector<Point> obstacles;
    int calls_1, calls_2;

    // future commands
    bool isNeedUp;
    bool isLowHeight;
    bool isNeedRight;
    bool isObjectBuilded;
    bool isDestinationReached;
};

void Ship::CalculateRoute() {
        // set finish point
        cerr<<"CanlculateRoute()\n";
        SetLandDestination();
        // If going to the right, then final distanation always righter
        // That mean, *** right hand side must be route[0]
        Point rhs = {x,y-accuracy}, lhs = route[0];
        if (isNeedRight) p_swap(rhs,lhs);
        BuildCurve(lhs, rhs);
        // ***  lhs ^    ^ rhs
        // backpropagation 
        BuildRoute();

}

void Ship::Processing()  {

        if (!isObjectBuilded) {
            CalculateMovementDirection();
            CalculateObstacles();
            CalculateRoute();
            isObjectBuilded = true;
            route_point = route.size()-1;
            CalculateStrategy();
        }
        
        CurrentState();
        rotate = _prev_rotate;
        power = _prev_power;
        cerr << "rot: "<<rotate<<" pow: "<<power<<endl;
        CheckDestination();
}

void Ship::CorrectLastPoint() {
    if (abs(h_speed) > 2 * landing_h_speed_limit) { 
        // 6 - speed of ship rotation
        int tx = 6 + (abs(h_speed) - landing_h_speed_limit)/4;
        cerr<<"time to reach point "<<tx<<endl;
        int sign = (route[route_point-1].x - x < 0 ? 1 : -1);
        int rt = h_speed * tx + 2 * sign * 0.7 * tx * tx;
        rt = rt * 1.1;
        cerr<<"need new coord "<<rt<<endl;
        route[route_point-1].x += rt;
    }
}

void Ship::CurrentState()  {

        // Function calculate necessary parameters for the reaching point.
        // System has limited resources.
        // If ship can do nothing - do nothing.
        // If need correcting -> Correcting.
        // Let's go!
        cerr<<"CurrentState"<<endl;
        power = 0;
        int dx = route[route_point].x - x;
        int dy = route[route_point].y - y;
        
        cerr<<"Moving to "<<route[route_point]<<"RS:"<<route.size()<<endl;
        // Does ship reach the point? Call route module!
        if (!RouteModule()) {
        // If ship need a correct, making corrections, till
        // RouteModule() returns (false)
            RouteCorrection(dx, dy);
        }
        else {
            // do nothing, if nothing needs to do
            cerr<<"Do nothing"<<endl;
            _prev_power = (_prev_power == 3 ? 4 : 3);
            _prev_rotate = 0;
            cerr<<"dnrs:"<<route.size()<<endl;

            
            if (route.size() <= 1)
                Landing();
            // but if need to land, run landing program

        }   
        _prev_rotate = 75;
        _prev_power = 0;   
        if (x < 5400)
        {
            _prev_power = 4;
            _prev_rotate = -90;
        }  
        if (fuel < 565)
            _prev_rotate = -10;

}



void Ship::CalculateStrategy() {
    // For each route interval ship must build strategy. 
    // Strategy calculating border condition between neighbor 
    // Route points allow reach every point with conditions, 
    // That allow to reach next point
    for (int i = route_point; i > 0; --i) {
        
    }

}

void Ship::CorrectRoute() {

}

bool Ship::RouteModule() {
        // The way S defined by speed's vector's direction as:
        // Sx = S0 + Vxt, Sy = S0 + Vyt + g*t*t/2
        // Our engines need to correct speed's vector's. 
        // In That calculation (if that called from CurrentState() function) 
        // they are not needed
        
        int t;
        if (h_speed != 0)
            t = (route[route_point].x - x) / h_speed; 
        else 
            t = (route[route_point].y - y) / ( v_speed ? v_speed : 1);
        if (t < 0 && !isLowHeight) return false;
        int sx = x + h_speed * t;
        int sy = y+v_speed*t+(power-g)*t*t/2;

        double dev = sqrt((sx-route[route_point].x)*(sx-route[route_point].x) 
        +((sy-route[route_point].y)*(sy-route[route_point].y)));

        cerr<<"expected deviation: "<<dev<<endl;
        cerr<<"iSLowHeight = "<<isLowHeight<<endl;
        if (dev < accuracy || isLowHeight) {
            return true;
        }
        return false;
}

void Ship::RouteCorrection(int vx, int vy) {
        cerr<<"RouteCorrection()"<<endl;
        
        CorrectRoute();
        
        power = _prev_power;
        rotate = _prev_rotate;
}


int main()
{
    int surface_n; // the number of points used to draw the surface of Mars.
    Surface * surface = new Surface(surface_n);
    cin >> surface_n; cin.ignore();
    for (int i = 0; i < surface_n; i++) {
        int land_x; // X coordinate of a surface point. (0 to 6999)
        int land_y; // Y coordinate of a surface point. By linking all the points together in a sequential fashion, you form the surface of Mars.
        cin >> land_x >> land_y; cin.ignore();
        Point p(land_x, land_y);
        surface->add_point(p);
    }

    Ship * ship = new Ship(*surface);

    // game loop
    while (1) {
        int x;
        int y;
        int h_speed; // the horizontal speed (in m/s), can be negative.
        int v_speed; // the vertical speed (in m/s), can be negative.
        int fuel; // the quantity of remaining fuel in liters.
        int rotate; // the rotation angle in degrees (-90 to 90).
        int power; // the thrust power (0 to 4).
        cin >> x >> y >> h_speed >> v_speed >> fuel >> rotate >> power; cin.ignore();
        
        ship->Systems( x,  y,  h_speed,  v_speed,  fuel,  rotate,  power);
        
        // 2 integers: rotate power. rotate is the desired rotation angle (should be 0 for level 1), power is the desired thrust power (0 to 4).
        //cout << (abs(h_speed) < 20 ? 0:-10) <<" "
        //<< (v_speed < 10 ? (v_speed > 0 ? 3 : 4) : 4)<<endl;
        cout<<ship->GetRotate()<< " "<< ship->GetPower() << endl;
        
    }

    delete surface;
    delete ship;
}


int Ship::GetPower() const {
    return power;
}

int Ship::GetRotate() const {
    return rotate;
}

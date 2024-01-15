#define _USE_MATH_DEFINES
#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <math.h>
#include <functional>
#include <iterator>
#include <chrono>
#include <iomanip>
// #include <corecrt_math_defines.h>
// #include "test"
// #include "matplotlibcpp.h"

typedef unsigned int uint;


using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;
// using std::chrono::milliseconds;
using std::to_string;
using std::string;
using std::stringstream;
using doubles = std::vector<double>;

struct Vector;
using Point = Vector;
using Velocity = Vector;
using Vectors = std::vector<Vector>;

class Particle;
class Treeparticle;
using Particles = std::vector<Particle>;

enum class NodeType { LEAF, NODE };
class Range;
class Tree;
class Leaf;
class Node;



// * my comment
// ! important
// ? question
// TODO: implement
// // commented out

const double G = 1;
const int printSize = 65;


Vector getForceRec(const Particle& p, const Tree* tree, double angle, double eps);

int clearMyPlot(const string& filename = "myplot");
int writePlotToMyPlot(const doubles& xs, const string& filename = "myplot");
int writePlotToMyPlot(const Vectors& xs, const string& filename = "myplot");
doubles linspace(double start, double end, int numPoints = 60, bool inclusive = true);
doubles logspace(double start, double end, int numPoints = 60, double base = 2);
doubles logspace0(double start, double end, int numIntervals = 60, double base = 2);
template <typename T1, typename T2>
std::vector<T2> map(const std::vector<T1>& xs, std::function<T2(T1)> f);
void runPython(string filename);
void boxedPrint(stringstream& ss); 
void boxedPrint(string s);
void paddedPrint(string s="", char fill='-', char end='|', int padLen=printSize);



/* MY TYPES */

/** @brief some operators for simple operations on doubles */
doubles operator+(const doubles& xs, double y) {
    doubles zs(xs.size());
    std::transform(xs.begin(), xs.end(), zs.begin(), [y](double x) {return x+y;});
    return zs;
}
doubles operator+(const doubles& xs, const doubles& ys) {
    const double mymin = std::min(xs.size(), ys.size());
    doubles zs(mymin);
    std::transform(xs.begin(), xs.begin() + mymin, ys.begin(), zs.begin(), std::plus<double>());
    return zs;
}
doubles operator/(const doubles& xs, double scalar) {
    doubles ys(xs.size());
    std::transform(xs.begin(), xs.end(), ys.begin(), [scalar](double x) {return x/scalar;});
    return ys;
}
std::ostream& operator<<(std::ostream& os, const doubles& xs) {
    os << "[";
    for (const auto& x : xs) {
        os << x << ", ";
    }
    os << "]";
    return os;
}
doubles addWithOffset(const doubles& xs, int offset) {
    double size = xs.size()-offset;
    doubles ys(size);
    std::transform(xs.begin(), xs.begin()+size, xs.begin()+offset, ys.begin(), std::plus<double>());
    return ys;
}

/** @brief a 3d-vector class with simple operator overloading, a print and distance function */
struct Vector {
    double x, y, z;

    double dist() const {
        // std::cout << "dist: " << sqrt(x*x + y*y + z*z) << std::endl;
        // std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
        // std::cout << "squared: " << x*x + y*y + z*z << std::endl;
        return sqrt(x*x + y*y + z*z);
    }
    Vector operator+(const Vector& other) const {
        return Vector{x+other.x, y+other.y, z+other.z};
    }
    Vector operator-(const Vector& other) const {
        return Vector{x-other.x, y-other.y, z-other.z};
    }
    Vector operator-() const {
        return Vector{-x, -y, -z};
    }
    Vector operator*(double scalar) const {
        return Vector{x*scalar, y*scalar, z*scalar};
    }
    friend Vector operator*(double scalar, const Vector& v) {
        return v*scalar;
    }
    /** dot product */
    double operator*(const Vector& other) const {
        return x*other.x + y*other.y + z*other.z;
    }
    Vector operator/(double scalar) const {
        return Vector{x/scalar, y/scalar, z/scalar};
    }
    Vector& operator+=(const Vector& other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }
    Vector& operator-=(const Vector& other) {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    
    double openingAngle(const Vector& other) const {
        return std::acos(*this*other / dist() / other.dist());
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector& v) {
        os << "V (" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }
};


class Particle {
protected:
    static double totalMass;
    static int idCount;
    static bool sorted;
    // // static Particles sortedParticles; // ! reset when particles moved

    int id;
    double mass;
    Point p;
    Velocity v;
    bool central;
public: 
    Particle(int id, double mass, Point p, Velocity v, bool central) : id(id), mass(mass), p(p), v(v), central(central) {
        if (id!=-1) { totalMass += mass; }
        else if (id==-1) { this->id = idCount; } // set unique id for tree particles
        idCount = std::max(idCount, this->id+1);
    }

    friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
        os << "Particle " << p.id << " with mass " << p.mass << " at position (" << p.p.x << ", " << p.p.y << ", " << p.p.z << ") and velocity (" << p.v.x << ", " << p.v.y << ", " << p.v.z << ")";
        return os;
    }

    bool operator==(const Particle& other) const {
        return id == other.id;
    }

    Vector getForceTo(const Particle& other, double eps) const {
        Vector r12 = -other.p + p;
        double r = r12.dist();
        double force = -G * mass * other.mass / (r*r + eps*eps);
        if (r!=0) return (force/r) * r12;
        else { std::cerr << "r=0 for particles " << *this << " and " << other << std::endl; return Vector{0,0,0}; }
    }

    void setApplyForce(const Vector& force, double dt) {
        sorted = false;

        // v += force/mass * dt;
        // p += v * dt;
        // or
        std::cout << "ID: " << id << "  Force: " << force << std::endl;
        Vector a = force/mass;
        p += v * dt + a * (dt*dt/2);
        v += a * dt;
        // or leapFrog
        // Vector f1{0,0,0}; // calculateForce(); -- can be omitted for speed
        // v += f1/mass * (dt/2);
        // p += v * dt;
        // Vector f2{0,0,0}; // calculateForce(); -- must be left here
        // v += f2/mass * (dt/2);


        throw std::runtime_error("Not implemented yet");
        
    }

    void setApplyLeapFrog1(const Vector& force, double dt) {
        v += force/mass * (dt/2);
        p += v * dt;
    }
    void setApplyLeapFrog2(const Vector& force, double dt) {
        v += force/mass * (dt/2);
    }

    double getDist() const { return p.dist(); }
    double getMass() const { return mass; }
    Point getP() const { return p; }
    int getId() const { return id; }


    /* some static methods */

    /** returns the total mass in O(1) */
    static double getTotalMass() {
        return totalMass;
    }

    /** assuming sorted particles */
    static double maxDistance(const Particles& particles) {
        if (!sorted) throw std::runtime_error("Particles not sorted");
        return particles.back().getDist();
        // // double max = 0;
        // // for (const auto& p : particles) {
        // //     max = std::max(max, p.getDist());
        // // }
        // // return max;
    }

    /** need to be sorted again after every iteration of a full gravity tree solver scheme */
    static void sortParticles(Particles& particles) {
        if (sorted) return;
        // sortedParticles = particles;
        std::sort(particles.begin(), particles.end(), [](const Particle& p1, const Particle& p2) {
            return p1.getDist() < p2.getDist();
        });
        sorted = true;
    }

    static bool isSorted() {
        return sorted;
    }

};
double Particle::totalMass = 0;
int Particle::idCount = 0;
bool Particle::sorted = false;





class Treeparticle : public Particle {
public: 
    Treeparticle() : Particle(-1, 0, Point{0, 0, 0}, Velocity{0, 0, 0}, false) {}
    Treeparticle(double mass, Vector p) : Particle(-1, mass, p, Velocity{0, 0, 0}, false) {}

    void add(const Treeparticle& other) {
        double newMass = mass + other.mass;
        if (newMass==0) return;
        Vector newP = (mass*p + other.mass*other.p) / newMass;
        mass = newMass;
        p = newP;
    }
};

//* DEFINING/DECLARING CLASSES
class Range {
private:
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;

public:
    Range(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), zmin(zmin), zmax(zmax) {}
    bool contains(const Vector& p) const { return xmin <= p.x && p.x <= xmax && ymin <= p.y && p.y <= ymax && zmin <= p.z && p.z <= zmax; }
    int getQuadrant(const Vector& p) const { return (p.x > (xmin+xmax)/2) + 2*(p.y > (ymin+ymax)/2) + 4*(p.z > (zmin+zmax)/2); }
    int getQuadrant(const Range& r) const { return (r.xmin != xmin) + 2*(r.ymin != ymin) + 4*(r.zmin != zmin); }
    Range getChildRange(int quadrant) const;
    double openingAngle(Vector p) const;

    static Range containingBox(const Particles& particles);
    friend std::ostream& operator<<(std::ostream& os, const Range& r) { os << "[(" << r.xmin << "," << r.xmax << ") x (" << r.ymin << "," << r.ymax << ") x (" << r.zmin << "," << r.zmax << ")]"; return os; }
};

class Tree {
protected:
    Range range;
    Node* parent;
public:
    Treeparticle treeParticle;
    Tree(const Range& range, Node* parent) : range(range), parent(parent) {}
    Range getRange() const { return range; }

    virtual ~Tree() = default;
    virtual bool add(const Particle& p) = 0;
    virtual NodeType getType() const = 0;
    virtual void print(std::ostream& os) const = 0;
    
    static Tree* makeRoot(const Range& range);
    friend std::ostream& operator<<(std::ostream&os, const Tree& tree) { tree.print(os); return os; };
};

class Leaf : public Tree {
private:
    std::vector<const Particle*> ps;
public:
    Leaf(const Range& range, Node* parent) : Tree(range, parent) {}

    bool add(const Particle& p) override;
    NodeType getType() const override { return NodeType::LEAF; }
    void print(std::ostream& os) const override { os << "Leaf with " << getNrParticles() << " particles in range " << range; }

    std::vector<const Particle*> getParticles() const { return ps; }
    int getNrParticles() const {return ps.size(); }

    Vector getForceFrom(const Particle& p, double eps) const {
        Vector force{ 0,0,0 };
        for (const auto& p2 : ps) {
            if (p == *p2) {
                continue;
            }
            force += p.getForceTo(*p2, eps);
        }
        return force;
    }
};

class Node : public Tree {
protected:
    Tree* children[8];
public:
    Node(const Range& range, Node* parent) : Tree(range, parent) {
        for (uint i=0; i<8; i++) children[i] = new Leaf(range.getChildRange(i), this);
    }
    
    ~Node() override { for (auto& child : children) delete child; }
    bool add(const Particle& p) override { 
        return children[range.getQuadrant(p.getP())]->add(p); 
    }
    NodeType getType() const override { return NodeType::NODE; }
    void print(std::ostream& os) const override { os << "Node with range " << range; }

    void setNewChild(Node* child) { children[range.getQuadrant(child->range)] = child; }
    Tree* getChild(int quadrant) const { return children[quadrant]; }

    Vector getForceFrom(const Particle& p, double angle, double eps) const {
        Vector force{ 0,0,0 };
        for (int i=0; i<8; i++) {
            force += getForceRec(p, children[i], angle, eps);
        }
        return force;
    }
};


Particles readIn(string filename, int maxParticles = 1e7) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    string line;
    getline(inputFile, line);           // skip first line

    Particles particles;
    int i=0;
    while (getline(inputFile, line)) {
        std::istringstream iss(line);
        int id; iss >> id;
        double mass; iss >> mass;
        // std::cout << mass << std::endl;
        Point p; iss >> p.x >> p.y >> p.z;
        Velocity v; iss >> v.x >> v.y >> v.z;

        particles.push_back(Particle(id, mass, p, v, id==0 && mass==1 && p.x==0 && p.y==0 && p.z==0 && v.x==0 && v.y==0 && v.z==0));

        i++;
        if (i >= maxParticles) break;
    }

    return particles;
}


//* TASK 1.1
doubles computeShells(const Particles& particles, const doubles& shellBoundaries0) {
    doubles shells(shellBoundaries0.size()-1);
    for (const auto& p : particles) {
        // SKIP CENTRAL PARTICLE???
        double dist = p.getDist();

        auto it = std::lower_bound(shellBoundaries0.begin(), shellBoundaries0.end(), dist);
        int shell = it - shellBoundaries0.begin() - 1;

        shells[shell] += p.getMass();
    }

    doubles densities(shells);
    for (uint i=0; i<shells.size(); i++) {
        double volume = 4./3. * M_PI * (pow(shellBoundaries0[i+1], 3) - pow(shellBoundaries0[i], 3));
        densities[i] /= volume;
    }
    return densities;
}

/** @brief compute `rho` according to Hernquist formula */
double rho(double r, double a=.1) {
    double M = Particle::getTotalMass();
    return M/(2*M_PI) * (a/r) * (1/(r+a)/(r+a)/(r+a));
}
/** @brief compute `rho`'s according to Hernquist formula */
doubles computeRhos(const doubles& xs) {
    return map<double, double>(xs, [](double x) {return rho(x);});
}
doubles computePoissonErrors(const doubles& shells, const doubles& lambdas) {
    /* in #of standard-deviations */
    doubles errors(shells.size());
    for (uint i=0; i<shells.size(); i++) {
        errors[i] = (shells[i] - lambdas[i]) / sqrt(lambdas[i]);
    }
    return errors;
}



//* TASK 1.2
/** @brief calculate number of particles within the half-mass radius */
int getHmNr(const Particles& particles) {
    if (!Particle::isSorted()) throw std::runtime_error("Particles not sorted");

    // Particles sorted(particles);

    // std::sort(sorted.begin(), sorted.end(), [](const Particle& p1, const Particle& p2) {
    //     return p1.getDist() < p2.getDist();
    // });

    // // TODO could store the sorted sequence!

    double totalMass = Particle::getTotalMass();


    double mass = 0;
    int i=0;
    while (mass < totalMass/2 && i < (int) particles.size()) {
        mass += particles[i++].getMass();
    }

    // std::cout << "Mass: " << mass << std::endl;
    // std::cout << "i: " << i << std::endl;

    /* Particles 0..i-1 have mass < totalMass/2, Particles 0..i have mass >= totalMass/2 */
    return i-1;
}

/** @brief calculates the half-mass radius `R_hm` */
double getRhm(const Particles& particles) {
    int half = getHmNr(particles);
    // std::cout << "Half: " << half << std::endl;
    // std::cout << "Dist: " << particles[half].getDist() << std::endl;
    return particles[half].getDist();

    // // TODO could also store after it is calculated once!
    // // doubles dists;
    // // for (const auto& p : particles) {
    // //     dists.push_back(p.getDist());
    // // }
    // // nth_element(dists.begin(), dists.begin()+half-1, dists.end());
    // // return dists[half];
}

/** @brief calculates the softening `d = (V/N)^(1/3)` */
double getSoftening(const Particles& particles) {
    double Rhm = getRhm(particles);
    double V = 4./3 * M_PI * Rhm*Rhm*Rhm;
    int N = getHmNr(particles);
    // std::cout << "Rhm: " << Rhm << std::endl;
    // std::cout << "N: " << N << std::endl;

    if (N==0) return 0;
    return pow(V/N, 1./3);
}

/** @brief calculates the direct forces in `O(n^2)`*/
Vectors getDirectForces(const Particles& particles, double softeningMultiplier=1, bool print=true) {
    double softening = getSoftening(particles)*softeningMultiplier;
    // std::cout << softening;

    if (print) {
        std::cout << "|";
        std::cout.flush();
    }
    int printed = 0; int max = printSize-2;
    Vectors forces(particles.size());
    for (uint i=0; i<particles.size(); i++) {
        if (print) {
            if (i*max > printed*particles.size()) {
                std::cout << "=";
                std::cout.flush();
                printed++;
            }
        }
        Particle p1 = particles[i];
        for (uint j=i+1; j<particles.size(); j++) {
            Particle p2 = particles[j];

            Vector force = p1.getForceTo(p2, softening);



            forces[i] += force;
            forces[j] -= force;
        }
    }
    if (print) {
        std::cout << std::setw(max-printed) << std::setfill('=') << std::right << "|" << std::endl;
        std::cout.flush();
    }
    return forces;
}

/** @brief calculates the analytical forces in `O(n)`
  * @note Newton II theorem for spherical potentials
  * @note `vec F(r) = -GM(r)/r^3 * vec r = -GM(r)/r^2 * hat e_r` with `M(r)=mass` within radius `r` */
int getAnalyticalForce(const Particles& particles, const double prec, doubles& shellBoundaries, doubles& analytical, doubles& shellMidpoints) {
    if (!Particle::isSorted()) throw std::runtime_error("Particles not sorted");

    double maxDist = Particle::maxDistance(particles);
    shellBoundaries = logspace0(0.3, maxDist, maxDist/prec, 1.2);
    shellMidpoints = addWithOffset(shellBoundaries, 1)/2;
    shellBoundaries.erase(shellBoundaries.begin());
    // shellBoundaries = linspace(prec, maxDist, maxDist/prec, false);
    // shellMidpoints = shellBoundaries+(-prec/2);

    int pI=0;
    int shellI=0;
    double prefixMass=0;
    while (pI < (int) particles.size()) {

        while (pI < (int) particles.size() && particles[pI].getDist() <= shellMidpoints[shellI]) {
            prefixMass += particles[pI].getMass();
            pI++;
        }
        analytical.push_back(-G*prefixMass / (shellMidpoints[shellI] * shellMidpoints[shellI]));
        shellI++;
        
        if (shellI >= (int) shellMidpoints.size()) break;
    }
    return 0;
}

int getComparisonDirectAnalyticalForce(const Particles& particles, const Vectors& cd, Vectors& scatter, const doubles& shellBoundaries, doubles& direct, doubles& angleOff) {
    if (!Particle::isSorted()) throw std::runtime_error("Particles not sorted");


    int pI=0;
    int shellI=0;
    while (pI < (int) particles.size()) {
        if (particles[pI].getDist()==0) continue;

        double sumMass = 0;
        double sumForce = 0;
        double sumAngleOff = 0;
        int count = 0;
        while (pI < (int) particles.size() && particles[pI].getDist() <= shellBoundaries[shellI]) {
            pI++;
            double openingAngle = cd[pI].openingAngle(-particles[pI].getP());
            double mass = particles[pI].getMass();
            // consider either only radial component or total force
            // double force = -cd[pI].dist() * cos(openingAngle);
            double force = -cd[pI].dist();

            sumMass += mass;
            sumForce += force;
            sumAngleOff += openingAngle;
            count++;

            scatter.push_back(Vector{particles[pI].getDist(), force/mass, openingAngle});
        }
        direct.push_back(sumForce / sumMass);
        angleOff.push_back(count>0 ? sumAngleOff/count : 0);
        shellI++;

        if (shellI >= (int) shellBoundaries.size()) break;
    }
    return 0;
}

/** @brief calculates the crossing time `t_cross`
 *  @returns `t_cross = 2*R_hm / v_c` with `v_c = sqrt(G*M_hm/R_hm)` */
double getTCross(const Particles& particles) {
    double Rhm = getRhm(particles);
    double mass_halfMassRadius = Particle::getTotalMass() / 2;
    double v_c = sqrt(G*mass_halfMassRadius/Rhm);

    double t_cross = 2*Rhm / v_c;       // ? calculate for half-mass radius

    return t_cross;
}
/** @brief calculates the relaxation timescale `t_relax` 
 *  @returns `t_relax = N/(8 ln N) * t_cross` with `t_cross = 2*R_hm / v_c` and `v_c = sqrt(G*M_hm/R_hm)` */
double getRelaxationTimescale(const Particles& particles) {
    // // double Rhm = getRhm(particles);
    // // double mass_halfMassRadius = Particle::getTotalMass() / 2;
    // // double v_c = sqrt(G*mass_halfMassRadius/Rhm);

    // // double t_cross = 2*Rhm / v_c;       // ? calculate for half-mass radius ?

    int N = particles.size();           //  calculate with central obj? no
    double t_relax = N/(8*log(N)) * getTCross(particles);

    return t_relax;
}

/* TASK 2 -- tree-code */
Range Range::getChildRange(int quadrant) const {
    double xmid = (xmin + xmax) / 2;
    double ymid = (ymin + ymax) / 2;
    double zmid = (zmin + zmax) / 2;
    switch (quadrant) {
        case 0: return Range(xmin, xmid, ymin, ymid, zmin, zmid);
        case 1: return Range(xmid, xmax, ymin, ymid, zmin, zmid);
        case 2: return Range(xmin, xmid, ymid, ymax, zmin, zmid);
        case 3: return Range(xmid, xmax, ymid, ymax, zmin, zmid);
        case 4: return Range(xmin, xmid, ymin, ymid, zmid, zmax);
        case 5: return Range(xmid, xmax, ymin, ymid, zmid, zmax);
        case 6: return Range(xmin, xmid, ymid, ymax, zmid, zmax);
        case 7: return Range(xmid, xmax, ymid, ymax, zmid, zmax);
        default: throw std::runtime_error("Invalid quadrant");
    }
}

/** simpler version */
double Range::openingAngle(Vector p) const {
    if (contains(p)) return 2*M_PI;
    double dist = (Vector{(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2}-p).dist();
    double rangeSize = std::max({xmax-xmin, ymax-ymin, zmax-zmin})*std::sqrt(3);
    return 2*std::atan(rangeSize/2/dist);
}
    
// double Range::openingAngle(Vector p) const {
//    if (contains(p)) return 2*M_PI;

//    // approximative, but fast solution
//    Vector mid = Vector{(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2};
//    Vectors poss = {
//        Vector{xmin, ymin, zmin},
//        Vector{xmin, ymin, zmax},
//        Vector{xmin, ymax, zmin},
//        Vector{xmin, ymax, zmax},
//        Vector{xmax, ymin, zmin},
//        Vector{xmax, ymin, zmax},
//        Vector{xmax, ymax, zmin},
//        Vector{xmax, ymax, zmax}
//    };

//    double rad = 0;
//    for (const Vector& pos : poss) {
//        rad = std::max(rad, openingAngle(mid-p, pos-p));
//    }

//    return 2*rad; /* note that rad is half of the opening angle */
// }

// ? why not static
Range Range::containingBox(const Particles& particles) {
    double xmin=particles[0].getP().x, xmax=particles[0].getP().x;
    double ymin=particles[0].getP().y, ymax=particles[0].getP().y;
    double zmin=particles[0].getP().z, zmax=particles[0].getP().z;
    for (const auto& p : particles) {
        xmin = std::min(xmin, p.getP().x);
        xmax = std::max(xmax, p.getP().x);
        ymin = std::min(ymin, p.getP().y);
        ymax = std::max(ymax, p.getP().y);
        zmin = std::min(zmin, p.getP().z);
        zmax = std::max(zmax, p.getP().z);
    }
    return Range(xmin, xmax, ymin, ymax, zmin, zmax);
}


/*   
       _____________
      /__6__/__7__ /|
     /__2__/__3__/|7|
    |     |     |3|/|
    |  2  |  3  |/|5/
    |-----|-----|1|/
    |  0  |  1  | /
    |_____|_____|/
*/
Tree* Tree::makeRoot(const Range& range) {
    return new Node(range, nullptr);
}


bool Leaf::add(const Particle& p) {
    ps.push_back(&p);
    if (ps.size() > 8) {
        /* upgrade from Leaf to Node */
        Node* newNode = new Node(range, parent);
        parent->setNewChild(newNode);
        for (const auto& p : ps) {
            newNode->add(*p);
        }

        delete this; /* dangerous */
    }
    return true;
}

/* for seeing if the tree was built correctly */
int treeTest(const Tree* tree) {
    if (tree->getType() == NodeType::LEAF) {
        return dynamic_cast<const Leaf*>(tree)->getNrParticles();
    } else {
        int sum = 0;
        for (int i=0; i<8; i++) {
            sum += treeTest(dynamic_cast<const Node*>(tree)->getChild(i));
        }
        if (sum > 220) {
            std::cout << tree->treeParticle << std::endl;
        }
        return sum;
    }
}
Tree* getTree(const Particles& particles, bool print) {
    Range range = Range::containingBox(particles);
    for (const auto& p : particles) { // comment out for speed
        if (!range.contains(p.getP())) throw std::runtime_error("Particle not in range");
    }

    if (print) paddedPrint(" Creating tree...",'.');

    Tree* root = Tree::makeRoot(range);
    for (const auto& p : particles) {
        root->add(p);
    }
    if (print) paddedPrint(" Added all particles...",'.');
    return root;
}

Treeparticle computeTreeCodeMasses(Tree* tree) {                // could be made const again and put into tree class??
    if (tree->getType() == NodeType::LEAF) {
        Leaf* leaf = dynamic_cast<Leaf*>(tree);
        for (const auto& p : leaf->getParticles()) {
            Treeparticle tp = Treeparticle(p->getMass(), p->getP());
            leaf->treeParticle.add(tp);
        }
    } 
    else {
        Node* node = dynamic_cast<Node*>(tree);
        for (int i=0; i<8; i++) {
            Treeparticle tp = computeTreeCodeMasses(node->getChild(i));
            node->treeParticle.add(tp);
        }
    }
    return tree->treeParticle;
}

Vector getForceRec(const Particle& p, const Tree* tree, double angle, double eps) {
    if (tree->treeParticle.getMass() == 0) return Vector{0,0,0};

    if (tree->getRange().openingAngle(p.getP()) > angle) {
        // Vector force{0,0,0};
        if (tree->getType() == NodeType::LEAF) {
            const Leaf* leaf = dynamic_cast<const Leaf*>(tree);
            return leaf->getForceFrom(p, eps);
            // for (const auto& p2 : leaf->getParticles()) {
            //     if (p == *p2) {
            //         continue;
            //     }
            //     force += p.getForceTo(*p2, eps);
            // }
        } else {
            const Node* node = dynamic_cast<const Node*>(tree);
            return node->getForceFrom(p, angle, eps);
            // for (int i=0; i<8; i++) {
            //     force += getForceRec(p, node->getChild(i), angle, eps);
            // }
        }
        // return force;
    }
    else {
        return p.getForceTo(tree->treeParticle, eps);
    }

}

Vectors getTreeCodeForces(const Particles& particles, const Tree* tree, double softeningMultiplier=1, double angleMultiplier=1, bool print=true) {
    double softening = getSoftening(particles)*softeningMultiplier;
    double angle = M_PI/4 * angleMultiplier; /* 45 degrees */
    Vectors forces;

    if (print) {
        std::cout << "|";
        std::cout.flush();
    }
    int printed = 0; int max = printSize-2;

    for (const auto& p : particles) {

        if (print) {
            if (forces.size()*max > printed*particles.size()) {
                std::cout << "=";
                std::cout.flush();
                printed++;
            }
        }

        forces.push_back(getForceRec(p, tree, angle, softening));
    }

    if (print) {
        std::cout << std::setw(max-printed) << std::setfill('=') << std::right << "|" << std::endl;
        std::cout.flush();
    }
        
    return forces;
}

void compareForces(const Vectors& v1, const Vectors& v2, double soft, double angle, bool print) {
    if (v1.size() != v2.size()) throw std::runtime_error("Vectors have different sizes");
    double sum = 0;
    double sumDiff = 0;
    for (uint i=0; i<v1.size(); i++) {
        Vector diff = v1[i] - v2[i];
        double sumi = (v1[i].dist() + v2[i].dist()) / 2;
        if (print && diff.dist() > sumi * 1e-2 && diff.dist() > 1e-10) {
            std::cout << "Difference in force for particle " << i << ": " << diff << " bigger than 1e-2 of " << v1[i] << ", " << v2[i] << "         ";
            std::cout << "Namely: " << diff.dist() << " > " << sumi * 1e-2 << std::endl;
        }


        sum += sumi;
        sumDiff += diff.dist();
    }

    stringstream ss;
    ss << " Comparing forces for eps=" << soft << " and angle=" << angle << std::endl;
    ss << " Total forces: " << sum << std::endl; 
    ss << " Total difference in forces: " << sumDiff << std::endl;
    boxedPrint(ss);
}


void paddedPrint(string s, char fill, char end, int padLen) {
    std::cout << end << std::setw(padLen-2) << std::setfill(fill) << std::left << s.substr(0,padLen-3) << end << std::endl;
    if (s.length() > (uint) (padLen-3)) paddedPrint(" " + s.substr(padLen-3), fill, end, padLen);
}
void boxedPrint(stringstream& ss) {
    paddedPrint("", '_', ' ');
    string s;
    while (std::getline(ss, s)) paddedPrint(s, ' ');
    paddedPrint("",'_');
}
void boxedPrint(string s) {
    std::stringstream ss(s);
    boxedPrint(ss);
}


string formattedDouble(double d, int precision=3) {
    std::stringstream ss;
    ss << std::fixed << std::setw(precision+7) << std::setprecision(precision) << d;
    return ss.str();
}

template <class TimeitFunctionType, class ... Args>
auto timeit(const string& name, TimeitFunctionType f, bool print, Args... args) {
    if (print) paddedPrint("", '_', ' ');
    if (print) paddedPrint(" FUNCTION " + name + " ", ' ');
    auto t1 = high_resolution_clock::now();
    auto result = f(args...);
    auto t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000.0;
    if (print) paddedPrint(" Time taken: " + formattedDouble(duration) + " milliseconds ", ' ');
    if (print) paddedPrint("", '_');
    return result;
}

static void saveVectorsToFile(const Vectors& vectors, const string& filename) {
    std::ofstream outputFile("bin/"+filename, std::ios::binary);
    for (const Vector& vec : vectors) {
        outputFile.write(reinterpret_cast<const char*>(&vec), sizeof(vec));
    }
}
static void loadVectorsFromFile(Vectors& vectors, const string& filename) {
    std::ifstream inputFile("bin/"+filename, std::ios::binary);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    Vector vec;
    while (inputFile.read(reinterpret_cast<char*>(&vec), sizeof(vec))) {
        vectors.push_back(vec);
    }
}

bool contains(const char* s, char** begin, int count) {
    for (int i=0; i<count; i++) {
        if (std::strcmp(begin[i],s) == 0) return true;
    }
    return false;
}

double containsValue(const char* s, char** begin, int count, double defa) {
    for (int i=0; i<count; i++) {
        if (std::strncmp(begin[i],s, std::strlen(s)) == 0) return std::stod(begin[i]+(std::strlen(s)));
    }
    return defa; // default values
    // // if (std::strcmp(s, "s:") == 0) return 1;
    // // else if (std::strcmp(s, "a:") == 0) return 1;
    // // else if (std::strcmp(s, "t:") == 0) return 10;
    // // else if (std::strcmp(s, "i=") == 0) return 5;
    // // // else if (std::strcmp(s, "dt=") == 0) return .

    // // std::cerr << "No matching argument found for " << s << std::endl;
    // // return 0;
}
string containsValue(const char* s, char** begin, int count, string defa) {
    for (int i=0; i<count; i++) {
        if (std::strncmp(begin[i],s, std::strlen(s)) == 0) return begin[i]+(std::strlen(s));
    }
    return defa; // default values
}

string getFilenameDirect(const string& dataname, int maxParticles, double softeningMult) {
    std::stringstream sd;
    sd << "bin_vectors-directforce-" << dataname << (maxParticles<1e5 ? "_"+to_string(maxParticles) : "") << "-" << softeningMult << ".bin";
    return sd.str();
}
string getFilenameTree(const string& dataname, int maxParticles, double softeningMult, double angleMult) {
    std::stringstream st;
    st << "bin_vectors-treecodeforce-" << dataname << (maxParticles<1e5 ? "_"+to_string(maxParticles) : "") << "-" << softeningMult << "-" << angleMult << ".bin";
    return st.str();
}


Vectors getCDVectors(const Particles& particles, const string& dataname, int maxParticles, double softeningMult, bool loadcd, bool print) {
    Vectors cd;
    string filenameDirect = getFilenameDirect(dataname, maxParticles, softeningMult);
    if (loadcd) loadVectorsFromFile(cd, filenameDirect);
    else {cd = timeit(filenameDirect, getDirectForces, print, particles, softeningMult, print); saveVectorsToFile(cd, filenameDirect);}
    return cd;
}
Vectors getCTVectors(const Particles& particles, const string& dataname, int maxParticles, const Tree* root, double softeningMult, double angleMult, bool loadct, bool print) {
    Vectors ct;
    string filenameTree = getFilenameTree(dataname, maxParticles, softeningMult, angleMult);
    if (loadct) loadVectorsFromFile(ct, filenameTree);
    else {ct = timeit(filenameTree, getTreeCodeForces, print, particles, root, softeningMult, angleMult, print); saveVectorsToFile(ct, filenameTree);}
    return ct;
}

string str(double d) {
    std::stringstream ss;
    ss << d;
    return ss.str();
}



int main(int argc, char* argv[]) {
    /*
    - no central particle in dataset data
    - central mass is 1 in datasets data0, data1
    - total particle mass is 10 in datasets data0, data1
    */

    // string dataname = "dataEarth";
    // int maxParticles = 1e1;


    std::cout << "Start" << std::endl << std::endl;

    string dataname      = containsValue("d:", argv+1, argc-1, "data");
    string dataset       = dataname + ".txt";
    int maxParticles     = containsValue("n:", argv+1, argc-1, 1e5);
    Particles particles  = readIn(dataset, maxParticles);
    
    double softeningMult = containsValue("s:", argv+1, argc-1, 1);
    double angleMult     = containsValue("a:", argv+1, argc-1, 1);
    double timeMult      = containsValue("t:", argv+1, argc-1, 10);
    double solverImgs    = containsValue("i:", argv+1, argc-1, 5);
    // double dt            = containsValue("dt=", argv+1, argc-1);
    bool loadcd          = contains("loadcd", argv+1, argc-1);
    bool loadct          = contains("loadct", argv+1, argc-1);
    bool printComp       = contains("printComp", argv+1, argc-1);
    bool printFunc       = contains("printFunc", argv+1, argc-1);

    boxedPrint(" SORTING PARTICLES ");
    Particle::sortParticles(particles);

    //* TESTING DATA INTEGRITY
    if (contains("test", argv+1, argc-1)) {
        stringstream ss;
        ss << " Testing data for " << dataset << " with " << particles.size() << " particles " << std::endl;
        ss << " Total mass: " << Particle::getTotalMass() << std::endl;
        boxedPrint(ss);
    }

    //* TASK 1.1
    if (contains("task1.1", argv+1, argc-1)) {
        paddedPrint("", '_', ' ');
        paddedPrint(" TASK 1.1 ", ' ');
        paddedPrint(" Computing shells ", '.');
        double maxDist = Particle::maxDistance(particles);
        doubles shellBoundaries0 = (!(contains("linspace", argv+1, argc-1))) ? logspace0(0.1, maxDist, 30, 1.2) : linspace(0, maxDist, 31);
        doubles shellMidpoints = addWithOffset(shellBoundaries0, 1)/2;
        doubles shellComp = computeShells(particles, shellBoundaries0);
        doubles rhos = computeRhos(shellMidpoints);
        doubles poissonErrors = computePoissonErrors(shellComp, rhos);
        paddedPrint(" " + to_string(shellComp.size()) + " shells computed ", ' ');
        paddedPrint("", '_');

        string filename = "plot_task1_1";
        clearMyPlot(filename);
        writePlotToMyPlot(shellBoundaries0, filename);
        writePlotToMyPlot(shellMidpoints, filename);
        writePlotToMyPlot(shellComp, filename);
        writePlotToMyPlot(rhos, filename);
        writePlotToMyPlot(poissonErrors, filename);
        runPython(filename);
    }
    
    //* TASK 1.2
    if (contains("task1.2", argv+1, argc-1)) {
        // part 1 - softening (& direct forces, seen in task 2)
        boxedPrint(" TASK 1.2 ");
        double softening = timeit("getSoftening", getSoftening, printFunc, particles); /* order of the mean interparticle separation */
        paddedPrint(" SOFTENING:  " + formattedDouble(softening, 6) + " ", '_');


        // part 2 - analytical force stuff
        double prec = .1;
        doubles shellBoundaries, analytical, shellMidpoints, direct, angleOff;
        getAnalyticalForce(particles, prec, shellBoundaries, analytical, shellMidpoints);

        Vectors cd = getCDVectors(particles, dataname, maxParticles, softeningMult, loadcd, printFunc);
        Vectors scatter;

        
        getComparisonDirectAnalyticalForce(particles, cd, scatter, shellBoundaries, direct, angleOff);

        string filename = "plot_task1_2";
        clearMyPlot(filename);
        writePlotToMyPlot(shellMidpoints, filename);
        writePlotToMyPlot(analytical, filename);
        writePlotToMyPlot(direct, filename);
        writePlotToMyPlot(angleOff, filename);
        writePlotToMyPlot(scatter, filename);
        runPython(filename);


        // part 3 - relaxation time
        double relaxationTime = getRelaxationTimescale(particles);
        paddedPrint("", '_', ' ');
        paddedPrint(" RELAXATION TIME: " + formattedDouble(relaxationTime) + " ", ' ');
        string info = " For an increased softening, one can expect a longer relaxation time, "
        "as the particles are less affected by each other (and the relaxation time gives an "
        "estimate of how long it takes for the particles to change their velocities due to "
        "gravitational interactions with other particles). ";
        paddedPrint(info, ' ');
        paddedPrint("", '_');
    }




    //* TASK 2 -- tree-code
    if (contains("task2", argv+1, argc-1)) {
        paddedPrint("", '_', ' ');
        paddedPrint(" TASK 2 ", ' ');
        Tree* root = getTree(particles, printFunc);             // root needs to be deleted later
        computeTreeCodeMasses(root);
        if (printFunc) paddedPrint(" Computed all treeparticle masses...", '.');
        if (contains("test", argv+1, argc-1)) treeTest(root);   // some test to see if the tree was built correctly
        paddedPrint("", '_');

        // // Vectors cd, ct;                                         // load or compute direct and tree-code forces
        // // string filenameDirect = getFilenameDirect(dataname, maxParticles, softeningMult);
        // // string filenameTree = getFilenameTree(dataname, maxParticles, softeningMult, angleMult);
        // // if (loadcd) loadVectorsFromFile(cd, filenameDirect);
        // // else {cd = timeit(filenameDirect, getDirectForces, particles, softeningMult, true); saveVectorsToFile(cd, filenameDirect);}

        // // if (loadct) loadVectorsFromFile(ct, filenameTree);
        // // else {ct = timeit(filenameTree, getTreeCodeForces, particles, root, softeningMult, angleMult, true); saveVectorsToFile(ct, filenameTree);}
        Vectors cd = getCDVectors(particles, dataname, maxParticles, softeningMult, loadcd, printFunc);
        Vectors ct = getCTVectors(particles, dataname, maxParticles, root, softeningMult, angleMult, loadct, printFunc);

        compareForces(cd, ct, softeningMult, angleMult, printComp); // compare direct and tree-code forces that were computed or loaded


        if (contains("compareMulti", argv+1, argc-1)) {         // compare for different softening and angle values
            boxedPrint(" COMPARING MULTIPLE SOFTENING AND ANGLE VALUES ");
            for (double soft : doubles{0,.1,1,10}) {
                Vectors cd = getCDVectors(particles, dataname, maxParticles, soft, loadcd, printFunc);
                for (double angle : doubles{.1,.3,1,3}) {

                    Vectors ct = getCTVectors(particles, dataname, maxParticles, root, soft, angle, loadct, printFunc);
                    compareForces(cd, ct, soft, angle, printComp);

                }
            }
        }

        delete root;
    }

    if (contains("gravityTreeSolver", argv+1, argc-1)) {
        boxedPrint(" GRAVITY TREE SOLVER ");
        // t_cross = .1167 for data0, .1182 for data1, .0001 for data, 0 for dataEarth

        double timescale = getTCross(particles) * timeMult;
        if (dataname == "dataEarth") timescale = 1./92.299; // one earth year
        double dt = timescale / solverImgs;

        string filename = "plot_task2";
        clearMyPlot(filename);


        Particle::sortParticles(particles);
        Vectors forces = getCDVectors(particles, dataname+"-tLFstart"+str(0)+"-"+str(dt), maxParticles, softeningMult, loadcd, printFunc);
        for (double t=0; t<timescale; t+=dt) {

            writePlotToMyPlot(map<Particle,Vector>(particles, [](const Particle& p) {return p.getP();}), filename);

            for (Particle& p : particles) p.setApplyLeapFrog1(forces[p.getId()], dt);

            // Tree* root = getTree(particles, printFunc);
            // computeTreeCodeMasses(root);
            // Vectors forces = getCTVectors(particles, dataname+"-t"+str(t), maxParticles, root, softeningMult, angleMult, loadct, printFunc);
            Particle::sortParticles(particles);
            forces = getCDVectors(particles, dataname+"-tLF"+str(t)+"-"+str(dt), maxParticles, softeningMult, loadcd, printFunc);
            for (Particle& p : particles) {
                p.setApplyLeapFrog2(forces[p.getId()], dt);
            }
        }
        // for (double t=0; t<timescale; t+=dt) {

        //     Particle::sortParticles(particles);
        //     writePlotToMyPlot(map<Particle,Vector>(particles, [](const Particle& p) {return p.getP();}), filename);
        //     // Tree* root = getTree(particles, printFunc);
        //     // computeTreeCodeMasses(root);
        //     // Vectors forces = getCTVectors(particles, dataname+"-t"+str(t), maxParticles, root, softeningMult, angleMult, loadct, printFunc);
        //     Vectors forces = getCDVectors(particles, dataname+"-t"+str(t), maxParticles, softeningMult, loadcd, printFunc);
        //     for (Particle& p : particles) {
        //         p.setApplyForce(forces[p.getId()], dt);
        //     }
        // }
        writePlotToMyPlot(map<Particle,Vector>(particles, [](const Particle& p) {return p.getP();}), filename);

        runPython(filename);
    }













    std::cout << std::endl << "Done!" << std::endl;
    return 0;
}

int clearMyPlot(const string& filename) {
    std::ofstream clearFile(filename+".txt");
    if (!clearFile.is_open()) throw std::runtime_error("Could not open file: " + filename+".txt");
    clearFile << "";

    return 0;
}

int writePlotToMyPlot(const Vectors& vs, const string& filename) {
    std::ofstream appendFile(filename+".txt", std::ios_base::app);
    if (!appendFile.is_open()) throw std::runtime_error("Could not open file: " + filename+".txt");
    for (uint i=0; i<vs.size(); i++) {
        if (i==0) appendFile << "[";
        else appendFile << ", ";
        appendFile << "[" << vs[i].x << ", " << vs[i].y << ", " << vs[i].z << "]";
        if (i==vs.size()-1) appendFile << "]";
    }
    appendFile << std::endl;
    appendFile.close();

    return 0;
}

int writePlotToMyPlot(const doubles& xs, const string& filename) {
    std::ofstream appendFile(filename+".txt", std::ios_base::app);
    if (!appendFile.is_open()) throw std::runtime_error("Could not open file: " + filename+".txt");
    for (double x : xs) appendFile << x << " ";
    appendFile << std::endl;
    appendFile.close();

    return 0;
}

/** @brief returns evenly spaced points between start and end, inclusive*/
doubles linspace(double start, double end, int numPoints, bool inclusive) {
    doubles xs(numPoints);
    int inclusiveInt = inclusive ? 1 : 0;
    double step = (end-start)/(numPoints-inclusiveInt);
    for (int i=0; i<numPoints; i++) xs[i] = start + i*step;
    return xs;
}

/** @brief in contrast to np.logspace, start and end are actually start and end */
doubles logspace(double start, double end, int numPoints, double base) {
    double logStart = log(start) / log(base);
    double logEnd = log(end) / log(base);
    doubles xs(numPoints);
    return map<double, double>(linspace(logStart, logEnd, numPoints), [base](double x) {return pow(base, x);});
}
/** @brief returns numIntervals many intervals, that is, numIntervals+1 many points ([0, logspace]) */
doubles logspace0(double start, double end, int numIntervals, double base) {
    doubles xs = logspace(start, end, numIntervals, base);
    xs.insert(xs.begin(), 0);
    return xs;
}


template <typename T1, typename T2>
std::vector<T2> map(const std::vector<T1>& xs, std::function<T2(T1)> f) {
    std::vector<T2> ys;
    std::transform(xs.begin(), xs.end(), std::back_inserter(ys), f);
    return ys;
}

void runPython(string filename) {
    std::string command = "python " + filename + ".py";
    int result = system(command.c_str());
    if (result != 0) throw std::runtime_error("Python script failed");
}
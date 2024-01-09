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


Vector getForceRec(const Particle& p, const Tree* tree, double angle, double eps);

int clearMyPlot(const string& filename = "myplot");
int writePlotToMyPlot(const doubles& xs, const string& filename = "myplot");
doubles linspace(double start, double end, int numPoints = 60, bool inclusive = true);
doubles logspace(double start, double end, int numPoints = 60, double base = 2);
doubles logspace0(double start, double end, int numIntervals = 60, double base = 2);
template <typename T1, typename T2>
std::vector<T2> map(const std::vector<T1>& xs, std::function<T2(T1)> f);
void runPython(string filename);


const double G = 1;

/* MY TYPES */

/** @brief some operators for simple operations on doubles */
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
        return sqrt(x*x + y*y + z*z);
    }
    Vector operator+(const Vector& other) const {
        return Vector{x+other.x, y+other.y, z+other.z};
    }
    Vector operator-(const Vector& other) const {
        return Vector{x-other.x, y-other.y, z-other.z};
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

    friend std::ostream& operator<<(std::ostream& os, const Vector& v) {
        os << "V (" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }
};


class Particle {
protected:
    static double totalMass;
    static int idCount;
    static Particles sortedParticles; // ! reset when particles moved

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
        Vector r12 = other.p - p;
        double r = r12.dist();
        double force = G * mass * other.mass / (r*r + eps*eps);
        if (r!=0) return (force/r) * r12;
        else { std::cerr << "r=0 for particles " << *this << " and " << other << std::endl; return Vector{0,0,0}; }
    }

    double getDist() const {
        return p.dist();
    }
    double getMass() const {
        return mass;
    }
    Point getP() const {
        return p;
    }

    int getId() const {
        return id;
    }

    static double getTotalMass() {
        return totalMass;
    }

    static double maxDistance(const Particles& particles) {
        double max = 0;
        for (const auto& p : particles) {
            max = std::max(max, p.getDist());
        }
        return max;
    }

};
double Particle::totalMass = 0;
int Particle::idCount = 0;





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
    double openingAngle(Vector x, Vector y) const {
        return std::acos(x*y / x.dist() / y.dist());
    }

public:
    Range(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), zmin(zmin), zmax(zmax) {}
    bool contains(const Vector& p) const { return xmin <= p.x && p.x <= xmax && ymin <= p.y && p.y <= ymax && zmin <= p.z && p.z <= zmax; }
    int getQuadrant(const Vector& p) const { return (p.x > (xmin+xmax)/2) + 2*(p.y > (ymin+ymax)/2) + 4*(p.z > (zmin+zmax)/2); }
    int getQuadrant(const Range& r) const { return (r.xmin != xmin) + 2*(r.ymin != ymin) + 4*(r.zmin != zmin); }
    Range getChildRange(int quadrant) const;
    double openingAngle(Vector p) const;

    static Range containingBox(const Particles& particles);
    friend std::ostream& operator<<(std::ostream& os, const Range& r);
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

    Particles sorted(particles);

    std::sort(sorted.begin(), sorted.end(), [](const Particle& p1, const Particle& p2) {
        return p1.getDist() < p2.getDist();
    });

    // TODO could store the sorted sequence!

    double totalMass = Particle::getTotalMass();


    double mass = 0;
    int i=0;
    while (mass < totalMass/2) {
        // std::cout << i << std::endl;
        mass += sorted[i++].getMass();
    }

    /* Particles 0..i-1 have mass < totalMass/2, Particles 0..i have mass >= totalMass/2 */
    return i-1;
}

/** @brief calculates the half-mass radius `R_hm` */
double getRhm(const Particles& particles) {
    int half = getHmNr(particles);

    // TODO could also store after it is calculated once!
    doubles dists;
    for (const auto& p : particles) {
        dists.push_back(p.getDist());
    }
    nth_element(dists.begin(), dists.begin()+half-1, dists.end());
    return dists[half];
}

/** @brief calculates the softening `d = (V/N)^(1/3)` */
double getSoftening(const Particles& particles) {
    double Rhm = getRhm(particles);
    double V = 4./3 * M_PI * Rhm*Rhm*Rhm;
    int N = getHmNr(particles);

    return pow(V/N, 1./3);
}

/** @brief calculates the direct forces in `O(n^2)`*/
Vectors getDirectForces(const Particles& particles, double softeningMultiplier=1, bool print=true) {
    double softening = getSoftening(particles)*softeningMultiplier;

    if (print) {
        std::cout << "|";
        std::cout.flush();
    }
    int printed = 0; int max = 58;
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

/* applies Newton II theorem for spherical potentials */
int getAnalyticalForce(const Particles& particles) {
    /* \vec F(r) = -GM(r)/r^3 * \vec r with M(r)=mass within radius r */
    return 0;
}

/** @brief calculates the relaxation timescale `t_relax` 
 *  @returns `t_relax = N/(8 ln N) * t_cross` with `t_cross = 2*R_hm / v_c` and `v_c = sqrt(G*M_hm/R_hm)`
 */
double getRelaxationTimescale(const Particles& particles) {
    double Rhm = getRhm(particles);
    double mass_halfMassRadius = Particle::getTotalMass() / 2;
    double v_c = sqrt(G*mass_halfMassRadius/Rhm);

    int N = particles.size();           // calculate with central obj.?
    double t_cross = 2*Rhm / v_c;       // calculate for half-mass radius ?

    double t_relax = N/(8*log(N)) * t_cross;

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

// ? why not friend and Range::
std::ostream& operator<<(std::ostream& os, const Range& r) {
    os << "[(" << r.xmin << "," << r.xmax << ") x (" << r.ymin << "," << r.ymax << ") x (" << r.zmin << "," << r.zmax << ")]";
    return os;
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
            // std::cout << "Treenode with " << sum << " particles: " << *tree << std::endl;
            std::cout << tree->treeParticle << std::endl;
        }
        return sum;
    }
}
Tree* getTree(const Particles& particles) {
    Range range = Range::containingBox(particles);
    for (const auto& p : particles) { // comment out for speed
        if (!range.contains(p.getP())) throw std::runtime_error("Particle not in range");
    }

    std::cout << "Creating tree..." << std::endl;

    Tree* root = Tree::makeRoot(range);
    for (const auto& p : particles) {
        root->add(p);
    }
    std::cout << "Added all particles..." << std::endl;
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
    int printed = 0; int max = 58;

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

void compareForces(const Vectors& v1, const Vectors& v2) {
    if (v1.size() != v2.size()) throw std::runtime_error("Vectors have different sizes");
    double sum = 0;
    double sumDiff = 0;
    for (uint i=0; i<v1.size(); i++) {
        Vector diff = v1[i] - v2[i];
        double sumi = (v1[i].dist() + v2[i].dist()) / 2;
        if (diff.dist() > sumi * 1e-2 && diff.dist() > 1e-10) {
            std::cout << "Difference in force for particle " << i << ": " << diff << " bigger than 1e-2 of " << v1[i] << ", " << v2[i] << "         ";
            std::cout << "Namely: " << diff.dist() << " > " << sumi * 1e-2 << std::endl;
        }


        sum += sumi;
        sumDiff += diff.dist();
    }
    std::cout << "Total forces: " << sum << std::endl; 
    std::cout << "Total difference in forces: " << sumDiff << std::endl;
}


void paddedPrint(string s = "", char fill = '-', char end = '|', int padLen = 60) {
    std::cout << end << std::setw(padLen-2) << std::setfill(fill) << std::left << s.substr(0,padLen-3) << end << std::endl;
    if (s.length() > (uint) (padLen-3)) paddedPrint(" " + s.substr(padLen-3), fill, end, padLen);
}
string formattedDouble(double d, int precision=3) {
    std::stringstream ss;
    ss << std::fixed << std::setw(precision+7) << std::setprecision(precision) << d;
    return ss.str();
}

template <class TimeitFunctionType, class ... Args>
auto timeit(const string& name, TimeitFunctionType f, Args... args) {
    paddedPrint("", '_', ' ');
    paddedPrint(" FUNCTION " + name + " ", ' ');
    auto t1 = high_resolution_clock::now();
    auto result = f(args...);
    auto t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000.0;
    paddedPrint(" Time taken: " + formattedDouble(duration) + " milliseconds ", ' ');
    paddedPrint("", '_');
    return result;
}

static void saveVectorsToFile(const Vectors& vectors, const string& filename) {
    std::ofstream outputFile(filename, std::ios::binary);
    for (const Vector& vec : vectors) {
        outputFile.write(reinterpret_cast<const char*>(&vec), sizeof(vec));
    }
}
static void loadVectorsFromFile(Vectors& vectors, const string& filename) {
    std::ifstream inputFile(filename, std::ios::binary);
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

double containsValue(const char* s, char** begin, int count) {
    for (int i=0; i<count; i++) {
        if (std::strcmp(begin[i],s) == 0) return std::stod(begin[i]+(std::strlen(s)));
    }
    // default values
    if (std::strcmp(s, "s=") == 0) return 1;
    else if (std::strcmp(s, "a=") == 0) return 1;

    std::cerr << "No matching argument found" << std::endl;
    return 0;
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



// std::stringstream sd; sd << "bin_vectors-directforce-" << dataname << (maxParticles<1e5 ? "_"+to_string(maxParticles) : "") << "-" << softeningMult << ".bin"; string filenameDirect = sd.str();
//     std::stringstream st; st << "bin_vectors-treecodeforce-" << dataname << (maxParticles<1e5 ? "_"+to_string(maxParticles) : "") << "-" << softeningMult << "-" << angleMult << ".bin"; string filenameTree = st.str();





int main(int argc, char* argv[]) {
    /*
    - no central particle in dataset data
    - central mass is 1 in datasets data0, data1
    - total particle mass is 10 in datasets data0, data1
    */

    string dataname = "data";
    string dataset = dataname + ".txt";
    int maxParticles = 1e4;
    Particles particles = readIn(dataset, maxParticles);


    std::cout << "Start" << std::endl << std::endl;
    
    double softeningMult = containsValue("s=", argv+1, argc-1);
    double angleMult = containsValue("a=", argv+1, argc-1);
    bool loadcd = contains("loadcd", argv+1, argc-1);
    bool loadct = contains("loadct", argv+1, argc-1);

    //* TESTING DATA INTEGRITY
    if (contains("test", argv+1, argc-1)) {
        paddedPrint("", '_', ' ');
        paddedPrint(" TESTING DATA INTEGRITY ", ' ');
        paddedPrint(" Total mass: " + formattedDouble(Particle::getTotalMass()) + " ", ' ');
        paddedPrint("",'_');
    }

    //* TASK 1.1
    if (contains("task1.1", argv+1, argc-1)) {
        paddedPrint("", '_', ' ');
        paddedPrint(" TASK 1.1 ", ' ');
        paddedPrint(" Computing shells ", ' ');
        double maxDist = Particle::maxDistance(particles);
        doubles shellBoundaries0 = true ? logspace0(0.1, maxDist, 30, 1.2) : linspace(0, maxDist, 31);
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
        paddedPrint("", '_', ' ');
        paddedPrint(" TASK 1.2 ", ' ');
        paddedPrint("", '_');
        double softening = timeit("getSoftening", getSoftening, particles); /* order of the mean interparticle separation */
        paddedPrint(" SOFTENING:  " + formattedDouble(softening, 6) + " ", '_');

        // part 2 - analytical force stuff
        // TODO IMPLEMENT ANALYTICAL FORCE STUFF




    }



    // double relaxationTime = getRelaxationTimescale(particles);
    // paddedPrint("", '_', ' ');
    // paddedPrint(" RELAXATION TIME: " + formattedDouble(relaxationTime) + " ", ' ');
    // string info = " For an increased softening, one can expect a longer relaxation time, "
    // "as the particles are less affected by each other (and the relaxation time gives an "
    // "estimate of how long it takes for the particles to change their velocities due to "
    // "gravitational interactions with other particles). ";
    // paddedPrint(info, ' ');
    // paddedPrint("", '_');

    //* TASK 2 -- tree-code
    if (contains("task2", argv+1, argc-1)) {
        Tree* root = getTree(particles);                        // root needs to be deleted later
        computeTreeCodeMasses(root);
        if (contains("test", argv+1, argc-1)) treeTest(root);   // some test to see if the tree was built correctly

        Vectors cd, ct;                                         // load or compute direct and tree-code forces
        string filenameDirect = getFilenameDirect(dataname, maxParticles, softeningMult);
        string filenameTree = getFilenameTree(dataname, maxParticles, softeningMult, angleMult);
        if (loadcd) loadVectorsFromFile(cd, filenameDirect);
        else {cd = timeit(filenameDirect, getDirectForces, particles, softeningMult, true); saveVectorsToFile(cd, filenameDirect);}

        if (loadct) loadVectorsFromFile(ct, filenameTree);
        else {ct = timeit(filenameTree, getTreeCodeForces, particles, root, softeningMult, angleMult, true); saveVectorsToFile(ct, filenameTree);}


        std::cout << "Comparing forces for eps=" << softeningMult << ",ang=" << angleMult << " " << std::endl; compareForces(cd, ct);


        if (contains("compareMulti", argv+1, argc-1)) {
            // TODO compare for soft=0, .1, 1, 10 and angle=.1, .3, 1, 3
        }

        delete root;
    }













    std::cout << std::endl << "Done!" << std::endl;
    return 0;
}

int clearMyPlot(const string& filename = "myplot") {
    std::ofstream clearFile(filename+".txt");
    if (!clearFile.is_open()) throw std::runtime_error("Could not open file: " + filename+".txt");
    clearFile << "";

    return 0;
}

int writePlotToMyPlot(const doubles& xs, const string& filename = "myplot") {
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
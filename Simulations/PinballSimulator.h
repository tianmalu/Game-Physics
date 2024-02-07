#ifndef PINBALLSIMULATOR_h
#define PINBALLSIMULATOR_h

#include "Simulator.h"
#include "collisionDetect.h"

class PinballSimulator: public Simulator{
private:
    /* data */
    //Attribues
    float damping;

    // UI Attributes
	Vec3 m_externalForce;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;
    // Extra UI Attributes
    //TODO: To detect the state of the mouse

    struct Pinball
    {
        Vec3 position;
        Vec3 linearVelocity;
        Vec3 pinballSize;
        float mass;

        Quat orientation;
        //Vec3 angularMomentum;
        //Vec3 forceLocation;
        Vec3 angularVelocity;
        //Vec3 externalForce;

        Mat4 worldMatrix;
    };

    struct Spring
    {
        float stiffness;
        float length;
    };

    
    struct RigidWall
    {
        Vec3 position;
        Quat orientation;
        Vec3 size;
        Mat4 worldMatrix;
        float direction;
    };

    std::vector<RigidWall> rigidwallVector;
    std::vector<Spring> springVector;
    std::vector<Pinball> pinballVector;

    void explicitEuler(float timeStep);
    void testSimulator(float timeStep);
    
public:
    //Constructor
    PinballSimulator(/* args */);

    //Functions
    void setDampingFactor(float damping);
    void setWorldMatrix(Mat4 worldMatrix);
    void setPositionOfPinball(Vec3 position);
    void setVelocityOfPinball(Vec3 velocity);
    void setDirectionOfPinball(float direction);

    Vec3 getPositionOfPinball();
	Vec3 getVelocityOfPinball();
    Mat4 getWorldMatrixOfPinball();
    Mat4 getWorldMatrixOfRigidwall(int index);
    
    void addPinball();
	void addSpring(float initialLength, float stiffness);
    void addRigidWall(Vec3 position, Vec3 orientation, Vec3 size);
    void addRigidWall2(Vec3 position, Vec3 orientation, Vec3 size);
    void addRigidWall_rotation(Vec3 position, Quat orientation, Vec3 size, float angle);
    void addPlayboard();
    
    void applyExternalForce(Vec3 force);
    void springForce(float timeStep);
    void forceEular(float timeStep);
    std::vector<std::vector<float>> calculateOuterProduct(const Vec3& w);

    void forceEular2(float timeStep);

    std::vector<std::vector<float>> calculateOuterProduct2(const Vec3& w);


    // UI Functions
	const char * getTestCasesStr();
	void initUI(DrawingUtilitiesClass * DUC);
	void reset();
	void drawFrame(ID3D11DeviceContext* pd3dImmediateContext);
	void notifyCaseChanged(int testCase);
    void generateRandomWalls(int numWalls, float wallHeight, float minSize, float maxSize);

	void externalForcesCalculations(float timeElapsed);
	void simulateTimestep(float timeStep);
	void onClick(int x, int y);
	void onMouse(int x, int y);
    void longPressOnMouse(float deltaTime);

};

#endif
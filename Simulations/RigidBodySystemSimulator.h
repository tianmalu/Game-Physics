#ifndef RIGIDBODYSYSTEMSIMULATOR_h
#define RIGIDBODYSYSTEMSIMULATOR_h
#include "Simulator.h"
#include "collisionDetect.h"
//add your header for your rigid body system, for e.g.,
//#include "rigidBodySystem.h" 

#define TESTCASEUSEDTORUNTEST 5

class RigidBodySystemSimulator:public Simulator{

private:
	// Attributes
	// add your RigidBodySystem data members, for e.g.,
	// RigidBodySystem * m_pRigidBodySystem; 
	struct Rigidbox {
		Vec3 position;
		Vec3 size;
		int mass;
		Quat orientation;
		Vec3 linearVelocity;
		Vec3 angularVelocity;
		Vec3 m_externalForce;
		Vec3 forceLocation;
		Vec3 angularMomentum;
		Mat4 inertiaTensor;
		Mat4 inertiaTensorInv;
		Mat4 worldMatrix;
		Vec3 torque;
	};
	std::vector<Rigidbox> rigidboxVector;
	float m_fGravity = 0.0f;
	float m_fAddVelocity_x = 0.0f;
	float m_fAddVelocity_y = 0.0f;
	float m_fAddVelocity_z= 0.0f;
	float inputScale = 0.000005f;
	// UI Attributes
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;

public:
	// Construtors
	RigidBodySystemSimulator();
	
	// Functions
	const char * getTestCasesStr();
	void initUI(DrawingUtilitiesClass * DUC);
	void reset();
	void drawFrame(ID3D11DeviceContext* pd3dImmediateContext);
	void notifyCaseChanged(int testCase);
	void externalForcesCalculations(float timeElapsed);
	void simulateTimestep(float timeStep);
	void onClick(int x, int y);
	void onMouse(int x, int y);

	// ExtraFunctions
	int getNumberOfRigidBodies();
	Vec3 getPositionOfRigidBody(int i);
	Vec3 getLinearVelocityOfRigidBody(int i);
	Vec3 getAngularVelocityOfRigidBody(int i);
	void applyForceOnBody(int i, Vec3 loc, Vec3 force);
	void addRigidBody(Vec3 position, Vec3 size, int mass);
	void setOrientationOf(int i,Quat orientation);
	void setVelocityOf(int i, Vec3 velocity);

	//AdditionFunctions
	void force_eular(float timeStep);
	Vec3 CalculateTorque(const Vec3& force, const Vec3& forceLocation, const Vec3& position);
	void IntegrateOrientation(Rigidbox& rb, float timeStep);
	void IntegrateAngularMomentum(Rigidbox& rb, const Vec3& torque, float timeStep);
	void UpdateInverseInertiaTensor(Rigidbox& rb);
	void UpdateAngularVelocity(Rigidbox& rb);
	Quat QuaternionMultiply(const Quat& q1, const Quat& q2);
	void force_Impulse(CollisionInfo info, Rigidbox& rb1, Rigidbox& rb2);
};
#endif
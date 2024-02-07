#include "RigidBodySystemSimulator.h"

// Construtors
RigidBodySystemSimulator::RigidBodySystemSimulator()
{
    // UI Attributes
}

// Functions
const char* RigidBodySystemSimulator::getTestCasesStr()
{
    return " Demo 1: One-step test,\
				Demo 2: Single-body simulation,\
				Demo 3: Two-rigid-body collision scene ,\
				Demo 4: Complex simulation";

}

void RigidBodySystemSimulator::initUI(DrawingUtilitiesClass* DUC)
{
	this->DUC = DUC;
	switch (m_iTestCase)
	{
	case 0:
		break;
	case 1:
		TwAddVarRW(DUC->g_pTweakBar, "Gravity", TW_TYPE_FLOAT, &m_fGravity, "min=0 max=10 step=0.001");
		break;
	case 2:
		TwAddVarRW(DUC->g_pTweakBar, "Gravity", TW_TYPE_FLOAT, &m_fGravity, "min=0 max=10 step=0.001");
		break;
	case 3:
		TwAddVarRW(DUC->g_pTweakBar, "Gravity", TW_TYPE_FLOAT, &m_fGravity, "min=0 max=10 step=0.001");
		TwAddVarRW(DUC->g_pTweakBar, "AddVelocity_X", TW_TYPE_FLOAT, &m_fAddVelocity_x, "min=-1 max=1 step=0.01");
		TwAddVarRW(DUC->g_pTweakBar, "AddVelocity_Y", TW_TYPE_FLOAT, &m_fAddVelocity_y, "min=-1 max=1 step=0.01");
		TwAddVarRW(DUC->g_pTweakBar, "AddVelocity_Z", TW_TYPE_FLOAT, &m_fAddVelocity_z, "min=-1 max=1 step=0.01");
		break;
	default:break;
	}

}

void RigidBodySystemSimulator::reset()
{
	rigidboxVector.clear();
    m_mouse.x = m_mouse.y = 0;
    m_trackmouse.x = m_trackmouse.y = 0;
    m_oldtrackmouse.x = m_oldtrackmouse.y = 0;

}

void RigidBodySystemSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext)
{
	DUC->setUpLighting(Vec3(0, 0, 0), 0.4 * Vec3(1, 1, 1), 2000.0, Vec3(0.5, 0.5, 0.5));
	for (Rigidbox& rb : rigidboxVector) {
		Mat4 position_M4 = Mat4(0.0);
		position_M4.initTranslation(rb.position.x, rb.position.y, rb.position.z);
		Mat4 size_M4 = Mat4(0.0);
		size_M4.initScaling(rb.size.x, rb.size.y, rb.size.z);
		rb.worldMatrix = size_M4 * rb.orientation.getRotMat() * position_M4;
		DUC->drawRigidBody(rb.worldMatrix);
	}
}

void RigidBodySystemSimulator::notifyCaseChanged(int testCase)
{
	m_iTestCase = testCase;
	reset();
	switch (m_iTestCase)
	{
	case 0:
	{
		cout << "Demo 1 : One - step test";
		addRigidBody(Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.6f, 0.5f), 2);
		setOrientationOf(0, Quat(Vec3(0.0f, 0.0f, 1.0f), (float)(M_PI) * 0.5f));
		applyForceOnBody(0, Vec3(0.3f, 0.5f, 0.25f), Vec3(1.0f, 1.0f, 0.0f));
		Vec3 point(0.3, 0.5, 0.25);
		force_eular(2);
		//print zone
		cout << "\nLinearVelocity: " << rigidboxVector[0].linearVelocity << "\n";
		cout << "AngularVelocity: " << rigidboxVector[0].angularVelocity << "\n";
		cout << "point-Velocity: " << rigidboxVector[0].linearVelocity + cross(rigidboxVector[0].angularVelocity,point) << "\n";
		inputScale = 0.000005f;
		break;
	}
	case 1:
		cout << "Demo 2: Single - body simulation\n";
		addRigidBody(Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.6f, 0.5f), 2);
		setOrientationOf(0, Quat(Vec3(0.0f, 0.0f, 1.0f), (float)(M_PI) * 0.5f));
		applyForceOnBody(0, Vec3(0.3f, 0.5f, 0.25f), Vec3(1.0f, 1.0f, 0.0f));
		inputScale = 0.000005f;
		break;
	case 2:
		cout << "Demo 3: Two - rigid - body collision scene.\n";
		addRigidBody(Vec3(2.0f, 0.5f, 0.0f), Vec3(0.7f, 1.0f, 0.5f), 2);
		addRigidBody(Vec3(-2.0f, 0.5f, 0.0f), Vec3(0.7f, 1.0f, 0.5f), 2);
		setVelocityOf(0, Vec3(-0.2f, 0.0f, 0.0f));
		setVelocityOf(1, Vec3(0.2f, -0.0f, 0.0f));
		applyForceOnBody(0, getPositionOfRigidBody(0)-Vec3(0.3,0.1,0.0), getLinearVelocityOfRigidBody(0));
		applyForceOnBody(1, getPositionOfRigidBody(1) - Vec3(0.0, 0.2, 0.0), getLinearVelocityOfRigidBody(1));
		setOrientationOf(0, Quat(Vec3(1.0f, 0.0f, 1.0f), (float)(M_PI) * 0.35f));
		inputScale = 0.005f;
		break;
	case 3:
		cout << "Demo 4: Complex simulation.\n";
		addRigidBody(Vec3(2.0f, 0.5f, 0.0f), Vec3(0.7f, 1.0f, 0.5f), 2);
		addRigidBody(Vec3(-2.0f, 0.5f, 0.0f), Vec3(0.7f, 1.0f, 0.5f), 2);
		addRigidBody(Vec3(2.0f, 3.0f, 2.0f), Vec3(0.7f, 1.0f, 0.5f), 2);
		addRigidBody(Vec3(-2.0f, 3.0f, -2.0f), Vec3(0.7f, 1.0f, 0.5f), 2);
		setVelocityOf(0, Vec3(-0.2f, 0.0f, 0.0f));
		setVelocityOf(1, Vec3(0.2f, -0.0f, 0.0f));
		setVelocityOf(2, Vec3(-0.3f, 0.2f, 0.1f));
		setVelocityOf(3, Vec3(0.1f, -0.2f, 0.6f));
		applyForceOnBody(0, getPositionOfRigidBody(0) - Vec3(0.3, 0.1, 0.0), getLinearVelocityOfRigidBody(0));
		applyForceOnBody(1, getPositionOfRigidBody(1) - Vec3(0.0, 0.2, 0.0), getLinearVelocityOfRigidBody(1));
		applyForceOnBody(2, getPositionOfRigidBody(2) - Vec3(0.1, 0.2, 0.0), getLinearVelocityOfRigidBody(2));
		applyForceOnBody(3, getPositionOfRigidBody(3) - Vec3(0.0, 0.3, 0.2), getLinearVelocityOfRigidBody(3));
		setOrientationOf(0, Quat(Vec3(1.0f, 0.0f, 1.0f), (float)(M_PI) * 0.35f));
		setOrientationOf(2, Quat(Vec3(1.0f, 0.0f, 1.0f), (float)(M_PI) * 0.25f));
		setOrientationOf(3, Quat(Vec3(0.0f, 1.0f, 1.0f), (float)(M_PI) * 0.45f));
		inputScale = 0.005f;
		break;
	default:
		cout << "Empty Test!\n";
		break;
	}
}

void RigidBodySystemSimulator::externalForcesCalculations(float timeElapsed)
{
	// Apply the mouse deltas to g_vfMovableObjectPos (move along cameras view plane)
	Point2D mouseDiff;
	mouseDiff.x = m_trackmouse.x - m_oldtrackmouse.x;
	mouseDiff.y = m_trackmouse.y - m_oldtrackmouse.y;
	if (mouseDiff.x != 0 || mouseDiff.y != 0)
	{
		Mat4 worldViewInv = Mat4(DUC->g_camera.GetWorldMatrix() * DUC->g_camera.GetViewMatrix());
		worldViewInv = worldViewInv.inverse();
		Vec3 inputView = Vec3((float)mouseDiff.x, (float)-mouseDiff.y, 0);
		Vec3 inputWorld = worldViewInv.transformVectorNormal(inputView);
		// find a proper scale!
		//float inputScale = 0.000005f; //7000fps+ :( so fast
		inputWorld = inputWorld * inputScale;
		for (size_t i = 0; i < getNumberOfRigidBodies(); i++)
		{
			//not sure is work fine, will push the thing back
			/*if(inputWorld.x < rigidboxVector[i].position.x + rigidboxVector[i].size.x/2 &&
				inputWorld.x > rigidboxVector[i].position.x - rigidboxVector[i].size.x / 2 &&
				inputWorld.y < rigidboxVector[i].position.y - rigidboxVector[i].size.y / 2 &&
				inputWorld.y > rigidboxVector[i].position.y - rigidboxVector[i].size.y / 2 )*/
			applyForceOnBody(i, rigidboxVector[i].position, inputWorld);
			applyForceOnBody(i, Vec3(inputWorld.x, inputWorld.y, rigidboxVector[i].position.z),
				Vec3(inputWorld.x, inputWorld.y, inputScale));
		}
		
	}

}

void RigidBodySystemSimulator::simulateTimestep(float timeStep)
{
	switch (m_iTestCase)
	{
		// handling different cases
	case 0:
		break;
	case 1:
		force_eular(0.01);
		break;
	case 2:
		for (int i = 0; i < rigidboxVector.size(); i++) {
			applyForceOnBody(i, getPositionOfRigidBody(i), Vec3(0.0f, -m_fGravity, 0.0f));
		}
		force_eular(timeStep);
		break;
	case 3:
		for (int i = 0; i < rigidboxVector.size(); i++) {
			applyForceOnBody(i, getPositionOfRigidBody(i), Vec3(0.0f, -m_fGravity, 0.0f));
			applyForceOnBody(i, getPositionOfRigidBody(i), Vec3(m_fAddVelocity_x, 0.0f, 0.0f));
			applyForceOnBody(i, getPositionOfRigidBody(i), Vec3( 0.0f, m_fAddVelocity_y, 0.0f));
			applyForceOnBody(i, getPositionOfRigidBody(i), Vec3( 0.0f, 0.0f, m_fAddVelocity_z));
		}
		force_eular(timeStep);
		break;
	default:
		force_eular(timeStep);
		break;
	}

}

void RigidBodySystemSimulator::onClick(int x, int y)
{
	m_trackmouse.x = x;
	m_trackmouse.y = y;

}

void RigidBodySystemSimulator::onMouse(int x, int y)
{
	m_oldtrackmouse.x = x;
	m_oldtrackmouse.y = y;
	m_trackmouse.x = x;
	m_trackmouse.y = y;

}

// ExtraFunctions

int RigidBodySystemSimulator::getNumberOfRigidBodies()
{
    return rigidboxVector.size();
}

Vec3 RigidBodySystemSimulator::getPositionOfRigidBody(int i)
{
    return rigidboxVector.at(i).position;
}

Vec3 RigidBodySystemSimulator::getLinearVelocityOfRigidBody(int i)
{
	return rigidboxVector.at(i).linearVelocity;
}

Vec3 RigidBodySystemSimulator::getAngularVelocityOfRigidBody(int i)
{
	return rigidboxVector.at(i).angularVelocity;
}

void RigidBodySystemSimulator::applyForceOnBody(int i, Vec3 loc, Vec3 force)
{
	Rigidbox& rb = rigidboxVector.at(i);
	rb.m_externalForce += force;

	rb.forceLocation = loc;
}

void RigidBodySystemSimulator::addRigidBody(Vec3 position, Vec3 size, int mass)
{
    Rigidbox rb;
	rb.position = position;
	rb.size = size;
	rb.mass = mass;
	rb.orientation = Quat(Vec3(0.0f, 0.0f, 0.0f), (float)(M_PI) * 0.0f);
	rb.linearVelocity = Vec3(0.0f, 0.0f, 0.0f);
	rb.angularVelocity = Vec3(0.0f, 0.0f, 0.0f);
	rb.m_externalForce = Vec3(0.0f, 0.0f, 0.0f);
	rb.forceLocation = Vec3(0.0f, 0.0f, 0.0f);
	rb.angularMomentum = Vec3(0.0f, 0.0f, 0.0f);
	// inertiaTensor
	float ixx = (rb.mass / 12.0f) * (rb.size.y * rb.size.y + rb.size.z * rb.size.z);
	float iyy = (rb.mass / 12.0f) * (rb.size.x * rb.size.x + rb.size.z * rb.size.z);
	float izz = (rb.mass / 12.0f) * (rb.size.x * rb.size.x + rb.size.y * rb.size.y);
	rb.inertiaTensor.initId();
	rb.inertiaTensor.value[0][0] = ixx;
	rb.inertiaTensor.value[1][1] = iyy;
	rb.inertiaTensor.value[2][2] = izz;
	rb.inertiaTensor.value[3][3] = 1;
	rb.inertiaTensorInv = rb.inertiaTensor.inverse();
	rigidboxVector.push_back(rb);
}

void RigidBodySystemSimulator::setOrientationOf(int i, Quat orientation)
{
	rigidboxVector.at(i).orientation = orientation;
}

void RigidBodySystemSimulator::setVelocityOf(int i, Vec3 velocity)
{
	rigidboxVector.at(i).linearVelocity = velocity;
}

void RigidBodySystemSimulator::force_eular(float timeStep) {

	for (Rigidbox& rb : rigidboxVector) {
		//for test
		//Rigidbox& rb = rigidboxVector.at(0);
		rb.torque += CalculateTorque(rb.m_externalForce, rb.forceLocation, rb.position);
		
		//Update the positions
		rb.position += timeStep * rb.linearVelocity;

		// linearVelocity
		rb.linearVelocity += (rb.m_externalForce * timeStep / (float)rb.mass);

		IntegrateOrientation(rb, timeStep);
		IntegrateAngularMomentum(rb, rb.torque, timeStep);
		UpdateInverseInertiaTensor(rb);
		UpdateAngularVelocity(rb);
		rb.m_externalForce = Vec3(0.0f, 0.0f, 0.0f);
		rb.torque = Vec3(0.0f, 0.0f, 0.0f);

		//Update the world space positions
		Mat4 position_M4;
		position_M4.initId();
		position_M4.initTranslation(rb.position.x, rb.position.y, rb.position.z);
		Mat4 size_M4;
		size_M4.initId();
		size_M4.initScaling(rb.size.x, rb.size.y, rb.size.z);
		rb.worldMatrix = size_M4 * rb.orientation.getRotMat() * position_M4;

		//// torque= xi ¡Á F
		//Vec3 xi = rb.forceLocation - rb.position;
		//Vec3 torque = cross(xi, rb.m_externalForce);

		//// r' = r + (0 w) * r  * t / 2

		//// default angularMomentum L+= h * q
		//rb.angularMomentum = timeStep * torque;

		////  Inertia Tensor I = rot r * i0 -1 * T rot r
		////
		////Update angular velocity w using I and L w(t+h) = I * L(t+h)
		//rb.angularVelocity.x = rb.inertiaTensorInv.x * rb.angularMomentum.x;
		//rb.angularVelocity.y = rb.inertiaTensorInv.y * rb.angularMomentum.y;
		//rb.angularVelocity.z = rb.inertiaTensorInv.z * rb.angularMomentum.z;
		
		//wall and floor
		//hit wall and floor reduce 1% Velocity and L
		Vec3 wallLimits(5.0f, 7.0f, 5.0f);
		float reduce = 0.99f;
		if (rb.position.x >= wallLimits.x || rb.position.x <= -wallLimits.x) {
			rb.linearVelocity.x = -rb.linearVelocity.x;
			rb.linearVelocity = rb.linearVelocity * reduce;
			rb.angularMomentum = -rb.angularMomentum * reduce;
		}
		if (rb.position.y >= wallLimits.y || rb.position.y <= -0.2) {
			rb.linearVelocity.y = -rb.linearVelocity.y;
			rb.linearVelocity = rb.linearVelocity * reduce;
			rb.angularMomentum = -rb.angularMomentum * reduce;
		}
		if (rb.position.z >= wallLimits.z || rb.position.z <= -wallLimits.z) {
			rb.linearVelocity.z = -rb.linearVelocity.z ;
			rb.linearVelocity = rb.linearVelocity * reduce;
			rb.angularMomentum = -rb.angularMomentum * reduce;
		}
		
	}
	for (int i = 0; i < rigidboxVector.size(); i++) {
		for (int j = i + 1; j < rigidboxVector.size(); j++) {
			CollisionInfo  cosInfo = checkCollisionSAT(rigidboxVector.at(i).worldMatrix, rigidboxVector[j].worldMatrix);
			if (cosInfo.isValid) {
				//cout << cosInfo.depth<<"\n";
				force_Impulse(cosInfo, rigidboxVector.at(i), rigidboxVector.at(j));
			}
		}
	}
	
}

Vec3 RigidBodySystemSimulator::CalculateTorque(const Vec3& force, const Vec3& forceLocation, const Vec3& position) {
	return cross(forceLocation - position, force);
}

void RigidBodySystemSimulator::IntegrateOrientation(Rigidbox& rb, float timeStep) {
	Quat quat;
	quat.w = 0;
	quat.x = rb.angularVelocity.x;
	quat.y = rb.angularVelocity.y;
	quat.z = rb.angularVelocity.z;
	Quat tmp1 = QuaternionMultiply(quat,rb.orientation);
	rb.orientation += tmp1.operator * (timeStep / 2);
	rb.orientation = rb.orientation.unit();
}

void  RigidBodySystemSimulator::IntegrateAngularMomentum(Rigidbox& rb, const Vec3& torque, float timeStep) {
	rb.angularMomentum += timeStep * torque;
}

void  RigidBodySystemSimulator::UpdateInverseInertiaTensor(Rigidbox& rb) {

	//rb.orientation.getRotMat().inverse().transformVector(rb.orientation.getRotMat().transformVector(rb.inertiaTensorInv))
	Mat4 tranRot = rb.orientation.getRotMat();
	tranRot.transpose();
	rb.inertiaTensorInv = rb.orientation.getRotMat() * rb.inertiaTensor.inverse() * tranRot;

}

void  RigidBodySystemSimulator::UpdateAngularVelocity(Rigidbox& rb) {
	rb.angularVelocity = rb.inertiaTensorInv.transformVector(rb.angularMomentum);
}

Quat RigidBodySystemSimulator::QuaternionMultiply(const Quat& q1, const Quat& q2) {
	Quat result;
	result.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
	result.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
	result.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
	result.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
	return result;
}

void RigidBodySystemSimulator::force_Impulse(CollisionInfo cosInfo, Rigidbox& rb1, Rigidbox& rb2) {
	float  c = 0.8f; // bounciness_coefficient
	Vec3 n = cosInfo.normalWorld; // normal of the collision

	// Collision point relative to the center of mass of rigid body A / B
	Vec3 x_a = cosInfo.collisionPointWorld - rb1.position;
	Vec3 x_b = cosInfo.collisionPointWorld - rb2.position;
	Vec3 v_a = (rb1.linearVelocity + cross(rb1.angularVelocity, x_a));
	Vec3 v_b = (rb2.linearVelocity + cross(rb2.angularVelocity, x_b));
	Vec3 v_rel = v_a - v_b;
	if (dot(v_rel,n)>0) {
		// separating
		return;
	}
	// -(1+c)v_rel . n
	double top = dot(-(1 + c) *  v_rel, n);
	// 1/Ma + 1/Mb
	double bot_mass = (1.0f / rb1.mass + 1.0f / rb2.mass);
	// I-1*(x_a * n) * x_a --------- same b
	Vec3 bot_cross_a = cross(rb1.inertiaTensor.inverse() * (cross(x_a, n)), x_a);
	Vec3 bot_cross_b = cross(rb2.inertiaTensor.inverse() * (cross(x_b, n)), x_b);
	double bottom = bot_mass + dot(bot_cross_a + bot_cross_b, n);
	double J = top / bottom;

	rb1.linearVelocity += J * n / rb1.mass;
	rb2.linearVelocity -= J * n / rb2.mass;

	rb1.angularMomentum += cross(x_a, J * n);
	rb2.angularMomentum -= cross(x_b, J* n);
}

// Wei
#include "MassSpringSystemSimulator.h"

MassSpringSystemSimulator::MassSpringSystemSimulator()
{
	reset();
	m_fMass = 10.0f;
	m_fStiffness = 40.0f;
	m_fDamping = 0.0f;
	m_externalForce = Vec3();
	massPointSize = 0.05f;
	lineColor = Vec3(255, 0, 0);
	defalut_initlength = 1.0f;
	m_fGravity = 0.0f;
	m_fTimestep = 0.005f;
	m_numberOfSpring = 10;
}


// UI Functions
const char* MassSpringSystemSimulator::getTestCasesStr()
{
	// TODO:
	return " Demo 1: One-step test,\
				Demo 2: Euler simulation,\
				Demo 3: Midpoint simulation,\
				Demo 4: Compare the stability of Euler and Midpoint method,\
				Demo 5: Leap-Frog,\
				Demo 4: another one";
}

void MassSpringSystemSimulator::initUI(DrawingUtilitiesClass* DUC)
{
	this->DUC = DUC;
	switch (m_iTestCase)
	{
	case 0:
	{
		TwAddSeparator(DUC->g_pTweakBar, "", NULL);
		TwEnumVal IntegrateEV[] = { {0, "EULER"}, {1, "LEAPFROG"}, {2, "MIDPOINT"} };
		TwType IntegrateTwType = TwDefineEnum("Integrate", IntegrateEV, 3);
		setIntegrator(EULER);
		TwAddVarRW(DUC->g_pTweakBar, "Integrate", IntegrateTwType, &m_iIntegrator, NULL);
		break;
	}
	case 1: {	
		TwRemoveVar(DUC->g_pTweakBar, "Timestep");
		break;
	}
	case 2: {
		TwRemoveVar(DUC->g_pTweakBar, "Timestep");
		break;
	}
	case 3:
	{
		TwAddSeparator(DUC->g_pTweakBar, "", NULL);
		TwEnumVal IntegrateEV[] = { {0, "EULER"}, {1, "LEAPFROG"}, {2, "MIDPOINT"} };
		TwType IntegrateTwType = TwDefineEnum("Integrate", IntegrateEV, 3);
		TwAddVarRW(DUC->g_pTweakBar, "Integrate", IntegrateTwType, &m_iIntegrator, NULL);
		TwAddVarRW(DUC->g_pTweakBar, "(reset required) Number of spring ", TW_TYPE_INT8, &m_numberOfSpring, "min=10 max=25");
		TwAddVarRW(DUC->g_pTweakBar, "Gravity", TW_TYPE_FLOAT, &m_fGravity, "min=0 max=5 step=0.01");
		TwAddVarRW(DUC->g_pTweakBar, "Damping", TW_TYPE_FLOAT, &m_fDamping, "min=0 max=3 step=0.01");
		break;
	}
	case 4:break;
	case 5: {
		TwAddSeparator(DUC->g_pTweakBar, "", NULL);
		TwEnumVal IntegrateEV[] = { {0, "EULER"}, {1, "LEAPFROG"}, {2, "MIDPOINT"} };
		TwType IntegrateTwType = TwDefineEnum("Integrate", IntegrateEV, 3);
		TwAddVarRW(DUC->g_pTweakBar, "Integrate", IntegrateTwType, &m_iIntegrator, NULL);
		TwAddVarRW(DUC->g_pTweakBar, "(reset required) Number of spring ", TW_TYPE_INT8, &m_numberOfSpring, "min=10 max=25");
		TwAddVarRW(DUC->g_pTweakBar, "Gravity", TW_TYPE_FLOAT, &m_fGravity, "min=0 max=5 step=0.01");
		TwAddVarRW(DUC->g_pTweakBar, "Damping", TW_TYPE_FLOAT, &m_fDamping, "min=0 max=3 step=0.01");
		break;
	}
	default:break;
	}
}

void MassSpringSystemSimulator::reset()
{
	masspointVector.clear();
	springVector.clear();
	m_mouse.x = m_mouse.y = 0;
	m_trackmouse.x = m_trackmouse.y = 0;
	m_oldtrackmouse.x = m_oldtrackmouse.y = 0;
	m_fDamping = 0.0f;
	m_fGravity = 0.0f;
}

void MassSpringSystemSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext)
{
		for (Spring s : springVector) {
			DUC->drawSphere(getPositionOfMassPoint(s.masspoint1), Vec3(massPointSize, massPointSize, massPointSize));
			DUC->drawSphere(getPositionOfMassPoint(s.masspoint2), Vec3(massPointSize, massPointSize, massPointSize));
			DUC->beginLine();
			DUC->drawLine(getPositionOfMassPoint(s.masspoint1), lineColor, getPositionOfMassPoint(s.masspoint2), lineColor);
			DUC->endLine();
		}

}

void MassSpringSystemSimulator::notifyCaseChanged(int testCase)
{
	m_iTestCase = testCase;
	reset();
	switch (m_iTestCase)
	{
	case 0:
	{
		cout << "Demo 1: One-step test.\n";
		int index_mp1 = addMassPoint(Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), false);
		int index_mp2 = addMassPoint(Vec3(0.0f, 2.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), false);
		addSpring(index_mp1, index_mp2, defalut_initlength);
		break;
	}
	case 1:
	{
		cout << "Demo 2: Euler simulation.\n";
		int index_mp1 = addMassPoint(Vec3(0.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0f), false);
		int index_mp2 = addMassPoint(Vec3(0.0f, 2.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f), false);
		addSpring(index_mp1, index_mp2, defalut_initlength);
		m_iIntegrator = EULER;
		m_fTimestep = 0.005f;
		break;
	}
	case 2:
	{
		cout << "Demo 3: Midpoint simulation.\n";
		int index_mp1 = addMassPoint(Vec3(0.0f, 0.0f, 0.0f), Vec3(-1.0, 0.0f, 0.0f), false);
		int index_mp2 = addMassPoint(Vec3(0.0f, 2.0f, 0.0f), Vec3(1.0, 0.0f, 0.0f), false);
		addSpring(index_mp1, index_mp2, defalut_initlength);
		m_iIntegrator = MIDPOINT;
		m_fTimestep = 0.005f;
		break;
	}
	case 3:
	{
		cout << "Demo 4: Compare the stability of Euler and Midpoint method.\n";
		
		for (int i = 0; i < m_numberOfSpring; i++)
		{
			float shift = i * 0.5f;
			int index_mp1 = addMassPoint(Vec3(-2.0f + shift, 0.0f + shift/2, -5.0f + shift), Vec3(-1.0, 0.0f, 0), false);
		}
		for (int i = 1; i < (int)masspointVector.size()/2 ; i++) {
			addSpring(0, i , defalut_initlength);
		}
		for (int i = (int)masspointVector.size() / 2; i < masspointVector.size()-1; i++) {
			addSpring(i, i+1, defalut_initlength);
		}

		break;
	}
	case 4:
	{
		cout << "Demo 5: Leap-Frog.\n";
		int index_mp1 = addMassPoint(Vec3(0.0, 0.0f, 0), Vec3(-1.0, 0.0f, 0), false);
		int index_mp2 = addMassPoint(Vec3(0.0, 2.0f, 0), Vec3(1.0, 0.0f, 0), false);
		addSpring(index_mp1, index_mp2, defalut_initlength);
		m_iIntegrator = LEAPFROG;
		break;
	}
	case 5:
	{
		cout << "Demo 4: Another Compare the stability of Euler and Midpoint method.\n";

		for (int i = 0; i < m_numberOfSpring; i++)
		{
			float shift = i * 0.5f;
			int index_mp1 = addMassPoint(Vec3(-5.0f + shift, 0.0f, -5.0f + shift), Vec3(-1.0, 0.0f, 0), false);
			int index_mp2 = addMassPoint(Vec3(-5.0f + shift, 2.0f, -5.0f + shift), Vec3(1.0, 0.0f, 0), false);
			addSpring(index_mp1, index_mp2, defalut_initlength);
		}

		break;
	}
	default:
		cout << "Empty Test!\n";
		break;
	}
}

void MassSpringSystemSimulator::externalForcesCalculations(float timeElapsed)
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
		float inputScale = 0.0000001f;
		inputWorld = inputWorld * inputScale;
		for (Masspoint& mp : masspointVector) {
			mp.position = mp.position + inputWorld;
		}
	}
}

void MassSpringSystemSimulator::simulateTimestep(float timeStep)
{	
	// update current setup for each frame
	switch (m_iTestCase)
	{
		// handling different cases
	case 0:
		//Demo 1: One-step test.
		if (m_iIntegrator == EULER) {
			explicitEuler(timeStep);
		}
		else if (m_iIntegrator == MIDPOINT) {
			midpoint(timeStep);
		}
		else if (m_iIntegrator == LEAPFROG) {
			leapfrog(timeStep);
		}
		break;
	case 1:
		explicitEuler(m_fTimestep);
		break;
	case 2:
		midpoint(m_fTimestep);
		break;
	case 3:
	{
		if (m_iIntegrator == EULER) {
			explicitEuler(timeStep);
		}
		else if (m_iIntegrator == MIDPOINT) {
			midpoint(timeStep);
		}
		else if (m_iIntegrator == LEAPFROG) {
			leapfrog(timeStep);
		}
		break;
	}
	case 4:
		leapfrog(timeStep);
		break;
	default:
		if (m_iIntegrator == EULER) {
			explicitEuler(timeStep);
		}
		else if (m_iIntegrator == MIDPOINT) {
			midpoint(timeStep);
		}
		else if (m_iIntegrator == LEAPFROG) {
			leapfrog(timeStep);
		}
		break;
	}
}

void MassSpringSystemSimulator::onClick(int x, int y)
{
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

void MassSpringSystemSimulator::onMouse(int x, int y)
{
	m_oldtrackmouse.x = x;
	m_oldtrackmouse.y = y;
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

// Specific Functions
void MassSpringSystemSimulator::setMass(float mass)
{
	m_fMass = mass;
}

void MassSpringSystemSimulator::setStiffness(float stiffness)
{
	m_fStiffness = stiffness;
}

void MassSpringSystemSimulator::setDampingFactor(float damping)
{
	m_fDamping = damping;
}

int MassSpringSystemSimulator::addMassPoint(Vec3 position, Vec3 Velocity, bool isFixed)
{	
	Masspoint mp{};
	mp.position = position;
	mp.Velocity = Velocity;
	mp.isFixed = isFixed;
	masspointVector.push_back(mp);
	return (int)masspointVector.size()-1;
}

void MassSpringSystemSimulator::addSpring(int masspoint1, int masspoint2, float initialLength)
{
	Spring sp{};
	sp.masspoint1 = masspoint1;
	sp.masspoint2 = masspoint2;
	sp.initialLength = initialLength;
	springVector.push_back(sp);
}

int MassSpringSystemSimulator::getNumberOfMassPoints()
{
	return (int)masspointVector.size();
}

int MassSpringSystemSimulator::getNumberOfSprings()
{
	return (int)springVector.size();
}

Vec3 MassSpringSystemSimulator::getPositionOfMassPoint(int index)
{
	return masspointVector.at(index).position;
}

Vec3 MassSpringSystemSimulator::getVelocityOfMassPoint(int index)
{
	return masspointVector.at(index).Velocity;
}

void MassSpringSystemSimulator::applyExternalForce(Vec3 force)
{
}

// Explicit Euler Methode Demo 2
void MassSpringSystemSimulator::explicitEuler(float timeStep)
{
	for (Spring& sp : springVector) {
		Masspoint& mp1 = masspointVector.at(sp.masspoint1);
		Masspoint& mp2 = masspointVector.at(sp.masspoint2);
		Masspoint tmpP1{};
		Masspoint tmpP2{};
		tmpP1.position = mp1.position;
		tmpP2.position = mp2.position;
		tmpP1.Velocity = mp1.Velocity;
		tmpP2.Velocity = mp2.Velocity;
		tmpP1.isFixed = false;
		tmpP2.isFixed = false;
		// calc acc first or will replace position
		std::tuple <Vec3, Vec3> acc_tuple = acceleration(tmpP1, tmpP2);
		// x(i+1) = x(i) + v(i) * dt
		mp1.position = tmpP1.position + tmpP1.Velocity * timeStep;
		mp2.position = tmpP2.position + tmpP2.Velocity * timeStep;
		// v(i+1) = v(i) + a(i) * dt
		mp1.Velocity = tmpP1.Velocity + std::get<0>(acc_tuple) * timeStep;
		mp2.Velocity = tmpP2.Velocity + std::get<1>(acc_tuple) * timeStep;
		wallCollision(mp1);
		wallCollision(mp2);
	}
}

// Midpoint Methode  Demo 3
void MassSpringSystemSimulator::midpoint(float timeStep) {
	for (Spring& sp : springVector) {
		Masspoint& mp1 = masspointVector.at(sp.masspoint1);
		Masspoint& mp2 = masspointVector.at(sp.masspoint2);
		// calc acc first or will replace position
		std::tuple <Vec3, Vec3> acc_tuple = acceleration(mp1, mp2);
		// x(i+1) = x(i) + v(i) * dt
		Vec3 midPo1 = mp1.position + mp1.Velocity * timeStep / 2;
		Vec3 midPo2 = mp2.position + mp2.Velocity * timeStep / 2;
		// v(i+1) = v(i) + a(i) * dt
		Vec3 midVel1 = mp1.Velocity + std::get<0>(acc_tuple) * timeStep / 2;
		Vec3 midVel2 = mp2.Velocity + std::get<1>(acc_tuple) * timeStep / 2;
		Masspoint tmpP1{};	
		Masspoint tmpP2{};
		tmpP1.position = midPo1;
		tmpP2.position = midPo2;
		tmpP1.Velocity = midVel1;
		tmpP2.Velocity = midVel2;
		tmpP1.isFixed = false;
		tmpP2.isFixed = false;
		acc_tuple = acceleration(tmpP1, tmpP2);
		// x(i+1) = x(i) + v(i) * dt
		mp1.position = mp1.position + midVel1 * timeStep;
		mp2.position = mp2.position + midVel2 * timeStep;
		// v(i+1) = v(i) + a(i) * dt
		mp1.Velocity = mp1.Velocity + std::get<0>(acc_tuple) * timeStep;
		mp2.Velocity = mp2.Velocity + std::get<1>(acc_tuple) * timeStep;
		wallCollision(mp1);
		wallCollision(mp2);
	}
}

// Leapfrog Methode  Demo 5
void MassSpringSystemSimulator::leapfrog(float timeStep) {
	for (Spring& sp : springVector) {
		Masspoint& mp1 = masspointVector.at(sp.masspoint1);
		Masspoint& mp2 = masspointVector.at(sp.masspoint2);

		std::tuple<Vec3, Vec3> acc_tuple = acceleration(mp1, mp2);
		mp1.Velocity = mp1.Velocity + std::get<0>(acc_tuple) * timeStep;
		mp2.Velocity = mp2.Velocity + std::get<1>(acc_tuple) * timeStep;

		mp1.position = mp1.position + mp1.Velocity * (timeStep / 2.0f);
		mp2.position = mp2.position + mp2.Velocity * (timeStep / 2.0f);

		wallCollision(mp1);
		wallCollision(mp2);
	}
}


std::tuple <Vec3, Vec3>  MassSpringSystemSimulator::acceleration(Masspoint mp1, Masspoint mp2)
{	
	// a=F/m, F=kx => a=kx/m, x is distance change.
	Vec3 dist = mp1.position - mp2.position;
	double distLength = sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z);
	Vec3 distDirection = Vec3(dist.x / distLength, dist.y / distLength, dist.z / distLength);
	Vec3  F_spring_mp1 = -m_fStiffness * (distLength - defalut_initlength) * distDirection;
	Vec3  F_spring_mp2 = -F_spring_mp1;
	Vec3  F_damping_mp1 = -m_fDamping * mp1.Velocity;
	Vec3  F_damping_mp2 = -m_fDamping * mp2.Velocity;
	Vec3  F_gravity = -m_fMass * m_fGravity;
	Vec3 F_total_mp1 = F_spring_mp1 + F_damping_mp1 + F_gravity;
	Vec3 F_total_mp2 = F_spring_mp2 + F_damping_mp2 + F_gravity;
	Vec3 acc_mp1 = F_total_mp1 / m_fMass;
	Vec3 acc_mp2 = F_total_mp2 / m_fMass;

	return { acc_mp1, acc_mp2 };
}

Vec3 MassSpringSystemSimulator::direction(Vec3 p1, Vec3 p2) {
	auto x = p1.x - p2.x;
	auto y = p1.y - p2.y;
	auto z = p1.z - p2.z;
	return Vec3(x, y, z);
}

void MassSpringSystemSimulator::wallCollision(Masspoint& mp) {
	Vec3 wallLimits(5.0f, 10.0f, 5.0f);
	if (mp.position.x >= wallLimits.x || mp.position.x <= -wallLimits.x) {
		mp.Velocity.x = -mp.Velocity.x; 
	}
	if (mp.position.y >= wallLimits.y || mp.position.y <= -1) {
		mp.Velocity.y = -mp.Velocity.y; 
	}
	if (mp.position.z >= wallLimits.z || mp.position.z <= -wallLimits.z) {
		mp.Velocity.z = -mp.Velocity.z; 
	}
}



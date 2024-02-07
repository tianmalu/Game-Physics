#include "PinballSimulator.h"

// Construtors
PinballSimulator:: PinballSimulator()
{
    reset();
    damping = 0.0f;
}

// UI Functions
const char* PinballSimulator::getTestCasesStr()
{
	return " Demo1: Pinball 1,\
            Demo 2: Pinball 2,\
			Demo 3: Pinball 3,\
			Demo 4: Pinball 4,";

}

void PinballSimulator:: initUI(DrawingUtilitiesClass* DUC)
{
    this->DUC=DUC;
    TwAddVarRW(DUC->g_pTweakBar, "Damping", TW_TYPE_FLOAT, &damping, "min=0 max=0.1 step=0.01");

}

void PinballSimulator::reset()
{
    m_mouse.x = m_mouse.y = 0;
    m_trackmouse.x = m_trackmouse.y = 0;
    m_oldtrackmouse.x = m_oldtrackmouse.y = 0;
    damping = 0.0f;
    rigidwallVector.clear();
    pinballVector.clear();
    springVector.clear();
}

void PinballSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext)
{
	DUC->setUpLighting(Vec3(0, 0, 0), 0.4 * Vec3(2, 2, 2), 2000.0, Vec3(0.5, 0.5, 0.5));
    //Pinball
    DUC->drawSphere(getPositionOfPinball(), pinballVector.at(0).pinballSize);
    //Spring
    DUC->beginLine();
    //Vec3=(length,height,width)
	if(getPositionOfPinball().x==0 && getPositionOfPinball().z <= springVector.at(0).length-0.55f) DUC->drawLine(getPositionOfPinball(), Vec3(255, 0, 0), Vec3(0.0f, 0.0f, -0.55f), Vec3(255, 0, 0));
    else if (getPositionOfPinball().x == -0.85f && getPositionOfPinball().y <= springVector.at(0).length - 0.55f)
    {
        DUC->drawLine(getPositionOfPinball(), Vec3(255, 0, 0), Vec3(-0.85f, 0.0f, -0.125f), Vec3(255, 0, 0));
    }
    else DUC->drawLine(Vec3(0.0f, 0.0f, springVector.at(0).length-0.55f), Vec3(255, 0, 0), Vec3(0.0f, 0.0f, -0.55f), Vec3(255, 0, 0));
	DUC->endLine();
    //Playboard
    RigidWall playboard;
    Mat4 position_M4 = Mat4(0.0);
	position_M4.initTranslation(playboard.position.x, playboard.position.y, playboard.position.z);
	Mat4 size_M4 = Mat4(0.0);
	size_M4.initScaling(playboard.size.x, playboard.size.y, playboard.size.z);
    Quat pbOrientation=Quat(Vec3(0.0f, 0.0f, 0.0f),(float)(M_PI) * 0.5f);
	Mat4 worldMatrix = size_M4 * pbOrientation.getRotMat() * position_M4;
	DUC->drawRigidBody(worldMatrix);
    //Wall
    for (RigidWall& rb : rigidwallVector) {
		Mat4 position_M4 = Mat4(0.0);
		position_M4.initTranslation(rb.position.x, rb.position.y, rb.position.z);
		Mat4 size_M4 = Mat4(0.0);
		size_M4.initScaling(rb.size.x, rb.size.y, rb.size.z);
		rb.worldMatrix = size_M4 * rb.orientation.getRotMat() * position_M4;
		DUC->drawRigidBody(rb.worldMatrix);
	}
	
}

void PinballSimulator::notifyCaseChanged(int testCase){
    m_iTestCase = testCase;
    reset();
    switch (testCase)
    {
    case 0:
        cout << "Demo: Pinball 1.\n";
        addPlayboard();
        // Basic Wall objects
        addRigidWall(Vec3(-0.05f, 0.06f, 0.45f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 2.0f));
        addRigidWall(Vec3(0.95f, 0.06f, -0.55f), Vec3(0.0f, 1.85f, 0.0f), Vec3(0.01f, 0.25f, 1.7f));
        addRigidWall(Vec3(1.95f, 0.06f, 0.40f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 1.9f));
        addRigidWall(Vec3(0.06f, 0.06f, 0.45f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 1.6f));
        addRigidWall(Vec3(0.5f, 0.06f, 2.00f), Vec3(0.0f, 1.0f, 0.0f), Vec3(0.01f, 0.25f, 1.55f));
        addRigidWall(Vec3(1.5f, 0.06f, 1.8f), Vec3(0.0f, -1.0f, 0.0f), Vec3(0.01f, 0.25f, 1.7f));
        //You can add any etxra WALL here.
        //addRigidWall(Vec3(1.2f, 0.06f, 0.4f), Vec3(0.0f, -0.8f, 0.0f), Vec3(0.01f, 0.25f, 2.2f));
        //Pinball
        addPinball();
        //Spring
        addSpring(0.85f, 4.0f);
        break;
    case 1:
        cout << "Demo: Pinball 2.\n";
        addPlayboard();
        // Basic Wall objects
        addRigidWall(Vec3(-0.05f, 0.06f, 0.45f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 2.0f));
        addRigidWall(Vec3(0.95f, 0.06f, -0.55f), Vec3(0.0f, 1.85f, 0.0f), Vec3(0.01f, 0.25f, 1.7f));
        addRigidWall(Vec3(1.95f, 0.06f, 0.40f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 1.9f));
        addRigidWall(Vec3(0.06f, 0.06f, 0.45f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 1.6f));
        addRigidWall(Vec3(0.5f, 0.06f, 2.00f), Vec3(0.0f, 1.0f, 0.0f), Vec3(0.01f, 0.25f, 1.55f));
        addRigidWall(Vec3(1.5f, 0.06f, 1.8f), Vec3(0.0f, -1.0f, 0.0f), Vec3(0.01f, 0.25f, 1.7f));
        //You can add any etxra WALL here.
        addRigidWall(Vec3(1.2f, 0.06f, 0.4f), Vec3(0.0f, -0.8f, 0.0f), Vec3(0.01f, 0.25f, 2.2f));
        //Pinball
        addPinball();
        //Spring
        addSpring(0.85f, 4.0f);
        break;
    case 2:
        cout << "Demo: Pinball 3.\n";
        addPlayboard();
        // Basic Wall objects
        addRigidWall(Vec3(-0.05f, 0.06f, 0.45f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 2.0f));
        addRigidWall(Vec3(0.95f, 0.06f, -0.55f), Vec3(0.0f, 1.85f, 0.0f), Vec3(0.01f, 0.25f, 1.7f));
        addRigidWall(Vec3(1.95f, 0.06f, 0.40f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 1.9f));
        addRigidWall(Vec3(0.06f, 0.06f, 0.45f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 1.6f));
        addRigidWall(Vec3(0.5f, 0.06f, 2.00f), Vec3(0.0f, 1.0f, 0.0f), Vec3(0.01f, 0.25f, 1.55f));
        addRigidWall(Vec3(1.5f, 0.06f, 1.8f), Vec3(0.0f, -1.0f, 0.0f), Vec3(0.01f, 0.25f, 1.7f));
        //You can add any etxra WALL here.
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 3; j++) {
                addRigidWall(Vec3(0.25f + i * 0.3f, 0.06f, 0.2f + j * 0.3f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 0.1f));
                addRigidWall(Vec3(0.25f + i * 0.3f, 0.06f, 0.2f + j * 0.3f), Vec3(0.0f, 0.7f, 0.0f), Vec3(0.01f, 0.25f, 0.1f));
                addRigidWall(Vec3(0.25f + i * 0.3f, 0.06f, 0.2f + j * 0.3f), Vec3(0.0f, -1.0f, 0.0f), Vec3(0.01f, 0.25f, 0.1f));
            }
        }

        //Pinball
        addPinball();


        //Spring
        addSpring(0.85f, 4.0f);
        break;
    case 3:
        {
        cout << "Demo4:  Pinball 4.\n";
        RigidWall playboard{};
        playboard.position = Vec3(0.0f, 1.5f, 0.0f);
        playboard.size = Vec3(2.0f, 0.01f, 3.0f);
        playboard.orientation = Quat(Vec3(-1.0f, 0.0f, 0.0f), (float)(M_PI) * 0.5f);
        rigidwallVector.push_back(playboard);
        // Basic Wall objects
        addRigidWall2(Vec3(-1.0, 1.0f, -0.125f), Vec3(-1.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 2.0f));
        float initialX = -1.0f;
        float initialY = 2.0f;
        float initialZ = -0.125f;
        float angleIncrement = (float)(M_PI) * 1.0f / 36.0f; 

        float minSizeZ = 0.005f;
        float maxSizeZ = 0.03f;

        for (int i = 0; i <37; ++i) {
            Quat rotation(Vec3(0.0f, 0.0f, 1.0f), angleIncrement * i); 

            Vec3 size(0.01f, 1.0f*(float)(M_PI) / 36.0f, 0.25f);
            float angle = angleIncrement * i;
            float radius = 1.0f;
            float x = radius * cos(angle);
            float y = radius * sin(angle);
            Vec3 position;
            if(i == 0)
            {
                position = Vec3(x, initialY, initialZ);
            }else
            {
                position = Vec3(x, initialY + y, initialZ);
            }

            addRigidWall_rotation(position, rotation, size, -angle);
        }
        addRigidWall2(Vec3(1.0, 1.0f, -0.125f), Vec3(-1.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 2.0f));
        addRigidWall2(Vec3(-0.7, 1.0f, -0.125f), Vec3(-1.0f, 0.0f, 0.0f), Vec3(0.01f, 0.25f, 1.8f));
        RigidWall rw{};
        rw.position = Vec3(0.0, 0.0f, -0.125f);
        rw.orientation = Quat(Vec3(0.0f, 0.0f, 1.0f), (float)(M_PI) * 0.5f);
        rw.size = Vec3(0.01f, 2.0f, 0.25f);
        rw.direction =  (float)(M_PI) * 0.5f;
        rigidwallVector.push_back(rw);
        //Pinball
        Pinball pb{};
        pb.position = Vec3(-0.85f, 0.11f, -0.125f);
        pb.linearVelocity = Vec3(0.0f, 0.0f, 0.0f);
        pb.angularVelocity = Vec3(0.0f, 0.0f, 0.0f);
        pb.pinballSize = Vec3(0.05f, 0.05f, 0.05f);
        pb.mass = 1.0f;
        pb.orientation = Quat(Vec3(0.0f, 0.0f, 0.0f), (float)(M_PI) * 0.25f);
        pinballVector.push_back(pb);
        //Spring
        addSpring(0.85f, 4.0f);
        break;
    }
    default:
    {
        break;
    }
    }
    

}

void PinballSimulator::addPinball()
{
    Pinball pb{};
    pb.position = Vec3 (0.0f, 0.0f, -0.35f);
    pb.linearVelocity = Vec3 (0.0f, 0.0f, 0.0f);
    pb.angularVelocity = Vec3 (0.0f, 0.0f, 0.0f);
    pb.pinballSize = Vec3 (0.05f, 0.05f, 0.05f);
    pb.mass = 3.0f;
    pb.orientation = Quat(Vec3(0.0f, 0.0f, 0.0f),(float)(M_PI) * 0.25f);
    pinballVector.push_back(pb);
};
void PinballSimulator::addSpring(float initialLength, float stiffness){
    Spring sp{};
	sp.length = initialLength;
    sp.stiffness = stiffness;
    springVector.push_back(sp);

};
void PinballSimulator::addRigidWall(Vec3 position, Vec3 orientation, Vec3 size){
    RigidWall rw{};
    rw.position = position;
    rw.orientation = Quat(orientation, (float)(M_PI) * 0.25f);
    rw.size = size;
    rw.direction = orientation.y *(float)(M_PI) * 0.25f;
    rigidwallVector.push_back(rw);
};

void PinballSimulator::addRigidWall2(Vec3 position, Vec3 orientation, Vec3 size) {
    RigidWall rw{};
    rw.position = position;
    rw.orientation = Quat(orientation, (float)(M_PI) * 0.5f);
    rw.size = size;
    rw.direction = orientation.z * (float)(M_PI) * 0.25f;
    rigidwallVector.push_back(rw);
};

void PinballSimulator::addRigidWall_rotation(Vec3 position, Quat orientation, Vec3 size, float angle) {
    RigidWall rw{};
    rw.position = position;
    rw.size = size;
    rw.orientation = orientation;
    rw.direction = angle;
    rigidwallVector.push_back(rw);
};

void PinballSimulator::addPlayboard(){
    RigidWall playboard{};
    playboard.position = Vec3 (0.95f, -0.06f, 0.95f);
    playboard.size = Vec3(2.0f, 0.01f, 3.0f);
    rigidwallVector.push_back(playboard);
};

void PinballSimulator::onClick(int x, int y)
{
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

void PinballSimulator::onMouse(int x, int y)
{
	m_oldtrackmouse.x = x;
	m_oldtrackmouse.y = y;
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

void PinballSimulator::longPressOnMouse(float deltaTime){

}

void PinballSimulator::externalForcesCalculations(float timeElapsed)
{
    //Apply the mouse long press to spring to get an initial velocity.
}

void PinballSimulator:: simulateTimestep(float timeStep){
    switch (m_iTestCase)
    {
        // handling different cases
    case 0:
        springForce(timeStep);
        forceEular(timeStep);
        break;
    case 3:
        springForce(timeStep);
        forceEular2(timeStep);
        break;
    default:
        springForce(timeStep);
        forceEular(timeStep);
        break;
    }
}

void PinballSimulator::setWorldMatrix(Mat4 worldMatrix){

}

Mat4 PinballSimulator::getWorldMatrixOfPinball(){
    Pinball pb = pinballVector.at(0);
    Mat4 position_M4 = Mat4(0.0);
	position_M4.initTranslation(pb.position.x, pb.position.y, pb.position.z);
	Mat4 size_M4 = Mat4(0.0);
	size_M4.initScaling(pb.pinballSize.x, pb.pinballSize.y, pb.pinballSize.z);
    return size_M4 * pb.orientation.getRotMat() * position_M4;
}

Mat4 PinballSimulator::getWorldMatrixOfRigidwall(int index){
    RigidWall rw = rigidwallVector.at(index);
    Mat4 position_M4 = Mat4(0.0);
	position_M4.initTranslation(rw.position.x, rw.position.y, rw.position.z);
	Mat4 size_M4 = Mat4(0.0);
	size_M4.initScaling(rw.size.x, rw.size.y, rw.size.z);
    return size_M4 * rw.orientation.getRotMat() * position_M4;
};

void PinballSimulator::explicitEuler(float timeStep)
{
    if (getPositionOfPinball().x==0 && getPositionOfPinball().z <= springVector.at(0).length-0.55f){
        for (Spring& sp : springVector) {
            Vec3 endPositionOfSpring = getPositionOfPinball();
            // calc acc first or will replace position
            Vec3 acceleration = Vec3(0.0f, 0.0f, 0.0f);
            acceleration.z = (-springVector.at(0).stiffness *(endPositionOfSpring.z + 0.55f - springVector.at(0).length)) / pinballVector.at(0).mass;
            cout<<"acc: "<<acceleration<<"\n";
            // x(i+1) = x(i) + v(i) * dt
            setPositionOfPinball(getPositionOfPinball() + getVelocityOfPinball() * timeStep);
            // v(i+1) = v(i) + a(i) * dt
            setVelocityOfPinball(getVelocityOfPinball() + acceleration * timeStep);
            cout<<"velocity: "<<getVelocityOfPinball()<<"\n";
            cout<< "***************************************************** "<<"\n";

	    }
    }else  if (getPositionOfPinball().x == -0.85f && getPositionOfPinball().y <= springVector.at(0).length - 0.55f) {
        for (Spring& sp : springVector) {
            Vec3 endPositionOfSpring = getPositionOfPinball();
            // calc acc first or will replace position
            Vec3 acceleration = Vec3(0.0f, 0.0f, 0.0f);
            acceleration.y = (-springVector.at(0).stiffness * (endPositionOfSpring.y + 0.55f - springVector.at(0).length)) / pinballVector.at(0).mass;
            cout << "acc: " << acceleration << "\n";
            // x(i+1) = x(i) + v(i) * dt
            setPositionOfPinball(getPositionOfPinball() + getVelocityOfPinball() * timeStep);
            // v(i+1) = v(i) + a(i) * dt
            setVelocityOfPinball(getVelocityOfPinball() + acceleration * timeStep);
            cout << "velocity: " << getVelocityOfPinball() << "\n";
            cout << "***************************************************** " << "\n";

        }
    }
	
}

Vec3 PinballSimulator::getPositionOfPinball(){
    return pinballVector.at(0).position;
}

void PinballSimulator::setPositionOfPinball(Vec3 position){
    pinballVector.at(0).position = position;
}

void PinballSimulator::setVelocityOfPinball(Vec3 velocity){
    pinballVector.at(0).linearVelocity = velocity;
}

Vec3 PinballSimulator::getVelocityOfPinball(){
    return pinballVector.at(0).linearVelocity;
}

void PinballSimulator::springForce(float timeStep){
    explicitEuler(timeStep);
}


void PinballSimulator::forceEular(float timeStep){
        for (int i = 0; i < rigidwallVector.size(); i++) {
			CollisionInfo  cosInfo = checkCollisionSAT(getWorldMatrixOfRigidwall(i), getWorldMatrixOfPinball());
			if (cosInfo.isValid) {
                cout<<"detect!! "<<cosInfo.collisionPointWorld <<"\n";
                Vec3 curVelocity = getVelocityOfPinball();
                curVelocity.x = curVelocity.x * (1-damping);
                curVelocity.z = curVelocity.z * (1-damping);
                Vec3 vectorOfRigidWall = Vec3(sin(rigidwallVector.at(i).direction),0, cos(rigidwallVector.at(i).direction));
                if ( vectorOfRigidWall.x*curVelocity.x + vectorOfRigidWall.z*curVelocity.z < 0) vectorOfRigidWall = Vec3(-sin(rigidwallVector.at(i).direction),0, -cos(rigidwallVector.at(i).direction));
                cout<< "V: "<<curVelocity <<"\n";
                cout<< "w: "<<vectorOfRigidWall <<"\n";
                //w*w^t
                std::vector<std::vector<float>> outerProductResult = calculateOuterProduct(vectorOfRigidWall);
                //I-2*w*w^t
                std::vector<std::vector<float>> tempResult(3, std::vector<float>(3, 0.0));
                std::vector<std::vector<float>> I(3, std::vector<float>(3, 0.0));
                for (int i = 0; i < 3; ++i) {
                    I[i][i] = 1.0;
                }
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        tempResult[i][j] = - I[i][j] + 2*outerProductResult[i][j];
                    }
                }
                cout<< "I-2wwt: "<<"\n";
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        std::cout << tempResult[i][j] << " ";
                    }
                    std::cout << std::endl;  
                }

                Vec3 newVelocity = Vec3 (0.0f, 0.0f, 0.0f);
                newVelocity.z = tempResult[2][0]*curVelocity.x+tempResult[2][2]*curVelocity.z;
                newVelocity.x = tempResult[0][0]*curVelocity.x+tempResult[0][2]*curVelocity.z;                    
                   
                cout<< "V_new: "<< newVelocity<< "\n";
                setVelocityOfPinball(newVelocity);
                cout<< "***************************************************** "<<"\n";
                break;
			}
	    }
    
    Vec3 curPosition = getPositionOfPinball();
    curPosition.z = curPosition.z + timeStep*getVelocityOfPinball().z;
    curPosition.x = curPosition.x + timeStep*getVelocityOfPinball().x;
    //cout<< "position: "<<curPosition <<"\n";
    setPositionOfPinball(curPosition);

}

std::vector<std::vector<float>> PinballSimulator::calculateOuterProduct(const Vec3& w) {
    // calculate ww^T
    std::vector<std::vector<float>> result(3, std::vector<float>(3, 0.0));
    result[0][0] = w.x * w.x;
    result[2][2] = w.z * w.z;
    result[0][2] = w.x * w.z;
    result[2][0] = w.z * w.x;
    return result;
}


void PinballSimulator::forceEular2(float timeStep) {
    for (int i = 0; i < rigidwallVector.size(); i++) {
        CollisionInfo  cosInfo = checkCollisionSAT(getWorldMatrixOfRigidwall(i), getWorldMatrixOfPinball());
        if (cosInfo.isValid) {
            cout << "detect!! " << cosInfo.collisionPointWorld << "\n";
            Vec3 curVelocity = getVelocityOfPinball();
            curVelocity.x = curVelocity.x * (1 - damping);
            curVelocity.y = curVelocity.y * (1 - damping);
            Vec3 vectorOfRigidWall = Vec3(sin(rigidwallVector.at(i).direction), cos(rigidwallVector.at(i).direction), 0);
            if (vectorOfRigidWall.x * curVelocity.x + vectorOfRigidWall.y * curVelocity.y < 0) vectorOfRigidWall = Vec3(-sin(rigidwallVector.at(i).direction), -cos(rigidwallVector.at(i).direction), 0);
            cout << "V: " << curVelocity << "\n";
            cout << "w: " << vectorOfRigidWall << "\n";
            //w*w^t
            std::vector<std::vector<float>> outerProductResult = calculateOuterProduct2(vectorOfRigidWall);
            //I-2*w*w^t
            std::vector<std::vector<float>> tempResult(3, std::vector<float>(3, 0.0));
            std::vector<std::vector<float>> I(3, std::vector<float>(3, 0.0));
            for (int i = 0; i < 3; ++i) {
                I[i][i] = 1.0;
            }
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    tempResult[i][j] = -I[i][j] + 2 * outerProductResult[i][j];
                }
            }
            cout << "I-2wwt: " << "\n";
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    std::cout << tempResult[i][j] << " ";
                }
                std::cout << std::endl;
            }

            Vec3 newVelocity = Vec3(0.0f, 0.0f, 0.0f);
            newVelocity.y = tempResult[1][0] * curVelocity.x + tempResult[1][1] * curVelocity.y;
            newVelocity.x = tempResult[0][0] * curVelocity.x + tempResult[0][1] * curVelocity.y;

            cout << "V_new: " << newVelocity << "\n";
            setVelocityOfPinball(newVelocity);
            cout << "***************************************************** " << "\n";
            break;
        }
    }

    Vec3 curPosition = getPositionOfPinball();
    curPosition.y = curPosition.y + timeStep * getVelocityOfPinball().y;
    curPosition.x = curPosition.x + timeStep * getVelocityOfPinball().x;
    //cout<< "position: "<<curPosition <<"\n";
    setPositionOfPinball(curPosition);

}

std::vector<std::vector<float>> PinballSimulator::calculateOuterProduct2(const Vec3& w) {
    // calculate ww^T
    std::vector<std::vector<float>> result(3, std::vector<float>(3, 0.0));
    result[0][0] = w.x * w.x;
    result[1][1] = w.y * w.y;
    result[0][1] = w.x * w.y;
    result[1][0] = w.y * w.x;
    result[2][2] = w.z * w.z;
    return result;
}


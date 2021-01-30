#include <cmath>
#include "aIKController.h"
#include "GL/glut.h"

#include "aActor.h"
#include "aMatrix.h"

#pragma warning (disable : 4018)

int IKController::gIKmaxIterations = 5;
double IKController::gIKEpsilon = 0.1;

// AIKchain class functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////
AIKchain::AIKchain()
{
	mWeight0 = 0.1;
}

AIKchain::~AIKchain()
{

}

AJoint* AIKchain::getJoint(int index) 
{ 
	return mChain[index]; 
}

void AIKchain::setJoint(int index, AJoint* pJoint) 
{ 
	mChain[index] = pJoint; 
}

double AIKchain::getWeight(int index) 
{ 
	return mWeights[index]; 
}

void AIKchain::setWeight(int index, double weight) 
{ 
	mWeights[index] = weight; 
}

int AIKchain::getSize() 
{ 
	return mChain.size(); 
}

std::vector<AJoint*>& AIKchain::getChain() 
{ 
	return mChain; 
}

std::vector<double>& AIKchain::getWeights() 
{ 
	return mWeights; 
}

void AIKchain::setChain(std::vector<AJoint*> chain) 
{
	mChain = chain; 
}

void AIKchain::setWeights(std::vector<double> weights) 
{ 
	mWeights = weights; 
}

// AIKController class functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////

IKController::IKController()
{
	m_pActor = NULL;
	m_pSkeleton = NULL;
	mvalidLimbIKchains = false;
	mvalidCCDIKchains = false;

	// Limb IK
	m_pEndJoint = NULL;
	m_pMiddleJoint = NULL;
	m_pBaseJoint = NULL;
	m_rotationAxis = vec3(0.0, 1.0, 0.0);

	ATransform desiredTarget = ATransform();
	mTarget0.setLocal2Parent(desiredTarget);  // target associated with end joint
	mTarget1.setLocal2Parent(desiredTarget);  // optional target associated with middle joint - used to specify rotation of middle joint about end/base axis
	mTarget0.setLocal2Global(desiredTarget);
	mTarget1.setLocal2Global(desiredTarget);

	//CCD IK
	mWeight0 = 0.1;  // default joint rotation weight value

}

IKController::~IKController()
{
}

ASkeleton* IKController::getSkeleton()
{
	return m_pSkeleton;
}

const ASkeleton* IKController::getSkeleton() const
{
	return m_pSkeleton;
}

ASkeleton* IKController::getIKSkeleton()
{
	return &mIKSkeleton;
}

const ASkeleton* IKController::getIKSkeleton() const
{
	return &mIKSkeleton;
}

AActor* IKController::getActor()
{
	return m_pActor;
}

void IKController::setActor(AActor* actor)

{
	m_pActor = actor;
	m_pSkeleton = m_pActor->getSkeleton();
}


AIKchain IKController::createIKchain(int endJointID, int desiredChainSize, ASkeleton* pSkeleton)
{
	// TODO: given the end joint ID and the desired size (i.e. length) of the IK chain, 
	// 1. add the corresponding skeleton joint pointers to the AIKChain "chain" vector data member starting with the end joint
	// 2. also add weight values to the associated AIKChain "weights" vector data member for use in the CCD IK implemention
	// Note: desiredChainSize = -1 should create an IK chain of maximum length (i.e. where the last chain joint is the joint before the root joint)
	bool getMaxSize = false;

	int EndJointID = endJointID;
	std::vector<AJoint*> chain;
	std::vector<double> weights;

	chain.clear();
	weights.clear();
	if (desiredChainSize == -1)
		getMaxSize = true;

	if ((EndJointID >= 0) && (EndJointID < pSkeleton->getNumJoints()))
	{
		AJoint* pJoint = pSkeleton->getJointByID(endJointID);

		// TODO: add code here to generate chain of desired size or terminate at the joint before root joint, so that root will not change during IK	
		// also add weight values to corresponding weights vector  (default value = 0.1)
		if (getMaxSize) 
		{
			while (pJoint->getParent()->getParent() != NULL)
			{
				chain.push_back(pJoint);
				pJoint = pJoint->getParent();
				weights.push_back(0.1f);
			}
		}
		else 
		{
			for (int i = 0; i < desiredChainSize; i++) 
			{
				chain.push_back(pJoint);
				pJoint = pJoint->getParent();
				weights.push_back(0.1f);
				if (pJoint == NULL) 
				{
					break;
				}
			}
		}
		
	}
	AIKchain result;
	result.setChain(chain);
	result.setWeights(weights);

	return result;
}



bool IKController::IKSolver_Limb(int endJointID, const ATarget& target)
{
	// Implements the analytic/geometric IK method assuming a three joint limb  

	if (!mvalidLimbIKchains)
	{
		mvalidLimbIKchains = createLimbIKchains();
		//assert(mvalidLimbIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;

	switch (endJointID)
	{
	case mLhandID:
		mLhandTarget = target;
		computeLimbIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
		break;
	case mRhandID:
		mRhandTarget = target;
		computeLimbIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
		break;
	case mLfootID:
		mLfootTarget = target;
		computeLimbIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
		break;
	case mRfootID:
		mRfootTarget = target;
		computeLimbIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
		break;
	case mRootID:
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeLimbIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
		computeLimbIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
		computeLimbIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
		computeLimbIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
		break;
	default:
		mIKchain = createIKchain(endJointID, 3, &mIKSkeleton);
		computeLimbIK(target, mIKchain, axisY, &mIKSkeleton);
		break;
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}



int IKController::createLimbIKchains()
{
	bool validChains = false;
	int desiredChainSize = 3;

	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);
	
	if (mLhandIKchain.getSize() == 3 && mRhandIKchain.getSize() == 3 && mLfootIKchain.getSize() == 3 && mRfootIKchain.getSize() == 3)
	{
		validChains = true;
		
		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}




int IKController::computeLimbIK(ATarget target, AIKchain& IKchain, const vec3 midJointAxis, ASkeleton* pIKSkeleton)
{
	// TODO: Implement the analytic/geometric IK method assuming a three joint limb  
	// The actual position of the end joint should match the target position within some episilon error 
	// the variable "midJointAxis" contains the rotation axis for the middle joint
	
	bool result = false;
	int endJointID;
	mTarget0 = target;

	if (IKchain.getSize() > 0)
		 endJointID = IKchain.getJoint(0)->getID();
	else endJointID = -1;

	int rotIndex = -1;
	int orientIndex = -1;
	int lastIndex = -1;
	vec3 orientAxis;
	int orientFlag = 0;
	vec3 rotOrient;
	int rotFlag = 0;


	if (!(fabs(midJointAxis[0] - 0.0) <= DBL_EPSILON))
	{
		rotIndex = 0;
	}

	if (!(fabs(midJointAxis[1] - 0.0) <= DBL_EPSILON))
	{
		rotIndex = 1;
	}

	if (!(fabs(midJointAxis[2] - 0.0) <= DBL_EPSILON))
	{
		rotIndex = 2;
	}



	if ((endJointID >= 0) && (endJointID < pIKSkeleton->getNumJoints()))
	{
		m_pEndJoint = IKchain.getJoint(0);
		m_pMiddleJoint = IKchain.getJoint(1);
		m_pBaseJoint = IKchain.getJoint(2);

		vec3 rDesire = target.getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
		vec3 orient = m_pBaseJoint->getLocalTranslation() - m_pMiddleJoint->getLocalTranslation();
		double orientEntry = fmax(abs(orient[0]), fmax(abs(orient[1]), abs(orient[2])));
		if (abs(abs(orient[0]) - orientEntry) <= DBL_EPSILON) 
		{
			orientIndex = 0;
			orientAxis = vec3(orient[0], 0.0, 0.0).Normalize();
			

			if (fabs(orient[0] - 0.0) > DBL_EPSILON && orient[0] > 0.0) 
			{
				orientFlag = 1;
			}
			else if (fabs(orient[0] - 0.0) > DBL_EPSILON&& orient[0] < 0.0) 
			{
				orientFlag = -1;
			}
			if (rotIndex == 1) 
			{
				lastIndex = 2;
				rotOrient = vec3(rDesire[0], 0.0, 0.0).Cross(vec3(0.0, 0.0, rDesire[2])).Normalize();
			}
			else 
			{
				lastIndex = 1;
				rotOrient = vec3(rDesire[0], 0.0, 0.0).Cross(vec3(0.0, rDesire[1], 0.0)).Normalize();
			}
		}

		if (abs(abs(orient[1]) - orientEntry) <= DBL_EPSILON)
		{
			orientIndex = 1;
			orientAxis = vec3(0.0, orient[1], 0.0).Normalize();
			if (fabs(orient[1] - 0.0) > DBL_EPSILON&& orient[1] > 0.0)
			{
				orientFlag = 1;
			}
			else if (fabs(orient[1] - 0.0) > DBL_EPSILON&& orient[1] < 0.0)
			{
				orientFlag = -1;
			}
			if (rotIndex == 0)
			{
				
				lastIndex = 2;
				rotOrient = vec3(0.0, rDesire[1], 0.0).Cross(vec3(0.0, 0.0, rDesire[2])).Normalize();
			}
			else
			{
				lastIndex = 0;
				rotOrient = vec3(rDesire[0], 0.0, 0.0).Cross(vec3(0.0, rDesire[1], 0.0)).Normalize();
			}
		}

		if (abs(abs(orient[2]) - orientEntry) <= DBL_EPSILON)
		{
			orientIndex = 2;
			orientAxis = vec3(0.0, 0.0, orient[2]).Normalize();
			if (fabs(orient[2] - 0.0) > DBL_EPSILON && orient[2] > 0.0)
			{
				orientFlag = 1;
			}
			else if (fabs(orient[2] - 0.0) > DBL_EPSILON && orient[2] < 0.0)
			{
				orientFlag = -1;
			}
			if (rotIndex == 1)
			{
				lastIndex = 0;
				rotOrient = vec3(rDesire[0], 0.0, 0.0).Cross(vec3(0.0, 0.0, rDesire[2])).Normalize();
			}
			else
			{
				lastIndex = 1;
				rotOrient = vec3(0.0, rDesire[1], 0.0).Cross(vec3(0.0, 0.0, rDesire[2])).Normalize();
			}
		}

		if ((fabs(rotOrient[0] - 0.0) > DBL_EPSILON&& rotOrient[0] > 0.0)
			|| (fabs(rotOrient[1] - 0.0) > DBL_EPSILON&& rotOrient[1] > 0.0)
			|| (fabs(rotOrient[2] - 0.0) > DBL_EPSILON&& rotOrient[2] > 0.0)) 
		{
			rotFlag = 1;
		}

		if ((fabs(rotOrient[0] - 0.0) > DBL_EPSILON&& rotOrient[0] < 0.0)
			|| (fabs(rotOrient[1] - 0.0) > DBL_EPSILON&& rotOrient[1] < 0.0)
			|| (fabs(rotOrient[2] - 0.0) > DBL_EPSILON&& rotOrient[2] < 0.0))
		{
			rotFlag = -1;
		}

		//TODO:
		// 1. compute error vector between target and end joint
		vec3 error = target.getGlobalTranslation() - m_pEndJoint->getGlobalTranslation();
		// 2. compute vector between end Joint and base joint
		vec3 r02 = m_pEndJoint->getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
		// 3. compute vector between target and base joint
		
		// 4. Compute desired angle for middle joint 
		vec3 l1 = m_pMiddleJoint->getLocalTranslation();
		vec3 l2 = m_pEndJoint->getLocalTranslation();
		float cosFi = (pow(l1.Length(), 2) + pow(l2.Length(), 2) - pow(rDesire.Length(), 2))
						/ (2 * l1.Length() * l2.Length());
		float fi = M_PI - acos(cosFi);
		if (isnan(fi)) 
		{
			return result;
		}
		// 5. given desired angle and midJointAxis, compute new local middle joint rotation matrix and update joint transform
		quat middleRotQ = quat();
		middleRotQ.FromAxisAngle(midJointAxis, fi);
		mat3 middleRotMat = middleRotQ.ToRotation();
		m_pMiddleJoint->setLocalRotation(middleRotMat);
		// 6. compute vector between target and base joint

		// 7. Compute base joint rotation axis (in global coords) and desired angle
		float thetaBase = asin(l2.Length() * sin(fi) / rDesire.Length());
		float thetaGamma = atan(rDesire[rotIndex] / rDesire[orientIndex]);
		float thetaBeta = atan(fabs(rDesire[lastIndex]) * double(rotFlag) * orientFlag / fabs(rDesire[orientIndex]));
		// 8. transform base joint rotation axis to local coordinates
		quat thetaBaseQ = quat();
		thetaBaseQ.FromAxisAngle(midJointAxis, -thetaBase);
		mat3 rotMatBase = thetaBaseQ.ToRotation();
		vec3 gammaAxis = midJointAxis.Cross(orientAxis);
		quat thetaGammaQ = quat();
		thetaGammaQ.FromAxisAngle(gammaAxis, -thetaGamma);
		mat3 rotMatGamma = thetaGammaQ.ToRotation(); 
		quat thetaBetaQ = quat();
		thetaBetaQ.FromAxisAngle(midJointAxis, thetaBeta);
		mat3 rotMatBeta = thetaBetaQ.ToRotation();
	

		// 9. given desired angle and local rotation axis, compute new local rotation matrix and update base joint transform
		m_pBaseJoint->setLocalRotation(rotMatBeta * rotMatGamma * rotMatBase);
		m_pBaseJoint->updateTransform();
		result = true;
	}
	return result;

}

bool IKController::IKSolver_CCD(int endJointID, const ATarget& target)
{
	// Implements the CCD IK method assuming a three joint limb 

	bool validChains = false;

	if (!mvalidCCDIKchains)
	{
		mvalidCCDIKchains = createCCDIKchains();
		//assert(mvalidCCDIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;

	switch (endJointID)
	{
	case mLhandID:
		mLhandTarget = target;
		computeCCDIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		break;
	case mRhandID:
		mRhandTarget = target;
		computeCCDIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		break;
	case mLfootID:
		mLfootTarget = target;
		computeCCDIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		break;
	case mRfootID:
		mRfootTarget = target;
		computeCCDIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
		break;
	case mRootID:
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeCCDIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		computeCCDIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		computeCCDIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		computeCCDIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
		break;
	default:
		mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
		computeCCDIK(target, mIKchain, &mIKSkeleton);
		break;
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}

int IKController::createCCDIKchains()
{
	bool validChains = false;

	int desiredChainSize = -1;  // default of -1 creates IK chain of maximum length from end joint to child joint of root


	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);

	if (mLhandIKchain.getSize() > 1 && mRhandIKchain.getSize() > 1 && mLfootIKchain.getSize() > 1 && mRfootIKchain.getSize() > 1)
	{
		validChains = true;

		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}


int IKController::computeCCDIK(ATarget target, AIKchain& IKchain, ASkeleton* pIKSkeleton)
{

	// TODO: Implement CCD IK  
	// The actual position of the end joint should match the desiredEndPos within some episilon error 

	bool result = false;

	mTarget0 = target;
	vec3 desiredEndPos = mTarget0.getGlobalTranslation();  // Get desired position of EndJoint

	int chainSize = IKchain.getSize();
	if (chainSize == 0) // There are no joints in the IK chain for manipulation
		return false;

	double epsilon = gIKEpsilon;
	int maxIterations = gIKmaxIterations;
	int numIterations = 0;

	m_pEndJoint = IKchain.getJoint(0);
	int endJointID = m_pEndJoint->getID();
	m_pBaseJoint = IKchain.getJoint(chainSize - 1);

	pIKSkeleton->copyTransforms(m_pSkeleton);

	if ((endJointID >= 0) && (endJointID < pIKSkeleton->getNumJoints()))
	{
		//TODO:
		vector<vec3> axes;
		vector<float> thetas;
		// 1. compute axis and angle for each joint in the IK chain (distal to proximal) in global coordinates
		for (int j = 0; j < 4; j++) 
		{
			for (int i = 1; i < chainSize; i++)
			{
				vec3 error = target.getGlobalTranslation() - m_pEndJoint->getGlobalTranslation();
				vec3 rDesire = target.getGlobalTranslation() - IKchain.getJoint(i)->getGlobalTranslation();
				vec3 rEndJoint = m_pEndJoint->getGlobalTranslation() - IKchain.getJoint(i)->getGlobalTranslation();
				vec3 crossAxis = rEndJoint.Cross(error).Length();
				float theta = rEndJoint.Cross(error).Length() / (rEndJoint * rEndJoint + rEndJoint * error);
				vec3 axis = rEndJoint.Cross(error).Normalize();

				// 2. once you have the desired axis and angle, convert axis to local joint coords 
				axis = IKchain.getJoint(i)->getGlobalRotation().Inverse() * axis;

				// 3. multiply angle by corresponding joint weight value
				theta = IKchain.getWeight(i) * theta;

				// 4. compute new local joint rotation matrix
				quat thetaQ = quat();
				thetaQ.FromAxisAngle(axis, theta);
				mat3 thetaRotMat = thetaQ.ToRotation();

				// 5. update joint transform
				IKchain.getJoint(i)->setLocalRotation(IKchain.getJoint(i)->getLocalRotation() * thetaRotMat);
				IKchain.getJoint(i)->updateTransform();
			}
		}
		
		// 6. repeat same operations above for each joint in the IKchain from end to base joint

	}
	return result;

}


bool IKController::IKSolver_PseudoInv(int endJointID, const ATarget& target)
{
	bool result = false;
	mIKSkeleton.copyTransforms(m_pSkeleton);
	int desiredChainSize = -1;

	// TODO: Implement Pseudo Inverse-based IK  
	// The actual position of the end joint should match the target position after the skeleton is updated with the new joint angles
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, 3, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, 3, &mIKSkeleton);

	switch (endJointID)
	{
	case mLhandID:
		mLhandTarget = target;
		mIKchain = mLhandIKchain;
		break;
	case mRhandID:
		mRhandTarget = target;
		mIKchain = mRhandIKchain;
		break;
	case mLfootID:
		mLfootTarget = target;
		mIKchain = mLfootIKchain;
		break;
	case mRfootID:
		mRfootTarget = target;
		mIKchain = mRfootIKchain;
		break;
	default:
		mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
		break;
	}

	//mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
	AJoint* endJoint = mIKchain.getJoint(0);
	vector<mat3> jEntries;
	matrix<double> thetaOrigin = matrix<double>(3 * (mIKchain.getSize() - 1), 1);

	vec3 desireDifference = target.getGlobalTranslation() - endJoint->getGlobalTranslation();
	matrix<double> desireDifferenceM = matrix<double>(3, 1);
	desireDifferenceM(0, 0) = desireDifference[0];
	desireDifferenceM(1, 0) = desireDifference[1];
	desireDifferenceM(2, 0) = desireDifference[2];

	for (int i = 1; i < mIKchain.getSize(); i++) 
	{
		vec3 radius = endJoint->getGlobalTranslation() - mIKchain.getJoint(i)->getGlobalTranslation();
		mat3 rotation = mIKchain.getJoint(i)->getGlobalRotation();
		//Set all rotation order Z->Y->X
		vec3 angularRad = vec3();
		rotation.ToEulerAngles(mat3::ZYX, angularRad);
		thetaOrigin(3 * (i - 1), 0) = angularRad[0];
		thetaOrigin(3 * (i - 1) + 1, 0) = angularRad[1];
		thetaOrigin(3 * (i - 1) + 2, 0) = angularRad[2];
		mat3 lMatrix = mat3(
			vec3(1.0, 0.0, -sin(angularRad[1])),
			vec3(0.0, cos(angularRad[0]), sin(angularRad[0]) * cos(angularRad[1])),
			vec3(0.0, -sin(angularRad[0]), cos(angularRad[0]) * cos(angularRad[1]))
		);
		vec3 ax = rotation.GetCol(0);
		vec3 ay = rotation.GetCol(1);
		vec3 az = rotation.GetCol(2);

		vec3 b1 = ax.Cross(radius);
		vec3 b2 = ay.Cross(radius);
		vec3 b3 = az.Cross(radius);

		mat3 bMatrix = mat3(
			vec3(b1[0], b2[0], b3[0]),
			vec3(b1[1], b2[1], b3[1]),
			vec3(b1[2], b2[2], b3[2])
		);

		mat3 jEntry = bMatrix * lMatrix;
		jEntries.push_back(jEntry);
	}
	matrix<double> jMatrix = matrix<double>(3, 3 * (mIKchain.getSize() - 1));
	for (int i = 0; i < mIKchain.getSize() - 1; i++) 
	{
		jMatrix(0, 3 * i) = jEntries[i][0][0];
		jMatrix(0, 3 * i + 1) = jEntries[i][0][1];
		jMatrix(0, 3 * i + 2) = jEntries[i][0][2];
		jMatrix(1, 3 * i) = jEntries[i][1][0];
		jMatrix(1, 3 * i + 1) = jEntries[i][1][1];
		jMatrix(1, 3 * i + 2) = jEntries[i][1][2];
		jMatrix(2, 3 * i) = jEntries[i][2][0];
		jMatrix(2, 3 * i + 1) = jEntries[i][2][1];
		jMatrix(2, 3 * i + 2) = jEntries[i][2][2];
	}
	matrix<double> jMatrixPlus;
	
	
	int col = 3 * (mIKchain.getSize() - 1);

	// Pseudo inverse
	if (col < 3) 
	{
		// Left inverse
		matrix<double> singualarIncase = matrix<double>(3 * (mIKchain.getSize() - 1), 3 * (mIKchain.getSize() - 1));
		singualarIncase.Unit();
		jMatrixPlus = ((~jMatrix * jMatrix + 0.001 * singualarIncase).Inv()) * (~jMatrix);
	}
	else if (col > 3) 
	{
		// Right inverse
		matrix<double> singualarIncase = matrix<double>(3, 3);
		singualarIncase.Unit();
		jMatrixPlus = (~jMatrix) * ((jMatrix * ~jMatrix + 0.00001 * singualarIncase).Inv());
	}
	else 
	{
		matrix<double> singualarIncase = matrix<double>(3, 3);
		singualarIncase.Unit();
		jMatrixPlus = (jMatrix + 0.00001 * singualarIncase).Inv();
	}

	matrix<double> targetRot = thetaOrigin + jMatrixPlus * desireDifferenceM;

	for (int i = 1; i < mIKchain.getSize(); i++) 
	{
		vec3 eulerAngle = vec3();
		eulerAngle[0] = targetRot(3 * (i - 1), 0);
		eulerAngle[1] = targetRot(3 * (i - 1) + 1, 0);
		eulerAngle[2] = targetRot(3 * (i - 1) + 2, 0);

		mat3 globalRot = mat3();
		globalRot.FromEulerAngles(mat3::ZYX, eulerAngle);
		globalRot = mIKchain.getJoint(i)->getParent()->getGlobalRotation().Inverse() * globalRot;
		mIKchain.getJoint(i)->setLocalRotation(globalRot);
	}
	mIKSkeleton.update();
	m_pSkeleton->copyTransforms(&mIKSkeleton);
	result = true;


	return result;
}

bool IKController::IKSolver_Other(int endJointID, const ATarget& target)
{
	
	bool result = false;
	
	// TODO: Put Optional IK implementation or enhancements here
	 
	return result;
}
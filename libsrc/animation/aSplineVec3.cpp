#include "aSplineVec3.h"
#include <algorithm>
#include <Eigen/Dense>
#include <vector>

#pragma warning(disable:4018)
#pragma warning(disable:4244)

ASplineVec3::ASplineVec3() : mInterpolator(new ALinearInterpolatorVec3())
{
}

ASplineVec3::~ASplineVec3()
{
  delete mInterpolator;
}

void ASplineVec3::setFramerate(double fps)
{
  mInterpolator->setFramerate(fps);
}

double ASplineVec3::getFramerate() const
{
  return mInterpolator->getFramerate();
}


void ASplineVec3::setLooping(bool loop)
{
  mLooping = loop;
}

bool ASplineVec3::getLooping() const
{
  return mLooping;
}

void ASplineVec3::setInterpolationType(ASplineVec3::InterpolationType type)
{
  double fps = getFramerate();

  delete mInterpolator;
  switch (type)
  {
  case LINEAR: mInterpolator = new ALinearInterpolatorVec3();
    break;
  case CUBIC_BERNSTEIN: mInterpolator = new ABernsteinInterpolatorVec3();
    break;
  case CUBIC_CASTELJAU: mInterpolator = new ACasteljauInterpolatorVec3();
    break;
  case CUBIC_MATRIX: mInterpolator = new AMatrixInterpolatorVec3();
    break;
  case CUBIC_HERMITE: mInterpolator = new AHermiteInterpolatorVec3();
    break;
  case CUBIC_BSPLINE: mInterpolator = new ABSplineInterpolatorVec3();
    break;
  case STAIR_CASE: mInterpolator = new AStaircaseInterpolatorVec3();
	break;
  };

  mInterpolator->setFramerate(fps);
  computeControlPoints();
  cacheCurve();
}

ASplineVec3::InterpolationType ASplineVec3::getInterpolationType() const
{
  return mInterpolator->getType();
}

void ASplineVec3::editKey(int keyID, const vec3& value)
{
  assert(keyID >= 0 && keyID < mKeys.size());
  mKeys[keyID].second = value;
  computeControlPoints();
  cacheCurve();
}

void ASplineVec3::editControlPoint(int ID, const vec3& value)
{
  assert(ID >= 0 && ID < mCtrlPoints.size()+2);
  if (ID == 0)
  {
    mStartPoint = value;
    computeControlPoints();
  }
  else if (ID == mCtrlPoints.size() + 1)
  {
    mEndPoint = value;
    computeControlPoints();
  }
  else mCtrlPoints[ID - 1] = value;
  cacheCurve();
}

void ASplineVec3::appendKey(double time, const vec3& value, bool updateCurve)
{
  mKeys.push_back(Key(time, value));

  if (mKeys.size() >= 2)
  {
    int totalPoints = mKeys.size();

    //If there are more than 1 interpolation point, set up the 2 end points to help determine the curve.
    //They lie on the tangent of the first and last interpolation points.
    vec3 tmp = mKeys[0].second - mKeys[1].second;
    double n = tmp.Length();
    mStartPoint = mKeys[0].second + (tmp / n) * n * 0.25;
    // distance to endpoint is 25% of distance between first 2 points

    tmp = mKeys[totalPoints - 1].second - mKeys[totalPoints - 2].second;
    n = tmp.Length();
    mEndPoint = mKeys[totalPoints - 1].second + (tmp / n) * n * 0.25;
  }

  if (updateCurve)
  {
    computeControlPoints();
    cacheCurve();
  }
}

void ASplineVec3::appendKey(const vec3& value, bool updateCurve)
{
  if (mKeys.size() == 0)
  {
    appendKey(0, value, updateCurve);
  }
  else
  {
    double lastT = mKeys[mKeys.size() - 1].first;
    appendKey(lastT + 1, value, updateCurve);
  }
}

void ASplineVec3::deleteKey(int keyID)
{
  assert(keyID >= 0 && keyID < mKeys.size());

  //Delete debug
  for (int i = keyID + 1; i < mKeys.size(); i++) 
  {
	  mKeys[i].first -= 1.0;
  }

  mKeys.erase(mKeys.begin() + keyID);
  computeControlPoints();
  cacheCurve();
}

vec3 ASplineVec3::getKey(int keyID)
{
  assert(keyID >= 0 && keyID < mKeys.size());
  return mKeys[keyID].second;
}

int ASplineVec3::getNumKeys() const
{
  return mKeys.size();
}

vec3 ASplineVec3::getControlPoint(int ID)
{
  assert(ID >= 0 && ID < mCtrlPoints.size()+2);
  if (ID == 0) return mStartPoint;
  else if (ID == mCtrlPoints.size() + 1) return mEndPoint;
  else return mCtrlPoints[ID - 1];
}

int ASplineVec3::getNumControlPoints() const
{
  return mCtrlPoints.size() + 2; // include endpoints
}

void ASplineVec3::clear()
{
  mKeys.clear();
}

double ASplineVec3::getDuration() const
{
  return mKeys[mKeys.size() - 1].first;
}

double ASplineVec3::getNormalizedTime(double t) const
{
  return (t / getDuration());
}

vec3 ASplineVec3::getValue(double t)
{
  if (mCachedCurve.size() == 0) return vec3();

  double dt = mInterpolator->getDeltaTime();
  int rawi = (int)(t / dt); // assumes uniform spacing
  int i = rawi % mCachedCurve.size();
  double frac = t - rawi * dt;
  int inext = i + 1;
  if (!mLooping) inext = std::min<int>(inext, mCachedCurve.size() - 1);
  else inext = inext % mCachedCurve.size();

  vec3 v1 = mCachedCurve[i];
  vec3 v2 = mCachedCurve[inext];
  vec3 v = v1 * (1 - frac) + v2 * frac;
  return v;
}

void ASplineVec3::cacheCurve()
{
  mInterpolator->interpolate(mKeys, mCtrlPoints, mCachedCurve);
}

void ASplineVec3::computeControlPoints()
{
  mInterpolator->computeControlPoints(mKeys, mCtrlPoints, mStartPoint, mEndPoint);
}

int ASplineVec3::getNumCurveSegments() const
{
  return mCachedCurve.size();
}

vec3 ASplineVec3::getCurvePoint(int i) const
{
  return mCachedCurve[i];
}

vec3 ASplineVec3::getCachedPoint(double t) 
{
	if (mCachedCurve.size() == 0) return vec3();
	double mDt = mInterpolator->getDeltaTime();
	int numFrames = (int)(t / mDt);
	int size = mCachedCurve.size();
	int i = numFrames % size;

	return mCachedCurve[i];
}

//---------------------------------------------------------------------
AInterpolatorVec3::AInterpolatorVec3(ASplineVec3::InterpolationType t) : mDt(1.0 / 120.0), mType(t)
{
}

void AInterpolatorVec3::setFramerate(double fps)
{
  mDt = 1.0 / fps;
}

double AInterpolatorVec3::getFramerate() const
{
  return 1.0 / mDt;
}

double AInterpolatorVec3::getDeltaTime() const
{
  return mDt;
}

// Compare two floats
int floatCompare(double x, double y)
{
	if (x > y && (x - y) > DBL_EPSILON)
	{
		return 1;
	}
	else if (x < y && (y - x) > DBL_EPSILON)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

void AInterpolatorVec3::interpolate(const std::vector<ASplineVec3::Key>& keys,
                                    const std::vector<vec3>& ctrlPoints, std::vector<vec3>& curve)
{
  vec3 val = 0.0;
  double u = 0.0;

  curve.clear();

  int numSegments = keys.size() - 1;
  for (int segment = 0; segment < numSegments; segment++)
  {
    for (double t = keys[segment].first; t < keys[segment + 1].first - FLT_EPSILON; t += mDt)
    {
      // TODO: Compute u, fraction of duration between segment and segmentnext, for example,
      // u = 0.0 when t = keys[segment-1].first  
      // u = 1.0 when t = keys[segment].first
	  u += mDt;
      val = interpolateSegment(keys, ctrlPoints, segment, u);
      curve.push_back(val);
    }
	u = 0.0f;
  }
  // add last point
  if (keys.size() > 1)
  {
    u = 1.0;
    val = interpolateSegment(keys, ctrlPoints, numSegments - 1, u);
    curve.push_back(val);
  }
}


vec3 ALinearInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double u)
{
  vec3 curveValue(0, 0, 0);
  vec3 key0 = keys[segment].second;
  vec3 key1 = keys[segment + 1].second;

  // TODO: 
  //Step 1: Create a Lerp helper function
  //Step 2: Linear interpolate between key0 and key1 so that u = 0 returns key0 and u = 1 returns key1

  // Linear interpolate
  curveValue = (1 - u) * key0 + u * key1;

  return curveValue;
}

vec3 ABernsteinInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double t)
{
  vec3 b0;
  vec3 b1;
  vec3 b2;
  vec3 b3;
  vec3 curveValue(0, 0, 0);
  // TODO: 
  // Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
  // Step2: Compute the interpolated value f(u) point using  Bernstein polynomials
 

  b0 = ctrlPoints[4 * segment];
  b1 = ctrlPoints[4 * segment + 1];
  b2 = ctrlPoints[4 * segment + 2];
  b3 = ctrlPoints[4 * segment + 3];
  curveValue = pow(1 - t, 3) * b0 + 3 * t * pow(1 - t, 2) * b1 + 3 * pow(t, 2) * (1 - t) * b2 + pow(t, 3) * b3;
  
  return curveValue;
}


vec3 AStaircaseInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double t)
{
	vec3 b0;
	vec3 center;
	vec3 b1;
	
	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  Bernstein polynomials
	b0 = keys[segment].second;
	b1 = keys[segment + 1].second;
	center = ctrlPoints[segment];

	double x = (1 - t) * b0[0] + t * b1[0];

	if (floatCompare(t, 0.5) == -1) {
		return vec3(x, b0[1], 0.0);
	}
	else if (floatCompare(t, 0.5) == 1) {
		return vec3(x, b1[1], 0.0);
	}
	else {
		return center;
	}
}

vec3 ACasteljauInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double t)
{
  vec3 b0;
  vec3 b1;
  vec3 b2;
  vec3 b3;
  vec3 curveValue(0, 0, 0);

  // TODO: 
  // Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
  // Step2: Compute the interpolated value f(u) point using  deCsteljau alogithm
  b0 = ctrlPoints[4 * segment];
  b1 = ctrlPoints[4 * segment + 1];
  b2 = ctrlPoints[4 * segment + 2];
  b3 = ctrlPoints[4 * segment + 3];

  vec3 b01 = (1 - t) * b0 + t * b1;
  vec3 b11 = (1 - t) * b1 + t * b2;
  vec3 b21 = (1 - t) * b2 + t * b3;

  vec3 b02 = (1 - t) * b01 + t * b11;
  vec3 b12 = (1 - t) * b11 + t * b21;

  curveValue = (1 - t) * b02 + t * b12;

  return curveValue;
}

vec3 AMatrixInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double t)
{
  vec3 b0;
  vec3 b1;
  vec3 b2;
  vec3 b3;
  vec3 curveValue(0, 0, 0);

  // TODO: 
  // Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
  // Step2: Compute the interpolated value f(u) point using  matrix method f(u) = GMU
  // Hint: Using Eigen::MatrixXd data representations for a matrix operations
  b0 = ctrlPoints[4 * segment];
  b1 = ctrlPoints[4 * segment + 1];
  b2 = ctrlPoints[4 * segment + 2];
  b3 = ctrlPoints[4 * segment + 3];
  
  Eigen::Matrix4d Fmon;
  Eigen::Matrix4d Mmon;
  Fmon << b0.n[0], b1.n[0], b2.n[0], b3.n[0],
	  b0.n[1], b1.n[1], b2.n[1], b3.n[1],
	  b0.n[2], b1.n[2], b2.n[2], b3.n[2],
	  1.0f, 1.0f, 1.0f, 1.0f;
 
  Mmon << 1.0f, -3.0f, 3.0f, -1.0f,
	  0.0f, 3.0f, -6.0f, 3.0f,
	  0.0f, 0.0f, 3.0f, -3.0f,
	  0.0f, 0.0f, 0.0f, 1.0f;

  Eigen::Vector4d u;
  u << 1.0f, t, pow(t, 2), pow(t, 3);
  Eigen::Vector4d re = Fmon * Mmon * u;
  curveValue = vec3(re(0), re(1), re(3));
  return curveValue;
}

vec3 AHermiteInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double t)
{
  vec3 p0 = keys[segment].second;
  vec3 p1 = keys[segment + 1].second;
  vec3 q0 = ctrlPoints[segment]; // slope at p0
  vec3 q1 = ctrlPoints[segment + 1]; // slope at p1
  vec3 curveValue(0, 0, 0);

  // TODO: Compute the interpolated value h(u) using a cubic Hermite polynomial  
  curveValue = (2 * pow(t, 3) - 3 * pow(t, 2) + 1) * p0 
	  + (3 * pow(t, 2) - 2 * pow(t, 3)) * p1 
	  + (pow(t, 3) - 2 * pow(t, 2) + t) * q0 
	  + (pow(t, 3) - pow(t, 2)) * q1;

  return curveValue;
}

double N(std::vector<int>& knots, int n, int j, double t) {
	if (n == 0) {
		if (t >= knots[j] && t < knots[j+1])
			return 1.0f;
		else
			return 0.0f;
	}

	double c1 = (t - knots[j]) / (knots[j + n] - knots[j]);
	double c2 = (knots[j + n + 1] - t) / (knots[j + n + 1] - knots[j + 1]);

	return c1 * N(knots, n - 1, j, t) + c2 * N(knots, n - 1, j + 1, t);

}

vec3 ABSplineInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double t)
{
  vec3 curveValue(0, 0, 0);

  // Hint: Create a recursive helper function N(knots,n,j,t) to calculate BSpline basis function values at t, where
  //     knots = knot array
  //	   n = degree of the spline curves (n =3 for cubic)
  //     j = curve interval on knot vector in which to interpolate
  //     t = time value	

  // Step 1: determine the index j
  // Step 2: compute the n nonzero Bspline Basis functions N given j
  // Step 3: get the corresponding control points from the ctrlPoints vector
  // Step 4: compute the Bspline curveValue at time t
  int degree = 3;
  int j = segment + degree;
  std::vector<int> knots;

  for (int i = 0; i < keys.size() + 6; i++) 
  {
	  knots.push_back(i - 3);
  }
  
  double nJn3 = N(knots, degree, j - 3, t + segment);
  double nJn2 = N(knots, degree, j - 2, t + segment);
  double nJn1 = N(knots, degree, j - 1, t + segment);
  double nJn0 = N(knots, degree, j, t + segment);

  curveValue = nJn3 * ctrlPoints[j - 3] + nJn2 * ctrlPoints[j - 2] + nJn1 * ctrlPoints[j - 1] + nJn0 * ctrlPoints[j];
 
  return curveValue;
}

void ACubicInterpolatorVec3::computeControlPoints(
  const std::vector<ASplineVec3::Key>& keys,
  std::vector<vec3>& ctrlPoints,
  vec3& startPoint, vec3& endPoint)
{
  ctrlPoints.clear();
  if (keys.size() <= 1) return;

  for (int i = 1; i < keys.size(); i++)
  {
    vec3 b0, b1, b2, b3;
    // TODO: compute b0, b1, b2, b3
	b0 = keys[i-1].second;
	if (i == 1) {
		vec3 keyM1 = keys[i - 1].second - (keys[i].second - keys[i - 1].second);
		b1 = keys[i - 1].second + (keys[i].second - keyM1) / 6;
	}
	else 
	{
		b1 = keys[i - 1].second + (keys[i].second - keys[i - 2].second) / 6;
	}
	if (i == keys.size() - 1) {
		vec3 keyP1 = keys[i].second + (keys[i].second - keys[i - 1].second);
		b2 = keys[i].second - (keyP1 - keys[i - 1].second) / 6;
	}
	else 
	{
		b2 = keys[i].second - (keys[i + 1].second - keys[i - 1].second) / 6;
	}
	b3 = keys[i].second;

    ctrlPoints.push_back(b0);
    ctrlPoints.push_back(b1);
    ctrlPoints.push_back(b2);
    ctrlPoints.push_back(b3);
  }
}

void AStaircaseInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints,
	vec3& startPoint, vec3& endPoint)
{
	ctrlPoints.clear();
	if (keys.size() <= 1) return;

	for (int i = 1; i < keys.size(); i++)
	{
		vec3 avg = (keys[i - 1].second + keys[i].second) / 2;
		ctrlPoints.push_back(avg);
	}
}

void AHermiteInterpolatorVec3::computeControlPoints(
  const std::vector<ASplineVec3::Key>& keys,
  std::vector<vec3>& ctrlPoints,
  vec3& startPoint, vec3& endPoint)
{
  ctrlPoints.clear();
  if (keys.size() <= 1) return;

  int numKeys = keys.size();


  // TODO: 
  // For each key point pi, compute the corresonding value of the slope pi_prime.
  // Hints: Using Eigen::MatrixXd for a matrix data structures, 
  // this can be accomplished by solving the system of equations AC=D for C.
  // Don't forget to save the values computed for C in ctrlPoints
  // For clamped endpoint conditions, set 1st derivative at first and last points (p0 and pm) to s0 and s1, respectively
  // For natural endpoints, set 2nd derivative at first and last points (p0 and pm) equal to 0

  // Step 1: Initialize A
  // Step 2: Initialize D
  // Step 3: Solve AC=D for C
  // Step 4: Save control points in ctrlPoints
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(keys.size(),keys.size());
  for (int i = 0; i < keys.size(); i++) 
  {
	  if (i == 0) {
		  A(0, 0) = 2.0;
		  A(0, 1) = 1.0;
	  }
	  else if (i == keys.size() - 1) {
		  A(keys.size() - 1, keys.size() - 2) = 1.0;
		  A(keys.size() - 1, keys.size() - 1) = 2.0;
	  }
	  else {
		  A(i, i - 1) = 1;
		  A(i, i) = 4;
		  A(i, i + 1) = 1;
	  }
  }

  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(keys.size(),3);
  for (int i = 0; i < keys.size(); i++) {
	  if (i == 0) {
		  vec3 re = keys[1].second - keys[0].second;
		  D(0, 0) = re.n[0];
		  D(0, 1) = re.n[1];
		  D(0, 2) = re.n[2];
	  }
	  else if (i == keys.size() - 1) {
		  vec3 re = keys[keys.size()-1].second - keys[keys.size()-2].second;
		  D(keys.size() - 1, 0) = re.n[0];
		  D(keys.size() - 1, 1) = re.n[1];
		  D(keys.size() - 1, 2) = re.n[2];
	  }
	  else {
		  vec3 re = 3 * (keys[i + 1].second - keys[i - 1].second);
		  D(i, 0) = re.n[0];
		  D(i, 1) = re.n[1];
		  D(i, 2) = re.n[2];
	  }
  }

  Eigen::MatrixXd C = A.inverse() * D;
  for (int i = 0; i < keys.size(); i++) {
	  vec3 point1Prime(C(i,0), C(i,1), C(i,2));
	  ctrlPoints.push_back(point1Prime);
  }
}



double dN(std::vector<int> &knots, int n, int j, double t, int l)
{
	if (l == 0) {
		return N(knots, n, j, t);
	}

	double c1 = 1.0 / (knots[j + n] - knots[j]);
	double c2 = 1.0 / (knots[j + n + 1] - knots[j + 1]);

	return n * (c1 * dN(knots, n - 1, j, t, l -1) - c2 * dN(knots, n - 1, j + 1, t, l - 1));
}


void ABSplineInterpolatorVec3::computeControlPoints(
  const std::vector<ASplineVec3::Key>& keys,
  std::vector<vec3>& ctrlPoints,
  vec3& startPt, vec3& endPt)
{
  ctrlPoints.clear();
  if (keys.size() <= 1) return;


  // TODO: c
  // Hints: 
  // 1. use Eigen::MatrixXd to calculate the control points by solving the system of equations AC=D for C

  // 2. Create a recursive helper function dN(knots,n,t,l) to calculate derivative BSpline values at t, where
  //     knots = knot array
  //	   n = degree of the spline curves (n =3 for cubic)
  //     j = interval on knot vector in which to interpolate
  //     t = time value
  //     l = derivative (l = 1 => 1st derivative)

  // Step 1: Calculate knot vector using a uniform BSpline
  //         (assume knots are evenly spaced 1 apart and the start knot is at time = 0.0)

  // Step 2: Calculate A matrix  for a natural BSpline
  //         (Set 2nd derivative at t0 and tm to zero, where tm is the last point knot; m = #segments)

  // Step 3: Calculate  D matrix composed of our target points to interpolate

  // Step 4: Solve AC=D for C 

  // Step 5: save control points in ctrlPoints
  std::vector<int> knots;

  for (int i = 0; i < keys.size() + 6; i++)
  {
	  knots.push_back(i - 3);
  }

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(keys.size() + 2, keys.size() + 2);

  int m = keys.size() - 1;
  for (int i = 0; i < keys.size() + 2; i++)
  {
	  if (i == 0) {
		  A(0, 0) = dN(knots, 3, 0, 0, 2);
		  A(0, 1) = dN(knots, 3, 1, 0, 2);
		  A(0, 2) = dN(knots, 3, 2, 0, 2);
		  A(0, 3) = dN(knots, 3, 3, 0, 2);
	  }
	  else if (i == keys.size() + 1) {
		  A(i, i - 3) = dN(knots, 3,  m  - 1, m, 2);
		  A(i, i - 2) = dN(knots, 3, m, m, 2);
		  A(i, i - 1) = dN(knots, 3, m + 1, m, 2);
		  A(i, i) = dN(knots, 3, m + 2, m, 2);
	  }
	  else if(i == keys.size()){
		  A(i, i - 2) = N(knots, 3, m - 1, m);
		  A(i, i - 1) = N(knots, 3, m, m);
		  A(i, i) = N(knots, 3, m + 1, m);
		  A(i, i + 1) = N(knots, 3, m + 2, m);
	  }
	  else {
		  A(i, i - 1) = N(knots, 3, i - 1, i - 1);
		  A(i, i) = N(knots, 3, i, i - 1);
		  A(i, i + 1) = N(knots, 3, i + 1, i - 1);
		  A(i, i + 2) = N(knots, 3, i + 2, i - 1);
	  }
  }

  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(keys.size() + 2, 3);
  for (int i = 0; i < keys.size() + 2; i++) 
  {
	  if (i == 0 || i == keys.size() + 1) {
		  D(i, 0) = 0.0;
		  D(i, 1) = 0.0;
		  D(i, 2) = 0.0;
	  }
	  else {
		  D(i, 0) = keys[i - 1].second[0];
		  D(i, 1) = keys[i - 1].second[1];
		  D(i, 2) = keys[i - 1].second[2];
		
	  }
  }

  Eigen::MatrixXd C = A.inverse() * D;
  for (int i = 0; i < keys.size() + 2; i++) {
	  vec3 point1Prime(C(i, 0), C(i, 1), C(i, 2));
	  ctrlPoints.push_back(point1Prime);
  }
}

// Processing code by Etienne JACOB
// motion blur template by beesandbombs

// kind of implemented towards future use of more than 2 rubiks levels
// linear complexity during rendering

////////////////////////////
// Utils :

int[][] result;
float t, c;

float c01(float x)
{
  return constrain(x,0,1);
}

float ease(float p) {
  return 3*p*p - 2*p*p*p;
}

float ease(float p, float g) {
  if (p < 0.5) 
    return 0.5 * pow(2*p, g);
  else
    return 1 - 0.5 * pow(2*(1 - p), g);
}

float map(float x, float a, float b, float c, float d, boolean constr)
{
  return constr ? constrain(map(x,a,b,c,d),min(c,d),max(c,d)) : map(x,a,b,c,d);
}

float mp01(float x, float a, float b)
{
  return map(x,a,b,0,1,true);
}

float pow_(float p,float g)
{
  return 1-pow(1-p,g);
}

float tanh(float x)
{
  return (float)Math.tanh(x);
}

float softplus(float q,float p){
  float qq = q+p;
  if(qq<=0){
    return 0;
  }
  if(qq>=2*p){
    return qq-p;
  }
  return 1/(4*p)*qq*qq;
}

float mn = .5*sqrt(3), ia = atan(sqrt(.5));

void push() {
  pushMatrix();
  pushStyle();
}

void pop() {
  popStyle();
  popMatrix();
}

void draw() {

  if (!recording) {
    t = (mouseX*1.3/width)%1;
    c = mouseY*1.0/height;
    if (mousePressed)
      println(c);
    draw_();
  } else {
    for (int i=0; i<width*height; i++)
      for (int a=0; a<3; a++)
        result[i][a] = 0;

    c = 0;
    for (int sa=0; sa<samplesPerFrame; sa++) {
      t = map(frameCount-1 + sa*shutterAngle/samplesPerFrame, 0, numFrames, 0, 1);
      draw_();
      loadPixels();
      for (int i=0; i<pixels.length; i++) {
        result[i][0] += pixels[i] >> 16 & 0xff;
        result[i][1] += pixels[i] >> 8 & 0xff;
        result[i][2] += pixels[i] & 0xff;
      }
    }

    loadPixels();
    for (int i=0; i<pixels.length; i++)
      pixels[i] = 0xff << 24 | 
        int(result[i][0]*1.0/samplesPerFrame) << 16 | 
        int(result[i][1]*1.0/samplesPerFrame) << 8 | 
        int(result[i][2]*1.0/samplesPerFrame);
    updatePixels();
    
    if (frameCount<=numFrames)
    {
      saveFrame("fr###.gif");
      println(frameCount,"/",numFrames);
    }
    
    if (frameCount==numFrames)
      stop();
  }
}

//////////////////////////////////////////////////////////////////////////////
// Main code :

int samplesPerFrame = 5;
int numFrames = 390;        
float shutterAngle = .6;

boolean recording = true;

float L = 85; // global size
int N = 6; // number of main events (big rotations)

PVector rotX(PVector u,float angle)
{
  float x = u.x;
  float y = u.y*cos(angle) - u.z*sin(angle);
  float z = u.y*sin(angle) + u.z*cos(angle);
  return new PVector(x,y,z);
}

PVector rotY(PVector u,float angle)
{
  float y = u.y;
  float x = u.x*cos(angle) - u.z*sin(angle);
  float z = u.x*sin(angle) + u.z*cos(angle);
  return new PVector(x,y,z);
}

PVector rotZ(PVector u,float angle)
{
  float z = u.z;
  float y = u.y*cos(angle) - u.x*sin(angle);
  float x = u.y*sin(angle) + u.x*cos(angle);
  return new PVector(x,y,z);
}

float eps = 0.01;

float sgn(float x)
{
  return (x>0?1:(x<0?-1:0));
}

class RubiksChange
{
  int type;
  float rotationAngle;
  
  RubiksChange(int type_, float angle)
  {
    type = type_;
    rotationAngle = angle;
  }
  
  PVector angleChanges(PVector u, float q)
  {
    PVector diff = new PVector(0,0,0);
    if((type==0 && u.z>eps) || (type==1 && u.z<-eps) || (type==2 && abs(u.z)>eps))
    {
        diff.z = q*rotationAngle*sgn(u.z);
    }
    if((type==3 && u.x>eps) || (type==4 && u.x<-eps) || (type==5 && abs(u.x)>eps))
    {
        diff.x = q*rotationAngle*sgn(u.x);
    }
    if((type==6 && u.y>eps) || (type==7 && u.y<-eps) || (type==8 && abs(u.y)>eps))
    {
        diff.y = q*rotationAngle*sgn(u.y);
    }
    return diff;
  }
  
  PVector applyRotation(PVector u, float q)
  {
    PVector angles = angleChanges(u,q);
    if((type==0 && u.z>eps) || (type==1 && u.z<-eps) || (type==2 && abs(u.z)>eps))
    {
        return rotZ(u,angles.z);
    }
    if((type==3 && u.x>eps) || (type==4 && u.x<-eps) || (type==5 && abs(u.x)>eps))
    {
        return rotX(u,angles.x);
    }
    if((type==6 && u.y>eps) || (type==7 && u.y<-eps) || (type==8 && abs(u.y)>eps))
    {
        return rotY(u,angles.y);
    }
    return u;
  }
}

ArrayList<RubiksChange> generate3x3x3RChangeSequence(int nChanges)
{
  ArrayList<RubiksChange> res = new ArrayList<RubiksChange>();
  res.add(new RubiksChange(2,HALF_PI));
  res.add(new RubiksChange(8,HALF_PI));
  res.add(new RubiksChange(3,HALF_PI));
  res.add(new RubiksChange(1,HALF_PI));
  res.add(new RubiksChange(5,HALF_PI));
  res.add(new RubiksChange(6,HALF_PI));
  return res;
}

ArrayList<RubiksChange> generate2x2x2RandomRChangeSequence(int nChanges)
{
  ArrayList<RubiksChange> res = new ArrayList<RubiksChange>();
  for(int i=0;i<nChanges;i++)
  {

    RubiksChange ch = new RubiksChange(0,(random(1)<0.5?-HALF_PI:HALF_PI));
    boolean ok = false;
    
    // while loop to avoid successive similar changes (find a good change when compared to previous one)
    while(!ok)
    {
      // decision to avoid double rotations for second level changes, not sure why
      int typeChoice = floor(random(6));
      int category = typeChoice/2;
      int orientation = typeChoice%2;
      ch = new RubiksChange(3*category+orientation,(random(1)<0.5?-HALF_PI:HALF_PI));
      if(i==0)
      {
        ok = true;
      }
      else
      {
        RubiksChange prev = res.get(i-1);
        if(ch.type/3!=prev.type/3)
        {
          ok = true;
        }
      }
    }
    res.add(ch);
  }
  return res;
}

PVector applyDiffRot(PVector u,PVector diffRot)
{
  u = rotX(u,diffRot.x);
  u = rotY(u,diffRot.y);
  u = rotZ(u,diffRot.z);
  return u;
}

class SquareStyle
{
  // some of these are unused, keeping all of it to have same results
  float dummy = random(1); // just to play on random values while keeping the main seed
  float dummy2 = random(1); // just to play on random values while keeping the main seed
  float strW = random(0.75,1.15);
  float fsquare = random(0.1,0.7);
  float strokeColor = random(120,280)*0+255;
  boolean qp = (random(1)<0.5?true:false);
  boolean fill = (random(1)<0.5?true:false);
  boolean edges = (random(1)<0.5?true:false);
}

class RubiksCube
{
  ArrayList<PVector> mainPositions = new ArrayList<PVector>();
  ArrayList<PVector> rotations = new ArrayList<PVector>();
  ArrayList<PVector> base1Checkpoint = new ArrayList<PVector>();
  ArrayList<PVector> base2Checkpoint = new ArrayList<PVector>();
  ArrayList<PVector> base3Checkpoint = new ArrayList<PVector>();
  
  ArrayList<RubiksChange> changes;
  int cn; // number of changes
  
  int startX,startY,startZ,endX,endY,endZ;
  
  // size and spacing parameters
  float l;
  float d;
  
  int depth;
  
  RubiksCube [] subCubeArray = new RubiksCube[8]; // for first level
  RubiksCube parent; // for second level
  
  SquareStyle ss;
  boolean squareStyleIsAssigned = false;
  
  
  RubiksCube(int X,int Y,int Z,ArrayList<RubiksChange> changes_,float l_,float d_,int depth_,RubiksCube parent_)
  {
    depth = depth_;
    parent = parent_;
    
    mainPositions.add(new PVector(X,Y,Z));
    base1Checkpoint.add(new PVector(1,0,0));
    base2Checkpoint.add(new PVector(0,1,0));
    base3Checkpoint.add(new PVector(0,0,1));
    
    startX = X;
    startY = Y;
    startZ = Z;
    
    changes = changes_;
    
    cn = changes.size();
    
    // precomputing the steps
    for(int i=0;i<cn;i++)
    {
      PVector pos = mainPositions.get(i).copy();
      PVector angleDiff  = changes.get(i).angleChanges(pos,1).copy();
      rotations.add(angleDiff);
      PVector pos2 = changes.get(i).applyRotation(pos,1);
      mainPositions.add(pos2);
      
      base1Checkpoint.add(applyDiffRot(base1Checkpoint.get(i),angleDiff));
      base2Checkpoint.add(applyDiffRot(base2Checkpoint.get(i),angleDiff));
      base3Checkpoint.add(applyDiffRot(base3Checkpoint.get(i),angleDiff));
    }
    
    if(depth==1)
    {
      // we want to find real positions at the end in canonical basis, not the rotated one ...
      
      // rotated basis at the end :
      PVector base1 = parent.base1Checkpoint.get(parent.cn);
      PVector base2 = parent.base2Checkpoint.get(parent.cn);
      PVector base3 = parent.base3Checkpoint.get(parent.cn);
      
      PVector vend = mainPositions.get(cn); // position in previous basis
      
      // basis change ...
      endX = (int)round(base1.x*vend.x+base2.x*vend.y+base3.x*vend.z);
      endY = (int)round(base1.y*vend.x+base2.y*vend.y+base3.y*vend.z);
      endZ = (int)round(base1.z*vend.x+base2.z*vend.y+base3.z*vend.z);
    }
    else
    {
      endX = (int)round(mainPositions.get(cn).x);
      endY = (int)round(mainPositions.get(cn).y);
      endZ = (int)round(mainPositions.get(cn).z);
    }
    
    l = l_;
    d = d_;
    
    ss = new SquareStyle();
  }
  
  void makeSubCubes()
  {
    ArrayList<RubiksChange> seq = generate2x2x2RandomRChangeSequence(floor(random(1.0*N,3.0*N)));
    
    int a = 0;
    for(int i=-1;i<=1;i+=2)
    {
      for(int j=-1;j<=1;j+=2)
      {
        for(int k=-1;k<=1;k+=2)
        {
          subCubeArray[a] = new RubiksCube(i,j,k,seq,l/2.2,d/5.5,depth+1,this);
          a++;
        }
      }
    }
  }
  
  // recursive square style assignment to match the permutation at the end of the loop
  void assign(SquareStyle ss_)
  {
    ss = ss_;
    squareStyleIsAssigned = true;
    
    // brute force (only done at setup) to match the cubes, because I haven't learned data structures in Java, I got away with that so far :)
    for(int i=0;i<n;i++)
    {
      for(int a=0;a<8;a++)
      {
        RubiksCube testCube = array[i].subCubeArray[a];
        if(!testCube.squareStyleIsAssigned && array[i].endX==parent.startX && array[i].endY==parent.startY && array[i].endZ==parent.startZ && testCube.endX==startX && testCube.endY==startY && testCube.endZ==startZ)
        {
          array[i].subCubeArray[a].squareStyleIsAssigned = true;
          array[i].subCubeArray[a].assign(ss_);
        }
      }
    }
  }
  
  void drawFace(float l)
  {
    stroke(255);
    strokeWeight(ss.edges?0:1.9);
    fill(0);
    rectMode(CENTER);
    rect(0,0,l,l);
    
    push();
    translate(0,0,1.5);
    strokeWeight(ss.strW);
    stroke(ss.strokeColor);
    fill(ss.fill?255:0);
    
    float f2 = (ss.qp?0.6:1)*0+(ss.fill?0.7:1);
    float lx = l*ss.fsquare*f2;
    float ly = l*ss.fsquare*f2;
    
    rect(0,0,lx,ly);
    
    pop();
  }
  
  void showSmallCube(float l)
  {
    fill(0);
  
    push();
    translate(0,0,l/2);
    drawFace(l);
    pop();
    
    push();
    rotateY(PI);
    translate(0,0,l/2);
    drawFace(l);
    pop();
    
    push();
    rotateX(HALF_PI);
    push();
    translate(0,0,l/2);
    drawFace(l);
    pop();
    
    push();
    rotateY(PI);
    translate(0,0,l/2);
    drawFace(l);
    pop();
    pop();
    
    push();
    rotateY(HALF_PI);
    push();
    translate(0,0,l/2);
    drawFace(l);
    pop();
    
    push();
    rotateY(PI);
    translate(0,0,l/2);
    drawFace(l);
    pop();
    pop();
  }
  
  float pToFloatIndex(float p)
  {
    return cn*p;
  }
  
  // https://easings.net/#easeOutElastic
  float easeOutElastic(float x)
  {
    float c4 = (2*PI)/3;
    if(x<=0) return 0;
    if(x>=1) return 1;
    return pow(2, -9 * x) * sin((x * 7 - 0.75) * c4) + 1;
  }
  
  void show(float p)
  {
    push();
    float fInd = pToFloatIndex(ease(p,1.0)); // no easing but it's a possibillity
    int ind = floor(fInd);
    float frac = fInd - ind;
    float easedFrac = easeOutElastic(pow(frac,1.9));
    PVector pos = mainPositions.get(ind);
    PVector pixelPos = changes.get(ind).applyRotation(pos,easedFrac).copy().mult(d);
    
    translate(pixelPos.x,pixelPos.y,pixelPos.z);
    
    // using "checkpoints" to avoid recomputing the base rotations from the start, just the last one is computed
    PVector base1 = base1Checkpoint.get(ind);
    PVector base2 = base2Checkpoint.get(ind);
    PVector base3 = base3Checkpoint.get(ind);
    
    PVector diffRot = changes.get(ind).angleChanges(pos,easedFrac);
    base1 = applyDiffRot(base1,diffRot);
    base2 = applyDiffRot(base2,diffRot);
    base3 = applyDiffRot(base3,diffRot);
    
    // looking for 2D rotations from orthonormal base change ...
    // https://math.stackexchange.com/questions/3690075/computing-a-3d-rotation-matrix-aligning-1-orthonormal-basis-to-another
    // https://nghiaho.com/?page_id=846
    
    float thetaZ = atan2(base1.y,base1.x);
    float thetaY = atan2(-base1.z,sqrt(base2.z*base2.z+base3.z*base3.z));
    float thetaX = atan2(base2.z,base3.z);
    
    rotateZ(thetaZ);
    rotateY(thetaY);
    rotateX(thetaX);
    
    if(depth==0)
    {
      for(int i=0;i<8;i++)
      {
        subCubeArray[i].show(p);
      }
    }
    else
    {
      showSmallCube(l/1.5);
    }
    
    pop();
  }
}

int n = (int)pow(3,3); // number of main cubes
RubiksCube [] array = new RubiksCube[n];

ArrayList<RubiksChange> mainRChangesSequence;

void setup(){
  size(600,600,P3D);
  result = new int[width*height][3];
  
  randomSeed(12351);
  
  mainRChangesSequence = generate3x3x3RChangeSequence(N);
  
  int a = 0;
  
  for(int i=-1;i<=1;i++)
  {
    for(int j=-1;j<=1;j++)
    {
      for(int k=-1;k<=1;k++)
      {
        array[a] = new RubiksCube(i,j,k,mainRChangesSequence,L,L,0,null); // main cubes don't have a parent
        a++;
      }
    }
  }
  
  for(int i=0;i<n;i++)
  {
    array[i].makeSubCubes();
  }
  
  ortho();
  
  smooth(4);
  
  for(int i=0;i<n;i++)
  {
    for(int b=0;b<8;b++)
    {
      if(!array[i].subCubeArray[b].squareStyleIsAssigned)
      {
        array[i].subCubeArray[b].assign(new SquareStyle());
      }
    }
  }
}


void draw_(){
  background(0);
  push();
  translate(width/2,height/2);
  
  rotateX(-0.393*HALF_PI); // (lol too lazy to do the maths)
  rotateY(0.5*HALF_PI);
  
  for(int i=0;i<n;i++)
  {
    array[i].show(ease(t,1.0)); // no easing but it's a possibility
  }

  pop();
}

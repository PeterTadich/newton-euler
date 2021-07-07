# newton-euler
Newton-Euler formulation

## Dependencies

There is 1 dependency 'matrix-computations'.

```bash
https://github.com/PeterTadich/matrix-computations
```

## Installation

### Node.js

```bash
npm install https://github.com/PeterTadich/newton-euler
```

### Google Chrome Web browser

No installation required for the Google Chrome Web browser.

## How to use

### Node.js

```js
import * as mcnef from 'newton-euler';
```

### Google Chrome Web browser

```js
import * as mcnef from './mcnef.mjs';
```

## Examples

### Node.js (server side)

Copy the following code to index.mjs

```js
//Two link manipulator - joint torque calculation
import * as mcnef from 'newton-euler';
import * as ttvm from 'trajectories';

//To do:
//   - fix I1zz calc
//   - mm1, mm2 are not being used

//parameters
//Denavit-Hartenberg parameters
var a1 = 1.0; //m
var a2 = 1.0; //m
//position of centre of mass from frame {i}
var l1 = 0.5; //m
var l2 = 0.5; //m
//inertia of the link
var Il1 = 10.0; //kg.m^2
var Il2 = 10.0; //kg.m^2
//mass of the link
var ml1 = 50.0; //kg
var ml2 = 50.0; //kg
//var Im1 = 0.01; //kg.m^2
//var Im2 = 0.01; //kg.m^2
//inertia of the rotor
var Im1 = 0.0; //kg.m^2 (IMPORTANT: mua.Ir = [0.0,0.01,0.01])
var Im2 = 0.0; //kg.m^2
//mass of the rotor
var mm1 = 5.0; //kg
var mm2 = 5.0; //kg
//Gear reduction ratio. IMPORTANT: it is not squared
var kr1 = 100.0;
var kr2 = 100.0;
//gravity
var g = 9.81; //ref: page 289.

var lc1 = l1-a1;
var lc2 = l2-a2;

//inertia
var I = [];
var I1zz = Il1+ml1*Math.pow((l1-a1),2)+Im2-ml1*Math.pow(lc1,2); //See page 292. (IMPORTANT: Im2 = 0.0)
I[0] = [
    [0.0,0.0, 0.0],
    [0.0,0.0, 0.0],
    [0.0,0.0,I1zz],
];
var I2zz = Il2+ml2*Math.pow((l2-a2),2)-ml2*Math.pow(lc2,2); //See page 292. 
I[1] = [
    [0.0,0.0, 0.0],
    [0.0,0.0, 0.0],
    [0.0,0.0,I2zz],
];

//joint position - initial condition
var v1 = -1.0*Math.PI/2.0; //rad
var v2 = Math.PI; //rad
var v1d = 0.0; //rad/s
var v2d = 0.0; //rad/s
var v1dd = 26.0; //rad.s^-2
var v2dd = 26.0; //rad.s^-2

var mua = {
    //Denavit-Hartenberg parameters
    DH: [ //[ai, alpha_i,  di, vi].
            [a1,     0.0, 0.0, v1], //'ai' in m, 'vi' in radians being joint position
            [a2,     0.0, 0.0, v2]
        ],
    //joint type
    jt:['R','R'],
    //joint kinematics
    v: [v1,v2], //joint position (same as q, qd, qdd (generalized coordinates) for revolute joint)
    vd: [v1d,v2d], //joint velocity
    vdd: [v1dd,v2dd], //joint acceleration
    //link kinematics
    w: [[[0.0],[0.0],[0.0]]], //link 0 (frame {0}) angular velocity (initial condition)
    wd: [[[0.0],[0.0],[0.0]]], //link 0 (frame {0}) angular acceleration (initial condition)
    wdr: [], //rotor
    pd: [], //link velocity
    pdd: [[[0.0],[0.0],[0.0]]], //link acceleration (initial condition)
    pcdd: [], //CoM acceleration
    R: [[
        [1.0,0.0,0.0],
        [0.0,1.0,0.0],
        [0.0,0.0,1.0]
    ]],
    //14 dynamic parameters per joint (page 262, 263):
    m: [ml1,ml2], //(1) mass of the links
    mr: [], //(1) mass of the rotor
    rc: [[[lc1],[0.0],[0.0]],[[lc2],[0.0],[0.0]]], //(3) components of the first moment of inertia (7.72)
    I: I, //(6) components of the inertia tensor in (7.73),
    Ir: [0.01,0.01], //(1) the moment of inertia of the rotor (IMPORTANT: need to check this)
    Fvi: [], //(1) viscous friction coefficient Fvi
    Fsi: [], //(1) Coulomb friction coefficient Fsi
    //link velocities derivation (page 108)
    r: [[[a1],[0.0],[0.0]],[[a2],[0.0],[0.0]]], //position of {i} w.r.t {i-1}
    //others
    g: [[0.0],[-9.81],[0.0]], //   - define 'g' (gravity) direction.
    //g: [[0.0],[0.0],[-9.81]], //   - define 'g' (gravity) direction.
    kr: [kr1,kr2], //   - gear reduction ratio. kr.
    Bm: [0.0,0.0],
    Tc: [[0.0,-0.0],[0.0,-0.0]]
};

//setup the time-step
var tf = 0.5; //total time (seconds).
var nstep = 100; //number of time steps.

//Joint 1.
var v1i = -1.0*Math.acos(0.1/1.0); //initial joint position (revolute).
var v1f = v1i + 0.5*Math.PI; //final joint position (revolute).
var pos_j1 = ttvm.lspb(v1i,v1f,tf,nstep); //joint angular position per time step.
var vel_j1 = ttvm.gradient(pos_j1); //joint angular velocity per time step.
var acc_j1 = ttvm.gradient(vel_j1); //joint angular acceleration per time step.

//Joint 2.
var v2i = 2.0*Math.abs(v1i); //initial joint position (revolute).
var v2f = v2i + 0.5*Math.PI; //final joint position (revolute).
var pos_j2 = ttvm.lspb(v2i,v2f,tf,nstep); //joint angular position per time step.
var vel_j2 = ttvm.gradient(pos_j2); //joint angular velocity per time step.
var acc_j2 = ttvm.gradient(vel_j2); //joint angular acceleration per time step.

//Joint torques.
var trq_j1 = [];
var trq_j2 = [];
var trq; //torque. [[torge joint 1],[torge joint 2]]
var t; //time.
for(var k=0;k<pos_j1.length;k=k+1){
    //adjust the manipular object
    mua.DH = [ //[ai, alpha_i,  di,           vi].
                 [a1,     0.0, 0.0, pos_j1[k][1]], //'ai' in m, 'vi' in radians being joint position
                 [a2,     0.0, 0.0, pos_j2[k][1]]
            ];
    mua.v = [pos_j1[k][1],pos_j2[k][1]];
    mua.vd = [vel_j1[k][1],vel_j2[k][1]];
    mua.vdd = [acc_j1[k][1],acc_j2[k][1]];
    
    //Lagrange Formulation 'twoLinkplanarArmLF()'
    //trq = mcnef.twoLinkplanarArmLF([pos_j1[k][1],pos_j2[k][1]],[vel_j1[k][1],vel_j2[k][1]],[acc_j1[k][1],acc_j2[k][1]]);
    
    //Newton-Euler (symbolic) implementation 'twoLinkplanarArmNE()'
    //trq = mcnef.twoLinkplanarArmNE([pos_j1[k][1],pos_j2[k][1]],[vel_j1[k][1],vel_j2[k][1]],[acc_j1[k][1],acc_j2[k][1]]);
    
    /*
    //Newton-Euler (recursive) - current frame
    //   - forward recursion
    kin = mcnef.linkAccelerationsCF(mua);
    //   - backward recursion
    trq = mcnef.linkForcesCF(mua);
    */
    //or
    trq = mcnef.NewtonEulerRecursion(mua,"CF");
    t = pos_j1[k][0];
    trq_j1.push([t,trq[1]]);
    trq_j2.push([t,trq[2]]);
    console.log("torques: T1 = " + trq[1].toFixed(4) + ", T2 = " + trq[2].toFixed(4));
}
```

Then run:

```bash
npm init -y
npm install https://github.com/PeterTadich/newton-euler https://github.com/PeterTadich/trajectories
node index.mjs
```

If the above does not work modify the package.json file as follows:
Helpful ref: [https://stackoverflow.com/questions/45854169/how-can-i-use-an-es6-import-in-node-js](https://stackoverflow.com/questions/45854169/how-can-i-use-an-es6-import-in-node-js)

```js
"scripts": {
    "test": "echo \"Error: no test specified\" && exit 1",
    "start": "node --experimental-modules index.mjs"
  },
"type": "module",
```

```bash
npm start
```

Result:

```js
torques: T1 = 1907.5984, T2 = 1538.7923
torques: T1 = 2812.1891, T2 = 2295.7922
torques: T1 = 3716.1144, T2 = 3052.5626
torques: T1 = 3714.7592, T2 = 3051.9862
torques: T1 = 3712.8995, T2 = 3051.1711
torques: T1 = 3710.5730, T2 = 3050.1094
torques: T1 = 3707.8283, T2 = 3048.7907
torques: T1 = 3704.7247, T2 = 3047.2026
torques: T1 = 3701.3328, T2 = 3045.3302
...
```
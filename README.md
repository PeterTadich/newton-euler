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
    //Lagrange Formulation 'twoLinkplanarArmLF()'
    //trq = mcnef.twoLinkplanarArmLF([pos_j1[k][1],pos_j2[k][1]],[vel_j1[k][1],vel_j2[k][1]],[acc_j1[k][1],acc_j2[k][1]]);
    //Newton-Euler (symbolic) implementation 'twoLinkplanarArmNE()'
    //trq = mcnef.twoLinkplanarArmNE([pos_j1[k][1],pos_j2[k][1]],[vel_j1[k][1],vel_j2[k][1]],[acc_j1[k][1],acc_j2[k][1]]);
    //Newton-Euler (recursive) implementation 'linkAccelerationsCF()' ref: newtonEuler.js
    trq = mcnef.linkAccelerationsCF([pos_j1[k][1],pos_j2[k][1]],[vel_j1[k][1],vel_j2[k][1]],[acc_j1[k][1],acc_j2[k][1]]);
    t = pos_j1[k][0];
    trq_j1.push([t,trq[0][0]]);
    trq_j2.push([t,trq[1][0]]);
    console.log("torques: T1 = " + trq[0][0].toFixed(4) + ", T2 = " + trq[0][0].toFixed(4));
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
torques: T1 = 1988.0273, T2 = 1988.0273
torques: T1 = 2930.3952, T2 = 2930.3952
torques: T1 = 3872.1285, T2 = 3872.1285
torques: T1 = 3870.8500, T2 = 3870.8500
torques: T1 = 3869.0975, T2 = 3869.0975
torques: T1 = 3866.9090, T2 = 3866.9090
torques: T1 = 3864.3327, T2 = 3864.3327
...
```
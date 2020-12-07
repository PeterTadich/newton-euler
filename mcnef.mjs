// mcnef = matrix computations Newton-Euler formulation

// ECMAScript module

// Newton-Euler formulation

import * as hlao from 'matrix-computations';

// Recursion Algorithm, REF: Robotics Modelling, Planning and Control, Page 286
function NewtonEulerRecursion(){
    //input:
    //   - q,           positions
    //   - qd,          velocities
    //   - qdd,         accelerations
    //   - w[0],        angular velocity
    //   - pdd[0] - g0, linear acceleration
    //   - wd[0],       angular acceleration
    
    //forward recursion: compute the link velocities and accelerations
    /*
    linkAccelerationsBF(); //base frame
    linkAccelerationsCF(); //current frame
    */
    
    //backward recursion: compute the forces and moments acting on each link
    /*
    linkForcesBF(); //base frame
    linkForcesCF(); //current frame
    */
}

// FORWARD RECURSION
// Link Accelerations, REF: Robotics Modelling, Planning and Control, Page 285
//    - ref. Base Frame
/*
function linkAccelerationsBF(){
    //i = 1,...,n
    
    //qd
    //qdd
    
    var w = [];
    var wd = [];
    var pd = [];
    var pdd = [];
    var pddc = [];
    var wdm = [];
    var vd = [];
    var z = [];
    var r = [];
    var dd = [];
    var ddd = [];
    
    w[0] = 0.0;
    
    //revolute joint:
    
    //   - link 'i':
    //      - angular velocity
    w[i] = hlao.matrix_arithmetic(
                w[i-1],
                hlao.matrix_multiplication_scalar(
                    z[i-1],
                    vd[i]
                ),
                '+'
            ); //equ. (7.93)
    //      - linear velocity
    pd[i] = hlao.matrix_arithmetic(
                pd[i-1],
                hlao.vector_cross(w[i],r[i]),
                '+'
            ); //equ. (7.94)
    //      - angular acceleration
    wd[i] = hlao.matrix_arithmetic(
                hlao.matrix_arithmetic(
                    wd[i-1],
                    hlao.matrix_multiplication_scalar(
                        vdd[i],
                        z[i-1]
                    )
                    '+'
                ),
                hlao.matrix_multiplication_scalar(
                    hlao.vector_cross(
                        w[i-1],
                        z[i-1]
                    ),
                    vd[i]
                )
            ); //equ. (7.96)
    //      - linear acceleration
    ////pdd[i] = 
    ////            pdd[i-1] + 
    ////            hlao.matrix_multiplication_scalar(z[i-1],ddd[i]) +
    ////            hlao.matrix_multiplication_scalar(w[i-1],dd[i]) +
    ////            hlao.vector_cross(
    ////                    hlao.matrix_multiplication(
    ////                        wd[i],
    ////                        r[i]
    ////                    ),
    ////                    z[i-1]
    ////                ) + 
    ////            hlao.vector_cross(w[i],hlao.matrix_multiplication_scalar(z[i-1],dd[i])) +
    ////            hlao.vector_cross(w[i],hlao.vector_cross(w[i],r[i])); //equ. (7.97)
    
    //      - linear acceleration
	pdd[i] = hlao.matrix_arithmetic(
                hlao.matrix_arithmetic(
                    pdd[i-1],
                    hlao.vector_cross(wd[i],r[i]),
                    '+'
                ),
                hlao.vector_cross(w[i],hlao.vector_cross(w[i],r[i])),
                '+'
            ); //equ. (7.99)
    
    //   - centre of mass of link 'i':
    //      - linear acceleration
    pddc[i] = hlao.matrix_arithmetic(
                    hlao.matrix_arithmetic(
                        pdd[i],
                        hlao.vector_cross(wd[i],rc[i]),
                        '+'
                    ),
                    hlao.vector_cross(w[i],hlao.vector_cross(w[i],rc[i])),
                    '+'
                ); //equ. (7.102)
    //      - angular acceleration of the rotor
    wdm[i] = hlao.matrix_arithmetic(
                    hlao.matrix_arithmetic(
                        wd[i-1],
                        hlao.matrix_multiplication_scalar(zm[i],kr[i]*qdd[i]),
                        '+'
                    ),
                    hlao.vector_cross(
                        hlao.matrix_multiplication_scalar(w[i-1],kr[i]*qd[i]),
                        zm[i]
                    ),
                    '+'
                ); //equ. (7.103)
}
*/

// BACKWARD RECURSION
// Link forces, REF: Robotics Modelling, Planning and Control, Page 287
//    - ref. Base Frame
/*
function linkForcesBF(){
    //i = n,...,1
    
    var f = [];
    var u = [];
    var T = [];
    
    //revolute joint:
    f[i] = f[i+1] + hlao.matrix_multiplication(m[i],pcdd[i]); //equ. (7.104)
    u[i] = -1.0*hlao.matrix_multiplication(f[i],(r[i]+rc[i])) + u[i+1] + hlao.matrix_multiplication(f[i+1],rc[i]) + I[i]*wd[i] +
               hlao.matrix_multiplication(w[i],(I[i]*w[i])) + kr[i+1]*qdd[i+1]*Im[i+1]*zm[i+1] +
               hlao.matrix_multiplication(kr[i+1]*qd[i+1]*Im[i+1]*w[i],zm[i+1]); //equ. (7.105)
    T[i] = matrix_transpose(u[i])*z[i-1] + kr[i]*Im[i]*matrix_transpose(wmd[i])*zm[i] +
               Fv[i]*vd[i] + Fs[i]*sgn(vd[i]); //equ. (7.106)
}
*/

// FORWARD RECURSION
// Link Accelerations, REF: Robotics Modelling, Planning and Control, Page 287
//    - ref. Current Frame
function linkAccelerationsCF(v,vd,vdd){
    //IMPORTANT:
    //   - check qdd[] is an array of vdd???

    //13 dynamic parameters per joint (page 262, 263):
    //   - mass (1)
    //   - the three components of the first moment of inertia in (7.72),
    //   - the six components of the inertia tensor in (7.73),
    //   - the moment of inertia of the rotor (1)
    //   - viscous friction coefficient Fvi
    //   - Coulomb friction coefficient Fsi
    
    //Required:
    //   - define 'g' (gravity) direction.
    //   - Joint position, velocity, acceleration. v, vd, vdd (the same as q, qd, qdd (generalized coordinates) for revolute joint).
    //   - Denavit-Hartenberg parameters. [ai, alpha_i, di, vi].
    //   - gear reduction ratio. kr.
    //   - masses. link (ml) and rotor (mm) mass.
    //   - moment of inertia of rotor (with respect to rotor axis - motors are located on joint axes). Im.
    //   - moment of inertial of link (with respect to centre of mass of link axis). Il.
    //   - distance of the centre of mass of the link from respective joint. 'l'. ref: page 265. IMPORTANT: this definition may not be correct.
    
    //var v1 = -1.0*Math.PI/2.0; //rad
    //var v2 = Math.PI; //rad
    //var v1d = 0.0; //rad/s
    //var v2d = 0.0; //rad/s
    //var v1dd = 26.0; //rad.s^-2
    //var v2dd = 26.0; //rad.s^-2
    var v1 = v[0]; //rad. Joint position.
    var v2 = v[1]; //rad
    var v1d = vd[0]; //rad/s. Joint velocity.
    var v2d = vd[1]; //rad/s
    var v1dd = vdd[0]; //rad.s^-2. Joint acceleration.
    var v2dd = vdd[1]; //rad.s^-2
    
    var a1 = 1.0; //m
    var a2 = 1.0; //m
    var l1 = 0.5; //m
    var l2 = 0.5; //m
    var Il1 = 10.0; //kg.m^2
    var Il2 = 10.0; //kg.m^2
    var ml1 = 50.0; //kg
    var ml2 = 50.0; //kg
    var Im1 = 0.01; //kg.m^2
    var Im2 = 0.01; //kg.m^2
    var mm1 = 5.0; //kg
    var mm2 = 5.0; //kg
    var kr1 = 100.0; //Gear reduction ratio.
    var kr2 = 100.0;
    var g = 9.81; //ref: page 289.
    
    //others
    var m1 = ml1 + mm2;
    var m2 = ml2;
    var lc1 = ml1*(l1-a1)/m1;
    var lc2 = ml2*(l2-a2)/m2;
    //console.log(lc1,lc2); //lc1 and lc2 should both be negative
    
    //IMPORTANT: required data, q[], qd[], qdd[]
    var qd  = [[0.0], [v1d], [v2d]]; //two joints -  qd[0][0] is dummy data
    var qdd = [[0.0],[v1dd],[v2dd]]; //two joints - qdd[0][0] is dummy data
    
    var w = []; var wd = []; var wmd = [];
    var pdd = []; var pcdd = [];
    var R = []; //coordinate transform from Frame {i} --> Frame {i-1}
    var vd = []; var vdd = [];
    var z = [];
    var r = [];
    var rc = [];
    var kr = [];
    var zm = [];
    
    w[0] = [[0.0],[0.0],[0.0]];
    wd[0] = [[0.0],[0.0],[0.0]];
    
    vd[1] = v1d; vd[2] = v2d;
    vdd[1] = v1dd; vdd[2] = v2dd;
    z[0] = [[0.0],[0.0],[1.0]];
    r[1] = [[a1],[0.0],[0.0]];
    r[2] = [[a2],[0.0],[0.0]];
    rc[1] = [[lc1],[0.0],[0.0]];
    rc[2] = [[lc2],[0.0],[0.0]];
    pdd[0] = [[0.0],[g],[0.0]];
    kr[1] = kr1; kr[2] = kr2;
    
    // Li = [ai, alpha_i, di, vi]
    var Li = [
        [a1, 0.0, 0.0, v1],
        [a2, 0.0, 0.0, v2]
    ];
    
    var n = Li.length; //number of links
    
    for(var i=1;i<=n;i=i+1){ //coordinate transforms between Frame {i} and Frame {i-1}
        //homogeneous transformation matrix:
        //note: Denavit-Hartenberg parameters for Link 1 are indexed at '0'
        var Aip = [ //Ai',i-1
            [Math.cos(Li[i-1][3]),-1.0*Math.sin(Li[i-1][3]), 0.0,         0.0],
            [Math.sin(Li[i-1][3]),     Math.cos(Li[i-1][3]), 0.0,         0.0],
            [                 0.0,                      0.0, 1.0,  Li[i-1][2]],
            [                 0.0,                      0.0, 0.0,         1.0]
        ];

        var Ai = [ //Ai,i'
            [      1.0,                  0.0,                      0.0,  Li[i-1][0]],
            [      0.0, Math.cos(Li[i-1][1]),-1.0*Math.sin(Li[i-1][1]),         0.0],
            [      0.0, Math.sin(Li[i-1][1]),     Math.cos(Li[i-1][1]),         0.0],
            [      0.0,                  0.0,                      0.0,         1.0]
        ];
        
        var Ai_ineg1 = hlao.matrix_multiplication(Aip,Ai);
        
        //rotation matrix
        R[i] = [
            [Ai_ineg1[0][0],Ai_ineg1[0][1],Ai_ineg1[0][2]],
            [Ai_ineg1[1][0],Ai_ineg1[1][1],Ai_ineg1[1][2]],
            [Ai_ineg1[2][0],Ai_ineg1[2][1],Ai_ineg1[2][2]]
        ];
    }
    
    //initialize
    for(var i=1;i<=n;i=i+1){ //i = 1,...,n
        zm[i] = z[0]; //zmi,i-1 (axis of rotation of rotors) coincide with z[0] (joint axis), page 289. zmi = zi-1 see page 266.
    }
    
    //forward recursion
    for(var i=1;i<=n;i=i+1){ //i = 1,...,n
        //revolute joint:
        
        //   - angular velocity
        w[i] = hlao.matrix_multiplication(matrix_transpose(R[i]),hlao.matrix_arithmetic(w[i-1],hlao.matrix_multiplication_scalar(z[0],vd[i]),'+')); //equ. (7.107)
        console.log(w[i]);
        
        //   - angular acceleration
        wd[i] = hlao.matrix_multiplication(
                   matrix_transpose(R[i]),
                   hlao.matrix_arithmetic(
                      hlao.matrix_arithmetic(wd[i-1],hlao.matrix_multiplication_scalar(z[0],vdd[i]),'+'), 
                      vector_cross(hlao.matrix_multiplication_scalar(w[i-1],vd[i]),z[0]),
                   '+')
                ); //equ. (7.108)
        //console.log(wd[i]);
        
        //   - linear acceleration, link
        pdd[i] = hlao.matrix_arithmetic(
                    hlao.matrix_arithmetic(
                       hlao.matrix_multiplication(matrix_transpose(R[i]),pdd[i-1]),
                       vector_cross(wd[i],r[i]),
                    '+'),
                     vector_cross(w[i],vector_cross(w[i],r[i])),
                 '+'); //equ. (7.109)
        //console.log(pdd[i]);
        
        //   - linear acceleration, centre of mass of link 'i'
        pcdd[i] = hlao.matrix_arithmetic(
                     hlao.matrix_arithmetic(
                        pdd[i],
                        vector_cross(wd[i],rc[i]),
                     '+'),
                     vector_cross(w[i],vector_cross(w[i],rc[i])),
                  '+'); //equ. (7.110)
        //console.log(pcdd[i]);
        
        //   - angular acceleration of the rotor
        wmd[i] = hlao.matrix_arithmetic(
                    hlao.matrix_arithmetic(
                       wd[i-1],
                       hlao.matrix_multiplication_scalar(zm[i],(kr[i]*qdd[i][0])),
                    '+'),
                    vector_cross(hlao.matrix_multiplication_scalar(w[i-1],(kr[i]*qd[i][0])),zm[i]),
                 '+'); //equ. (7.111)
        //console.log(wmd[i]);
    }
    
    //IMPORTANT: move backward recursion into function 'linkForcesCF()'
    var f = [];
    var u = [];
    var T = [];
    var m = [];
    var I = [];
    
    m[0] = 0.0; /*dummy*/ m[1] = m1; m[2] = m2;

    //coordinate transforms between Frame {i+1} and Frame {i}
    R[n+1] = identity_matrix(3);
    
    //forces:
    f[n+1] = [[0.0],[0.0],[0.0]];
    //and
    u[n+1] = [[0.0],[0.0],[0.0]];
    
    //inertia
    var I1zz = Il1+ml1*Math.pow((l1-a1),2)+Im2-m1*Math.pow(lc1,2); //See page 292. 
    I[1] = [
        [0.0,0.0, 0.0],
        [0.0,0.0, 0.0],
        [0.0,0.0,I1zz],
    ];
    var I2zz = Il2+ml2*Math.pow((l2-a2),2)-m2*Math.pow(lc2,2); //See page 292. 
    I[2] = [
        [0.0,0.0, 0.0],
        [0.0,0.0, 0.0],
        [0.0,0.0,I2zz],
    ];
    var Im = [0.0,Im1,Im2]; //Im[0] is dummy
    
    //IMPORTANT check the following:
    kr[3] = 0.0;
    qd.push([0.0]);
    qdd.push([0.0]);
    zm[3] = [[0.0],[0.0],[0.0]];
    Im.push(0.0);
    
    //backward recursion
    for(var i=n;i>=1;i=i-1){//i = n,...,1
        //revolute joint:
        
        //   - force
        f[i] = hlao.matrix_arithmetic(
                  hlao.matrix_multiplication(R[i+1],f[i+1]),
                  hlao.matrix_multiplication_scalar(pcdd[i],m[i]),
               '+'); //equ. (7.112)
        console.log(f[i]);
        
        //   - u
        //NEED TO FIX - not working correctly
        //   - check zm[] is in the right reference frame
        u[i] = hlao.matrix_arithmetic(
                  hlao.matrix_multiplication_scalar(vector_cross(f[i],hlao.matrix_arithmetic(r[i],rc[i],'+')),-1.0),
                  hlao.matrix_arithmetic(
                     hlao.matrix_arithmetic(hlao.matrix_multiplication(R[i+1],u[i+1]), 
                        vector_cross(hlao.matrix_multiplication(R[i+1],f[i+1]),rc[i]),'+'),
                     hlao.matrix_arithmetic(
                        hlao.matrix_arithmetic(hlao.matrix_multiplication(I[i],wd[i]),
                           vector_cross(w[i],hlao.matrix_multiplication(I[i],w[i])),'+'),
                        hlao.matrix_arithmetic(hlao.matrix_multiplication_scalar(zm[i+1],kr[i+1]*qdd[i+1][0]*Im[i+1]),
                           vector_cross(hlao.matrix_multiplication_scalar(w[i],kr[i+1]*qd[i+1][0]*Im[i+1]),zm[i+1]),'+'),
                      '+'),
                   '+'),
                '+'); //equ. (7.113)
        console.log(u[i]);
        
        //   - torque
        T[i] = vector_dot(hlao.matrix_multiplication(matrix_transpose(u[i]),matrix_transpose(R[i])),z[0]) +
               vector_dot(vector_transpose(wmd[i]),zm[i])*kr[i]*Im[i]; //equ. (7.114)
        //T[i] = matrix_transpose(u[i])*matrix_transpose(R[i])*z[0] + kr[i]*Im[i]*matrix_transpose(wmd[i])*zm(i) +
        //           Fv[i]*vd[i] + Fs[i]*sgn(vd[i]); //equ. (7.114)
        console.log(T[i]);
    }
    
    return [[T[1]],[T[2]]];
}

// BACKWARD RECURSION
// Link forces, REF: Robotics Modelling, Planning and Control, Page 288
//    - ref. Current Frame
/*
function linkForcesCF(){
    var f = [];
    var u = [];
    var T = [];

    //backward recursion
    for(var i=n;i>=1;i=i-1){//i = n,...,1
        //revolute joint:
        
        //   - force
        f[i] = R[i+1]*f[i+1] + m[i]*pcdd[i]; //equ. (7.112)
        
        //   - u
        u[i] = -1.0*hlao.matrix_multiplication(f[i],(r[i]+rc[i])) + R[i+1]*u[i+1] + hlao.matrix_multiplication(R[i+1]*f[i+1],rc[i]) +
                   I[i]*wd[i] + hlao.matrix_multiplication(w[i],(I[i]*w[i])) +
                   hlao.matrix_multiplication(w[i],(I[i]*w[i])) + kr[i+1]*qdd[i+1]*Im[i+1]*zm[i+1] +
                   hlao.matrix_multiplication(kr[i+1]*qd[i+1]*Im[i+1]*w[i],zm[i+1]); //equ. (7.113)
        
        //   - torque
        T[i] = matrix_transpose(u[i])*matrix_transpose(R[i])*z[0] + kr[i]*Im[i]*matrix_transpose(wmd[i])*zm(i) +
                   Fv[i]*vd[i] + Fs[i]*sgn(vd[i]); //equ. (7.114)
    }
}
*/

// Example 7.5.3, REF: Robotics Modelling, Planning and Control, Page 289
// Newton-Euler formulation
function twoLinkplanarArmNE(v,vd,vdd){
    //var v1 = -1.0*Math.PI/2.0; //rad
    //var v2 = Math.PI; //rad
    //var v1d = 0.0; //rad/s
    //var v2d = 0.0; //rad/s
    //var v1dd = 26.0; //rad.s^-2
    //var v2dd = 26.0; //rad.s^-2
    var v1 = v[0]; //rad
    var v2 = v[1]; //rad
    var v1d = vd[0]; //rad/s
    var v2d = vd[1]; //rad/s
    var v1dd = vdd[0]; //rad.s^-2
    var v2dd = vdd[1]; //rad.s^-2
    
    var a1 = 1.0; //m
    var a2 = 1.0; //m
    var l1 = 0.5; //m
    var l2 = 0.5; //m
    var Il1 = 10.0; //kg.m^2
    var Il2 = 10.0; //kg.m^2
    var ml1 = 50.0; //kg
    var ml2 = 50.0; //kg
    var Im1 = 0.01; //kg.m^2
    var Im2 = 0.01; //kg.m^2
    var mm1 = 5.0; //kg
    var mm2 = 5.0; //kg
    var kr1 = 100.0;
    var kr2 = 100.0;
    var g = 9.81;
    
    //others
    var m1 = ml1 + mm2;
    var m2 = ml2;
    var lc1 = ml1*(l1-a1)/m1;
    var lc2 = ml2*(l2-a2)/m2;
    var I1zz = Il1+ml1*Math.pow((l1-a1),2)+Im2-m1*Math.pow(lc1,2);
    var I2zz = Il2+ml2*Math.pow((l2-a2),2)-m2*Math.pow(lc2,2);
    
    var c1 = Math.cos(v1);
    var c2 = Math.cos(v2);
    var c12 = Math.cos(v1+v2);
    var s1 = Math.sin(v1);
    var s2 = Math.sin(v2);
    var s12 = Math.sin(v1+v2);
    
    //forward link 1:
    var w11 = [[0.0],[0.0],[v1d]];
    var w11d = [[0.0],[0.0],[v1dd]];
    var p11dd = [
        [-1.0*a1*Math.pow(v1d,2)+g*s1],
        [                a1*v1dd+g*c1],
        [                         0.0]
    ];
    var p1c1dd = [
        [-1.0*(lc1+a1)*Math.pow(v1d,2)+g*s1],
        [                (lc1+a1)*v1dd+g*c1],
        [                               0.0]
    ];
    var w0m1d = [[0.0],[0.0],[kr1*v1dd]];
    //console.log([[w11],[w11d],[p11dd],[p1c1dd],[w0m1d]]);
    
    //forward link 2:
    var w22 = [[0.0],[0.0],[v1d+v2d]];
    var w22d = [[0.0],[0.0],[v1dd+v2dd]];
    var p22dd = [
        [a1*s2*v1dd-a1*c2*Math.pow(v1d,2)-a2*Math.pow((v1d+v2d),2)+g*s12],
        [          a1*c2*v1dd+a2*(v1dd+v2dd)+a1*s2*Math.pow(v1d,2)+g*c12],
        [                                                            0.0]
    ];
    var p2c2dd = [
        [  a1*s2*v1dd-a1*c2*Math.pow(v1d,2)-(lc2+a2)*Math.pow((v1d+v2d),2)+g*s12],
        [a1*c2*v1dd+(lc2+a2)*(v1dd+v2dd)+a1*s2*Math.pow(v1d,2)+g*c12],
        [                                                        0.0]
    ];
    var w1m2d = [[0.0],[0.0],[v1dd+kr2*v2dd]];
    //console.log([[w22],[w22d],[p22dd],[p2c2dd],[w1m2d]]);
    //console.log(w11,w22);
    //console.log(w11d,w22d);
    //console.log(p11dd,p22dd);
    //console.log(p1c1dd,p2c2dd);
    //console.log(w0m1d,w1m2d);
    
    //backward recursion: link 2
    var f22 = [
        [m2*(a1*s2*v1dd-a1*c2*Math.pow(v1d,2)-(lc2+a2)*Math.pow((v1d+v2d),2)+g*s12)],
        [          m2*(a1*c2*v1dd+(lc2+a2)*(v1dd+v2dd)+a1*s2*Math.pow(v1d,2)+g*c12)],
        [                                                                       0.0]
    ];
    var u22 = [
        [0.0], //IMPORTANT: has not been computed
        [0.0], //IMPORTANT: has not been computed
        [I2zz*(v1dd+v2dd)+m2*Math.pow((lc2+a2),2)*(v1dd+v2dd)+m2*a1*(lc2+a2)*c2*v1dd+m2*a1*(lc2+a2)*s2*Math.pow(v1d,2)+m2*(lc2+a2)*g*c12]
    ];
    var T2 = (I2zz+m2*(Math.pow((lc2+a2),2)+a1*(lc2+a2)*c2)+kr2*Im2)*v1dd+(I2zz+m2*Math.pow((lc2+a2),2)+Math.pow(kr2,2)*Im2)*v2dd+m2*a1*(lc2+a2)*s2*Math.pow(v1d,2)+m2*(lc2+a2)*g*c12;
    //console.log([[f22],[u22],[T2]]);
    
    //backward recursion: link 1
    var f11 = [
        [-1.0*m2*(lc2+a2)*s2*(v1dd+v2dd)-m1*(lc1+a1)*Math.pow(v1d,2)-m2*a1*Math.pow(v1d,2)-m2*(lc2+a2)*c2*Math.pow((v1d+v2d),2)+(m1+m2)*g*s1],
        [                           m1*(lc1+a1)*v1dd+m2*a1*v1dd+m2*(lc2+a2)*c2*(v1dd+v2dd)-m2*(lc2+a2)*s2*Math.pow((v1d+v2d),2)+(m1+m2)*g*c1],
        [                                                                                                                                0.0]
    ];
    var u11 = [
        [0.0], //IMPORTANT: has not been computed
        [0.0], //IMPORTANT: has not been computed
        [I1zz*v1dd+m2*Math.pow(a1,2)*v1dd+m1*Math.pow((lc1+a1),2)*v1dd+m2*a1*(lc2+a2)*c2*v1dd+I2zz*(v1dd+v2dd)+m2*a1*(lc2+a2)*c2*(v1dd+v2dd)+
         m2*Math.pow((lc2+a2),2)*(v1dd+v2dd)+kr2*Im2*v2dd+m2*a1*(lc2+a2)*s2*Math.pow(v1d,2)-m2*a1*(lc2+a2)*s2*Math.pow((v1d+v2d),2)+
         m1*(lc1+a1)*g*c1+m2*a1*g*c1+m2*(lc2+a2)*g*c12]
    ];
    var T1 = (I1zz+m1*Math.pow((lc1+a1),2)+Math.pow(kr1,2)*Im1+I2zz+m2*(Math.pow(a1,2)+Math.pow((lc2+a2),2)+2.0*a1*(lc2+a2)*c2))*v1dd+
             (I2zz+m2*(Math.pow((lc2+a2),2)+a1*(lc2+a2)*c2)+kr2*Im2)*v2dd-2.0*m2*a1*(lc2+a2)*s2*v1d*v2d-m2*a1*(lc2+a2)*s2*Math.pow(v2d,2)+
             (m1*(lc1+a1)+m2*a1)*g*c1+m2*(lc2+a2)*g*c12;
    //console.log([[f11],[u11],[T1]]);
    console.log(f11,f22);
    console.log(u11,u22);
    
    return [[T1],[T2]];
}

// 7.3.2, REF: Robotics Modelling, Planning and Control, Page 265
// Lagrange Formulation
function twoLinkplanarArmLF(v,vd,vdd){
    /*
    var v1 = -1.0*Math.PI/2.0; //rad
    var v2 = Math.PI; //rad
    var v1d = 0.0; //rad/s
    var v2d = 0.0; //rad/s
    var v1dd = 26.0; //rad.s^-2
    var v2dd = 26.0; //rad.s^-2
    //T = [
    //    3861.3
    //    3146.0
    //]
    */
    var v1 = v[0]; //rad
    var v2 = v[1]; //rad
    var v1d = vd[0]; //rad/s
    var v2d = vd[1]; //rad/s
    var v1dd = vdd[0]; //rad.s^-2
    var v2dd = vdd[1]; //rad.s^-2
    var a1 = 1.0; //m
    var a2 = 1.0; //m
    var l1 = 0.5; //m
    var l2 = 0.5; //m
    var Il1 = 10.0; //kg.m^2
    var Il2 = 10.0; //kg.m^2
    var ml1 = 50.0; //kg
    var ml2 = 50.0; //kg
    var Im1 = 0.01; //kg.m^2
    var Im2 = 0.01; //kg.m^2
    var mm1 = 5.0; //kg
    var mm2 = 5.0; //kg
    var kr1 = 100.0;
    var kr2 = 100.0;
    var g = 9.81;

    var c1  = Math.cos(v1);
    var c2  = Math.cos(v2);
    var c12 = Math.cos(v1 + v2);
    var s2  = Math.sin(v2);

    var T1 = (Il1+ml1*Math.pow(l1,2)+Math.pow(kr1,2)*Im1+Il2+ml2*(Math.pow(a1,2)+Math.pow(l2,2)+2.0*a1*l2*c2)+Im2+mm2*Math.pow(a1,2))*v1dd +
             (Il2+ml2*(Math.pow(l2,2)+a1*l2*c2)+kr2*Im2)*v2dd -
             2.0*ml2*a1*l2*s2*v1d*v2d - ml2*a1*l2*s2*Math.pow(v2d,2) +
             (ml1*l1+mm2*a1+ml2*a1)*g*c1 + ml2*l2*g*c12;
             
    var T2 = (Il2+ml2*(Math.pow(l2,2)+a1*l2*c2)+kr2*Im2)*v1dd+(Il2+ml2*Math.pow(l2,2)+Math.pow(kr2,2)*Im2)*v2dd +
             ml2*a1*l2*s2*Math.pow(v1d,2)+ml2*l2*g*c12;
             
    return [[T1],[T2]];
}

// linearity
// Example 7.2, REF: Robotics Modelling, Planning and Control, Page 269
function twoLinkplanarArm_parameterization(){
    var v1 = -1.0*Math.PI/2.0; //rad
    var v2 = Math.PI; //rad
    var v1d = 0.0; //rad/s
    var v2d = 0.0; //rad/s
    var v1dd = 26.0; //rad.s^-2
    var v2dd = 26.0; //rad.s^-2
    /*
    var v1 = v[0]; //rad
    var v2 = v[1]; //rad
    var v1d = vd[0]; //rad/s
    var v2d = vd[1]; //rad/s
    var v1dd = vdd[0]; //rad.s^-2
    var v2dd = vdd[1]; //rad.s^-2
    */
    var a1 = 1.0; //m
    var a2 = 1.0; //m
    var l1 = 0.5; //m
    var l2 = 0.5; //m
    var Il1 = 10.0; //kg.m^2
    var Il2 = 10.0; //kg.m^2
    var ml1 = 50.0; //kg
    var ml2 = 50.0; //kg
    var Im1 = 0.01; //kg.m^2
    var Im2 = 0.01; //kg.m^2
    //var mm1 = 5.0; //kg
    var mm2 = 5.0; //kg
    var kr1 = 100.0;
    var kr2 = 100.0;
    var g = 9.81;
    
    var c1  = Math.cos(v1);
    var c2  = Math.cos(v2);
    var c12 = Math.cos(v1 + v2);
    var s2  = Math.sin(v2);
    
    //Parameter vector eq. 7.83
    var pi1 = ml1 + mm2;
    var pi2 = ml1*(l1 - a1);
    var pi3 = Il1 + ml1*Math.pow((l1- a1),2) + Im2;
    var pi4 = Im1;
    var pi5 = ml2;
    var pi6 = ml2*(l2 - a2);
    var pi7 = Il2 + ml2*Math.pow((l2 - a2),2);
    var pi8 = Im2;
    var pi = [[pi1],[pi2],[pi3],[pi4],[pi5],[pi6],[pi7],[pi8]]; //eq. 7.83

    //Regressor 7.84
    var y11 = Math.pow(a1,2)*v1dd + a1*g*c1;
    var y12 = 2.0*a1*v1dd + g*c1;
    var y13 = v1dd;
    var y14 = Math.pow(kr1,2)*v1dd;
    var y15 = (Math.pow(a1,2) + 2.0*a1*a2*c2 + Math.pow(a2,2))*v1dd + (a1*a2*c2 + Math.pow(a2,2))*v2dd - 2.0*a1*a2*s2*v1d*v2d - a1*a2*s2*Math.pow(v2d,2) + a1*g*c1 + a2*g*c12;
    var y16 = (2.0*a1*c2 + 2.0*a2)*v1dd + (a1*c2 + 2.0*a2)*v2dd - 2.0*a1*s2*v1d*v2d - a1*s2*Math.pow(v2d,2) + g*c12;
    var y17 = v1dd + v2dd;
    var y18 = kr2*v2dd;
    var y21 = 0.0;
    var y22 = 0.0;
    var y23 = 0.0;
    var y24 = 0.0;
    var y25 = (a1*a2*c2 + Math.pow(a2,2))*v1dd + Math.pow(a2,2)*v2dd + a1*a2*s2*Math.pow(v1d,2) + a2*g*c12;
    var y26 = (a1*c2 + 2.0*a2)*v1dd + 2.0*a2*v2dd + a1*s2*Math.pow(v1d,2) + g*c12;
    var y27 = v1dd + v2dd;
    var y28 = kr2*v1dd + Math.pow(kr2,2)*v2dd;
    var Y = [
        [y11, y12, y13, y14, y15, y16, y17, y18],
        [y21, y22, y23, y24, y25, y26, y27, y28]
    ];

    //Relation 7.81
    var T = hlao.matrix_multiplication(Y,pi);
    console.log(T);
    //T = [
    //    3861.3
    //    3146.0
    //]
}

export {
    NewtonEulerRecursion,
    //linkAccelerationsBF,
    //linkForcesBF,
    linkAccelerationsCF,
    //linkForcesCF,
    twoLinkplanarArmNE,
    twoLinkplanarArmLF,
    twoLinkplanarArm_parameterization
};
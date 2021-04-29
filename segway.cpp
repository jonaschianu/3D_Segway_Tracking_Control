// File in C++ format
// Application title: Segway Regulation Control
//
//    copyright (C) 2009 Olivier VERLINDEN
//    Service de Mecanique rationnelle, Dynamique et Vibrations
//    Faculte Polytechnique de Mons
//    31, Bd Dolez, 7000 MONS (Belgium)
//    Olivier.Verlinden@fpms.ac.be
//
// This file is part of EasyDyn
//
// EasyDyn is free software; you can redistribute it and/or modify it under the
// terms of the GNU GeneralPublic License as published by the Free Software
// Foundation; either version 2, or (at your option) any later version.
//
// EasyDyn is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// EasyDyn; see the file COPYING.  If not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

//

#define EASYDYNMBSMAIN  // to declare the global variable
#define EASYDYNMBSADVANCED // This option implies that a specific
                           // version of ComputePartialVelocities() is provided
#include <stdio.h>
#include <math.h>
#include <EasyDyn/mbs.h>
#include <EasyDyn/visu.h>
#include <EasyDyn/vec.h>
#include <fstream>

using namespace std;
scene thescene;
ofstream VanFile;

// Global variables

/// Geometrical constants
double L=0.9;
double r=0.2;
double c=0.4;
double a=0.25;
double n=25.0;

/// Electrical constants
double Km =0.085;
double w_left, w_right;
double i_left, i_right;
double R = 0.35;

/// Control
double desired_wz;
double actual_wz;
double error_wz;

double u_steer;
double u_left;
double u_right;

double K_wz = 30.0;

double SetPoint;

/// Time parameters
double t_total = 15.0; // s
double dt=0.01; // s

/// Forward control constants
double CurrentPosition_Z;
double CurrentVelocity_Z;
double CurrentAcceleration_Z;
double Z_Axis_init = 0; // [m]

double TCP_Angular_velocity = 3.5; // [rad/s]
double Duration_Accel = 2.0; // [s] We want the desired velocity after half of a revolution
double TCP_velocity_Z=0.5;
double PostPosition_Z=0.0;
double Duration_Decel = M_PI/(0.5*TCP_Angular_velocity);
double Duration_CstVelocity = 3.0; // For the constant velocity
double Duration_Rest = 9.0;


///************************************///
///  Segment: Linear velocity profile  ///
///************************************///

//---------------------------------------------------

void TrajLinVel_Segment(bool Reverse,double &dx, double desired_vel_TCP, double t_accel, double t,
                        double &xd_end, double &xdd_end, double t_init, double x_init, double &x_final)
{
    // Compute position profile so that the velocity has trapezoidal shape
    // No blend between segments

    // dx              : computed position increment to be added to the initial position
    // desired_vel_TCP : desired TCP velocity in m/s
    // t_accel         : time duration during which acceleration is constant
    // t               : current time of EasyDdyn simulation
    // xd_end          : computed TCP velocity (theoretical)
    // xdd_end         : computed TCP acceleration (theoretical)
    // t_init          : initial time when calling this function for the first time
    // x_init          : initial position of the end effector

    // Time parameter
    double t_linear_acc = t_accel;

    // End effector velocity
    double FeedRate_SI = desired_vel_TCP;

    // Initial acceleration
    double a = desired_vel_TCP/t_accel;

    // Current variables to compute the kinematic profile
    double Acc_t;
    double Vel_t;
    double Position_t;
    double t_phase;

    // Precomputation of the constants (for linear velocity and linear deceleration)
    // constants for linear velocity
    double a0_phase1;
    if(Reverse == false)
    {
        a0_phase1=a;
    }
    else
    {
        a0_phase1=-a;
    }

    //double a0_phase1=a;

    double v0_phase1=0.0;
    double x0_phase1= x_init;

    // Computation of increments
    if( ((t-t_init)<t_accel) && ((t-t_init)>=0))
    {
        // phase 1: constant acceleration
        if( (t-t_init)<=t_linear_acc)
        {
            t_phase = t-t_init;
            Acc_t = a0_phase1;
            Vel_t = v0_phase1 + a0_phase1*t_phase;
            Position_t = x0_phase1 + v0_phase1*t_phase + (1.0/2.0)*a0_phase1*pow(t_phase,2);
        }
    }
    else if( (t-t_init)<0.0 )
    {
        Acc_t = xdd_end;
        Vel_t= xd_end;
        Position_t = dx;
    }
    else
    {
        Acc_t = 0.0;
        Vel_t= 0.0;
        Position_t = 0.0;
    }

    // Final position after duration of constant jerk
    x_final = x0_phase1 + v0_phase1*t_accel + (1.0/2.0)*a0_phase1*pow(t_accel,2);

    // save current acceleration velocity and position
    xdd_end = Acc_t;
    xd_end = Vel_t;
    dx = Position_t;
}

///************************************///
///  Segment: Constant velocity profile  ///
///************************************///

//---------------------------------------------------

void TrajConstVel_Segment(double &dx, double desired_vel_TCP, double t_CstVel, double t,
                        double &xd_end, double &xdd_end, double t_init, double x_init, double &x_final)
{
    // Constant velocity over a particular duration t_CstVel

    // dx              : computed position increment of the TCP to be added to the initial position of the TCP
    // desired_vel_TCP : desired TCP velocity in m/s
    // t_CstVel         : time duration during which velocity of TCP is constant
    // t               : current time of EasyDdyn simulation
    // xd_end          : computed TCP velocity (theoretical)
    // xdd_end         : computed TCP acceleration (theoretical)
    // t_init          : initial time when calling this function for the first time
    // x_init          : initial position of the end effector

    // endtraj: last segment for rest, keep same values

    // End effector velocity
    double FeedRate_SI = desired_vel_TCP;

    // Current variables to compute the kinematic profile
    // double Jerk_t
    double Acc_t;
    double Vel_t;
    double Position_t;
    double t_phase;

    // phase 1: constant velocity
    if( ((t-t_init)<t_CstVel) && ((t-t_init)>=0)) // Constant velocity
    {

        t_phase = t-t_init;
        Acc_t = 0.0;
        Vel_t = FeedRate_SI;
        Position_t = x_init + FeedRate_SI*t_phase;

    }
    else if( (t-t_init)<0.0 )
    {
        Acc_t = xdd_end;
        Vel_t= xd_end;
        Position_t = dx;
    }
    else
    {
        Acc_t = 0.0;
        Vel_t= 0.0;
        Position_t = 0.0;
    }

    // Final position after duration of constant jerk
    x_final = x_init + FeedRate_SI*t_CstVel;

    // save current acceleration velocity and position
    xdd_end = Acc_t;
    xd_end = Vel_t;
    dx = Position_t;
}

//---------------------------------------------------

//---------------------------------------------------

void WriteDataHeader(ostream &OutFile)
{
WriteStateVariablesHeader(OutFile);
OutFile << endl;
}

//----------------------------------------------------

void SaveData(ostream &OutFile)
{
/// Forward velocity
/// Controller (same as in 2D)
u[0] = -(-67.9014*q[0] - 331.0229*q[4] - 86.1153*qd[0] - 120.3355*qd[4] + 26.0448*qd[8]); // LQR

SaveStateVariables(OutFile);
OutFile << endl;
thescene.WriteCoord(VanFile);
}

//----------------------------------------------------

void SetInertiaData()
{
body[0].mass=1.0e-10;
body[0].PhiG.Ixx=1.0e-10;
body[0].PhiG.Iyy=1.0e-10;
body[0].PhiG.Izz=1.0e-10;
body[0].PhiG.Ixy=0;
body[0].PhiG.Ixz=0;
body[0].PhiG.Iyz=0;
body[1].mass=90;
body[1].PhiG.Ixx=21.600000000000001;
body[1].PhiG.Iyy=21.600000000000001;
body[1].PhiG.Izz=2.5;
body[1].PhiG.Ixy=0;
body[1].PhiG.Ixz=0;
body[1].PhiG.Iyz=0;
body[2].mass=15;
body[2].PhiG.Ixx=0.20000000000000001;
body[2].PhiG.Iyy=0.29999999999999999;
body[2].PhiG.Izz=0.20000000000000001;
body[2].PhiG.Ixy=0;
body[2].PhiG.Ixz=0;
body[2].PhiG.Iyz=0;
body[3].mass=15;
body[3].PhiG.Ixx=0.20000000000000001;
body[3].PhiG.Iyy=0.29999999999999999;
body[3].PhiG.Izz=0.20000000000000001;
body[3].PhiG.Ixy=0;
body[3].PhiG.Ixz=0;
body[3].PhiG.Iyz=0;
body[4].mass=1;
body[4].PhiG.Ixx=0.00012;
body[4].PhiG.Iyy=0.0001;
body[4].PhiG.Izz=0.00012;
body[4].PhiG.Ixy=0;
body[4].PhiG.Ixz=0;
body[4].PhiG.Iyz=0;
body[5].mass=1;
body[5].PhiG.Ixx=0.00012;
body[5].PhiG.Iyy=0.0001;
body[5].PhiG.Izz=0.00012;
body[5].PhiG.Ixy=0;
body[5].PhiG.Ixz=0;
body[5].PhiG.Iyz=0;
}

//-----------------------------------------

void ComputeMotion()
{
//Homogenous transformation matrices of each body
//Insert kinematics generated by python/sympy

//Body[0]
body[0].T0G.R.r11=1;
body[0].T0G.R.r12=0;
body[0].T0G.R.r13=0;
body[0].T0G.R.r21=0;
body[0].T0G.R.r22=1;
body[0].T0G.R.r23=0;
body[0].T0G.R.r31=0;
body[0].T0G.R.r32=0;
body[0].T0G.R.r33=1;
body[0].T0G.e.x=0;
body[0].T0G.e.y=0;
body[0].T0G.e.z=0;
//Body[1]
body[1].T0G.R.r11=cos(q[4])*cos(q[5]);
body[1].T0G.R.r12=-sin(q[5])*cos(q[4]);
body[1].T0G.R.r13=sin(q[4]);
body[1].T0G.R.r21=sin(q[3])*sin(q[4])*cos(q[5]) + sin(q[5])*cos(q[3]);
body[1].T0G.R.r22=-sin(q[3])*sin(q[4])*sin(q[5]) + cos(q[3])*cos(q[5]);
body[1].T0G.R.r23=-sin(q[3])*cos(q[4]);
body[1].T0G.R.r31=sin(q[3])*sin(q[5]) - sin(q[4])*cos(q[3])*cos(q[5]);
body[1].T0G.R.r32=sin(q[3])*cos(q[5]) + sin(q[4])*sin(q[5])*cos(q[3]);
body[1].T0G.R.r33=cos(q[3])*cos(q[4]);
body[1].T0G.e.x=q[0] + 0.90000000000000002*sin(q[4]);
body[1].T0G.e.y=q[1] - 0.90000000000000002*sin(q[3])*cos(q[4]);
body[1].T0G.e.z=q[2] + 0.90000000000000002*cos(q[3])*cos(q[4]) + 0.20000000000000001;
body[1].vG.x=qd[0] + 0.90000000000000002*qd[4]*cos(q[4]);
body[1].vG.y=qd[1] - 0.90000000000000002*qd[3]*cos(q[3])*cos(q[4]) + 0.90000000000000002*qd[4]*sin(q[3])*sin(q[4]);
body[1].vG.z=qd[2] - 0.90000000000000002*qd[3]*sin(q[3])*cos(q[4]) - 0.90000000000000002*qd[4]*sin(q[4])*cos(q[3]);
body[1].aG.x=-0.90000000000000002*pow(qd[4], 2)*sin(q[4]) + qdd[0] + 0.90000000000000002*qdd[4]*cos(q[4]);
body[1].aG.y=0.90000000000000002*pow(qd[3], 2)*sin(q[3])*cos(q[4]) + 1.8*qd[3]*qd[4]*sin(q[4])*cos(q[3]) + 0.90000000000000002*pow(qd[4], 2)*sin(q[3])*cos(q[4]) + qdd[1] - 0.90000000000000002*qdd[3]*cos(q[3])*cos(q[4]) + 0.90000000000000002*qdd[4]*sin(q[3])*sin(q[4]);
body[1].aG.z=-0.90000000000000002*pow(qd[3], 2)*cos(q[3])*cos(q[4]) + 1.8*qd[3]*qd[4]*sin(q[3])*sin(q[4]) - 0.90000000000000002*pow(qd[4], 2)*cos(q[3])*cos(q[4]) + qdd[2] - 0.90000000000000002*qdd[3]*sin(q[3])*cos(q[4]) - 0.90000000000000002*qdd[4]*sin(q[4])*cos(q[3]);
body[1].omega.x=qd[3] + qd[5]*sin(q[4]);
body[1].omega.y=qd[4]*cos(q[3]) - qd[5]*sin(q[3])*cos(q[4]);
body[1].omega.z=qd[4]*sin(q[3]) + qd[5]*cos(q[3])*cos(q[4]);
body[1].omegad.x=qd[4]*qd[5]*cos(q[4]) + qdd[3] + qdd[5]*sin(q[4]);
body[1].omegad.y=-qd[3]*qd[4]*sin(q[3]) - qd[3]*qd[5]*cos(q[3])*cos(q[4]) + qd[4]*qd[5]*sin(q[3])*sin(q[4]) + qdd[4]*cos(q[3]) - qdd[5]*sin(q[3])*cos(q[4]);
body[1].omegad.z=qd[3]*qd[4]*cos(q[3]) - qd[3]*qd[5]*sin(q[3])*cos(q[4]) - qd[4]*qd[5]*sin(q[4])*cos(q[3]) + qdd[4]*sin(q[3]) + qdd[5]*cos(q[3])*cos(q[4]);
//Body[2]
body[2].TrefG.R.r11=cos(q[6]);
body[2].TrefG.R.r12=0;
body[2].TrefG.R.r13=sin(q[6]);
body[2].TrefG.R.r21=0;
body[2].TrefG.R.r22=1;
body[2].TrefG.R.r23=0;
body[2].TrefG.R.r31=-sin(q[6]);
body[2].TrefG.R.r32=0;
body[2].TrefG.R.r33=cos(q[6]);
body[2].TrefG.e.x=0;
body[2].TrefG.e.y=-0.40000000000000002;
body[2].TrefG.e.z=-0.90000000000000002;
body[2].omegarel.y=qd[6];
body[2].omegadrel.y=qdd[6];
ComposeMotion(2,1);
//Body[3]
body[3].TrefG.R.r11=cos(q[7]);
body[3].TrefG.R.r12=0;
body[3].TrefG.R.r13=sin(q[7]);
body[3].TrefG.R.r21=0;
body[3].TrefG.R.r22=1;
body[3].TrefG.R.r23=0;
body[3].TrefG.R.r31=-sin(q[7]);
body[3].TrefG.R.r32=0;
body[3].TrefG.R.r33=cos(q[7]);
body[3].TrefG.e.x=0;
body[3].TrefG.e.y=0.40000000000000002;
body[3].TrefG.e.z=-0.90000000000000002;
body[3].omegarel.y=qd[7];
body[3].omegadrel.y=qdd[7];
ComposeMotion(3,1);
//Body[4]
body[4].TrefG.R.r11=cos(25*q[4] - 25*q[6]);
body[4].TrefG.R.r12=0;
body[4].TrefG.R.r13=-sin(25*q[4] - 25*q[6]);
body[4].TrefG.R.r21=0;
body[4].TrefG.R.r22=1;
body[4].TrefG.R.r23=0;
body[4].TrefG.R.r31=sin(25*q[4] - 25*q[6]);
body[4].TrefG.R.r32=0;
body[4].TrefG.R.r33=cos(25*q[4] - 25*q[6]);
body[4].TrefG.e.x=0;
body[4].TrefG.e.y=-0.25;
body[4].TrefG.e.z=-0.90000000000000002;
body[4].omegarel.y=-25*qd[4] + 25*qd[6];
body[4].omegadrel.y=-25*qdd[4] + 25*qdd[6];
ComposeMotion(4,1);
//Body[5]
body[5].TrefG.R.r11=cos(25*q[4] - 25*q[7]);
body[5].TrefG.R.r12=0;
body[5].TrefG.R.r13=-sin(25*q[4] - 25*q[7]);
body[5].TrefG.R.r21=0;
body[5].TrefG.R.r22=1;
body[5].TrefG.R.r23=0;
body[5].TrefG.R.r31=sin(25*q[4] - 25*q[7]);
body[5].TrefG.R.r32=0;
body[5].TrefG.R.r33=cos(25*q[4] - 25*q[7]);
body[5].TrefG.e.x=0;
body[5].TrefG.e.y=0.25;
body[5].TrefG.e.z=-0.90000000000000002;
body[5].omegarel.y=-25*qd[4] + 25*qd[7];
body[5].omegadrel.y=-25*qdd[4] + 25*qdd[7];
ComposeMotion(5,1);
}

//-----------------------------------------

void ComputePartialVelocities()
{
//Body[0]
//Body[1]
body[1].vGpartial[0].x=1;
body[1].vGpartial[1].y=1;
body[1].vGpartial[2].z=1;
body[1].vGpartial[3].y=-0.90000000000000002*cos(q[3])*cos(q[4]);
body[1].vGpartial[3].z=-0.90000000000000002*sin(q[3])*cos(q[4]);
body[1].omegapartial[3].x=1;
body[1].vGpartial[4].x=0.90000000000000002*cos(q[4]);
body[1].vGpartial[4].y=0.90000000000000002*sin(q[3])*sin(q[4]);
body[1].vGpartial[4].z=-0.90000000000000002*sin(q[4])*cos(q[3]);
body[1].omegapartial[4].y=cos(q[3]);
body[1].omegapartial[4].z=sin(q[3]);
body[1].omegapartial[5].x=sin(q[4]);
body[1].omegapartial[5].y=-sin(q[3])*cos(q[4]);
body[1].omegapartial[5].z=cos(q[3])*cos(q[4]);
//Body[2]
body[2].omegarelpartial[6].y=1;
ComposePartialVelocities(2,1);
//Body[3]
body[3].omegarelpartial[7].y=1;
ComposePartialVelocities(3,1);
//Body[4]
body[4].omegarelpartial[4].y=-25;
body[4].omegarelpartial[6].y=25;
ComposePartialVelocities(4,1);
//Body[5]
body[5].omegarelpartial[4].y=-25;
body[5].omegarelpartial[7].y=25;
ComposePartialVelocities(5,1);
}

//-----------------------------------------

void AddAppliedEfforts()
{
//Contribution of gravity
vec gravity(0,0,-9.8100000000000005);
AddGravityForces(gravity);
//Contribution of user defined forces


/// Steering command

/// Desired steering
if( t> 4.5 && t<5.0) {
    desired_wz = 3.0; /// rad/s
} else if( t> 6.0 && t<6.5) {
    desired_wz = -3.0; /// rad/s
} else {
    desired_wz = 0.0;
}

/// Error on angular velocity

actual_wz = (a/r) * (qd[7]-qd[6]);
error_wz = desired_wz - actual_wz;

/// Input for each motor
u_steer= K_wz*error_wz;

u_left = u[0] - u_steer;
u_right = u[0] + u_steer;

/// Saturation
if(u_left>72.0) u[0] = 72.0;
if(u_left<-72.0) u[0] = -72.0;

if(u_right>72.0) u[0] = 72.0;
if(u_right<-72.0) u[0] = -72.0;

/// Motors

w_left=n*(qd[6]-qd[4]);
w_right=n*(qd[7]-qd[4]);

i_left = (u_left-Km*w_left)/R;
i_right = (u_right-Km*w_right)/R;

body[4].MG += (Km*i_left)*body[4].T0G.R.uy(); // Torque on Rotor 1
body[1].MG -= (Km*i_left)*body[4].T0G.R.uy();

body[5].MG += (Km*i_right)*body[5].T0G.R.uy(); // Torque on Rotor 2
body[1].MG -= (Km*i_right)*body[5].T0G.R.uy();

/// Tyre
structtyre tyre_segway;
tyre_segway.r1=0.2;
tyre_segway.r2=0.03;
tyre_segway.Kz=100000;
tyre_segway.Cz=600;
tyre_segway.Fznom=700;
tyre_segway.Clongnom=10000;
tyre_segway.nlong=0.1;
tyre_segway.Clatnom=10000;
tyre_segway.nlat=0.1;
tyre_segway.Ccambernom=1000;
tyre_segway.ncamber=0.1;
tyre_segway.fClbs=1.0;
tyre_segway.fClbd=0.8;

// Left Wheel
AddTyreEfforts(2,vcoord(0,1,0),tyre_segway);
// Right Wheel
AddTyreEfforts(3,vcoord(0,1,0),tyre_segway);

}

//-----------------------------------------

void ComputeResidual()
{
// Definition of Set point (trapezoidal profile for velocity)
TrajLinVel_Segment(false,CurrentPosition_Z, TCP_velocity_Z, Duration_Accel, t,
                   CurrentVelocity_Z, CurrentAcceleration_Z, 0.0, Z_Axis_init, PostPosition_Z);

TrajConstVel_Segment(CurrentPosition_Z, TCP_velocity_Z, Duration_CstVelocity, t,
                   CurrentVelocity_Z, CurrentAcceleration_Z, Duration_Accel, PostPosition_Z,PostPosition_Z);

TrajLinVel_Segment(true,CurrentPosition_Z, TCP_velocity_Z, Duration_Decel, t,
                   CurrentVelocity_Z, CurrentAcceleration_Z, Duration_Accel+Duration_CstVelocity, PostPosition_Z,PostPosition_Z);

TrajConstVel_Segment(CurrentPosition_Z, 0.0, Duration_Rest+dt, t,
                   CurrentVelocity_Z, CurrentAcceleration_Z, Duration_Accel+Duration_CstVelocity+Duration_Decel, PostPosition_Z,PostPosition_Z);

ComputeResidualmbs();

SetPoint = CurrentPosition_Z; /// Position profile
f[8]=qdd[8]-(SetPoint - ((q[6]+q[7])/2.0 * r) );
}

//-----------------------------------------

int main()
{
// Initialization and memory allocation
nbrdof=8+1;
nbrdep=0;
nbrbody=6;

nbrinput=1; // u[0]

application=new char[8];
strcpy(application,"segway");
InitEasyDynmbs();
// Let's create the shapes
shape *s0;
s0=new box(&(body[0].T0G),vcoord(-1,-1.5,-0.01),vcoord(3,1.5,0.0),0,6);
thescene.AddShape(s0);
shape *s1;
s1=new box(&(body[1].T0G),vcoord(-0.1,-0.1,-L-a/4),vcoord(0.1,0.1,L),9,9);
thescene.AddShape(s1);
shape *s12;
s12=new box(&(body[1].T0G),vcoord(-a/2,-a,-L-a/4),vcoord(a/2,a,-L-a/2),9,9);
thescene.AddShape(s12);
shape *s13;
s13=new sphere(&(body[1].T0G),vcoord(0,0,L),L/5,12,12,13,13);
thescene.AddShape(s13);
shape *s2;
s2=new frustum(&(body[2].T0G),vcoord(0,-r/4,0),vcoord(0,r/4,0),r,r,24,0,8);
thescene.AddShape(s2);
shape *s3;
s3=new frustum(&(body[3].T0G),vcoord(0,-r/4,0),vcoord(0,r/4,0),r,r,24,0,8);
thescene.AddShape(s3);
shape *s4;
s4=new frustum(&(body[4].T0G),vcoord(0,-a/2,0),vcoord(0,a/2,0),a/4,a/4,6,0,8);
thescene.AddShape(s4);
shape *s5;
s5=new frustum(&(body[5].T0G),vcoord(0,-a/2,0),vcoord(0,a/2,0),a/4,a/4,6,0,8);
thescene.AddShape(s5);
// Uncomment the following line if you want a moving observer
// thescene.SetVisuFrame(&Tref);
// Let's open an animation file
VanFile.open("segway.van");
// Initial configuration
// Let's save the structure of the scene
ComputeMotion();
thescene.CreateVolFile("segway.vol");
// Searching for equilibrium position
StaticEquilibrium();

// Let's calculate the poles
cout << "Eigen value analysis" << endl;
SaveLinearizedSystem();
ComputePoles();
CreateVmoFile(thescene);
cout << "Eigen values computed !" << endl;

// Let's perform the integration !
NewmarkIntegration(t_total,dt,1e-05);
// The clean way to finish
EndEasyDynmbs();
VanFile.close();

}
//-----------------------------------------

/*
 * ARDroneTTSim_data.c
 *
 * Code generation for model "ARDroneTTSim".
 *
 * Model version              : $Id: UAV_SIL.mdl 912 2012-10-10 23:11:42Z joh07594 $
 * Simulink Coder version : 9.1 (R2019a) 23-Nov-2018
 * C source code generated on : Wed Dec 21 00:02:56 2022
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Freescale->32-bit PowerPC
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "ARDroneTTSim.h"
#include "ARDroneTTSim_private.h"

/* Block parameters (default storage) */
P_ARDroneTTSim_T ARDroneTTSim_P = {
  /* Variable: K
   * Referenced by: '<S3>/Constant1'
   */
  { 1.4142135623730909, 9.3343423553306474E-17, 1.54158548861486E-17,
    2.6148651629130455E-17, 1.4142135623730905, 7.13622527615965E-17,
    -2.3112660784475402E-17, -4.0252913459767436E-17, 1.4142135623730909,
    4.7779103303375372, 2.5331769282861866E-16, -7.3874727532158141E-17,
    8.1645786651949584E-16, 4.777910330337539, -1.0803405761460532E-16,
    -2.4963973794448446E-16, -4.866579594391248E-16, 4.7779103303375408 },

  /* Variable: k_w
   * Referenced by: '<S3>/Constant'
   */
  1.5,

  /* Variable: timeDelay
   * Referenced by: '<S1>/time delay'
   */
  0.12,

  /* Expression: 1
   * Referenced by: '<Root>/commands enabled'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<Root>/commands disabled'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S15>/Gain2'
   */
  1.0,

  /* Expression: -0.5
   * Referenced by: '<S11>/Gain'
   */
  -0.5,

  /* Expression: 1.5
   * Referenced by: '<S12>/proportional control gain - yaw'
   */
  1.5,

  /* Expression: 1
   * Referenced by: '<S13>/proportional control gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S15>/Gain3'
   */
  1.0,

  /* Expression: 0.4
   * Referenced by: '<S14>/Gain1'
   */
  0.4,

  /* Computed Parameter: StateSpace_A
   * Referenced by: '<S1>/State-Space'
   */
  { -4.2683003477481094, -3.1716278132782731, 4.0 },

  /* Computed Parameter: StateSpace_B
   * Referenced by: '<S1>/State-Space'
   */
  2.0,

  /* Computed Parameter: StateSpace_C
   * Referenced by: '<S1>/State-Space'
   */
  { 0.74169067532173338, 0.44045524009948395 },

  /* Expression: 0
   * Referenced by: '<S1>/State-Space'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S1>/gain1'
   */
  1.0,

  /* Expression: 180/pi
   * Referenced by: '<S1>/deg 2 rad1'
   */
  57.295779513082323,

  /* Computed Parameter: StateSpace1_A
   * Referenced by: '<S1>/State-Space1'
   */
  { -3.9783762642064673, -2.979599735183454, 4.0 },

  /* Computed Parameter: StateSpace1_B
   * Referenced by: '<S1>/State-Space1'
   */
  2.0,

  /* Computed Parameter: StateSpace1_C
   * Referenced by: '<S1>/State-Space1'
   */
  { 1.2569330751058592, 0.60825727527713747 },

  /* Expression: 0
   * Referenced by: '<S1>/State-Space1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S1>/gain'
   */
  1.0,

  /* Expression: 180/pi
   * Referenced by: '<S1>/deg 2 rad'
   */
  57.295779513082323,

  /* Computed Parameter: StateSpace4_A
   * Referenced by: '<S1>/State-Space4'
   */
  -0.0058794721680327112,

  /* Computed Parameter: StateSpace4_B
   * Referenced by: '<S1>/State-Space4'
   */
  1.0,

  /* Computed Parameter: StateSpace4_C
   * Referenced by: '<S1>/State-Space4'
   */
  1.2652682582028261,

  /* Expression: 0
   * Referenced by: '<S1>/State-Space4'
   */
  0.0,

  /* Expression: 180/pi
   * Referenced by: '<S1>/deg 2 rad2'
   */
  57.295779513082323,

  /* Expression: pi/180
   * Referenced by: '<Root>/deg 2 rad'
   */
  0.017453292519943295,

  /* Computed Parameter: StateSpace3_A
   * Referenced by: '<S1>/State-Space3'
   */
  -0.665040528368003,

  /* Computed Parameter: StateSpace3_B
   * Referenced by: '<S1>/State-Space3'
   */
  2.0,

  /* Computed Parameter: StateSpace3_C
   * Referenced by: '<S1>/State-Space3'
   */
  -3.0771603323034462,

  /* Expression: 0
   * Referenced by: '<S1>/State-Space3'
   */
  0.0,

  /* Computed Parameter: StateSpace2_A
   * Referenced by: '<S1>/State-Space2'
   */
  -0.45955158912909533,

  /* Computed Parameter: StateSpace2_B
   * Referenced by: '<S1>/State-Space2'
   */
  2.0,

  /* Computed Parameter: StateSpace2_C
   * Referenced by: '<S1>/State-Space2'
   */
  2.386838950161883,

  /* Expression: 0
   * Referenced by: '<S1>/State-Space2'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S1>/The ARDrone sends zero for the vertical velocity.  '
   */
  0.0,

  /* Expression: [0 0]
   * Referenced by: '<S19>/Constant1'
   */
  { 0.0, 0.0 },

  /* Computed Parameter: StateSpace5_A
   * Referenced by: '<S1>/State-Space5'
   */
  { -5.81998044867891, -3.6045888405166024E-6, 3.814697265625E-6 },

  /* Computed Parameter: StateSpace5_B
   * Referenced by: '<S1>/State-Space5'
   */
  1024.0,

  /* Computed Parameter: StateSpace5_C
   * Referenced by: '<S1>/State-Space5'
   */
  { 0.00014907141964484375, 1319.1453493750244 },

  /* Expression: 0
   * Referenced by: '<S1>/State-Space5'
   */
  0.0,

  /* Expression: 180/pi
   * Referenced by: '<S5>/deg 2 rad1'
   */
  57.295779513082323,

  /* Expression: 0
   * Referenced by: '<Root>/RTrans4'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S1>/Saturation 1'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S1>/Saturation 1'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<Root>/RTrans3'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S1>/Saturation 2'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S1>/Saturation 2'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<Root>/RTrans2'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S1>/Saturation 3'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S1>/Saturation 3'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<Root>/RTrans1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S1>/Saturation 4'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S1>/Saturation 4'
   */
  -1.0,

  /* Expression: 10
   * Referenced by: '<S1>/total communication time delay'
   */
  10.0,

  /* Expression: [0 0 0 0]
   * Referenced by: '<S1>/total communication time delay'
   */
  { 0.0, 0.0, 0.0, 0.0 },

  /* Expression: 0
   * Referenced by: '<S8>/Constant'
   */
  0.0,

  /* Computed Parameter: DiscreteTimeIntegrator_gainval
   * Referenced by: '<S3>/Discrete-Time Integrator'
   */
  0.03,

  /* Expression: [0;0;0;0]
   * Referenced by: '<S3>/Discrete-Time Integrator'
   */
  { 0.0, 0.0, 0.0, 0.0 },

  /* Computed Parameter: ManualSwitch2_CurrentSetting
   * Referenced by: '<Root>/Manual Switch2'
   */
  0U,

  /* Computed Parameter: ManualSwitch3_CurrentSetting
   * Referenced by: '<S3>/Manual Switch3'
   */
  0U
};

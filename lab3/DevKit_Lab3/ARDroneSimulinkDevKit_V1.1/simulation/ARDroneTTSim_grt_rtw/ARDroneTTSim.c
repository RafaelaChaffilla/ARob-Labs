/*
 * ARDroneTTSim.c
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

/* Block signals (default storage) */
B_ARDroneTTSim_T ARDroneTTSim_B;

/* Continuous states */
X_ARDroneTTSim_T ARDroneTTSim_X;

/* Block states (default storage) */
DW_ARDroneTTSim_T ARDroneTTSim_DW;

/* Real-time model */
RT_MODEL_ARDroneTTSim_T ARDroneTTSim_M_;
RT_MODEL_ARDroneTTSim_T *const ARDroneTTSim_M = &ARDroneTTSim_M_;
static void rate_monotonic_scheduler(void);

/*
 * Time delay interpolation routine
 *
 * The linear interpolation is performed using the formula:
 *
 *          (t2 - tMinusDelay)         (tMinusDelay - t1)
 * u(t)  =  ----------------- * u1  +  ------------------- * u2
 *              (t2 - t1)                  (t2 - t1)
 */
real_T rt_TDelayInterpolate(
  real_T tMinusDelay,                 /* tMinusDelay = currentSimTime - delay */
  real_T tStart,
  real_T *tBuf,
  real_T *uBuf,
  int_T bufSz,
  int_T *lastIdx,
  int_T oldestIdx,
  int_T newIdx,
  real_T initOutput,
  boolean_T discrete,
  boolean_T minorStepAndTAtLastMajorOutput)
{
  int_T i;
  real_T yout, t1, t2, u1, u2;

  /*
   * If there is only one data point in the buffer, this data point must be
   * the t= 0 and tMinusDelay > t0, it ask for something unknown. The best
   * guess if initial output as well
   */
  if ((newIdx == 0) && (oldestIdx ==0 ) && (tMinusDelay > tStart))
    return initOutput;

  /*
   * If tMinusDelay is less than zero, should output initial value
   */
  if (tMinusDelay <= tStart)
    return initOutput;

  /* For fixed buffer extrapolation:
   * if tMinusDelay is small than the time at oldestIdx, if discrete, output
   * tailptr value,  else use tailptr and tailptr+1 value to extrapolate
   * It is also for fixed buffer. Note: The same condition can happen for transport delay block where
   * use tStart and and t[tail] other than using t[tail] and t[tail+1].
   * See below
   */
  if ((tMinusDelay <= tBuf[oldestIdx] ) ) {
    if (discrete) {
      return(uBuf[oldestIdx]);
    } else {
      int_T tempIdx= oldestIdx + 1;
      if (oldestIdx == bufSz-1)
        tempIdx = 0;
      t1= tBuf[oldestIdx];
      t2= tBuf[tempIdx];
      u1= uBuf[oldestIdx];
      u2= uBuf[tempIdx];
      if (t2 == t1) {
        if (tMinusDelay >= t2) {
          yout = u2;
        } else {
          yout = u1;
        }
      } else {
        real_T f1 = (t2-tMinusDelay) / (t2-t1);
        real_T f2 = 1.0 - f1;

        /*
         * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
         */
        yout = f1*u1 + f2*u2;
      }

      return yout;
    }
  }

  /*
   * When block does not have direct feedthrough, we use the table of
   * values to extrapolate off the end of the table for delays that are less
   * than 0 (less then step size).  This is not completely accurate.  The
   * chain of events is as follows for a given time t.  Major output - look
   * in table.  Update - add entry to table.  Now, if we call the output at
   * time t again, there is a new entry in the table. For very small delays,
   * this means that we will have a different answer from the previous call
   * to the output fcn at the same time t.  The following code prevents this
   * from happening.
   */
  if (minorStepAndTAtLastMajorOutput) {
    /* pretend that the new entry has not been added to table */
    if (newIdx != 0) {
      if (*lastIdx == newIdx) {
        (*lastIdx)--;
      }

      newIdx--;
    } else {
      if (*lastIdx == newIdx) {
        *lastIdx = bufSz-1;
      }

      newIdx = bufSz - 1;
    }
  }

  i = *lastIdx;
  if (tBuf[i] < tMinusDelay) {
    /* Look forward starting at last index */
    while (tBuf[i] < tMinusDelay) {
      /* May occur if the delay is less than step-size - extrapolate */
      if (i == newIdx)
        break;
      i = ( i < (bufSz-1) ) ? (i+1) : 0;/* move through buffer */
    }
  } else {
    /*
     * Look backwards starting at last index which can happen when the
     * delay time increases.
     */
    while (tBuf[i] >= tMinusDelay) {
      /*
       * Due to the entry condition at top of function, we
       * should never hit the end.
       */
      i = (i > 0) ? i-1 : (bufSz-1);   /* move through buffer */
    }

    i = ( i < (bufSz-1) ) ? (i+1) : 0;
  }

  *lastIdx = i;
  if (discrete) {
    /*
     * tempEps = 128 * eps;
     * localEps = max(tempEps, tempEps*fabs(tBuf[i]))/2;
     */
    double tempEps = (DBL_EPSILON) * 128.0;
    double localEps = tempEps * fabs(tBuf[i]);
    if (tempEps > localEps) {
      localEps = tempEps;
    }

    localEps = localEps / 2.0;
    if (tMinusDelay >= (tBuf[i] - localEps)) {
      yout = uBuf[i];
    } else {
      if (i == 0) {
        yout = uBuf[bufSz-1];
      } else {
        yout = uBuf[i-1];
      }
    }
  } else {
    if (i == 0) {
      t1 = tBuf[bufSz-1];
      u1 = uBuf[bufSz-1];
    } else {
      t1 = tBuf[i-1];
      u1 = uBuf[i-1];
    }

    t2 = tBuf[i];
    u2 = uBuf[i];
    if (t2 == t1) {
      if (tMinusDelay >= t2) {
        yout = u2;
      } else {
        yout = u1;
      }
    } else {
      real_T f1 = (t2-tMinusDelay) / (t2-t1);
      real_T f2 = 1.0 - f1;

      /*
       * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
       */
      yout = f1*u1 + f2*u2;
    }
  }

  return(yout);
}

time_T rt_SimUpdateDiscreteEvents(
  int_T rtmNumSampTimes, void *rtmTimingData, int_T *rtmSampleHitPtr, int_T
  *rtmPerTaskSampleHits )
{
  rtmSampleHitPtr[1] = rtmStepTask(ARDroneTTSim_M, 1);
  rtmSampleHitPtr[2] = rtmStepTask(ARDroneTTSim_M, 2);
  UNUSED_PARAMETER(rtmNumSampTimes);
  UNUSED_PARAMETER(rtmTimingData);
  UNUSED_PARAMETER(rtmPerTaskSampleHits);
  return(-1);
}

/*
 *   This function updates active task flag for each subrate
 * and rate transition flags for tasks that exchange data.
 * The function assumes rate-monotonic multitasking scheduler.
 * The function must be called at model base rate so that
 * the generated code self-manages all its subrates and rate
 * transition flags.
 */
static void rate_monotonic_scheduler(void)
{
  /* To ensure a deterministic data transfer between two rates,
   * data is transferred at the priority of a fast task and the frequency
   * of the slow task.  The following flags indicate when the data transfer
   * happens.  That is, a rate interaction flag is set true when both rates
   * will run, and false otherwise.
   */

  /* tid 1 shares data with slower tid rate: 2 */
  if (ARDroneTTSim_M->Timing.TaskCounters.TID[1] == 0) {
    ARDroneTTSim_M->Timing.RateInteraction.TID1_2 =
      (ARDroneTTSim_M->Timing.TaskCounters.TID[2] == 0);

    /* update PerTaskSampleHits matrix for non-inline sfcn */
    ARDroneTTSim_M->Timing.perTaskSampleHits[5] =
      ARDroneTTSim_M->Timing.RateInteraction.TID1_2;
  }

  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (ARDroneTTSim_M->Timing.TaskCounters.TID[2])++;
  if ((ARDroneTTSim_M->Timing.TaskCounters.TID[2]) > 5) {/* Sample time: [0.03s, 0.0s] */
    ARDroneTTSim_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 11;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  ARDroneTTSim_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * Output and update for atomic system:
 *    '<S1>/normalize angle  between -pi and pi radians'
 *    '<S12>/normalize angle  between -pi and pi radians'
 *    '<S9>/normalize angle  between -pi and pi radians'
 */
void normalizeanglebetweenpiandp(real_T rtu_angleIn,
  B_normalizeanglebetweenpiandp_T *localB)
{
  boolean_T rEQ0;
  real_T q;

  /* MATLAB Function 'ARDrone Simulation Block/normalize angle  between -pi and pi radians': '<S6>:1' */
  /* '<S6>:1:4' */
  if (rtIsNaN(rtu_angleIn) || rtIsInf(rtu_angleIn)) {
    localB->angleOut = (rtNaN);
  } else if (rtu_angleIn == 0.0) {
    localB->angleOut = 0.0;
  } else {
    localB->angleOut = fmod(rtu_angleIn, 6.2831853071795862);
    rEQ0 = (localB->angleOut == 0.0);
    if (!rEQ0) {
      q = fabs(rtu_angleIn / 6.2831853071795862);
      rEQ0 = (fabs(q - floor(q + 0.5)) <= 2.2204460492503131E-16 * q);
    }

    if (rEQ0) {
      localB->angleOut = 0.0;
    } else {
      if (rtu_angleIn < 0.0) {
        localB->angleOut += 6.2831853071795862;
      }
    }
  }

  if (localB->angleOut > 3.1415926535897931) {
    /* '<S6>:1:5' */
    /* '<S6>:1:6' */
    localB->angleOut += -6.2831853071795862;
  }
}

/* Model output function for TID0 */
void ARDroneTTSim_output0(void)        /* Sample time: [0.0s, 0.0s] */
{
  /* local block i/o variables */
  real_T rtb_StateSpace4;
  real_T Vel_xy_tmp;
  real_T Vel_xy_tmp_0;
  if (rtmIsMajorTimeStep(ARDroneTTSim_M)) {
    /* set solver stop time */
    if (!(ARDroneTTSim_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&ARDroneTTSim_M->solverInfo,
                            ((ARDroneTTSim_M->Timing.clockTickH0 + 1) *
        ARDroneTTSim_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&ARDroneTTSim_M->solverInfo,
                            ((ARDroneTTSim_M->Timing.clockTick0 + 1) *
        ARDroneTTSim_M->Timing.stepSize0 + ARDroneTTSim_M->Timing.clockTickH0 *
        ARDroneTTSim_M->Timing.stepSize0 * 4294967296.0));
    }

    {                                  /* Sample time: [0.0s, 0.0s] */
      rate_monotonic_scheduler();
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(ARDroneTTSim_M)) {
    ARDroneTTSim_M->Timing.t[0] = rtsiGetT(&ARDroneTTSim_M->solverInfo);
  }

  /* Gain: '<S1>/gain1' incorporates:
   *  StateSpace: '<S1>/State-Space'
   */
  ARDroneTTSim_B.gain1 = (ARDroneTTSim_P.StateSpace_C[0] *
    ARDroneTTSim_X.StateSpace_CSTATE[0] + ARDroneTTSim_P.StateSpace_C[1] *
    ARDroneTTSim_X.StateSpace_CSTATE[1]) * ARDroneTTSim_P.gain1_Gain;

  /* Gain: '<S1>/deg 2 rad1' */
  ARDroneTTSim_B.deg2rad1 = ARDroneTTSim_P.deg2rad1_Gain * ARDroneTTSim_B.gain1;

  /* Gain: '<S1>/gain' incorporates:
   *  StateSpace: '<S1>/State-Space1'
   */
  ARDroneTTSim_B.gain = (ARDroneTTSim_P.StateSpace1_C[0] *
    ARDroneTTSim_X.StateSpace1_CSTATE[0] + ARDroneTTSim_P.StateSpace1_C[1] *
    ARDroneTTSim_X.StateSpace1_CSTATE[1]) * ARDroneTTSim_P.gain_Gain;

  /* Gain: '<S1>/deg 2 rad' */
  ARDroneTTSim_B.deg2rad = ARDroneTTSim_P.deg2rad_Gain * ARDroneTTSim_B.gain;

  /* StateSpace: '<S1>/State-Space4' */
  rtb_StateSpace4 = 0.0;
  rtb_StateSpace4 += ARDroneTTSim_P.StateSpace4_C *
    ARDroneTTSim_X.StateSpace4_CSTATE;

  /* MATLAB Function: '<S1>/normalize angle  between -pi and pi radians' */
  normalizeanglebetweenpiandp(rtb_StateSpace4,
    &ARDroneTTSim_B.sf_normalizeanglebetweenpiandpi);

  /* Gain: '<S1>/deg 2 rad2' */
  ARDroneTTSim_B.deg2rad2 = ARDroneTTSim_P.deg2rad2_Gain *
    ARDroneTTSim_B.sf_normalizeanglebetweenpiandpi.angleOut;
  if (rtmIsMajorTimeStep(ARDroneTTSim_M)) {
  }

  /* Gain: '<Root>/deg 2 rad' */
  ARDroneTTSim_B.deg2rad_p[0] = ARDroneTTSim_P.deg2rad_Gain_a *
    ARDroneTTSim_B.deg2rad1;
  ARDroneTTSim_B.deg2rad_p[1] = ARDroneTTSim_P.deg2rad_Gain_a *
    ARDroneTTSim_B.deg2rad;
  ARDroneTTSim_B.deg2rad_p[2] = ARDroneTTSim_P.deg2rad_Gain_a *
    ARDroneTTSim_B.deg2rad2;

  /* RateTransition: '<Root>/RTrans5' */
  if (rtmIsMajorTimeStep(ARDroneTTSim_M) &&
      ARDroneTTSim_M->Timing.RateInteraction.TID1_2) {
    ARDroneTTSim_DW.RTrans5_1_Buffer = ARDroneTTSim_B.deg2rad_p[0];
    ARDroneTTSim_DW.RTrans5_2_Buffer = ARDroneTTSim_B.deg2rad_p[1];
    ARDroneTTSim_DW.RTrans5_3_Buffer = ARDroneTTSim_B.deg2rad_p[2];
  }

  /* StateSpace: '<S1>/State-Space3' */
  ARDroneTTSim_B.StateSpace3 = 0.0;
  ARDroneTTSim_B.StateSpace3 += ARDroneTTSim_P.StateSpace3_C *
    ARDroneTTSim_X.StateSpace3_CSTATE;

  /* RateTransition: '<Root>/RTrans5' */
  if (rtmIsMajorTimeStep(ARDroneTTSim_M) &&
      ARDroneTTSim_M->Timing.RateInteraction.TID1_2) {
    ARDroneTTSim_DW.RTrans5_4_Buffer = ARDroneTTSim_B.StateSpace3;
  }

  /* StateSpace: '<S1>/State-Space2' */
  ARDroneTTSim_B.StateSpace2 = 0.0;
  ARDroneTTSim_B.StateSpace2 += ARDroneTTSim_P.StateSpace2_C *
    ARDroneTTSim_X.StateSpace2_CSTATE;

  /* RateTransition: '<Root>/RTrans5' incorporates:
   *  Constant: '<S19>/Constant1'
   *  Constant: '<S1>/The ARDrone sends zero for the vertical velocity.  '
   */
  if (rtmIsMajorTimeStep(ARDroneTTSim_M)) {
    if (ARDroneTTSim_M->Timing.RateInteraction.TID1_2) {
      ARDroneTTSim_DW.RTrans5_5_Buffer = ARDroneTTSim_B.StateSpace2;
      ARDroneTTSim_DW.RTrans5_6_Buffer =
        ARDroneTTSim_P.TheARDronesendszeroforthevertic;
    }

    ARDroneTTSim_B.Constant1[0] = ARDroneTTSim_P.Constant1_Value[0];
    ARDroneTTSim_B.Constant1[1] = ARDroneTTSim_P.Constant1_Value[1];
  }

  /* Integrator: '<S19>/Integrator' */
  if (ARDroneTTSim_DW.Integrator_IWORK != 0) {
    ARDroneTTSim_X.Integrator_CSTATE[0] = ARDroneTTSim_B.Constant1[0];
    ARDroneTTSim_X.Integrator_CSTATE[1] = ARDroneTTSim_B.Constant1[1];
  }

  ARDroneTTSim_B.Integrator[0] = ARDroneTTSim_X.Integrator_CSTATE[0];
  ARDroneTTSim_B.Integrator[1] = ARDroneTTSim_X.Integrator_CSTATE[1];

  /* End of Integrator: '<S19>/Integrator' */

  /* RateTransition: '<Root>/RTrans5' */
  if (rtmIsMajorTimeStep(ARDroneTTSim_M) &&
      ARDroneTTSim_M->Timing.RateInteraction.TID1_2) {
    ARDroneTTSim_DW.RTrans5_7_Buffer = ARDroneTTSim_B.Integrator[0];
    ARDroneTTSim_DW.RTrans5_8_Buffer = ARDroneTTSim_B.Integrator[1];
  }

  /* StateSpace: '<S1>/State-Space5' */
  ARDroneTTSim_B.StateSpace5 = 0.0;
  ARDroneTTSim_B.StateSpace5 += ARDroneTTSim_P.StateSpace5_C[0] *
    ARDroneTTSim_X.StateSpace5_CSTATE[0];
  ARDroneTTSim_B.StateSpace5 += ARDroneTTSim_P.StateSpace5_C[1] *
    ARDroneTTSim_X.StateSpace5_CSTATE[1];

  /* RateTransition: '<Root>/RTrans5' */
  if (rtmIsMajorTimeStep(ARDroneTTSim_M) &&
      ARDroneTTSim_M->Timing.RateInteraction.TID1_2) {
    ARDroneTTSim_DW.RTrans5_9_Buffer = ARDroneTTSim_B.StateSpace5;
  }

  /* Clock: '<Root>/Time' */
  ARDroneTTSim_B.Time = ARDroneTTSim_M->Timing.t[0];

  /* RateTransition: '<Root>/RTrans7' */
  if (rtmIsMajorTimeStep(ARDroneTTSim_M)) {
    if (ARDroneTTSim_M->Timing.RateInteraction.TID1_2) {
      ARDroneTTSim_DW.RTrans7_Buffer = ARDroneTTSim_B.Time;
    }

    /* ToWorkspace: '<S5>/To Workspace' */
    {
      double locTime = ARDroneTTSim_M->Timing.t[1];
      ;
      if (rtmIsMajorTimeStep(ARDroneTTSim_M)) {
        rt_UpdateStructLogVar((StructLogVar *)
                              ARDroneTTSim_DW.ToWorkspace_PWORK_d.LoggedData,
                              &locTime, &ARDroneTTSim_B.StateSpace5);
      }
    }
  }

  /* End of RateTransition: '<Root>/RTrans7' */

  /* Gain: '<S5>/deg 2 rad1' */
  ARDroneTTSim_B.deg2rad1_l = ARDroneTTSim_P.deg2rad1_Gain_n *
    ARDroneTTSim_B.deg2rad_p[2];
  if (rtmIsMajorTimeStep(ARDroneTTSim_M)) {
    /* ToWorkspace: '<S5>/To Workspace1' */
    {
      real_T u[1];

      {
        int32_T i;
        for (i = 0; i < 1; i++) {
          u[i] = 0.0;
        }
      }

      {
        double locTime = ARDroneTTSim_M->Timing.t[1];
        ;
        if (rtmIsMajorTimeStep(ARDroneTTSim_M)) {
          rt_UpdateStructLogVar((StructLogVar *)
                                ARDroneTTSim_DW.ToWorkspace1_PWORK.LoggedData,
                                &locTime, u);
        }
      }
    }

    /* RateTransition: '<Root>/RTrans4' incorporates:
     *  RateTransition: '<Root>/RTrans3'
     */
    if (ARDroneTTSim_M->Timing.RateInteraction.TID1_2) {
      ARDroneTTSim_B.RTrans4 = ARDroneTTSim_DW.RTrans4_Buffer0;
      ARDroneTTSim_B.RTrans3 = ARDroneTTSim_DW.RTrans3_Buffer0;
    }

    /* End of RateTransition: '<Root>/RTrans4' */

    /* Saturate: '<S1>/Saturation 1' */
    if (ARDroneTTSim_B.RTrans4 > ARDroneTTSim_P.Saturation1_UpperSat) {
      ARDroneTTSim_B.Saturation1 = ARDroneTTSim_P.Saturation1_UpperSat;
    } else if (ARDroneTTSim_B.RTrans4 < ARDroneTTSim_P.Saturation1_LowerSat) {
      ARDroneTTSim_B.Saturation1 = ARDroneTTSim_P.Saturation1_LowerSat;
    } else {
      ARDroneTTSim_B.Saturation1 = ARDroneTTSim_B.RTrans4;
    }

    /* End of Saturate: '<S1>/Saturation 1' */

    /* Saturate: '<S1>/Saturation 2' */
    if (ARDroneTTSim_B.RTrans3 > ARDroneTTSim_P.Saturation2_UpperSat) {
      ARDroneTTSim_B.Saturation2 = ARDroneTTSim_P.Saturation2_UpperSat;
    } else if (ARDroneTTSim_B.RTrans3 < ARDroneTTSim_P.Saturation2_LowerSat) {
      ARDroneTTSim_B.Saturation2 = ARDroneTTSim_P.Saturation2_LowerSat;
    } else {
      ARDroneTTSim_B.Saturation2 = ARDroneTTSim_B.RTrans3;
    }

    /* End of Saturate: '<S1>/Saturation 2' */

    /* RateTransition: '<Root>/RTrans2' incorporates:
     *  RateTransition: '<Root>/RTrans1'
     */
    if (ARDroneTTSim_M->Timing.RateInteraction.TID1_2) {
      ARDroneTTSim_B.RTrans2 = ARDroneTTSim_DW.RTrans2_Buffer0;
      ARDroneTTSim_B.RTrans1 = ARDroneTTSim_DW.RTrans1_Buffer0;
    }

    /* End of RateTransition: '<Root>/RTrans2' */

    /* Saturate: '<S1>/Saturation 3' */
    if (ARDroneTTSim_B.RTrans2 > ARDroneTTSim_P.Saturation3_UpperSat) {
      ARDroneTTSim_B.Saturation3 = ARDroneTTSim_P.Saturation3_UpperSat;
    } else if (ARDroneTTSim_B.RTrans2 < ARDroneTTSim_P.Saturation3_LowerSat) {
      ARDroneTTSim_B.Saturation3 = ARDroneTTSim_P.Saturation3_LowerSat;
    } else {
      ARDroneTTSim_B.Saturation3 = ARDroneTTSim_B.RTrans2;
    }

    /* End of Saturate: '<S1>/Saturation 3' */

    /* Saturate: '<S1>/Saturation 4' */
    if (ARDroneTTSim_B.RTrans1 > ARDroneTTSim_P.Saturation4_UpperSat) {
      ARDroneTTSim_B.Saturation4 = ARDroneTTSim_P.Saturation4_UpperSat;
    } else if (ARDroneTTSim_B.RTrans1 < ARDroneTTSim_P.Saturation4_LowerSat) {
      ARDroneTTSim_B.Saturation4 = ARDroneTTSim_P.Saturation4_LowerSat;
    } else {
      ARDroneTTSim_B.Saturation4 = ARDroneTTSim_B.RTrans1;
    }

    /* End of Saturate: '<S1>/Saturation 4' */
  }

  /* VariableTransportDelay: '<S1>/total communication time delay' incorporates:
   *  Constant: '<S1>/time delay'
   */
  {
    real_T **uBuffer = (real_T**)
      &ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[4];
    real_T simTime = ARDroneTTSim_M->Timing.t[0];
    real_T appliedDelay;
    appliedDelay = ARDroneTTSim_P.timeDelay;

    /* For variable time delay, output here */
    if (appliedDelay > ARDroneTTSim_P.totalcommunicationtimedelay_Max) {
      appliedDelay = ARDroneTTSim_P.totalcommunicationtimedelay_Max;
    }

    if (appliedDelay < 0.0) {
      /* negative delay is not supported
       *  set delay to 0
       */
      appliedDelay = 0.0;
    }

    ARDroneTTSim_B.totalcommunicationtimedelay[0] = rt_TDelayInterpolate(
      simTime - appliedDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[0],
      &ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Last[0],
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[0],
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[0],
      (ARDroneTTSim_P.totalcommunicationtimedelay_Ini[0]),
      1,
      0);
    tBuffer++;
    uBuffer++;
    appliedDelay = ARDroneTTSim_P.timeDelay;

    /* For variable time delay, output here */
    if (appliedDelay > ARDroneTTSim_P.totalcommunicationtimedelay_Max) {
      appliedDelay = ARDroneTTSim_P.totalcommunicationtimedelay_Max;
    }

    if (appliedDelay < 0.0) {
      /* negative delay is not supported
       *  set delay to 0
       */
      appliedDelay = 0.0;
    }

    ARDroneTTSim_B.totalcommunicationtimedelay[1] = rt_TDelayInterpolate(
      simTime - appliedDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[1],
      &ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Last[1],
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[1],
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[1],
      (ARDroneTTSim_P.totalcommunicationtimedelay_Ini[1]),
      1,
      0);
    tBuffer++;
    uBuffer++;
    appliedDelay = ARDroneTTSim_P.timeDelay;

    /* For variable time delay, output here */
    if (appliedDelay > ARDroneTTSim_P.totalcommunicationtimedelay_Max) {
      appliedDelay = ARDroneTTSim_P.totalcommunicationtimedelay_Max;
    }

    if (appliedDelay < 0.0) {
      /* negative delay is not supported
       *  set delay to 0
       */
      appliedDelay = 0.0;
    }

    ARDroneTTSim_B.totalcommunicationtimedelay[2] = rt_TDelayInterpolate(
      simTime - appliedDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[2],
      &ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Last[2],
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[2],
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[2],
      (ARDroneTTSim_P.totalcommunicationtimedelay_Ini[2]),
      1,
      0);
    tBuffer++;
    uBuffer++;
    appliedDelay = ARDroneTTSim_P.timeDelay;

    /* For variable time delay, output here */
    if (appliedDelay > ARDroneTTSim_P.totalcommunicationtimedelay_Max) {
      appliedDelay = ARDroneTTSim_P.totalcommunicationtimedelay_Max;
    }

    if (appliedDelay < 0.0) {
      /* negative delay is not supported
       *  set delay to 0
       */
      appliedDelay = 0.0;
    }

    ARDroneTTSim_B.totalcommunicationtimedelay[3] = rt_TDelayInterpolate(
      simTime - appliedDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[3],
      &ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Last[3],
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[3],
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[3],
      (ARDroneTTSim_P.totalcommunicationtimedelay_Ini[3]),
      1,
      0);
  }

  /* MATLAB Function: '<S19>/Velocity from vehicle body frame  to inertial NED  frame' incorporates:
   *  SignalConversion: '<S20>/TmpSignal ConversionAt SFunction Inport1'
   */
  /* MATLAB Function 'Position estimation Important:This block provides an  inaccurate extimation of position  based on  velocity information. /Position estimation/Velocity from vehicle body frame  to inertial NED  frame': '<S20>:1' */
  /* '<S20>:1:4' */
  /* '<S20>:1:5' */
  ARDroneTTSim_B.Vel_xy[0] = 0.0;
  Vel_xy_tmp_0 = cos(ARDroneTTSim_B.deg2rad_p[2]);
  ARDroneTTSim_B.Vel_xy[0] += Vel_xy_tmp_0 * ARDroneTTSim_B.StateSpace3;
  Vel_xy_tmp = sin(ARDroneTTSim_B.deg2rad_p[2]);
  ARDroneTTSim_B.Vel_xy[0] += -Vel_xy_tmp * ARDroneTTSim_B.StateSpace2;
  ARDroneTTSim_B.Vel_xy[1] = 0.0;
  ARDroneTTSim_B.Vel_xy[1] += Vel_xy_tmp * ARDroneTTSim_B.StateSpace3;
  ARDroneTTSim_B.Vel_xy[1] += Vel_xy_tmp_0 * ARDroneTTSim_B.StateSpace2;
}

/* Model update function for TID0 */
void ARDroneTTSim_update0(void)        /* Sample time: [0.0s, 0.0s] */
{
  /* Update for Integrator: '<S19>/Integrator' */
  ARDroneTTSim_DW.Integrator_IWORK = 0;

  /* Update for VariableTransportDelay: '<S1>/total communication time delay' incorporates:
   *  Constant: '<S1>/time delay'
   */
  {
    real_T **uBuffer = (real_T**)
      &ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[4];
    real_T simTime = ARDroneTTSim_M->Timing.t[0];
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[0] =
      ((ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[0] <
        (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[0]-1)) ?
       (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[0]+1) : 0);
    if (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[0] ==
        ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[0]) {
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[0] =
        ((ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[0] <
          (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[0]-1))
         ? (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[0]+1) : 0);
    }

    (*tBuffer++)[ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[0]] =
      simTime;
    (*uBuffer++)[ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[0]] =
      ARDroneTTSim_B.Saturation1;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[1] =
      ((ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[1] <
        (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[1]-1)) ?
       (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[1]+1) : 0);
    if (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[1] ==
        ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[1]) {
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[1] =
        ((ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[1] <
          (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[1]-1))
         ? (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[1]+1) : 0);
    }

    (*tBuffer++)[ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[1]] =
      simTime;
    (*uBuffer++)[ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[1]] =
      ARDroneTTSim_B.Saturation2;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[2] =
      ((ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[2] <
        (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[2]-1)) ?
       (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[2]+1) : 0);
    if (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[2] ==
        ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[2]) {
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[2] =
        ((ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[2] <
          (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[2]-1))
         ? (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[2]+1) : 0);
    }

    (*tBuffer++)[ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[2]] =
      simTime;
    (*uBuffer++)[ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[2]] =
      ARDroneTTSim_B.Saturation3;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[3] =
      ((ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[3] <
        (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[3]-1)) ?
       (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[3]+1) : 0);
    if (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[3] ==
        ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[3]) {
      ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[3] =
        ((ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[3] <
          (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[3]-1))
         ? (ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[3]+1) : 0);
    }

    (*tBuffer)[ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[3]] =
      simTime;
    (*uBuffer)[ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[3]] =
      ARDroneTTSim_B.Saturation4;

    /* when use fixed buffer, reset solver at when buffer is updated
     * to avoid output consistency fail.
     */
  }

  if (rtmIsMajorTimeStep(ARDroneTTSim_M)) {
    rt_ertODEUpdateContinuousStates(&ARDroneTTSim_M->solverInfo);
  }

  /* Update absolute time */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++ARDroneTTSim_M->Timing.clockTick0)) {
    ++ARDroneTTSim_M->Timing.clockTickH0;
  }

  ARDroneTTSim_M->Timing.t[0] = rtsiGetSolverStopTime
    (&ARDroneTTSim_M->solverInfo);

  /* Update absolute time */
  /* The "clockTick1" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick1"
   * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick1 and the high bits
   * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++ARDroneTTSim_M->Timing.clockTick1)) {
    ++ARDroneTTSim_M->Timing.clockTickH1;
  }

  ARDroneTTSim_M->Timing.t[1] = ARDroneTTSim_M->Timing.clockTick1 *
    ARDroneTTSim_M->Timing.stepSize1 + ARDroneTTSim_M->Timing.clockTickH1 *
    ARDroneTTSim_M->Timing.stepSize1 * 4294967296.0;
}

/* Derivatives for root system: '<Root>' */
void ARDroneTTSim_derivatives(void)
{
  XDot_ARDroneTTSim_T *_rtXdot;
  _rtXdot = ((XDot_ARDroneTTSim_T *) ARDroneTTSim_M->derivs);

  /* Derivatives for StateSpace: '<S1>/State-Space' */
  _rtXdot->StateSpace_CSTATE[0] = 0.0;

  /* Derivatives for StateSpace: '<S1>/State-Space1' */
  _rtXdot->StateSpace1_CSTATE[0] = 0.0;

  /* Derivatives for Integrator: '<S19>/Integrator' */
  _rtXdot->Integrator_CSTATE[0] = ARDroneTTSim_B.Vel_xy[0];

  /* Derivatives for StateSpace: '<S1>/State-Space' */
  _rtXdot->StateSpace_CSTATE[1] = 0.0;

  /* Derivatives for StateSpace: '<S1>/State-Space1' */
  _rtXdot->StateSpace1_CSTATE[1] = 0.0;

  /* Derivatives for Integrator: '<S19>/Integrator' */
  _rtXdot->Integrator_CSTATE[1] = ARDroneTTSim_B.Vel_xy[1];

  /* Derivatives for StateSpace: '<S1>/State-Space' */
  _rtXdot->StateSpace_CSTATE[0] += ARDroneTTSim_P.StateSpace_A[0] *
    ARDroneTTSim_X.StateSpace_CSTATE[0];
  _rtXdot->StateSpace_CSTATE[0] += ARDroneTTSim_P.StateSpace_A[1] *
    ARDroneTTSim_X.StateSpace_CSTATE[1];
  _rtXdot->StateSpace_CSTATE[1] += ARDroneTTSim_P.StateSpace_A[2] *
    ARDroneTTSim_X.StateSpace_CSTATE[0];
  _rtXdot->StateSpace_CSTATE[0] += ARDroneTTSim_P.StateSpace_B *
    ARDroneTTSim_B.totalcommunicationtimedelay[0];

  /* Derivatives for StateSpace: '<S1>/State-Space1' */
  _rtXdot->StateSpace1_CSTATE[0] += ARDroneTTSim_P.StateSpace1_A[0] *
    ARDroneTTSim_X.StateSpace1_CSTATE[0];
  _rtXdot->StateSpace1_CSTATE[0] += ARDroneTTSim_P.StateSpace1_A[1] *
    ARDroneTTSim_X.StateSpace1_CSTATE[1];
  _rtXdot->StateSpace1_CSTATE[1] += ARDroneTTSim_P.StateSpace1_A[2] *
    ARDroneTTSim_X.StateSpace1_CSTATE[0];
  _rtXdot->StateSpace1_CSTATE[0] += ARDroneTTSim_P.StateSpace1_B *
    ARDroneTTSim_B.totalcommunicationtimedelay[1];

  /* Derivatives for StateSpace: '<S1>/State-Space4' */
  _rtXdot->StateSpace4_CSTATE = 0.0;
  _rtXdot->StateSpace4_CSTATE += ARDroneTTSim_P.StateSpace4_A *
    ARDroneTTSim_X.StateSpace4_CSTATE;
  _rtXdot->StateSpace4_CSTATE += ARDroneTTSim_P.StateSpace4_B *
    ARDroneTTSim_B.totalcommunicationtimedelay[2];

  /* Derivatives for StateSpace: '<S1>/State-Space3' */
  _rtXdot->StateSpace3_CSTATE = 0.0;
  _rtXdot->StateSpace3_CSTATE += ARDroneTTSim_P.StateSpace3_A *
    ARDroneTTSim_X.StateSpace3_CSTATE;
  _rtXdot->StateSpace3_CSTATE += ARDroneTTSim_P.StateSpace3_B *
    ARDroneTTSim_B.gain;

  /* Derivatives for StateSpace: '<S1>/State-Space2' */
  _rtXdot->StateSpace2_CSTATE = 0.0;
  _rtXdot->StateSpace2_CSTATE += ARDroneTTSim_P.StateSpace2_A *
    ARDroneTTSim_X.StateSpace2_CSTATE;
  _rtXdot->StateSpace2_CSTATE += ARDroneTTSim_P.StateSpace2_B *
    ARDroneTTSim_B.gain1;

  /* Derivatives for StateSpace: '<S1>/State-Space5' */
  _rtXdot->StateSpace5_CSTATE[0] = 0.0;
  _rtXdot->StateSpace5_CSTATE[1] = 0.0;
  _rtXdot->StateSpace5_CSTATE[0] += ARDroneTTSim_P.StateSpace5_A[0] *
    ARDroneTTSim_X.StateSpace5_CSTATE[0];
  _rtXdot->StateSpace5_CSTATE[0] += ARDroneTTSim_P.StateSpace5_A[1] *
    ARDroneTTSim_X.StateSpace5_CSTATE[1];
  _rtXdot->StateSpace5_CSTATE[1] += ARDroneTTSim_P.StateSpace5_A[2] *
    ARDroneTTSim_X.StateSpace5_CSTATE[0];
  _rtXdot->StateSpace5_CSTATE[0] += ARDroneTTSim_P.StateSpace5_B *
    ARDroneTTSim_B.totalcommunicationtimedelay[3];

  /* Derivatives for VariableTransportDelay: '<S1>/total communication time delay' incorporates:
   *  Constant: '<S1>/time delay'
   */
  {
  }
}

/* Model output function for TID2 */
void ARDroneTTSim_output2(void)        /* Sample time: [0.03s, 0.0s] */
{
  /* local block i/o variables */
  real_T rtb_yawanglerad;
  real_T rtb_Add1[3];
  real_T rtb_Sum1;
  real_T rtb_pd[3];
  real_T rtb_dot_pd[3];
  real_T rtb_psi_d;
  real_T rtb_Add[3];
  real_T dt;
  real_T Rz[9];
  real_T dot_p[3];
  real_T u[3];
  real_T x[9];
  int32_T p1;
  int32_T p2;
  int32_T p3;
  static const real_T c[3] = { 0.0, 0.0, 9.81 };

  real_T rtb_ManualSwitch2;
  real_T p[6];
  real_T tmp;
  real_T Rz_tmp;
  real_T tmp_0;
  real_T rtb_ManualSwitch2_tmp;
  real_T p_tmp;
  real_T p_tmp_0;

  /* ManualSwitch: '<Root>/Manual Switch2' incorporates:
   *  Constant: '<Root>/commands disabled'
   *  Constant: '<Root>/commands enabled'
   */
  if (ARDroneTTSim_P.ManualSwitch2_CurrentSetting == 1) {
    rtb_ManualSwitch2 = ARDroneTTSim_P.commandsdisabled_Value;
  } else {
    rtb_ManualSwitch2 = ARDroneTTSim_P.commandsenabled_Value;
  }

  /* End of ManualSwitch: '<Root>/Manual Switch2' */

  /* RateTransition: '<Root>/RTrans5' */
  rtb_yawanglerad = ARDroneTTSim_DW.RTrans5_3_Buffer;

  /* MATLAB Function: '<Root>/MATLAB Function1' incorporates:
   *  RateTransition: '<Root>/RTrans5'
   *  RateTransition: '<Root>/RTrans7'
   */
  /* MATLAB Function 'MATLAB Function1': '<S2>:1' */
  /* '<S2>:1:21' */
  /* '<S2>:1:22' */
  if (!ARDroneTTSim_DW.previous_status_not_empty) {
    /* '<S2>:1:24' */
    ARDroneTTSim_DW.previous_status_not_empty = true;

    /* '<S2>:1:26' */
    ARDroneTTSim_DW.p0[0] = ARDroneTTSim_DW.RTrans5_7_Buffer;
    ARDroneTTSim_DW.p0[1] = ARDroneTTSim_DW.RTrans5_8_Buffer;
    ARDroneTTSim_DW.p0[2] = -1.0;
  }

  if (rtb_ManualSwitch2 == 0.0) {
    /* '<S2>:1:29' */
    if (ARDroneTTSim_DW.previous_status == 1.0) {
      /* '<S2>:1:30' */
      /* '<S2>:1:32' */
      ARDroneTTSim_DW.p0[0] = ARDroneTTSim_DW.RTrans5_7_Buffer;
      ARDroneTTSim_DW.p0[1] = ARDroneTTSim_DW.RTrans5_8_Buffer;
      ARDroneTTSim_DW.p0[2] = -1.0;
    }

    /* '<S2>:1:34' */
    /* '<S2>:1:35' */
    /* '<S2>:1:36' */
    rtb_pd[0] = ARDroneTTSim_DW.p0[0];
    rtb_dot_pd[0] = 0.0;
    rtb_pd[1] = ARDroneTTSim_DW.p0[1];
    rtb_dot_pd[1] = 0.0;
    rtb_pd[2] = ARDroneTTSim_DW.p0[2];
    rtb_dot_pd[2] = 0.0;

    /* '<S2>:1:37' */
    rtb_psi_d = 0.0;
  } else {
    if ((!ARDroneTTSim_DW.t0_not_empty) || (ARDroneTTSim_DW.previous_status ==
         0.0)) {
      /* '<S2>:1:39' */
      /* '<S2>:1:41' */
      ARDroneTTSim_DW.t0 = ARDroneTTSim_DW.RTrans7_Buffer;
      ARDroneTTSim_DW.t0_not_empty = true;

      /* '<S2>:1:42' */
      ARDroneTTSim_DW.p0[0] = ARDroneTTSim_DW.RTrans5_7_Buffer;
      ARDroneTTSim_DW.p0[1] = ARDroneTTSim_DW.RTrans5_8_Buffer;
      ARDroneTTSim_DW.p0[2] = -1.0;
    }

    /* '<S2>:1:44' */
    dt = ARDroneTTSim_DW.RTrans7_Buffer - ARDroneTTSim_DW.t0;

    /* '<S2>:1:63' */
    rtb_pd[0] = ARDroneTTSim_DW.p0[0];
    rtb_pd[1] = ARDroneTTSim_DW.p0[1];
    rtb_pd[2] = ARDroneTTSim_DW.p0[2];
    if (ARDroneTTSim_DW.RTrans7_Buffer > 20.0) {
      /* '<S2>:1:64' */
      /* '<S2>:1:65' */
      rtb_dot_pd[0] = 0.1;
      rtb_dot_pd[1] = 0.0;
      rtb_dot_pd[2] = 0.0;

      /* '<S2>:1:66' */
      rtb_pd[0] = 0.1 * dt + ARDroneTTSim_DW.p0[0];
    } else {
      /* '<S2>:1:68' */
      rtb_dot_pd[0] = 0.0;
      rtb_dot_pd[1] = 0.0;
      rtb_dot_pd[2] = 0.0;
    }

    if (ARDroneTTSim_DW.RTrans7_Buffer > 20.0 - dt) {
      /* '<S2>:1:71' */
      /* '<S2>:1:72' */
      rtb_dot_pd[0] = 0.1;
      rtb_dot_pd[1] = 0.0;
      rtb_dot_pd[2] = 0.0;
    }

    if (ARDroneTTSim_DW.RTrans7_Buffer > 10.0) {
      /* '<S2>:1:75' */
      /* '<S2>:1:76' */
      rtb_pd[2] = 0.75;
    } else {
      /* '<S2>:1:78' */
      rtb_pd[2] = 0.0;
    }

    /* '<S2>:1:80' */
    /* '<S2>:1:81' */
    rtb_psi_d = 0.0;
  }

  /* '<S2>:1:86' */
  ARDroneTTSim_DW.previous_status = rtb_ManualSwitch2;

  /* ToWorkspace: '<Root>/To Workspace' */
  {
    double locTime = ARDroneTTSim_M->Timing.t[2];
    ;
    if (rtmIsMajorTimeStep(ARDroneTTSim_M)) {
      rt_UpdateStructLogVar((StructLogVar *)
                            ARDroneTTSim_DW.ToWorkspace_PWORK.LoggedData,
                            &locTime, &rtb_pd[0]);
    }
  }

  /* RelationalOperator: '<S8>/Compare' incorporates:
   *  Constant: '<S8>/Constant'
   */
  ARDroneTTSim_B.Compare = (rtb_ManualSwitch2 <= ARDroneTTSim_P.Constant_Value);

  /* DiscreteIntegrator: '<S3>/Discrete-Time Integrator' */
  if (ARDroneTTSim_B.Compare || (ARDroneTTSim_DW.DiscreteTimeIntegrator_PrevRese
       != 0)) {
    ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[0] =
      ARDroneTTSim_P.DiscreteTimeIntegrator_IC[0];
    ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[1] =
      ARDroneTTSim_P.DiscreteTimeIntegrator_IC[1];
    ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[2] =
      ARDroneTTSim_P.DiscreteTimeIntegrator_IC[2];
    ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[3] =
      ARDroneTTSim_P.DiscreteTimeIntegrator_IC[3];
  }

  /* End of DiscreteIntegrator: '<S3>/Discrete-Time Integrator' */

  /* MATLAB Function: '<S3>/MATLAB Function' incorporates:
   *  Constant: '<S3>/Constant1'
   *  MATLAB Function: '<Root>/MATLAB Function1'
   *  MATLAB Function: '<S15>/Coordinate trnasformation  from inertial frame to body frame '
   *  RateTransition: '<Root>/RTrans5'
   *  SignalConversion: '<S2>/TmpSignal ConversionAt SFunction Inport2'
   *  Sum: '<S9>/Add'
   */
  /* MATLAB Function 'Outer-loop Controller/MATLAB Function': '<S10>:1' */
  /* '<S10>:1:24' */
  /* '<S10>:1:21' */
  /* '<S10>:1:22' */
  /* '<S10>:1:23' */
  /* '<S10>:1:61' */
  /* '<S10>:1:62' */
  tmp = sin(rtb_yawanglerad);
  Rz_tmp = cos(rtb_yawanglerad);
  Rz[0] = Rz_tmp;
  Rz[3] = -tmp;
  Rz[6] = 0.0;
  Rz[1] = tmp;
  Rz[4] = Rz_tmp;
  Rz[7] = 0.0;

  /* '<S10>:1:27' */
  Rz[2] = 0.0;
  Rz[5] = 0.0;
  Rz[8] = 1.0;

  /* '<S10>:1:28' */
  /* '<S10>:1:29' */
  /* '<S10>:1:30' */
  rtb_ManualSwitch2 = cos(ARDroneTTSim_DW.RTrans5_2_Buffer);
  x[0] = rtb_ManualSwitch2 * Rz_tmp;
  dt = cos(ARDroneTTSim_DW.RTrans5_1_Buffer);
  rtb_ManualSwitch2_tmp = sin(ARDroneTTSim_DW.RTrans5_2_Buffer);
  p_tmp = sin(ARDroneTTSim_DW.RTrans5_1_Buffer);
  p_tmp_0 = p_tmp * rtb_ManualSwitch2_tmp;
  x[3] = p_tmp_0 * Rz_tmp - dt * tmp;
  tmp_0 = dt * rtb_ManualSwitch2_tmp;
  x[6] = tmp_0 * Rz_tmp + p_tmp * tmp;
  x[1] = rtb_ManualSwitch2 * tmp;
  x[4] = p_tmp_0 * tmp + dt * Rz_tmp;
  x[7] = tmp_0 * tmp - p_tmp * Rz_tmp;
  x[2] = -rtb_ManualSwitch2_tmp;
  x[5] = p_tmp * rtb_ManualSwitch2;
  x[8] = dt * rtb_ManualSwitch2;
  for (p1 = 0; p1 < 3; p1++) {
    dot_p[p1] = x[p1 + 6] * ARDroneTTSim_DW.RTrans5_6_Buffer + (x[p1 + 3] *
      ARDroneTTSim_DW.RTrans5_5_Buffer + x[p1] *
      ARDroneTTSim_DW.RTrans5_4_Buffer);
  }

  /* '<S10>:1:31' */
  dot_p[2] = -dot_p[2];

  /* '<S10>:1:46' */
  memcpy(&x[0], &Rz[0], 9U * sizeof(real_T));
  p1 = 0;
  p2 = 3;
  p3 = 6;
  dt = fabs(tmp);
  if ((dt > fabs(Rz_tmp)) && (dt > 0.0)) {
    p1 = 3;
    p2 = 0;
    x[0] = tmp;
    x[1] = Rz_tmp;
    x[3] = Rz_tmp;
    x[4] = -sin(rtb_yawanglerad);
    x[6] = 0.0;
    x[7] = 0.0;
  }

  rtb_ManualSwitch2_tmp = x[1] / x[0];
  x[1] = rtb_ManualSwitch2_tmp;
  dt = x[2] / x[0];
  x[2] /= x[0];
  x[4] -= rtb_ManualSwitch2_tmp * x[3];
  x[5] -= dt * x[3];
  x[7] -= rtb_ManualSwitch2_tmp * x[6];
  x[8] -= dt * x[6];
  if (fabs(x[5]) > fabs(x[4])) {
    p3 = p2;
    p2 = 6;
    x[1] = dt;
    x[2] = rtb_ManualSwitch2_tmp;
    rtb_ManualSwitch2 = x[4];
    x[4] = x[5];
    x[5] = rtb_ManualSwitch2;
    rtb_ManualSwitch2 = x[7];
    x[7] = x[8];
    x[8] = rtb_ManualSwitch2;
  }

  rtb_ManualSwitch2_tmp = x[5] / x[4];
  x[8] -= rtb_ManualSwitch2_tmp * x[7];
  rtb_ManualSwitch2 = (rtb_ManualSwitch2_tmp * x[1] - x[2]) / x[8];
  dt = -(x[7] * rtb_ManualSwitch2 + x[1]) / x[4];
  Rz[p1] = ((1.0 - x[3] * dt) - x[6] * rtb_ManualSwitch2) / x[0];
  Rz[p1 + 1] = dt;
  Rz[p1 + 2] = rtb_ManualSwitch2;
  rtb_ManualSwitch2 = -rtb_ManualSwitch2_tmp / x[8];
  dt = (1.0 - x[7] * rtb_ManualSwitch2) / x[4];
  Rz[p2] = -(x[3] * dt + x[6] * rtb_ManualSwitch2) / x[0];
  Rz[p2 + 1] = dt;
  Rz[p2 + 2] = rtb_ManualSwitch2;
  rtb_ManualSwitch2 = 1.0 / x[8];
  dt = -x[7] * rtb_ManualSwitch2 / x[4];
  Rz[p3] = -(x[3] * dt + x[6] * rtb_ManualSwitch2) / x[0];
  Rz[p3 + 1] = dt;
  Rz[p3 + 2] = rtb_ManualSwitch2;
  p_tmp = ARDroneTTSim_DW.RTrans5_7_Buffer - rtb_pd[0];
  p[0] = p_tmp;
  p[3] = dot_p[0] - rtb_dot_pd[0];
  p_tmp_0 = ARDroneTTSim_DW.RTrans5_8_Buffer - rtb_pd[1];
  p[1] = p_tmp_0;
  p[4] = dot_p[1] - rtb_dot_pd[1];
  rtb_ManualSwitch2_tmp = -ARDroneTTSim_DW.RTrans5_9_Buffer - rtb_pd[2];
  p[2] = rtb_ManualSwitch2_tmp;
  p[5] = dot_p[2] - rtb_dot_pd[2];
  for (p1 = 0; p1 < 3; p1++) {
    rtb_ManualSwitch2 = 0.0;
    for (p2 = 0; p2 < 6; p2++) {
      rtb_ManualSwitch2 += ARDroneTTSim_P.K[3 * p2 + p1] * p[p2];
    }

    dot_p[p1] = rtb_ManualSwitch2 + c[p1];
  }

  for (p1 = 0; p1 < 3; p1++) {
    u[p1] = Rz[p1 + 6] * dot_p[2] + (Rz[p1 + 3] * dot_p[1] + Rz[p1] * dot_p[0]);
  }

  /* '<S10>:1:48' */
  /* '<S10>:1:49' */
  /* '<S10>:1:50' */
  /* '<S10>:1:51' */
  /* '<S10>:1:52' */
  rtb_ManualSwitch2 = atan(u[0] / u[2]);
  dt = atan(-u[1] / sqrt(u[0] * u[0] + u[2] * u[2]));
  ARDroneTTSim_B.dot_xi[0] = 0.0;
  ARDroneTTSim_B.dot_xi[1] = 0.0;
  ARDroneTTSim_B.dot_xi[2] = 0.0;
  ARDroneTTSim_B.dot_xi[3] = 0.0;

  /* Sum: '<S9>/Add1' incorporates:
   *  RateTransition: '<Root>/RTrans5'
   */
  rtb_Add1[0] = ARDroneTTSim_DW.RTrans5_1_Buffer - dt;
  rtb_Add1[1] = ARDroneTTSim_DW.RTrans5_2_Buffer - rtb_ManualSwitch2;
  rtb_Add1[2] = rtb_yawanglerad - rtb_psi_d;

  /* Sum: '<S9>/Add' incorporates:
   *  RateTransition: '<Root>/RTrans5'
   */
  rtb_Add[0] = p_tmp;
  rtb_Add[1] = p_tmp_0;
  rtb_Add[2] = ARDroneTTSim_DW.RTrans5_9_Buffer - rtb_pd[2];

  /* MATLAB Function: '<S9>/normalize angle  between -pi and pi radians' */
  normalizeanglebetweenpiandp(rtb_Add1[2],
    &ARDroneTTSim_B.sf_normalizeanglebetweenpiand_c);

  /* Scope: '<S9>/psi  psi_d' */
  if (rtmIsMajorTimeStep(ARDroneTTSim_M)) {
    StructLogVar *svar = (StructLogVar *)
      ARDroneTTSim_DW.psipsi_d_PWORK.LoggedData;
    LogVar *var = svar->signals.values;

    /* time */
    {
      double locTime = ARDroneTTSim_M->Timing.t[2];
      ;
      rt_UpdateLogVar((LogVar *)svar->time, &locTime, 0);
    }

    /* signals */
    {
      real_T up0[3];
      up0[0] = rtb_yawanglerad;
      up0[1] = rtb_psi_d;
      up0[2] = ARDroneTTSim_B.sf_normalizeanglebetweenpiand_c.angleOut;
      rt_UpdateLogVar((LogVar *)var, up0, 0);
    }
  }

  /* Sum: '<S15>/Sum4' incorporates:
   *  RateTransition: '<Root>/RTrans5'
   */
  /* MATLAB Function 'Outer-loop Controller/Baseline Controller/Position controller/Coordinate trnasformation  from inertial frame to body frame ': '<S17>:1' */
  /* '<S17>:1:3' */
  /* '<S17>:1:4' */
  p_tmp = rtb_pd[0] - ARDroneTTSim_DW.RTrans5_7_Buffer;
  p_tmp_0 = rtb_pd[1] - ARDroneTTSim_DW.RTrans5_8_Buffer;

  /* Sum: '<S12>/Sum1' */
  rtb_Sum1 = rtb_psi_d - rtb_yawanglerad;

  /* MATLAB Function: '<S12>/normalize angle  between -pi and pi radians' */
  normalizeanglebetweenpiandp(rtb_Sum1,
    &ARDroneTTSim_B.sf_normalizeanglebetweenpiand_e);

  /* ManualSwitch: '<S3>/Manual Switch3' incorporates:
   *  Constant: '<S3>/Constant'
   *  Gain: '<S11>/Gain'
   *  Gain: '<S12>/proportional control gain - yaw'
   *  Gain: '<S13>/proportional control gain'
   *  Gain: '<S14>/Gain1'
   *  Gain: '<S15>/Gain2'
   *  Gain: '<S15>/Gain3'
   *  MATLAB Function: '<S15>/Coordinate trnasformation  from inertial frame to body frame '
   *  MATLAB Function: '<S3>/MATLAB Function'
   *  RateTransition: '<Root>/RTrans5'
   *  Sum: '<S11>/Sum2'
   *  Sum: '<S13>/Sum3'
   *  Sum: '<S14>/Sum5'
   */
  if (ARDroneTTSim_P.ManualSwitch3_CurrentSetting == 1) {
    ARDroneTTSim_B.ManualSwitch3[0] = (rtb_pd[2] -
      ARDroneTTSim_DW.RTrans5_9_Buffer) *
      ARDroneTTSim_P.proportionalcontrolgain_Gain;
    ARDroneTTSim_B.ManualSwitch3[1] =
      ARDroneTTSim_P.proportionalcontrolgainyaw_Gain *
      ARDroneTTSim_B.sf_normalizeanglebetweenpiand_e.angleOut;
    ARDroneTTSim_B.ManualSwitch3[2] = ((Rz_tmp * p_tmp + tmp * p_tmp_0) *
      ARDroneTTSim_P.Gain2_Gain - ARDroneTTSim_DW.RTrans5_4_Buffer) *
      ARDroneTTSim_P.Gain_Gain;
    ARDroneTTSim_B.ManualSwitch3[3] = ((-sin(rtb_yawanglerad) * p_tmp + Rz_tmp *
      p_tmp_0) * ARDroneTTSim_P.Gain3_Gain - ARDroneTTSim_DW.RTrans5_5_Buffer) *
      ARDroneTTSim_P.Gain1_Gain;
  } else {
    ARDroneTTSim_B.ManualSwitch3[0] = rtb_ManualSwitch2_tmp * ARDroneTTSim_P.k_w;
    ARDroneTTSim_B.ManualSwitch3[1] = 0.0;
    ARDroneTTSim_B.ManualSwitch3[2] = rtb_ManualSwitch2;
    ARDroneTTSim_B.ManualSwitch3[3] = dt;
  }

  /* End of ManualSwitch: '<S3>/Manual Switch3' */
}

/* Model update function for TID2 */
void ARDroneTTSim_update2(void)        /* Sample time: [0.03s, 0.0s] */
{
  /* Update for RateTransition: '<Root>/RTrans4' */
  ARDroneTTSim_DW.RTrans4_Buffer0 = ARDroneTTSim_B.ManualSwitch3[3];

  /* Update for RateTransition: '<Root>/RTrans3' */
  ARDroneTTSim_DW.RTrans3_Buffer0 = ARDroneTTSim_B.ManualSwitch3[2];

  /* Update for RateTransition: '<Root>/RTrans2' */
  ARDroneTTSim_DW.RTrans2_Buffer0 = ARDroneTTSim_B.ManualSwitch3[1];

  /* Update for RateTransition: '<Root>/RTrans1' */
  ARDroneTTSim_DW.RTrans1_Buffer0 = ARDroneTTSim_B.ManualSwitch3[0];

  /* Update for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' */
  if (!ARDroneTTSim_B.Compare) {
    ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[0] +=
      ARDroneTTSim_P.DiscreteTimeIntegrator_gainval * ARDroneTTSim_B.dot_xi[0];
    ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[1] +=
      ARDroneTTSim_P.DiscreteTimeIntegrator_gainval * ARDroneTTSim_B.dot_xi[1];
    ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[2] +=
      ARDroneTTSim_P.DiscreteTimeIntegrator_gainval * ARDroneTTSim_B.dot_xi[2];
    ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[3] +=
      ARDroneTTSim_P.DiscreteTimeIntegrator_gainval * ARDroneTTSim_B.dot_xi[3];
  }

  ARDroneTTSim_DW.DiscreteTimeIntegrator_PrevRese = (int8_T)
    ARDroneTTSim_B.Compare;

  /* End of Update for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' */

  /* Update absolute time */
  /* The "clockTick2" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick2"
   * and "Timing.stepSize2". Size of "clockTick2" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick2 and the high bits
   * Timing.clockTickH2. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++ARDroneTTSim_M->Timing.clockTick2)) {
    ++ARDroneTTSim_M->Timing.clockTickH2;
  }

  ARDroneTTSim_M->Timing.t[2] = ARDroneTTSim_M->Timing.clockTick2 *
    ARDroneTTSim_M->Timing.stepSize2 + ARDroneTTSim_M->Timing.clockTickH2 *
    ARDroneTTSim_M->Timing.stepSize2 * 4294967296.0;
}

/* Model output wrapper function for compatibility with a static main program */
void ARDroneTTSim_output(int_T tid)
{
  switch (tid) {
   case 0 :
    ARDroneTTSim_output0();
    break;

   case 2 :
    ARDroneTTSim_output2();
    break;

   default :
    break;
  }
}

/* Model update wrapper function for compatibility with a static main program */
void ARDroneTTSim_update(int_T tid)
{
  switch (tid) {
   case 0 :
    ARDroneTTSim_update0();
    break;

   case 2 :
    ARDroneTTSim_update2();
    break;

   default :
    break;
  }
}

/* Model initialize function */
void ARDroneTTSim_initialize(void)
{
  /* SetupRuntimeResources for ToWorkspace: '<Root>/To Workspace' */
  {
    static int_T rt_ToWksWidths[] = { 3 };

    static int_T rt_ToWksNumDimensions[] = { 2 };

    static int_T rt_ToWksDimensions[] = { 3, 1 };

    static boolean_T rt_ToWksIsVarDims[] = { 0 };

    static void *rt_ToWksCurrSigDims[] = { (NULL), (NULL) };

    static int_T rt_ToWksCurrSigDimsSize[] = { 4, 4 };

    static BuiltInDTypeId rt_ToWksDataTypeIds[] = { SS_DOUBLE };

    static int_T rt_ToWksComplexSignals[] = { 0 };

    static int_T rt_ToWksFrameData[] = { 0 };

    static RTWPreprocessingFcnPtr rt_ToWksLoggingPreprocessingFcnPtrs[] = {
      (NULL)
    };

    static const char_T *rt_ToWksLabels[] = { "" };

    static RTWLogSignalInfo rt_ToWksSignalInfo = {
      1,
      rt_ToWksWidths,
      rt_ToWksNumDimensions,
      rt_ToWksDimensions,
      rt_ToWksIsVarDims,
      rt_ToWksCurrSigDims,
      rt_ToWksCurrSigDimsSize,
      rt_ToWksDataTypeIds,
      rt_ToWksComplexSignals,
      rt_ToWksFrameData,
      rt_ToWksLoggingPreprocessingFcnPtrs,

      { rt_ToWksLabels },
      (NULL),
      (NULL),
      (NULL),

      { (NULL) },

      { (NULL) },
      (NULL),
      (NULL)
    };

    static const char_T rt_ToWksBlockName[] = "ARDroneTTSim/To Workspace";
    ARDroneTTSim_DW.ToWorkspace_PWORK.LoggedData = rt_CreateStructLogVar(
      ARDroneTTSim_M->rtwLogInfo,
      0.0,
      rtmGetTFinal(ARDroneTTSim_M),
      ARDroneTTSim_M->Timing.stepSize0,
      (&rtmGetErrorStatus(ARDroneTTSim_M)),
      "pd",
      1,
      0,
      1,
      0.03,
      &rt_ToWksSignalInfo,
      rt_ToWksBlockName);
    if (ARDroneTTSim_DW.ToWorkspace_PWORK.LoggedData == (NULL))
      return;
  }

  /* SetupRuntimeResources for ToWorkspace: '<S5>/To Workspace' */
  {
    static int_T rt_ToWksWidths[] = { 1 };

    static int_T rt_ToWksNumDimensions[] = { 1 };

    static int_T rt_ToWksDimensions[] = { 1 };

    static boolean_T rt_ToWksIsVarDims[] = { 0 };

    static void *rt_ToWksCurrSigDims[] = { (NULL) };

    static int_T rt_ToWksCurrSigDimsSize[] = { 4 };

    static BuiltInDTypeId rt_ToWksDataTypeIds[] = { SS_DOUBLE };

    static int_T rt_ToWksComplexSignals[] = { 0 };

    static int_T rt_ToWksFrameData[] = { 0 };

    static RTWPreprocessingFcnPtr rt_ToWksLoggingPreprocessingFcnPtrs[] = {
      (NULL)
    };

    static const char_T *rt_ToWksLabels[] = { "<h (m)>" };

    static RTWLogSignalInfo rt_ToWksSignalInfo = {
      1,
      rt_ToWksWidths,
      rt_ToWksNumDimensions,
      rt_ToWksDimensions,
      rt_ToWksIsVarDims,
      rt_ToWksCurrSigDims,
      rt_ToWksCurrSigDimsSize,
      rt_ToWksDataTypeIds,
      rt_ToWksComplexSignals,
      rt_ToWksFrameData,
      rt_ToWksLoggingPreprocessingFcnPtrs,

      { rt_ToWksLabels },
      (NULL),
      (NULL),
      (NULL),

      { (NULL) },

      { (NULL) },
      (NULL),
      (NULL)
    };

    static const char_T rt_ToWksBlockName[] =
      "ARDroneTTSim/Visualization of Drone states/To Workspace";
    ARDroneTTSim_DW.ToWorkspace_PWORK_d.LoggedData = rt_CreateStructLogVar(
      ARDroneTTSim_M->rtwLogInfo,
      0.0,
      rtmGetTFinal(ARDroneTTSim_M),
      ARDroneTTSim_M->Timing.stepSize0,
      (&rtmGetErrorStatus(ARDroneTTSim_M)),
      "h_sim",
      1,
      0,
      1,
      0.005,
      &rt_ToWksSignalInfo,
      rt_ToWksBlockName);
    if (ARDroneTTSim_DW.ToWorkspace_PWORK_d.LoggedData == (NULL))
      return;
  }

  /* SetupRuntimeResources for ToWorkspace: '<S5>/To Workspace1' */
  {
    static int_T rt_ToWksWidths[] = { 1 };

    static int_T rt_ToWksNumDimensions[] = { 1 };

    static int_T rt_ToWksDimensions[] = { 1 };

    static boolean_T rt_ToWksIsVarDims[] = { 0 };

    static void *rt_ToWksCurrSigDims[] = { (NULL) };

    static int_T rt_ToWksCurrSigDimsSize[] = { 4 };

    static BuiltInDTypeId rt_ToWksDataTypeIds[] = { SS_DOUBLE };

    static int_T rt_ToWksComplexSignals[] = { 0 };

    static int_T rt_ToWksFrameData[] = { 0 };

    static RTWPreprocessingFcnPtr rt_ToWksLoggingPreprocessingFcnPtrs[] = {
      (NULL)
    };

    static const char_T *rt_ToWksLabels[] = { "" };

    static RTWLogSignalInfo rt_ToWksSignalInfo = {
      1,
      rt_ToWksWidths,
      rt_ToWksNumDimensions,
      rt_ToWksDimensions,
      rt_ToWksIsVarDims,
      rt_ToWksCurrSigDims,
      rt_ToWksCurrSigDimsSize,
      rt_ToWksDataTypeIds,
      rt_ToWksComplexSignals,
      rt_ToWksFrameData,
      rt_ToWksLoggingPreprocessingFcnPtrs,

      { rt_ToWksLabels },
      (NULL),
      (NULL),
      (NULL),

      { (NULL) },

      { (NULL) },
      (NULL),
      (NULL)
    };

    static const char_T rt_ToWksBlockName[] =
      "ARDroneTTSim/Visualization of Drone states/To Workspace1";
    ARDroneTTSim_DW.ToWorkspace1_PWORK.LoggedData = rt_CreateStructLogVar(
      ARDroneTTSim_M->rtwLogInfo,
      0.0,
      rtmGetTFinal(ARDroneTTSim_M),
      ARDroneTTSim_M->Timing.stepSize0,
      (&rtmGetErrorStatus(ARDroneTTSim_M)),
      "h_sim1",
      1,
      0,
      1,
      0.005,
      &rt_ToWksSignalInfo,
      rt_ToWksBlockName);
    if (ARDroneTTSim_DW.ToWorkspace1_PWORK.LoggedData == (NULL))
      return;
  }

  /* Start for Constant: '<S19>/Constant1' */
  ARDroneTTSim_B.Constant1[0] = ARDroneTTSim_P.Constant1_Value[0];
  ARDroneTTSim_B.Constant1[1] = ARDroneTTSim_P.Constant1_Value[1];

  /* Start for RateTransition: '<Root>/RTrans4' */
  ARDroneTTSim_B.RTrans4 = ARDroneTTSim_P.RTrans4_InitialCondition;

  /* Start for RateTransition: '<Root>/RTrans3' */
  ARDroneTTSim_B.RTrans3 = ARDroneTTSim_P.RTrans3_InitialCondition;

  /* Start for RateTransition: '<Root>/RTrans2' */
  ARDroneTTSim_B.RTrans2 = ARDroneTTSim_P.RTrans2_InitialCondition;

  /* Start for RateTransition: '<Root>/RTrans1' */
  ARDroneTTSim_B.RTrans1 = ARDroneTTSim_P.RTrans1_InitialCondition;

  /* Start for VariableTransportDelay: '<S1>/total communication time delay' incorporates:
   *  Constant: '<S1>/time delay'
   */
  {
    real_T *pBuffer =
      &ARDroneTTSim_DW.totalcommunicationtimedelay_RWO.TUbufferArea[0];
    int_T j;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[0] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[0] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Last[0] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[0] = 4096;
    for (j=0; j < 4096; j++) {
      pBuffer[j] = (ARDroneTTSim_P.totalcommunicationtimedelay_Ini[0]);
      pBuffer[4096 + j] = ARDroneTTSim_M->Timing.t[0];
    }

    ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[0] = (void *)
      &pBuffer[0];
    ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[4] = (void *)
      &pBuffer[4096];
    pBuffer += 8192;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[1] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[1] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Last[1] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[1] = 4096;
    for (j=0; j < 4096; j++) {
      pBuffer[j] = (ARDroneTTSim_P.totalcommunicationtimedelay_Ini[1]);
      pBuffer[4096 + j] = ARDroneTTSim_M->Timing.t[0];
    }

    ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[1] = (void *)
      &pBuffer[0];
    ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[5] = (void *)
      &pBuffer[4096];
    pBuffer += 8192;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[2] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[2] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Last[2] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[2] = 4096;
    for (j=0; j < 4096; j++) {
      pBuffer[j] = (ARDroneTTSim_P.totalcommunicationtimedelay_Ini[2]);
      pBuffer[4096 + j] = ARDroneTTSim_M->Timing.t[0];
    }

    ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[2] = (void *)
      &pBuffer[0];
    ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[6] = (void *)
      &pBuffer[4096];
    pBuffer += 8192;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Tail[3] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Head[3] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.Last[3] = 0;
    ARDroneTTSim_DW.totalcommunicationtimedelay_IWO.CircularBufSize[3] = 4096;
    for (j=0; j < 4096; j++) {
      pBuffer[j] = (ARDroneTTSim_P.totalcommunicationtimedelay_Ini[3]);
      pBuffer[4096 + j] = ARDroneTTSim_M->Timing.t[0];
    }

    ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[3] = (void *)
      &pBuffer[0];
    ARDroneTTSim_DW.totalcommunicationtimedelay_PWO.TUbufferPtrs[7] = (void *)
      &pBuffer[4096];
  }

  /* Start for Scope: '<S9>/psi  psi_d' */
  {
    RTWLogSignalInfo rt_ScopeSignalInfo;
    static int_T rt_ScopeSignalWidths[] = { 3 };

    static int_T rt_ScopeSignalNumDimensions[] = { 1 };

    static int_T rt_ScopeSignalDimensions[] = { 3 };

    static void *rt_ScopeCurrSigDims[] = { (NULL) };

    static int_T rt_ScopeCurrSigDimsSize[] = { 4 };

    static const char_T *rt_ScopeSignalLabels[] = { "" };

    static char_T rt_ScopeSignalTitles[] = "";
    static int_T rt_ScopeSignalTitleLengths[] = { 0 };

    static boolean_T rt_ScopeSignalIsVarDims[] = { 0 };

    static int_T rt_ScopeSignalPlotStyles[] = { 1, 1, 1 };

    BuiltInDTypeId dTypes[1] = { SS_DOUBLE };

    static char_T rt_ScopeBlockName[] =
      "ARDroneTTSim/Outer-loop Controller/Error Scopes/psi \npsi_d";
    static int_T rt_ScopeFrameData[] = { 0 };

    static RTWPreprocessingFcnPtr rt_ScopeSignalLoggingPreprocessingFcnPtrs[] =
      {
      (NULL)
    };

    rt_ScopeSignalInfo.numSignals = 1;
    rt_ScopeSignalInfo.numCols = rt_ScopeSignalWidths;
    rt_ScopeSignalInfo.numDims = rt_ScopeSignalNumDimensions;
    rt_ScopeSignalInfo.dims = rt_ScopeSignalDimensions;
    rt_ScopeSignalInfo.isVarDims = rt_ScopeSignalIsVarDims;
    rt_ScopeSignalInfo.currSigDims = rt_ScopeCurrSigDims;
    rt_ScopeSignalInfo.currSigDimsSize = rt_ScopeCurrSigDimsSize;
    rt_ScopeSignalInfo.dataTypes = dTypes;
    rt_ScopeSignalInfo.complexSignals = (NULL);
    rt_ScopeSignalInfo.frameData = rt_ScopeFrameData;
    rt_ScopeSignalInfo.preprocessingPtrs =
      rt_ScopeSignalLoggingPreprocessingFcnPtrs;
    rt_ScopeSignalInfo.labels.cptr = rt_ScopeSignalLabels;
    rt_ScopeSignalInfo.titles = rt_ScopeSignalTitles;
    rt_ScopeSignalInfo.titleLengths = rt_ScopeSignalTitleLengths;
    rt_ScopeSignalInfo.plotStyles = rt_ScopeSignalPlotStyles;
    rt_ScopeSignalInfo.blockNames.cptr = (NULL);
    rt_ScopeSignalInfo.stateNames.cptr = (NULL);
    rt_ScopeSignalInfo.crossMdlRef = (NULL);
    rt_ScopeSignalInfo.dataTypeConvert = (NULL);
    ARDroneTTSim_DW.psipsi_d_PWORK.LoggedData = rt_CreateStructLogVar(
      ARDroneTTSim_M->rtwLogInfo,
      0.0,
      rtmGetTFinal(ARDroneTTSim_M),
      ARDroneTTSim_M->Timing.stepSize0,
      (&rtmGetErrorStatus(ARDroneTTSim_M)),
      "psi",
      1,
      500000,
      1,
      0.03,
      &rt_ScopeSignalInfo,
      rt_ScopeBlockName);
    if (ARDroneTTSim_DW.psipsi_d_PWORK.LoggedData == (NULL))
      return;
  }

  /* InitializeConditions for StateSpace: '<S1>/State-Space' */
  ARDroneTTSim_X.StateSpace_CSTATE[0] =
    ARDroneTTSim_P.StateSpace_InitialCondition;

  /* InitializeConditions for StateSpace: '<S1>/State-Space1' */
  ARDroneTTSim_X.StateSpace1_CSTATE[0] =
    ARDroneTTSim_P.StateSpace1_InitialCondition;

  /* InitializeConditions for StateSpace: '<S1>/State-Space' */
  ARDroneTTSim_X.StateSpace_CSTATE[1] =
    ARDroneTTSim_P.StateSpace_InitialCondition;

  /* InitializeConditions for StateSpace: '<S1>/State-Space1' */
  ARDroneTTSim_X.StateSpace1_CSTATE[1] =
    ARDroneTTSim_P.StateSpace1_InitialCondition;

  /* InitializeConditions for StateSpace: '<S1>/State-Space4' */
  ARDroneTTSim_X.StateSpace4_CSTATE =
    ARDroneTTSim_P.StateSpace4_InitialCondition;

  /* InitializeConditions for StateSpace: '<S1>/State-Space3' */
  ARDroneTTSim_X.StateSpace3_CSTATE =
    ARDroneTTSim_P.StateSpace3_InitialCondition;

  /* InitializeConditions for StateSpace: '<S1>/State-Space2' */
  ARDroneTTSim_X.StateSpace2_CSTATE =
    ARDroneTTSim_P.StateSpace2_InitialCondition;

  /* InitializeConditions for Integrator: '<S19>/Integrator' */
  if (rtmIsFirstInitCond(ARDroneTTSim_M)) {
    ARDroneTTSim_X.Integrator_CSTATE[0] = 0.0;
    ARDroneTTSim_X.Integrator_CSTATE[1] = 0.0;
  }

  ARDroneTTSim_DW.Integrator_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S19>/Integrator' */

  /* InitializeConditions for StateSpace: '<S1>/State-Space5' */
  ARDroneTTSim_X.StateSpace5_CSTATE[0] =
    ARDroneTTSim_P.StateSpace5_InitialCondition;
  ARDroneTTSim_X.StateSpace5_CSTATE[1] =
    ARDroneTTSim_P.StateSpace5_InitialCondition;

  /* InitializeConditions for RateTransition: '<Root>/RTrans4' */
  ARDroneTTSim_DW.RTrans4_Buffer0 = ARDroneTTSim_P.RTrans4_InitialCondition;

  /* InitializeConditions for RateTransition: '<Root>/RTrans3' */
  ARDroneTTSim_DW.RTrans3_Buffer0 = ARDroneTTSim_P.RTrans3_InitialCondition;

  /* InitializeConditions for RateTransition: '<Root>/RTrans2' */
  ARDroneTTSim_DW.RTrans2_Buffer0 = ARDroneTTSim_P.RTrans2_InitialCondition;

  /* InitializeConditions for RateTransition: '<Root>/RTrans1' */
  ARDroneTTSim_DW.RTrans1_Buffer0 = ARDroneTTSim_P.RTrans1_InitialCondition;

  /* InitializeConditions for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' */
  ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[0] =
    ARDroneTTSim_P.DiscreteTimeIntegrator_IC[0];
  ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[1] =
    ARDroneTTSim_P.DiscreteTimeIntegrator_IC[1];
  ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[2] =
    ARDroneTTSim_P.DiscreteTimeIntegrator_IC[2];
  ARDroneTTSim_DW.DiscreteTimeIntegrator_DSTATE[3] =
    ARDroneTTSim_P.DiscreteTimeIntegrator_IC[3];
  ARDroneTTSim_DW.DiscreteTimeIntegrator_PrevRese = 0;

  /* SystemInitialize for MATLAB Function: '<Root>/MATLAB Function1' */
  ARDroneTTSim_DW.t0_not_empty = false;
  ARDroneTTSim_DW.previous_status_not_empty = false;
  ARDroneTTSim_DW.previous_status = 1.0;

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(ARDroneTTSim_M)) {
    rtmSetFirstInitCond(ARDroneTTSim_M, 0);
  }
}

/* Model terminate function */
void ARDroneTTSim_terminate(void)
{
  /* (no terminate code required) */
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  if (tid == 1)
    tid = 0;
  ARDroneTTSim_output(tid);
}

void MdlUpdate(int_T tid)
{
  if (tid == 1)
    tid = 0;
  ARDroneTTSim_update(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  ARDroneTTSim_initialize();
}

void MdlTerminate(void)
{
  ARDroneTTSim_terminate();
}

/* Registration function */
RT_MODEL_ARDroneTTSim_T *ARDroneTTSim(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)ARDroneTTSim_M, 0,
                sizeof(RT_MODEL_ARDroneTTSim_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&ARDroneTTSim_M->solverInfo,
                          &ARDroneTTSim_M->Timing.simTimeStep);
    rtsiSetTPtr(&ARDroneTTSim_M->solverInfo, &rtmGetTPtr(ARDroneTTSim_M));
    rtsiSetStepSizePtr(&ARDroneTTSim_M->solverInfo,
                       &ARDroneTTSim_M->Timing.stepSize0);
    rtsiSetdXPtr(&ARDroneTTSim_M->solverInfo, &ARDroneTTSim_M->derivs);
    rtsiSetContStatesPtr(&ARDroneTTSim_M->solverInfo, (real_T **)
                         &ARDroneTTSim_M->contStates);
    rtsiSetNumContStatesPtr(&ARDroneTTSim_M->solverInfo,
      &ARDroneTTSim_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&ARDroneTTSim_M->solverInfo,
      &ARDroneTTSim_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&ARDroneTTSim_M->solverInfo,
      &ARDroneTTSim_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&ARDroneTTSim_M->solverInfo,
      &ARDroneTTSim_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&ARDroneTTSim_M->solverInfo, (&rtmGetErrorStatus
      (ARDroneTTSim_M)));
    rtsiSetRTModelPtr(&ARDroneTTSim_M->solverInfo, ARDroneTTSim_M);
  }

  rtsiSetSimTimeStep(&ARDroneTTSim_M->solverInfo, MAJOR_TIME_STEP);
  ARDroneTTSim_M->intgData.f[0] = ARDroneTTSim_M->odeF[0];
  ARDroneTTSim_M->contStates = ((real_T *) &ARDroneTTSim_X);
  rtsiSetSolverData(&ARDroneTTSim_M->solverInfo, (void *)
                    &ARDroneTTSim_M->intgData);
  rtsiSetSolverName(&ARDroneTTSim_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = ARDroneTTSim_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    mdlTsMap[2] = 2;
    ARDroneTTSim_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    ARDroneTTSim_M->Timing.sampleTimes =
      (&ARDroneTTSim_M->Timing.sampleTimesArray[0]);
    ARDroneTTSim_M->Timing.offsetTimes =
      (&ARDroneTTSim_M->Timing.offsetTimesArray[0]);

    /* task periods */
    ARDroneTTSim_M->Timing.sampleTimes[0] = (0.0);
    ARDroneTTSim_M->Timing.sampleTimes[1] = (0.005);
    ARDroneTTSim_M->Timing.sampleTimes[2] = (0.03);

    /* task offsets */
    ARDroneTTSim_M->Timing.offsetTimes[0] = (0.0);
    ARDroneTTSim_M->Timing.offsetTimes[1] = (0.0);
    ARDroneTTSim_M->Timing.offsetTimes[2] = (0.0);
  }

  rtmSetTPtr(ARDroneTTSim_M, &ARDroneTTSim_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = ARDroneTTSim_M->Timing.sampleHitArray;
    int_T *mdlPerTaskSampleHits = ARDroneTTSim_M->Timing.perTaskSampleHitsArray;
    ARDroneTTSim_M->Timing.perTaskSampleHits = (&mdlPerTaskSampleHits[0]);
    mdlSampleHits[0] = 1;
    ARDroneTTSim_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(ARDroneTTSim_M, 200.0);
  ARDroneTTSim_M->Timing.stepSize0 = 0.005;
  ARDroneTTSim_M->Timing.stepSize1 = 0.005;
  ARDroneTTSim_M->Timing.stepSize2 = 0.03;
  rtmSetFirstInitCond(ARDroneTTSim_M, 1);

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = NULL;
    ARDroneTTSim_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(ARDroneTTSim_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(ARDroneTTSim_M->rtwLogInfo, (NULL));
    rtliSetLogT(ARDroneTTSim_M->rtwLogInfo, "tout");
    rtliSetLogX(ARDroneTTSim_M->rtwLogInfo, "");
    rtliSetLogXFinal(ARDroneTTSim_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(ARDroneTTSim_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(ARDroneTTSim_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(ARDroneTTSim_M->rtwLogInfo, 0);
    rtliSetLogDecimation(ARDroneTTSim_M->rtwLogInfo, 1);
    rtliSetLogY(ARDroneTTSim_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(ARDroneTTSim_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(ARDroneTTSim_M->rtwLogInfo, (NULL));
  }

  ARDroneTTSim_M->solverInfoPtr = (&ARDroneTTSim_M->solverInfo);
  ARDroneTTSim_M->Timing.stepSize = (0.005);
  rtsiSetFixedStepSize(&ARDroneTTSim_M->solverInfo, 0.005);
  rtsiSetSolverMode(&ARDroneTTSim_M->solverInfo, SOLVER_MODE_MULTITASKING);

  /* block I/O */
  ARDroneTTSim_M->blockIO = ((void *) &ARDroneTTSim_B);
  (void) memset(((void *) &ARDroneTTSim_B), 0,
                sizeof(B_ARDroneTTSim_T));

  /* parameters */
  ARDroneTTSim_M->defaultParam = ((real_T *)&ARDroneTTSim_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &ARDroneTTSim_X;
    ARDroneTTSim_M->contStates = (x);
    (void) memset((void *)&ARDroneTTSim_X, 0,
                  sizeof(X_ARDroneTTSim_T));
  }

  /* states (dwork) */
  ARDroneTTSim_M->dwork = ((void *) &ARDroneTTSim_DW);
  (void) memset((void *)&ARDroneTTSim_DW, 0,
                sizeof(DW_ARDroneTTSim_T));

  /* Initialize Sizes */
  ARDroneTTSim_M->Sizes.numContStates = (11);/* Number of continuous states */
  ARDroneTTSim_M->Sizes.numPeriodicContStates = (0);
                                      /* Number of periodic continuous states */
  ARDroneTTSim_M->Sizes.numY = (0);    /* Number of model outputs */
  ARDroneTTSim_M->Sizes.numU = (0);    /* Number of model inputs */
  ARDroneTTSim_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  ARDroneTTSim_M->Sizes.numSampTimes = (3);/* Number of sample times */
  ARDroneTTSim_M->Sizes.numBlocks = (93);/* Number of blocks */
  ARDroneTTSim_M->Sizes.numBlockIO = (29);/* Number of block outputs */
  ARDroneTTSim_M->Sizes.numBlockPrms = (96);/* Sum of parameter "widths" */
  return ARDroneTTSim_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/

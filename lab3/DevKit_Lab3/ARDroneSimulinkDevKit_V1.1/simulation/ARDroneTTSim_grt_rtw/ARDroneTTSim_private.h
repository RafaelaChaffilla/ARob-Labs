/*
 * ARDroneTTSim_private.h
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

#ifndef RTW_HEADER_ARDroneTTSim_private_h_
#define RTW_HEADER_ARDroneTTSim_private_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#include "zero_crossing_types.h"
#include "ARDroneTTSim.h"

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
  ;
extern void normalizeanglebetweenpiandp(real_T rtu_angleIn,
  B_normalizeanglebetweenpiandp_T *localB);
void ARDroneTTSim_output0(void);
void ARDroneTTSim_update0(void);
void ARDroneTTSim_output2(void);
void ARDroneTTSim_update2(void);       /* private model entry point functions */
extern void ARDroneTTSim_derivatives(void);

#endif                                 /* RTW_HEADER_ARDroneTTSim_private_h_ */
